from pathlib import Path
import argparse
from datetime import date
import time
import pandas as pd
import numpy as np
import sys
import re
import bioframe as bf

"""
Author     : GA
Description: filters out duplicate regions in catalog
             based on coordinate size, purity, frequency, & motif size
Outputs    : merged bed file based on format, 
             merged tsv extra columns containing merge information:
                merged_from : original locus IDs merged into this record
                other_motifs: motifs of merged loci not chosen as representative
             column-reduced input catalog with new column containing unique ids
                _orig_id: unique int identifier (match up with merged_from in merged tv)
Example    : python filter_duplicates <catalog_path> <output dir> --format atarva
"""

def mergeKey(row):
    return (
        -(row["end_1based"] - row["start_0based"]), # grab the longest loci
        -getSafe(row, "ReferenceRepeatPurity", 0), # then get the nighest purity
        -convertFreq(getSafe(row, "AlleleFrequenciesFromIllumina174k")), # higher population frequency
        getSafe(row, "MotifSize", 0), # smaller motif size
        row["_tiebreak"] # if all above are tied, then choose random
    )


def getSafe(row, col, default=""):
    val = row.get(col, default)
    return default if pd.isna(val) else val


def convertFreq(row):

    if row:
        freq_list = row.split(",")

        freq_dict = {}
        for freq in freq_list:
            fs = freq.split(":")
            freq_dict[fs[0]] = int(fs[1])

        chosen_freq = max(freq_dict, key=freq_dict.get)
        chosen_freq = int(re.sub(r"\D", "", chosen_freq)) # strip all non-digits and convert to int

    else:
        chosen_freq = 0

    return chosen_freq


def moveColToEnd(df, col):
    return df[[c for c in df.columns if c != col] + [col]]


def ToAtarva(cdf: pd.DataFrame, out_path: str | Path):
    cols = [*cdf.columns[:3], "ReferenceMotif", "MotifSize"]

    fdf = cdf[cols].copy()

    fdf.columns = ["#CHROM", "START", "END", "MOTIF", "MOTIF_LEN"]
    fdf.to_csv(out_path, index=False, sep="\t", compression="gzip")

    return fdf.shape


def getStats(cdf):

    bins = [1, 10, 30, 60, 90, 100]
    labels = [r"1% - 10%",
              r"10% - 30%",
              r"30% - 60%",
              r"60% - 90%",
              r"90% - 100%",
              r"100%"]
    cols = tuple(cdf.columns[:3]) # get explicit cdf column names


    # get the overlaps of the catalog with itself:
    overlaps = bf.overlap(cdf, cdf, cols1=cols, cols2=cols, suffixes=("_1", "_2"), return_index=True)

    # drop self-matches
    overlaps = overlaps[overlaps["index_1"] != overlaps["index_2"]]

    # drop duplicate direction (A-B and B-A are the same pair)
    overlaps = overlaps[overlaps["index_1"] < overlaps["index_2"]]


    # compute overlap length
    overlaps["overlap_len"] = (
        overlaps[["end_1based_1", "end_1based_2"]].min(axis=1) - overlaps[["start_0based_1", "start_0based_2"]].max(axis=1)
    )

    # percent overlap
    overlaps["locus_len"] = overlaps["end_1based_1"] - overlaps["start_0based_1"]  # get length to use as denominator
    overlaps["pct_overlap"] = 100 * overlaps["overlap_len"] / overlaps["locus_len"]

    # define bins
    exact_100 = overlaps["pct_overlap"] == 100
    overlaps["overlap_bin"] = pd.cut(
                        overlaps["pct_overlap"], 
                        bins=bins, 
                        labels=labels[:-1],
                        include_lowest=True, 
                        right=False
                    ).astype(object) 
    overlaps.loc[exact_100, "overlap_bin"] = "100%"

    st_dict = {label: 0 for label in labels}
    counts = overlaps["overlap_bin"].value_counts()
    for label in labels:
        st_dict[label] = int(counts.get(label, 0))

    return st_dict


def mergeOverlaps(cdf: pd.DataFrame, min_dist = 0, seed=42):
    # set seed for tiebeaker randomness
    rng = np.random.default_rng(seed=seed)  
    cdf["_tiebreak"] = rng.random(len(cdf))
    coord_cols = tuple(cdf.columns[:3])

    # run bioframe overlap clustering
    clustered = bf.cluster(cdf, cols=coord_cols, min_dist=min_dist)

    # collect all member ids per cluster
    merged_ids = (
        clustered.groupby("cluster")["_orig_id"]
        .apply(list)
        .rename("merged_from")
    )

    # apply mergeKey criteria to determine rank within each cluster
    clustered["_rank"] = clustered.apply(mergeKey, axis=1)
    clustered = clustered.sort_values(["cluster", "_rank"])
    rep_df = clustered.groupby("cluster", as_index=False).first()

    # add a column containing all of the OTHER motifs from the different merged regions
    winner_idx = clustered.groupby("cluster")["_rank"].idxmin()
    is_win = clustered.index.isin(winner_idx)
    merged_motifs = (
        clustered.loc[~is_win]
        .groupby("cluster")["ReferenceMotif"]
        .agg(list)
        .rename("other_motifs")
    )

    # attach the two additional merge info cols to the rep_df
    rep_df = rep_df.merge(merged_ids, on="cluster", how="left")
    rep_df = rep_df.merge(merged_motifs, on="cluster", how="left")

    # change start and end coords to cluster coords
    rep_df[coord_cols[1]] = rep_df["cluster_start"]
    rep_df[coord_cols[2]] = rep_df["cluster_end"]

    # remove tiebeaker column from og df
    cdf.drop(columns=["_tiebreak"], inplace=True)

    # prep rep_df for return
    rep_df = rep_df.drop(columns=["cluster_start", "cluster_end", "_rank"])
    rep_df = moveColToEnd(rep_df, "cluster")

    return rep_df


def removeDupes(cdf):
    clean_df = (
        cdf.assign(_rank=cdf.apply(mergeKey, axis=1))
            .sort_values(by="_rank", ascending=True)
            .drop_duplicates(subset="ReferenceRegion", keep="first")
            .sort_values("_orig_id")
            .drop(columns=["_orig_id", "_rank"])
    )

    # double check to see if any duplicates remain
    num_dupes_cleaned = clean_df[clean_df.duplicated(subset= ["ReferenceRegion"], keep=False)].shape[0]

    if num_dupes_cleaned > 0:
        print(f"WARNING {num_dupes_cleaned} duplicates still detected.")
        sys.exit(1)

    return clean_df


def main():
    stime = time.perf_counter()
    dt = date.today().strftime("%Y%m%d")

    parser = argparse.ArgumentParser()
    parser.add_argument("cat_path", type=Path)
    parser.add_argument("out_dir", type=Path)
    parser.add_argument("--stat_only", 
                        action="store_true", 
                        default=False)
    parser.add_argument("--format",
                        type=str,
                        default="same",
                        required=False)
    parser.add_argument("--rem_dupes", 
                        action="store_true", 
                        default=False)
    parser.add_argument("--no_merge", 
                    action="store_true", 
                    default=False)

    args = parser.parse_args()
    cat_path = args.cat_path
    stat_only = args.stat_only
    out_format = args.format
    rem_dupes = args.rem_dupes
    no_merge = args.no_merge
    out_dir = args.out_dir

    # format output path
    path_suffix = cat_path.suffixes[-2:]
    out_base = out_dir / Path(cat_path.stem).stem
    out_cpy_path = out_dir / (Path(cat_path.stem).stem + "_with-ids" + path_suffix[0] + path_suffix[1])


    # read in catalog data to pd dataframe
    print("Reading data...")
    cat_df = pd.read_csv(cat_path,
                          sep="\t",
                          usecols=["chrom", 
                                   "start_0based", 
                                   "end_1based", 
                                   "ReferenceRegion",
                                   "ReferenceMotif", 
                                   "MotifSize",
                                   "ReferenceRepeatPurity",
                                   "AlleleFrequenciesFromIllumina174k",
                                   ],
                            engine="python" # using because the C parser was hitting an error in AlleleFrequenciesFromIllumina174k
                          )
    in_shape = cat_df.shape
    cat_df["_orig_id"] = range(len(cat_df)) # make explicit temp col to keep og order/unique row identifiers

    # calculate overlaps
    print("Calculating stats...")
    stat_dict = getStats(cat_df)


    if not stat_only:
        
        
        # handle exact DUPLICATES (not overlaps where 100% is contained within another)
        num_dupes = cat_df[cat_df.duplicated(subset= ["ReferenceRegion"], keep=False)].shape[0]

        if rem_dupes:
            cat_df = removeDupes(cat_df) # remove exact duplicates

        if not no_merge:
            print("Running merge...")
            clean_df = mergeOverlaps(cat_df) # handle overlaps by merging into single extended row

            # output updated copy of input catalog with id column
            cat_df[[*clean_df.columns[:3], "ReferenceMotif", "MotifSize", "_orig_id"]].to_csv(out_cpy_path, index=False, sep="\t", compression="gzip")

        # output bed using desired format
        if out_format.lower() == "atarva":
            at_path = Path(str(out_base) + "_mrgd" + "_atarva" + ".bed" + path_suffix[1])
            out_shape = ToAtarva(clean_df, at_path)

        # output tsv containg extended merge info
        out_path = Path(str(out_base) + "_mrgd" + "_with-info" + path_suffix[0] + path_suffix[1])
        clean_df[[*clean_df.columns[:3], "ReferenceMotif", "MotifSize", "merged_from", "other_motifs"]].to_csv(out_path, index=False, sep="\t", compression="gzip")


    etime = time.perf_counter()
    rtime = etime - stime

    print("---Program Complete---")
    print(f"Runtime                 : {rtime:.3}s")
    print(f"Input rows              : {in_shape[0]}")
    print(f"Output rows             : {out_shape[0]}")
    print(f"Duplicate regions found : {num_dupes}")
    print("Overlap bins            :")
    for label, count in stat_dict.items():
        print(f"    {label:<12}: {count}")
    if not stat_only:
        print(f"Output Dir     : {out_base.parent}")



if __name__ == "__main__":
    main()
