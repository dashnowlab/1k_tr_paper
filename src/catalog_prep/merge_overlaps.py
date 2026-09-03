from pathlib import Path
import argparse
from datetime import date
import time
import pandas as pd
import numpy as np
import bioframe as bf

"""
Author: GA
Description: filters out duplicat regions in catalog
             based on frequency, purity, & motif size
Example: python filter_duplicates <catalog_path> --format atarva
python filter_duplicates.py /projects/garner1@xsede.org/dashnow_lab/1k_tr_paper/data/local_data/TR_catalog.5599658_loci.20260123_034640.tsv.gz --format atarva
"""

def merge_key(row):
    return (
        -(row["end_1based"] - row["start_0based"]), # grab the longest loci
        -row.get("ReferenceRepeatPurity", 0), # then get the nighest purity
        -row.get("MotifSize", 0), # larger motif size
        row["_tiebreak"]
    )


def ToAtarva(cdf: pd.DataFrame, out_path: str | Path):

    cols = [*cdf.columns[1:4], "ReferenceMotif", "MotifSize", "merged_from"]
    fdf = cdf[cols]
    fdf.columns = ["#CHROM", "START", "END", "MOTIF", "MOTIF_LEN", "merged_from"]

    print(fdf["END"])

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
    cols = tuple(cdf.columns[:3])

    clustered = bf.cluster(cdf, cols=cols, min_dist=min_dist)

    # collect all member ids per cluster
    merged_ids = (
        clustered.groupby("cluster")["_orig_id"]
        .apply(list)
        .rename("merged_from")
    )

    # apply merge_key criteria to determine rank within each cluster
    clustered["_rank"] = clustered.apply(merge_key, axis=1)
    clustered = clustered.sort_values(["cluster", "_rank"])
    representative = clustered.groupby("cluster", as_index=False).first()

    print(clustered["cluster"].isna().sum())

    # attach the full member list to each representative row
    representative = representative.merge(merged_ids, on="cluster", how="left")

    # change start and end coords to cluster coords
    representative[cols[1]] = representative["cluster_start"]
    representative[cols[2]] = representative["cluster_end"]

    print(representative.head)

    dropped = clustered[~clustered.index.isin(
        clustered.groupby("cluster")["_rank"].idxmin()
    )]

    cdf.drop(columns=["_tiebreak"], inplace=True)

    return representative.drop(columns=["cluster_start", "cluster_end", "_rank"]), dropped


def main():
    stime = time.perf_counter()
    dt = date.today().strftime("%Y%m%d")

    parser = argparse.ArgumentParser()
    parser.add_argument("cat_path", type=Path)
    #parser = argparse.ArgumentParser()
    #parser.add_argument("out_dir", type=Path)
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

    args = parser.parse_args()
    cat_path = args.cat_path
    stat_only = args.stat_only
    out_format = args.format
    rem_dupes = args.rem_dupes


    # format output path
    suffix = cat_path.suffixes[-2:]
    out_path = cat_path.parent / (Path(cat_path.stem).stem + "_ovlp-mrg" + suffix[0] + suffix[1])
    out_cpy_path = cat_path.parent / (Path(cat_path.stem).stem + "_with_ids" + suffix[0] + suffix[1])


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
                            engine="python" # had to use because the C parser was hitting an error in AlleleFrequenciesFromIllumina174k
                          )
    in_shape = cat_df.shape
    cat_df["_orig_id"] = range(len(cat_df)) # make explicit temp col to keep og order * unique row identifiers

    # calculate overlaps
    print("Calculating stats...")
    stat_dict = getStats(cat_df)


    if not stat_only:

        # handle exact DUPLICATES (not overlaps where 100% is contained within another)
        num_dupes = cat_df[cat_df.duplicated(subset= ["ReferenceRegion"], keep=False)].shape[0]

        if rem_dupes:
            cat_df = (
                cat_df.sort_values(
                    by=["AlleleFrequenciesFromIllumina174k", "ReferenceRepeatPurity", "MotifSize"], 
                    ascending=[False, False, True]
                    )
                    .drop_duplicates(subset="ReferenceRegion", keep="first") # keep highest order duplicate based on above values
                    .sort_values("_orig_id") # ensure orignal ordering is held
                    .drop(columns="_orig_id") # drop ordering column
            )

            # double check to see if any duplicates remain
            num_dupes_cleaned = cat_df[cat_df.duplicated(subset= ["ReferenceRegion"], keep=False)].shape[0]

            if num_dupes_cleaned > 0:
                print(f"WARNING {num_dupes_cleaned} duplicates still detected.")
                return 1

        # handle overlaps
        clean_df, drp_df = mergeOverlaps(cat_df)

        # output using desired format
        if out_format.lower() == "atarva":
            out_shape = ToAtarva(clean_df, out_path)
        else:
            out_shape = clean_df.shape
            clean_df.to_csv(out_path, index=False, sep="\t", compression="gzip")

        # output updated copy of input catalog with id column
        cat_df.to_csv(out_cpy_path, index=False, sep="\t", compression="gzip")


    etime = time.perf_counter()
    rtime = etime - stime

    print("---Program Complete---")
    print(f"Runtime                 : {rtime:.3}s")
    print(f"Input size              : {in_shape}")
    print(f"Output size             : {out_shape}")
    print(f"Duplicate regions found : {num_dupes}")
    print("Overlap bins            :")
    for label, count in stat_dict.items():
        print(f"    {label:<12}: {count}")
    if not stat_only:
        print(f"Merged output path      : {out_path}")
        print(f"Id updated output path  : {out_cpy_path}")



if __name__ == "__main__":
    main()
