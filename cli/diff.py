#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created at 2020.10.10 by Zhang Yiming
"""
import os
from glob import glob
from multiprocessing import cpu_count, Pool
from shutil import rmtree

import click
import pybedtools as pb

from src import diff_cluster, diff_segmentReadCounts, diff_ru, diff_test
from src.functions import sort_bedfile


def multiprocessing_diff_ru(data):

    for i in data["files"]:
        sort_bedfile(i, i, sort_by_bedtools=True, add_header = False)

    diff_ru.diff_ru(
        data["files"],
        segments=data["segments"], 
        condition=data["condition"], 
        input_cp=data["input_cp"], 
        output=data["output"], 
        min_segments=data["min_segments"], 
        verbose=data["verbose"]
    )
    pb.helpers.cleanup()


@click.command()
@click.option(
    "-i", "--input_dir", type=click.Path(exists=True),
    help="Path to climb output directory", required = True
)
@click.option("-o", "--output", type=click.Path(),help="Prefix of output file", required=True)
@click.option("--min-conditions", type=float, default=1, help="Minimum number of conditions for a gene to be clustered across conditions.", show_default=True)
@click.option(
    "--minpts", type=int, default=0, show_default=True, 
    help="List of space-delimited DBSCAN minPts values. These indicate the minimum # points for DBSCAN to consider a core point. The minimum of this list will be used to cluster across conditions."
)
@click.option("--min-fc", type=float, default=-1,help="Minimum fold change for change points.", show_default=True)
@click.option("--min-expn", type=float, default=0, help="Minimum expression in exons for a gene to be clustered.", show_default=True)
@click.option("--lm-flag", type=click.BOOL, default=False, show_default=True, help="Input are results from diff_cluster.")
@click.option("--ss-flag", type=click.BOOL, default=False, show_default=True, help="Flag: RNA-Seq is strand-specific.")
@click.option(
    "--eps", type=float, default=-1, show_default=True,
    help="Maximum distance between 2 points in a neighborhood. -1.0 indicates using the minimum optimal window size from mountain climber."
)
@click.option("--bgminus", type=click.Path(exists=True),help="List of space-delimited bedgraphs: minus strand. One file per line", show_default=True)
@click.option("--min-segments", type=int, default=3, help="Minimum number of segments required in the TU to calculate relative end usage", show_default=True)
@click.option("-p", "--processes", type=click.IntRange(1, cpu_count(), clamp=True), default=1, help="How many processes to use", show_default=True)
@click.option("--min-dstl-cov", type=float, default=5, help="Minimum average reads per bp in distal segment across samples in at least 1 condition.", show_default=True)
@click.option("--min-prxl-cov", type=float, default=0, help="Minimum average reads per bp in proximal segment across samples in all conditions.", show_default=True)
@click.option("--pmax", type=float, default=0.05, help="Maximum p-value.", show_default=True)
@click.option("--dtop-abs-dif-min", type=float, default=0.05, help="Minimum relative usage (RU) difference.", show_default=True)
@click.option("--root", type=str, help="The root to compare to.")
@click.option('--verbose', type=click.BOOL, default=False, show_default=True, help='Print progress.')
@click.option('--keep', type=click.BOOL, default=False, show_default=True, help='Keep temp file of diff_test.')
def diff(
    input_dir: str, output: str, min_conditions: int,
    minpts: int, min_fc: float, min_expn: float, 
    lm_flag: bool, eps: float, ss_flag: bool,
    bgminus: str, min_segments: int, processes: int,
    min_dstl_cov:float, min_prxl_cov:float, pmax: float, dtop_abs_dif_min: float,
    verbose: bool, keep: bool, root: str
):
    u"""
    differential analysis after climb pipline
    
    \f
    Output files include _cluster_totals.txt, _segments.bed, _cp.bed, and one _cp.bed file for each condition. _cp.bed name field = label_prioritized;condition_labels:gene:TUstart:TUend:chrom:strand:dbscan_epsilon:min_clustered_change_point:max_clustered_change_point:cluster_standard_deviation:total_clusters. _segments.bed name field = label_prioritized_cp1;condition_labels_cp1|label_prioritized_cp2;condition_labels_cp2:gene:TUstart:TUend:chrom:strand:dbscan_epsilon.
    """
    output = os.path.abspath(output)
    temp_dir = output + "_tmp"
    os.makedirs(temp_dir, exist_ok=True)
    pb.helpers.set_tempdir(temp_dir)
 
    # diff cluster
    cp_files = {}
    for x in glob(os.path.join(input_dir, "*_CP.bed")):
        c = os.path.basename(x).replace("_CP.bed", "").split("_")[1:]
        c = "_".join(c).strip("_")

        temp = cp_files.get(c, [])
        temp.append(x)
        cp_files[c] = temp

    # construct the pairs
    comparisions = []
    if root:
        for c in cp_files.keys():
            if c != root:
                comparisions.append([root, c])
    else:
        cs = sorted(cp_files.keys())
        for i in range(len(cs)):
            for j in range(i+1, len(cs)):
                comparisions.append([cs[i], cs[j]])

    for i in comparisions:

        temp_diff_cluster_output = os.path.join(output, "_vs_".join(i), "diff_cluster")
        temp_diff_segments_output = os.path.join(output, "_vs_".join(i), "segments_read_counts")
        temp_diff_ru_output = os.path.join(output, "_vs_".join(i), "diff_ru")
        temp_diff_test_output = os.path.join(output, "_vs_".join(i), "diff_test")
        
        os.makedirs(temp_diff_cluster_output, exist_ok=True)
        temp_files, conditions = [], []
        for c in i:
            temp_files += cp_files[c]
            conditions += [c for _ in range(len(cp_files[c]))]

        diff_cluster.cluster(
            temp_files, os.path.join(temp_diff_cluster_output, "diff"), 
            conditions, min_conditions, [minpts for _ in range(len(temp_files))], 
            min_fc, min_expn, lm_flag, eps, ss_flag, 
            verbose=verbose
        )
        pb.helpers.cleanup()

        # diff segment read count
        input_file, conditions = [], []
        for c in i:
            for x in glob(os.path.join(input_dir, f"*_{c}.bedgraph")):
                rc = os.path.basename(x).replace(".bedgraph", "").split("_")
                rc = "_".join(rc[1:])

                if rc == c:
                    input_file.append(x)
                    # conditions.append(os.path.basename(x).replace(".bedgraph", ""))
                    conditions.append("")

        bgminus_list = None
        if bgminus:
            with open(bgminus) as r:
                bgminus_list = [x.strip() for x in r if x and os.path.exists(x.strip())]

        seg = os.path.join(temp_diff_cluster_output, "diff_segments.bed")
        sort_bedfile(seg, seg, add_header = False, sort_by_bedtools = True)
        os.makedirs(temp_diff_segments_output, exist_ok=True)
        diff_segmentReadCounts.read_count(
            segments=seg, 
            conditions=conditions,
            bgplus=input_file, 
            bgminus=bgminus_list, 
            output=os.path.join(temp_diff_segments_output, "diff"),
            n_jobs=processes,
        )
        pb.helpers.cleanup()

        # diff_ru
        input_file = {}
        for x in glob(os.path.join(temp_diff_segments_output, "*_readCounts.bed")):
            c = os.path.basename(x).replace("_readCounts.bed", "").replace("diff_", "").split("_")
            c = "_".join(c[1:]).strip("_")
            temp = input_file.get(c, [])
            temp.append(x)
            input_file[c] = temp

        cmds = []
  
        for c, files in input_file.items():
            cmds.append({
                "files": files, "condition": c,
                "segments": os.path.join(temp_diff_cluster_output, "diff_segments.bed"), 
                "input_cp": os.path.join(temp_diff_cluster_output, f"diff_cp_{c}.bed"),
                "output": os.path.join(temp_diff_ru_output, "diff"), "min_segments": min_segments,
                "verbose": verbose
            })

        os.makedirs(temp_diff_ru_output, exist_ok=True)
        with Pool(min(processes, len(cmds))) as p:
            p.map(multiprocessing_diff_ru, cmds)

        cmds = {
            "input_file": [], "conditions_input": [],
            "ru_segments": set(), "conditions_ru_segments": set()
        }

        for c in i:
            files = input_file[c]
            for f in files:
                cmds["input_file"].append(f)
                cmds["conditions_input"].append(c)
            
            cmds["ru_segments"].add(os.path.join(temp_diff_ru_output, f"diff_ru_segments_{c}.bed"))
            cmds["conditions_ru_segments"].add(c)
            
        os.makedirs(temp_diff_test_output, exist_ok=True)
        diff_test.test(
            input_file=cmds["input_file"], ru_segments=list(cmds["ru_segments"]), 
            output=os.path.join(temp_diff_test_output, f"diff_{'_'.join(i)}"), 
            conditions_input=cmds["conditions_input"], conditions_ru_segments=list(cmds["conditions_ru_segments"]), 
            dtop_abs_dif_min=dtop_abs_dif_min, min_dstlCov=min_dstl_cov, 
            min_prxlCov=min_prxl_cov, pmax=pmax, 
            verbose=verbose, keep=keep
        )

    if os.path.exists(temp_dir):
        rmtree(temp_dir)
    # python src/diff_cluster.py -i tests/science/C000R7B4_neutrophil_CP.bed tests/science/S001QBB1_HSC_CP.bed -c nerutrophil HSC -o tests/diff_test1
    # python src/diff_segmentReadCounts.py -i tests/diff_test_segments.bed -p tests/science/C000R7B4_neutrophil.bedgraph tests/science/S001QBB1_HSC.bedgraph -c nerutrophil HSC -o tests/diff
    # python src/diff_ru.py -i tests/diff_C000R7B4_neutrophil_nerutrophil_readCounts.bed -s tests/diff_test_segments.bed -c neutrophil -o tests/diff_neutrophil -l tests/diff_test_cp_nerutrophil.bed