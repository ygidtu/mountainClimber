#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created at 2020.10.10 by Zhang Yiming
"""
import os
from glob import glob

import click

from src import diff_cluster, diff_segmentReadCounts, diff_ru


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
@click.option('--verbose', type=click.BOOL, default=False, show_default=True, help='Print progress.')
def diff(
    input_dir: str, output: str, min_conditions: int,
    minpts: int, min_fc: float, min_expn: float, 
    lm_flag: bool, eps: float, ss_flag: bool,
    bgminus: str, min_segments: int,
    verbose: bool
):
    u"""
    differential analysis after climb pipline
    
    \f
    Output files include _cluster_totals.txt, _segments.bed, _cp.bed, and one _cp.bed file for each condition. _cp.bed name field = label_prioritized;condition_labels:gene:TUstart:TUend:chrom:strand:dbscan_epsilon:min_clustered_change_point:max_clustered_change_point:cluster_standard_deviation:total_clusters. _segments.bed name field = label_prioritized_cp1;condition_labels_cp1|label_prioritized_cp2;condition_labels_cp2:gene:TUstart:TUend:chrom:strand:dbscan_epsilon.
    """
    # diff cluster
    input_file, conditions = [], []
    for x in glob(os.path.join(input_dir, "*_CP.bed")):
        input_file.append(x)
        conditions.append(os.path.basename(x).replace("_CP.bed", "").split("_")[1])

    temp_diff_cluster_output = os.path.join(output, "diff_cluster")
    os.makedirs(temp_diff_cluster_output, exist_ok=True)
    diff_cluster.cluster(
        input_file, os.path.join(temp_diff_cluster_output, "diff"), 
        conditions, min_conditions, [minpts for _ in range(len(input_file))], 
        min_fc, min_expn, lm_flag, eps, ss_flag, 
        verbose=verbose
    )

    # diff segment read count
    input_file, conditions = [], []
    for x in glob(os.path.join(input_dir, "*.bedgraph")):
        input_file.append(x)
        conditions.append(os.path.basename(x).replace(".bedgraph", ""))

    temp_diff_segments_output = os.path.join(output, "segments_read_counts")

    bgminus_list = None
    if bgminus:
        with open(bgminus) as r:
            bgminus_list = [x.strip() for x in r if x and os.path.exists(x.strip())]

    diff_segmentReadCounts.read_count(
        segments=os.path.join(temp_diff_cluster_output, "diff_segments.bed"), 
        conditions=conditions,
        bgplus=input_file, 
        bgminus=bgminus_list, 
        output=os.path.join(temp_diff_segments_output, "diff")
    )

    # diff_ru
    input_file = {}
    for x in glob(os.path.join(temp_diff_segments_output, "_readCounts.bed")):
        c = os.path.basename(x).replace("_readCounts.bed", "").split("_")[-1]
        temp = input_file.get(c, [])
        temp.append(x)
        input_file[c] = temp

    for c, files in input_file.items():
        diff_ru.diff_ru(
            files, 
            segments=os.path.join(temp_diff_cluster_output, "diff_segments.bed"), 
            condition=c, 
            input_cp=os.path.join(temp_diff_cluster_output, f"diff_cp_{c}.bed"), 
            output=output, 
            min_segments=min_segments, 
            verbose=verbose
        )

    # python src/diff_cluster.py -i tests/science/C000R7B4_neutrophil_CP.bed tests/science/S001QBB1_HSC_CP.bed -c nerutrophil HSC -o tests/diff_test1
    # python src/diff_segmentReadCounts.py -i tests/diff_test_segments.bed -p tests/science/C000R7B4_neutrophil.bedgraph tests/science/S001QBB1_HSC.bedgraph -c nerutrophil HSC -o tests/diff
    # python src/diff_ru.py -i tests/diff_C000R7B4_neutrophil_nerutrophil_readCounts.bed -s tests/diff_test_segments.bed -c neutrophil -o tests/diff_neutrophil -l tests/diff_test_cp_nerutrophil.bed