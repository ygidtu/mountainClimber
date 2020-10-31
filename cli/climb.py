#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created at 2020.09.04 by Zhang Yiming
"""
import os
from glob import glob

from multiprocessing import Pool, cpu_count
from subprocess import check_call

import click
import pybedtools as pb
import pysam

from loguru import logger

from src import get_junction_counts, mountainClimberTU, merge_tus, mountainClimberRU, mountainClimberCP
from src.functions import run_command, exists


def __single_process__(args):
    if not os.path.exists(args['input_file'] + ".bai"):
        pysam.index(args['input_file'])

    jxn = args['output'] + "_jxn.bed"
    if not exists(jxn):
        get_junction_counts.make_intron_bed(
            args['input_file'], jxn, 
            overhang = args['overhang'], min_intron = args['min_intron'], 
            max_intron = args['max_intron'], strandedness = args['strand']
        )
    else:
        logger.info("{} exists, pass", jxn)

    bedgraph = args['output'] + ".bedgraph"
    if not exists(bedgraph):
        check_call(f"bedtools genomecov -trackline -bg -split -ibam {args['input_file']} -g {args['output']} > {bedgraph}", shell=True, cwd = os.path.dirname(args['output']))
    else:
        logger.info("{} exists, pass", bedgraph)


    tu = args['output'] + "_tu.bed"
    if not exists(tu):
        mountainClimberTU.process(
            bedgraph, jxn, args['gsize'], tu, args['strand_in_int'].get(args['strand'], 0), 
            args['min_jxn_count'], args['min_percent'], 
            args['min_reads'], args['window_size']
        )
    else:
        logger.info("{} exists, pass", tu)

    args["jxn"] = jxn
    args["tu"] = tu
    args["bedgraph"] = bedgraph
    return args


def __single_process_cp_ru__(args):

    if not exists(args["output"] + "_CP.bed"):
        mountainClimberCP.run(
            input_bg=args["bedgraph"], input_bg_minus=args["strand"] == "fr-secondstrand",
            input_regions=args['merged_ref'], output=args["output"] + "_CP.bed",
            genome=args['genome'], plot=args['plot'], junc=args["jxn"], juncdist=args["juncdist"],
            minjxncount=args["min_jxn_count"], min_length=args['min_length'], min_expn=args['min_expn'],
            test_thresh=args['test_thresh'], winsize=args['winsize'], peak_min_dist=args['peak_min_dist'],
            peak_thresh=args['peak_thresh'], fcthresh=args['fcthresh'], max_end_ru=args['max_end_ru'], 
            verbose=args['verbose'], min_expn_distal=args["min_expn_distal"], # min_segments=args["min_segments"]
        )
    else:
        logger.info("{} exists, pass", args["output"] + "_CP.bed")

    if not exists(args["output"] + "_RU.bed"):
        mountainClimberRU.run(
            input_file=args["output"] + "_CP.bed", 
            output=args["output"] + "_RU.bed",
            min_segments=args['min_segments'], 
            verbose=args['verbose']
        )
    else:
        logger.info("{} exists, pass", args["output"] + "_RU.bed")


@click.command()
@click.option(
    "-g", "--genome", type=click.Path(exists=True),
    help="Path to genome fasta file", required = True
)
@click.option(
    "-r", "--reference", type=click.Path(exists=True),
    help="Path to reference gtf file", required = True
)
@click.option(
    "-o", "--output", type=click.Path(),
    help="Prefix of output file", required=True
)
@click.option(
    "-s", "--strand", type=click.Choice(['fr-firststrand', 'fr-secondstrand', 'fr-unstrand', 'single']),
    help="Strandedness. Options", default="fr-unstrand",  show_default=True
)
@click.option(
    "-e", "--overhang", type=int, default=8,
    help="Minimum number of base pairs in each exon", show_default=True
)
@click.option(
    '-m', '--min-intron', type=int, default=30,
    help="Minimum intron length", show_default=True
)
@click.option(
    '-a', '--max-intron', type=int, default=500000,
    help="Maximum intron length", show_default=True
)
@click.option(
    '-c', '--min-jxn-count', type=int, default=2,
    help="Minimum junction read count.", show_default=True
)
@click.option(
    '-w', '--window-size', type=int, default=1000,
    help="Window size.", show_default=True
)
@click.option(
    '-p', '--min_percent', type=float, default=1.0,
    help="Minimum percentage of the window covered.", show_default=True
)
@click.option(
    '-n', '--min-reads', type=int, default=10,
    help="Minimum number of reads per window.", show_default=True
)
@click.option(
    "-t", "--n_jobs", type=click.IntRange(1, cpu_count(), clamp=True), default=1,
    help="The number of processes to use", show_default=True
)
@click.option(
    '--peak_thresh', type=float, default=-1.0, 
	help='Normalized threshold (float between [0., 1.]). Only the peaks with amplitude higher than the threshold will be detected. -1.0 indicates optimizing between [0.01, 0.05, 0.1] for each TU.', show_default=True
)
@click.option(
    '--peak_min_dist',type=int, default=-1, show_default=True,
	help='Minimum distance between each detected peak. The peak with the highest amplitude is preferred to satisfy this constraint. -1 indicates optimizing between [10, 50] for each TU.'
)
@click.option(
    '--winsize', type=int, default=-1, show_default=True,
	help='Window size for de-noising and increasing precision. -1 indicates optimizing between [50, max(100, gene_length / 100) * 2].'
)
@click.option(
    '--test_thresh', type=float, default=0.001, show_default=True,
	help='Maximum p-value threshold for KS test and t-test.'
)
@click.option(
    '--min_length', type=int, default=1000, show_default=True,
	help='Minimum gene length for running mountain climber.'
)
@click.option(
    '--min_expn', type=int, default=10, show_default=True,
	help='Minimum expression level (average # reads per bp)'
)
@click.option(
    '--min_expn_distal', type=int, default=1, show_default=True,
	help='Minimum distal expression level (average # reads per bp).'
)
@click.option(
    '--fcthresh', type=float, default=1.5, show_default=True,
	help='Minimum fold change.'
)
@click.option(
    '--juncdist',type=int, default=10, show_default=True,
	help='Minimum distance to exon-intron junction.'
)
@click.option(
    '--max_end_ru', type=float, default=0.01, show_default=True,
	help='Maximum end relative usage = coverage of end / max segment coverage.'
)
@click.option(
    '--max_end_ru', type=float, default=0.01, show_default=True,
	help='Maximum end relative usage = coverage of end / max segment coverage.'
)
@click.option(
    '--plot', type=click.BOOL, default=False, show_default=True,
	help='Plot the cumulative read sum (CRS), the distance from CRS to line y=ax, and the coverage with predicted change points.'
)
@click.option(
    '--min_segments', type=int, default=3, show_default=True,
    help='Minimum number of segments required in the TU to calculate relative end usage.'
)
@click.option(
    '--pair-end', type=click.BOOL, default=False, show_default=True,
	help='Whether input is paired-end'
)
@click.option(
    '--verbose', type=click.BOOL, default=False, show_default=True,
	help='Print progress.'
)
@click.argument('input-files', type=click.Path(exists=True), nargs=-1)
def climb(
    input_files, genome: str, reference: str, strand: str, output: str, 
    overhang: int, min_intron: int, max_intron: int, n_jobs: int,
    min_jxn_count: int, window_size: int, min_percent: float, min_reads: int,
    peak_thresh: float, peak_min_dist: int, winsize: int, test_thresh: float,
    min_length: int, min_expn: int, min_expn_distal: int, fcthresh: float,
    juncdist: int, max_end_ru: float, plot: bool, min_segments: int, verbose: bool,
    pair_end: bool
):
    u""" Climb pipeline """
    if not os.path.exists(genome + ".fai"):
        pysam.faidx(genome)
    
    output = os.path.abspath(output)
    os.makedirs(output, exist_ok=True)
    temp_dir = output + "_tmp"
    os.makedirs(temp_dir, exist_ok=True)
    pb.helpers.set_tempdir(temp_dir)

    logger.info("Get genome size")
    gsize = os.path.join(output, "genome.chrom.sizes")
    check_call(f"cut -f1,2 {genome}.fai > {gsize}", shell=True, cwd=os.path.dirname(output))
    strand_in_int = {'fr-firststrand': 1, 'fr-secondstrand': -1}

    cmds = [{
        "input_file": os.path.abspath(i), "output": os.path.join(output, os.path.basename(i).split(".")[0]),
        "genome": os.path.abspath(genome), "reference": os.path.abspath(reference), "gsize": gsize,
        "strand": strand, "overhang": overhang,
        "min_intron": min_intron, "max_intron": max_intron,
        "min_jxn_count": min_jxn_count, "window_size": window_size, 
        "min_percent": min_percent, "min_reads": min_reads,
        "strand_in_int": strand_in_int, 
        "peak_thresh": peak_thresh, "peak_min_dist": peak_min_dist, 
        "winsize": winsize, "test_thresh": test_thresh,
        "min_length": min_length, "min_expn": min_expn, 
        "min_expn_distal": min_expn_distal, "fcthresh": fcthresh,
        "juncdist": juncdist, "max_end_ru": max_end_ru, "plot": plot,
        "verbose": verbose, "min_segments": min_segments
    } for i in input_files if os.path.exists(i)]
    
    with Pool(n_jobs) as p:
        res = list(p.map(__single_process__, cmds))
    
    p.close()
    p.join()

    logger.info("merge")
    tu_merge = os.path.join(output, "tus_merged")
    if not os.path.exists(tu_merge):
        merge_tus.merge(infiles=[x["tu"] for x in res if os.path.exists(x["tu"])], output=tu_merge, refgtf=reference, ss="y" if strand in strand_in_int.keys() else "n")

    # logger.info("run rsem")
    # temp_rsem = os.path.join(output, "RSEM")
    # temp_rsem_ref = os.path.join(temp_rsem, "ref")
    # os.makedirs(temp_rsem_ref, exist_ok=True)
    # tu_merge = glob(tu_merge + "*.gtf")
    # tu_merge = tu_merge[0]
    # run_command(f"rsem-prepare-reference -p {n_jobs} --gtf {tu_merge} --star {genome} {temp_rsem_ref}")

    # for i in input_files:
    #     paired_end = "--paired-end" if pair_end else ""
    #     transcriptome_bam = i.replace("Aligned.sortedByCoord.out.bam", "Aligned.toTranscriptome.out.bam")
    #     run_command(f"rsem-calculate-expression -p {n_jobs} {paired_end} --append-names --seed 0 --estimate-rspd --sampling-for-bam --output-genome-bam --alignments {transcriptome_bam} {temp_rsem_ref} {os.path.join(temp_rsem, os.path.basename(i).split('.')[0])}")

    merged_ref = glob(os.path.join(output, "*_singleGenes.bed"))
    if merged_ref:
        merged_ref = merged_ref[0]

        new_cmds = []
        for c in res:
            c["merged_ref"] = merged_ref
            new_cmds.append(c)

        with Pool(n_jobs) as p:
            res = list(p.map(__single_process_cp_ru__, new_cmds))
        
        p.close()
        p.join()

    pb.helpers.cleanup()
