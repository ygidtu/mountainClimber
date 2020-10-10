#!/user/bin/python -tt
"""
Calculate the average reads/bp for each segment after clustering.
"""

import os
import sys
import argparse
import numpy as np  # v1.10.4
import pybedtools as pb
from collections import defaultdict
from datetime import datetime
from loguru import logger

try:
	from bed6 import BedUtils, BedPair
except ImportError:
	from src.bed6 import BedUtils, BedPair


def bedgraph_per_gene_intersect(src, target, output=None, strandness = False, keep_all = False):
	if isinstance(src, str):
		src = BedUtils.read(src)

	if isinstance(target, str):
		target = BedUtils.read(target, required_regions = src)

	logger.info(f"{len(src)}, {len(target)}")

	res = BedUtils.overlap(src, target, strandness = strandness)
	if output:
		BedPair.save(res, output, keep_all = keep_all)
	else:
		return res


def get_seg2cov(intersect, cond, sample, outfile):
	"""Get coverage of each segment & write to outfile in bed format"""
	seg2cov = {}
	o = open(outfile, 'w')

	maxl = 0
	with open(intersect, 'r') as f:
		for l, line in enumerate(f):
			maxl = l

	# @2020.10.10 by Zhang Yiming - init variables
	prev_gene, prev_cov_array = "", [],
	prev_geneid, prev_chrom, prev_start, prev_end, prev_strand = "", "", -1, -1, "."
	with open(intersect, 'r') as f:
		for l, line in enumerate(f):
			if line != '':
				x = line.rstrip().split('\t')
				if len(x) == 11:
					# KI270438.1       9442    19867   Intron;HSC|Junction;HSC:novel64:9251:22706:KI270438.1:NA:50.0     2       .       KI270438.1      14517   14593   .       1       .       76
					(achrom, astart, aend, ageneid, ascore, astrand, bchrom, bstart, bend, bcov, overlap_len) = x
				elif len(x) == 10:
					(achrom, astart, aend, ageneid, ascore, bchrom, bstart, bend, bcov, overlap_len) = x
					astrand = 0
				else:
					logger.error('EXIT: do not recognize bedgraph intersect format')
					logger.error(line)
					sys.exit(1)

				astart = int(astart)
				aend = int(aend)
				bstart = int(bstart)
				bend = int(bend)
				bcov = float(bcov)
			else:
				continue

			if overlap_len == '0':
				continue

			if l == 0:  # first line
				prev_gene = ':'.join(x[:5]) if astrand == 0 else ':'.join(x[:6])
				this_start = max(astart, bstart)
				this_end = min(aend, bend)
				prev_cov_array = np.zeros(aend - astart)
				prev_cov_array[(this_start - astart):(this_end - astart)] += bcov

				# === next round  ===
				prev_start = astart
				prev_end = aend
				prev_geneid = ageneid
				prev_chrom = achrom
				prev_strand = astrand

				if l == maxl:
					cov_avg = np.mean(prev_cov_array)
					if (prev_geneid, prev_chrom, prev_start, prev_end, prev_strand, sample, cond) not in seg2cov:
						seg2cov[(prev_geneid, prev_chrom, prev_start, prev_end, prev_strand, sample, cond)] = cov_avg
						if prev_strand == 0:
							o.write('\t'.join(map(str, [prev_chrom, prev_start, prev_end, ':'.join([prev_geneid, sample, cond]), cov_avg])) + '\n')
						else:
							o.write('\t'.join(map(str, [prev_chrom, prev_start, prev_end, ':'.join([prev_geneid, sample, cond]), cov_avg, prev_strand])) + '\n')
					else:
						logger.error(' '.join(map(str, ['EXIT: seen seg2cov:', prev_geneid, prev_chrom, prev_start, prev_end, prev_strand, sample, cond])) + '\n')
						logger.error(' '.join(map(str, [cov_avg, seg2cov[(prev_geneid, prev_chrom, prev_start, prev_end, prev_strand, sample, cond)]])) + '\n')
						sys.exit(1)
			else:
				this_gene = ':'.join(x[:5]) if astrand == 0 else ':'.join(x[:6])
				if line == '' and this_gene == prev_gene and this_gene == '':  # EOF
					break
				elif this_gene == prev_gene and l != maxl:  # get coverage
					this_start = max(astart, bstart)
					this_end = min(aend, bend)
					prev_cov_array[(this_start - astart):(this_end - astart)] += bcov

					# === next round  ===
					prev_start = astart
					prev_end = aend
					prev_geneid = ageneid
					prev_chrom = achrom
					prev_strand = astrand
				else:  # finished reading all info for one gene
					# === get per-base coverage ===
					if l == maxl:
						# === last line of this gene & last line of the file ===
						this_start = max(astart, bstart)
						this_end = min(aend, bend)

						# === previous gene ===
						prev_cov_array[(this_start - astart):(this_end - astart)] += bcov
						cov_avg = np.mean(prev_cov_array)
						if (prev_geneid, prev_chrom, prev_start, prev_end, prev_strand, sample, cond) not in seg2cov:
							seg2cov[(prev_geneid, prev_chrom, prev_start, prev_end, prev_strand, sample, cond)] = cov_avg
							if prev_strand == 0:
								o.write('\t'.join(map(str, [prev_chrom, prev_start, prev_end, ':'.join([prev_geneid, sample, cond]), cov_avg])) + '\n')
							else:
								o.write('\t'.join(map(str, [prev_chrom, prev_start, prev_end, ':'.join([prev_geneid, sample, cond]), cov_avg, prev_strand])) + '\n')
						else:
							logger.error(' '.join(map(str, ['EXIT: seen seg2cov:', prev_geneid, prev_chrom, prev_start, prev_end, prev_strand, sample, cond])) + '\n')
							logger.error(' '.join(map(str, [cov_avg, seg2cov[(prev_geneid, prev_chrom, prev_start, prev_end, prev_strand, sample, cond)]])) + '\n')
							sys.exit(1)

						if this_gene != prev_gene:
							prev_start = astart
							prev_end = aend
							prev_geneid = ageneid
							prev_chrom = achrom
							prev_strand = astrand

							cov_avg = np.mean(prev_cov_array)
							if (prev_geneid, prev_chrom, prev_start, prev_end, prev_strand, sample, cond) not in seg2cov:
								seg2cov[(prev_geneid, prev_chrom, prev_start, prev_end, prev_strand, sample, cond)] = cov_avg
								if prev_strand == 0:
									o.write('\t'.join(map(str, [prev_chrom, prev_start, prev_end, ':'.join([prev_geneid, sample, cond]), cov_avg])) + '\n')
								else:
									o.write('\t'.join(map(str, [prev_chrom, prev_start, prev_end, ':'.join([prev_geneid, sample, cond]), cov_avg, prev_strand])) + '\n')
							else:
								logger.error(' '.join(map(str, ['EXIT: seen seg2cov:', prev_geneid, prev_chrom, prev_start, prev_end, prev_strand, sample, cond])) + '\n')
								logger.error(' '.join(map(str, [cov_avg, seg2cov[(prev_geneid, prev_chrom, prev_start, prev_end, prev_strand, sample, cond)]])) + '\n')
								sys.exit(1)

					else:
						cov_avg = np.mean(prev_cov_array)
						if (prev_geneid, prev_chrom, prev_start, prev_end, prev_strand, sample, cond) not in seg2cov:
							seg2cov[(prev_geneid, prev_chrom, prev_start, prev_end, prev_strand, sample, cond)] = cov_avg
							if prev_strand == 0:
								o.write('\t'.join(map(str, [prev_chrom, prev_start, prev_end, ':'.join([prev_geneid, sample, cond]), cov_avg])) + '\n')
							else:
								o.write('\t'.join(map(str, [prev_chrom, prev_start, prev_end, ':'.join([prev_geneid, sample, cond]), cov_avg, prev_strand])) + '\n')
						else:
							logger.error(' '.join(map(str, ['EXIT: seen seg2cov:', prev_geneid, prev_chrom, prev_start, prev_end, prev_strand, sample, cond])) + '\n')
							logger.error(' '.join(map(str, [cov_avg, seg2cov[(prev_geneid, prev_chrom, prev_start, prev_end, prev_strand, sample, cond)]])) + '\n')
							sys.exit(1)

						# === first line of the next gene ===
						this_start = max(astart, bstart)
						this_end = min(aend, bend)
						prev_cov_array = np.zeros(aend - astart)
						prev_cov_array[(this_start - astart):(this_end - astart)] += bcov

						# === next round ===
						prev_start = astart
						prev_end = aend
						prev_geneid = ageneid
						prev_chrom = achrom
						prev_strand = astrand
						prev_gene = this_gene
	o.close()


def read_count(segments, conditions, bgplus, bgminus, output):
	# --------------------------------------------------
	# main routine
	# --------------------------------------------------

	if not segments:
		logger.error('EXIT: Please provide --segments')
		sys.exit(1)

	if not conditions:
		logger.error('EXIT: Please provide --conditions')
		sys.exit(1)

	if bgplus:
		if len(conditions) != len(bgplus):
			logger.error('EXIT: number of samples don\'t match!')
			sys.exit(1)
	else:
		logger.error('EXIT: Please provide --bgplus')
		sys.exit(1)

	if bgminus:
		if len(conditions) != len(bgminus):
			logger.error('EXIT: number of samples don\'t match!')
			sys.exit(1)
			
	# === set temporary dir ===
	os.makedirs(output, exist_ok=True)

	# === get the input files for each condition ===
	cond2bgplus = defaultdict(list)
	cond2bgminus = defaultdict(list)
	for c, cond in enumerate(conditions):
		cond2bgplus[cond].append(bgplus[c])
		if bgminus:
			cond2bgminus[cond].append(bgminus[c])

	# === get bedgraph for each segment ===
	logger.info(f'getting read counts per bp for each gene {datetime.now().time()}')
	# condSeg2cov = {}
	for cond in cond2bgplus:
		for b, bgplus in enumerate(cond2bgplus[cond]):
			sample = os.path.splitext(os.path.basename(bgplus))[0]
			logger.info(f'> bedtools intersect: {cond} {sample} {datetime.now().time()}')

			# get bedgraph for each segment
			intersect = '_'.join([output, cond, sample, 'intersect.txt'])
			if bgminus:
				segments = BedUtils.read(segments)
				res = bedgraph_per_gene_intersect(segments, bgplus, strandness = True)
				res += bedgraph_per_gene_intersect(segments, cond2bgminus[cond][b], strandness = True)
				BedPair.save(res, intersect)
			else:
				bedgraph_per_gene_intersect(segments, bgplus, intersect, strandness = True)

			logger.info(f'avg coverage per bp for each segment {datetime.now().time()}')
			get_seg2cov(intersect, cond, sample, outfile='_'.join([output, sample, cond + '_readCounts.bed']))

			# remove temporary file
			os.remove(intersect)


def main(argv):
	# --------------------------------------------------
	# get args
	# --------------------------------------------------
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description='Calculate the average reads/bp for each segment after clustering.')
	group = parser.add_argument_group('Input')
	group.add_argument('-i', '--segments', dest='segments', type=str, metavar='', help='_segments.bed output from clusterCPs')
	group.add_argument('-p', '--bgplus', dest='bgplus', type=str, nargs='*', metavar='', help='List of space-delimited bedgraphs: non-strand-specific or plus strand.')
	group.add_argument('-m', '--bgminus', dest='bgminus', type=str, nargs='*', metavar='', help='List of space-delimited bedgraphs: minus strand.')
	group.add_argument('-c', '--conditions', dest='conditions', type=str, nargs='*', metavar='', help='List of space-delimited condition labels for each --bgplus file')
	group = parser.add_argument_group('Output')
	group.add_argument('-o', '--output', dest='output', type=str, metavar='',
		help='Output prefix. Outputs one _readCounts.bed file per sample.')
	args = parser.parse_args()
	logger.debug(args)
	logger.info(f'job starting: {str(datetime.now().time())}')

	read_count(
		segments=args.segments, 
		conditions=args.conditions, 
		bgplus=args.bgplus, 
		bgminus=args.bgminus, 
		output=args.output
	)
	logger.info(f'finished: {datetime.now().time()}')


# boilerplate
if __name__ == '__main__':
	main(sys.argv[1:])
