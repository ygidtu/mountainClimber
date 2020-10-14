#!/user/bin/python -tt

import os
import sys

from subprocess import Popen, PIPE

import pybedtools as pb

from loguru import logger


def run_command(cmd, stdin=0, stdoutfile=0, printflag=False):
	"""run a subprocess"""
	if printflag:
		print (' '.join(cmd))
	if stdin == 0 and stdoutfile == 0:
		p = Popen(cmd, stdout=PIPE, stderr=PIPE)
		stdout, stderr = p.communicate()
	elif stdoutfile == 0:
		p = Popen(cmd, stdout=PIPE, stderr=PIPE, stdin=PIPE)
		stdout, stderr = p.communicate(stdin)
	else:
		p = Popen(cmd, stdout=stdoutfile, stderr=PIPE)
		p.wait()
		stdoutfile.flush()
		stdout = None
		stderr = None

	if p.returncode != 0:
		if stdoutfile == 0:
			logger.error('command failed')
			logger.error(stdout)
			logger.error(stderr)
			sys.exit(1)
		else:
			logger.error('command failed')
			logger.error(stderr)
			sys.exit(1)

	return stdout, stderr


def sort_bedfile(infile, outfile, add_header: bool = True, sort_by_bedtools: bool = False):
	"""
	sort bed file

	@2020.10.10 by Zhang Yiming: several modifications
	1. if infile and outfile is same, use a temp file
	2. add parameter to contol the bed header
	3. using check_output to better handle the output from command line tools
	
	"""
	if infile == outfile:
		outfile = outfile + ".sorted"

	if sort_by_bedtools:
		pb.BedTool(infile).sort().saveas(outfile)
	else:
		with open(outfile, 'w') as o:
			if add_header:
				o.write('track name=\"' + os.path.basename(infile).replace('.bed', '') + '\"\n')
			cmd = ['sort', '-k1,1', '-k2,2n', '-T', os.path.dirname(outfile), infile]
			stdout, _ = run_command(cmd)
			if stdout:
				o.write(stdout.decode("utf-8"))
				o.close()

	if outfile.endswith(".sorted"):
		os.rename(outfile, outfile.replace(".sorted", ""))
