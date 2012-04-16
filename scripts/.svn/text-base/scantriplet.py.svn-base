#!/usr/bin/python2.6
import sys
from PftScan.Fasta import Fasta
from PftScan.core import scan_triplet
from PftScan.utils import *
from optparse import OptionParser,TitledHelpFormatter

VERSION = "1.02"



nreport = 1 
min_spacer = 0
max_spacer = 13
cutoff = 0

id = "mattch"
SCORE_FRACTION = 0.85

usage = "usage: %prog -i <FILE> [optional arguments]"
parser = OptionParser(version=VERSION, usage=usage, formatter=TitledHelpFormatter(max_help_position=40, short_first=1))
parser.add_option("-i", "--input", dest="inputfile", help="FASTA-formatted inputfile", metavar="FILE")
parser.add_option("-n", "--nreport", dest="nreport", help="report the N best matches", metavar="N")
parser.add_option("-p", "--pwm", dest="pwmfile", help="specify your PWM file", metavar="FILE")
parser.add_option("-s", "--spacer", dest="spacer", help="use spacer length L (either specify one spacer, or use the format min-max)", metavar="L")
parser.add_option("-c", "--cutoff", dest="cutoff", help="cutoff", metavar="", type="float", default="0")
parser.add_option("-b", "--bed", action="store_true", dest="bed", default=False, help="output bed format")

(options, args) = parser.parse_args()

if not options.inputfile:
	parser.print_help()
	sys.exit(0)

inputfile = options.inputfile

if options.nreport:
	nreport = int(options.nreport)

if options.spacer:
	if options.spacer.find("-") != -1:
		min_spacer,max_spacer = [int(x.strip()) for x in options.spacer.split("-")]
	else:
		min_spacer = int(options.spacer) 
		max_spacer = int(options.spacer)

if options.cutoff:
	cutoff = options.cutoff

if options.pwmfile:
	id, pwm_left,pwm_middle,pwm_right = read_pwm_triplet(options.pwmfile) 
	sys.stderr.write("PWM-file supplied, left site: %s, middle site: %s, right site: %s\n" % (len(pwm_left), len(pwm_middle), len(pwm_right)))
	if not options.cutoff:
		s_max = max_score(pwm_left) + max_score(pwm_right) + max_score(pwm_middle) 
		s_min = min_score(pwm_left) + min_score(pwm_right) + min_score(pwm_middle)
		sys.stderr.write("No cutoffs specified, using %s as a cutoff (%s of maxscore)\n" % (SCORE_FRACTION * (s_max - s_min) + s_min, SCORE_FRACTION))
		cutoff = SCORE_FRACTION * (s_max - s_min) + s_min 
	
bed = options.bed
f = Fasta(inputfile)

for (id,seq) in f.items():
	result = scan_triplet(seq.upper(), pwm_left, pwm_middle, pwm_right, min_spacer, max_spacer, cutoff, nreport)
	for (score, pos, spacer1, spacer2, strand) in result:
		if bed:
			first = id.split(" ")[0]	
			(chr,loc) = first.split(":")
			if loc:
				(start, end) = map(int, loc.split("-"))
				print "%s\t%s\t%s\t%s" % (chr, start + pos, start + pos + len(pwm_left) + spacer1 + len(pwm_middle) + spacer2 + len(pwm_right), score)
			else:
				print "%s\t%s\t%s\t%s" % (id, pos, pos + len(pwm_left) + spacer1 + len(pwm_middle) + spacer2 + len(pwm_right), score)
		else:
			print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s %s %s %s %s" % (
			id, 
			"p53scan", 
			"misc_feature", 
			pos, pos + len(pwm_left) + spacer1 + len(pwm_middle) + spacer2 + len(pwm_right), 
			score, 
			strand, 
			".", 
			seq[pos: pos + len(pwm_left)],
			seq[pos + len(pwm_left): pos + len(pwm_left) + spacer1],
			seq[pos + len(pwm_left) + spacer1: pos + len(pwm_left) + spacer1 + len(pwm_middle)],
			seq[pos + len(pwm_left) + spacer1 + len(pwm_middle): pos + len(pwm_left) + spacer1 + len(pwm_middle) + spacer2],
			seq[pos + len(pwm_left) + spacer1 + len(pwm_middle) + spacer2: pos + len(pwm_left) + spacer1 + len(pwm_middle) + spacer2 + len(pwm_right)],
		)
