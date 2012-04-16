#!/usr/bin/env python
import sys
from PftScan.Fasta import Fasta
from PftScan.core import scan
from PftScan.utils import *
from optparse import OptionParser,TitledHelpFormatter

VERSION = "1.05"

pwm_left = [
[0.3611,0.0769,0.4003,0.1664],
[0.2716,0.0283,0.5667,0.1381],
[0.6358,0.0016,0.3344,0.0330],
[0.0016,0.9859,0.0016,0.0157],
[0.8085,0.0063,0.0502,0.1397],
[0.2276,0.0157,0.0330,0.7284],
[0.0031,0.0016,0.9984,0.0016],
[0.0377,0.3799,0.0016,0.5856],
[0.0816,0.7096,0.0173,0.1962],
[0.1350,0.4035,0.0675,0.3987],
]

pwm_right = [
[0.4035,0.0534,0.4129,0.1350],
[0.1743,0.0126,0.7394,0.0785],
[0.6044,0.0016,0.3673,0.0314],
[0.0016,0.9984,0.0016,0.0031],
[0.6970,0.0298,0.0188,0.2590],
[0.1381,0.0471,0.0078,0.8116],
[0.0173,0.0031,0.9827,0.0016],
[0.0377,0.3501,0.0031,0.6138],
[0.1413,0.5385,0.0235,0.3014],
[0.1570,0.4066,0.0911,0.3501],
]

nreport = 1 
min_spacer = 0
max_spacer = 13
#cutoffs = [-100,13.7,14.7,15.9,14.8,13.8,15.2,15.6,13.4,14.3,13.0,14.7,13.4,15.5]
cutoffs = [4.393,13.7,14.7,15.9,14.8,13.8,15.2,15.6,13.4,14.3,13.0,14.7,13.4,15.5]

id = "p53bs"
SCORE_FRACTION = 0.85

usage = "usage: %prog -i <FILE> [optional arguments]"
parser = OptionParser(version=VERSION, usage=usage, formatter=TitledHelpFormatter(max_help_position=40, short_first=1))
parser.add_option("-i", "--input", dest="inputfile", help="FASTA-formatted inputfile", metavar="FILE")
parser.add_option("-n", "--nreport", dest="nreport", help="report the N best matches", metavar="N")
parser.add_option("-p", "--pwm", dest="pwmfile", help="specify your own PWM file", metavar="FILE")
parser.add_option("-s", "--spacer", dest="spacer", help="use spacer length L (either specify one spacer, or use the format min-max)", metavar="L")
parser.add_option("-c", "--cutoff", dest="cutoff", help="cutoff or comma-separated list of cutoffs, one for each spacer length", metavar="")
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
	cutoffs = map(float, options.cutoff.split(","))
	if len(cutoffs) != max_spacer - min_spacer + 1:
		print "Please specify same number of cutoffs as spacers, separated by commas"
		sys.exit(0)

if options.pwmfile:
	id, pwm_left,pwm_right = read_pwm_with_spacer(options.pwmfile) 
	sys.stderr.write("PWM-file supplied, left halfsite: %s, right halfsite: %s\n" % (len(pwm_left), len(pwm_right)))
	if not options.cutoff:
		
		s_max = max_score(pwm_left) + max_score(pwm_right) 
		s_min = min_score(pwm_left) + min_score(pwm_right)
		sys.stderr.write("No cutoffs specified, using %s as a cutoff (%s of maxscore)\n" % (SCORE_FRACTION * (s_max - s_min) + s_min, SCORE_FRACTION))
		cutoffs = [SCORE_FRACTION * (s_max - s_min) + s_min for i in range(max_spacer - min_spacer + 1)]
	
bed = options.bed
f = Fasta(inputfile)

for (id,seq) in f.items():
	result = scan(seq.upper(), pwm_left, pwm_right, min_spacer, max_spacer, cutoffs, nreport)
	for (score, pos, spacer, strand) in result:
		if bed:
			first = id.split(" ")[0]	
			(chr,loc) = first.split(":")
			if loc:
				(start, end) = map(int, loc.split("-"))
				print "%s\t%s\t%s\t%s" % (chr, start + pos, start + pos + len(pwm_left) + spacer + len(pwm_right), score)
			else:
				print "%s\t%s\t%s\t%s" % (id, pos, pos +  len(pwm_left) + spacer + len(pwm_right), score)
		else:
			print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s %s %s" % (
			id, 
			"p53scan", 
			"misc_feature", 
			pos, pos + len(pwm_left) + spacer + len(pwm_right), 
			score, 
			strand, 
			".", 
			seq[pos: pos + len(pwm_left)],
			seq[pos + len(pwm_left): pos + len(pwm_left) + spacer],
			seq[pos + len(pwm_left) + spacer: pos + len(pwm_left) + spacer + len(pwm_right)]
		)
