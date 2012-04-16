#!/usr/bin/env python
import sys
from PftScan.Fasta import Fasta
from PftScan.core import scan
from PftScan.utils import *
from optparse import OptionParser,TitledHelpFormatter

VERSION = "1.05"

pwm_left = [
[0.282487110278,0.0977411337588,0.374860098537,0.244911657426],
[0.406620302733,0.0295237937609,0.440169814207,0.123686089299],
[2.23663409829e-07,0.968910115044,0.00626279913863,0.0248268621545],
[0.692909467314,0.116752523594,0.0733618220874,0.116976187004],
[0.296354241687,0.0185642866792,0.150973025298,0.534108446336],
[2.23663409829e-07,2.23663409829e-07,0.99999932901,2.23663409829e-07],
[0.103332719005,0.382911981291,0.00559180890914,0.508163490795],
[0.164392829888,0.439498823978,0.0883472705459,0.307761075588],
[0.160814215331,0.333258704309,0.202191946149,0.303735134211]
]

pwm_right = [
[0.309550382867,0.199731648641,0.329009099522,0.16170886897],
[0.307761075588,0.0885709339558,0.443748428765,0.159919561691],
[0.511742105353,0.00671012595829,0.380675347193,0.100872421496],
[2.23663409829e-07,0.99999932901,2.23663409829e-07,2.23663409829e-07],
[0.532990129286,0.151420352118,0.0172223062203,0.298367212376],
[0.116081533365,0.0740328123169,0.115857869955,0.694027784363],
[0.0259451792036,0.00648646254846,0.967568134585,2.23663409829e-07],
[0.123909752709,0.43659119965,0.029300130351,0.41019891729],
[0.237307101492,0.379109703324,0.0975174703489,0.286065724835],
[0.199284321821,0.302840480572,0.158801244642,0.339073952964]
]

nreport = 1 
min_spacer = 0
max_spacer = 0 
#cutoffs = [-100,13.7,14.7,15.9,14.8,13.8,15.2,15.6,13.4,14.3,13.0,14.7,13.4,15.5]
cutoffs = [4.737]

id = "p63bs"
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
			"p63scan", 
			"misc_feature", 
			pos, pos + len(pwm_left) + spacer + len(pwm_right), 
			score, 
			strand, 
			".", 
			seq[pos: pos + len(pwm_left)],
			seq[pos + len(pwm_left): pos + len(pwm_left) + spacer],
			seq[pos + len(pwm_left) + spacer: pos + len(pwm_left) + spacer + len(pwm_right)]
		)


