from math import log
import re
import sys

def max_score(pwm):
	score = 0
	for row in pwm:
		score += log(max(row) / 0.25 + 0.01)
	return score

def min_score(pwm):
	score = 0
	for row in pwm:
		score += log(min(row) / 0.25 + 0.01)
	return score


class Motif: 
	def __init__(self, name): 
		self.pwm = [] 
		self.name = name 

	def append(self, a): 
		self.pwm.append(a) 

	def __repr__(self): 
		return ">%s\n" % (self.name) + "\n".join(["\t".join(map(str, row)) for row in self.pwm]) 

class MDModuleResult: 
	def __init__(self, file): 
		p = re.compile(r'^\d+\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)') 
		pmotifnum = re.compile(r'Motif\.(\d+)\.(\d+)') 
		motifnum = 0 
		self.motifs = [] 
		for line in open(file): 
			m = pmotifnum.search(line) 
			if m: 
				motiflength = int(m.group(1))
				motifnum = int(m.group(2)) 
				self.motifs.append(Motif("Motif%s.%s" % (motiflength,motifnum))) 
			m = p.search(line) 
			if m: 
				self.motifs[-1].append(map(lambda x: float(x) / 100, (m.group(1), m.group(2),m.group(3),m.group(4)))) 

def read_pwm_with_spacer(file):
	f = open(file)
	line = f.readline()
	if (not(line.startswith(">"))):
		sys.stderr.write("Error: PWM file not in the right format!")
		sys.exit(1)
	id = line.strip()[1:]
	pwm_left = []
	line = f.readline()
	while line and not line.startswith(">"):
		try:
			pwm_left.append(map(float, line.strip().split("\t")))
		except:
			sys.stderr.write("Error reading PWM file, first halfsite! %s\n" % line.strip())
			sys.exit()
		line = f.readline()
	if not line:	
		sys.stderr.write("Error: found only one halfsite! Two halfsites should be present. Use pwmscan.py for motifs without a spacer!")
		sys.exit()
	pwm_right = []
	line = f.readline()
	while line:
		try:
			pwm_right.append(map(float, line.strip().split("\t")))
		except:
			sys.stderr.write("Error reading PWM file, second halfsite! %s\n" % line.strip())
			sys.exit()
		line = f.readline()
	f.close()
	return id,pwm_left,pwm_right

def read_pwm_triplet(file):
	f = open(file)
	line = f.readline()
	if (not(line.startswith(">"))):
		sys.stderr.write("Error: PWM file not in the right format!")
		sys.exit(1)
	id = line.strip()[1:]
	pwm_left = []
	line = f.readline()
	while line and not line.startswith(">"):
		try:
			pwm_left.append(map(float, line.strip().split("\t")))
		except:
			sys.stderr.write("Error reading PWM file, first PWM! %s\n" % line.strip())
			sys.exit()
		line = f.readline()
	if not line:	
		sys.stderr.write("Error: found only one PWM! Three PWMs should be present. Use pwmscan.py for motifs without a spacer!")
		sys.exit()
	
	pwm_middle = []
	line = f.readline()
	while line and not line.startswith(">"):
		try:
			pwm_middle.append(map(float, line.strip().split("\t")))
		except:
			sys.stderr.write("Error reading PWM file, second PWM! %s\n" % line.strip())
			sys.exit()
		line = f.readline()
	
	pwm_right = []
	line = f.readline()
	while line and not line.startswith(">"):
		try:
			pwm_right.append(map(float, line.strip().split("\t")))
		except:
			sys.stderr.write("Error reading PWM file, third PWM! %s\n" % line.strip())
			sys.exit()
		line = f.readline()

	f.close()
	return id,pwm_left,pwm_middle,pwm_right
