#!/usr/bin/env python

import math
import numpy as np
import itertools
import sys
import json
import argparse

class ValidationError(Exception):
	def __init__(self, value):
		self.value = value
	def __str__(self):
		return repr(self.value)

###subC1A
def validate1A(data):
	data = data.split('\n')
	data = filter(None,data)
	if len(data) < 1:
		raise ValidationError("Input file contains zero lines")
	if len(data) > 1:
		raise ValidationError("Input file contains more than one line")
	data = data[0].strip()
	try:
		numeric = float(data)
	except ValueError:
		raise ValidationError("Data could not be converted to float: %s" % data)
	if math.isinf(numeric):
		raise ValidationError("Non-finite Cellularity")
	if math.isnan(numeric):
		raise ValidationError("Cellularity is NaN")
	if numeric < 0:
		raise ValidationError("Cellularity was < 0: %f" % numeric)
	if numeric > 1:
		raise ValidationError("Cellularity was > 1: %f" % numeric)

	return numeric

def calculate1A(pred,truth):
	return 1 - abs(truth - pred)

### subC1B
def validate1B(data):
	data = data.split('\n')
	data = filter(None,data)
	if len(data) != 1:
		if len(data) == 0:
			raise ValidationError("Input file contains zero lines")
		else:
			raise ValidationError("Input file contains more than one line")
	data = data[0].strip()
	try:
		numeric = int(data)
	except ValueError:
		raise ValidationError("Data could not be converted to int: %s" % data)
	if numeric < 1:
		raise ValidationError("Number of lineages was less than 1: %d" % numeric)
	if numeric > 20:
		raise ValidationError("Number of lineages was greater than 20: %d" % numeric)
	return numeric

def calculate1B(pred,truth):
	return abs(pred-truth) / float(truth)

### subC1C
def validate1C(data):
	data = data.split('\n')
	data = filter(None,data)
	data = [x.strip() for x in data]
	#if len(data) != K:
	#	raise ValidationError("Number of lines (%d) is different than number of clusters (%d)" % (len(data),K))
	data2 = [x.split('\t') for x in data]
	for i in range(len(data)):
		if len(data2[i]) != 3:
			raise ValidationError("Number of tab separated columns in line %d is not 3" % (i+1))
		try:
			id = int(data2[i][0])
			if id != i+1:
				raise ValidationError("Cluster ID in line %d is not %d" % (i+1,i+1))
		except ValueError:
			raise ValidationError("Cluster ID in line %d can not be cast as an integer: %s" % (i+1,data2[i][0]))
		try:
			nm = int(data2[i][1])
			if nm < 1:
				raise ValidationError("Number of mutations in line %d is less than 1." % (i+1))
		except ValueError:
			raise ValidationError("Number of mutations in line %d can not be cast as an integer: %s" % (i+1,data2[i][1]))
		try:
			cf = float(data2[i][2])
			if math.isinf(cf):
				raise ValidationError("Cellular Frequency for cluster %d is non-finite" % (i+1))
			if math.isnan(cf):
				raise ValidationError("Cellular Frequency for cluster %d is NaN" % (i+1))
			if cf < 0:
				raise ValidationError("Cellular Frequency for cluster %d is negative: %f" % (i+1,cf))
			if cf > 1:
				raise ValidationError("Cellular Frequency for cluster %d is > 1: %f" % (i+1,cf))

		except ValueError:
			raise ValidationError("Cellular Frequency for cluster %d can not be cast as a float: %s" % (i+1,data2[i][2]))
	#reported_nssms = sum([int(x[1]) for x in data2])
	#if reported_nssms != nssms:
	#	raise ValidationError("Total number of reported mutations is %d.  Should be %d" % (reported_nssms,nssms))
	return zip([int(x[1]) for x in data2], [float(x[2]) for x in data2])

def calculate1C(pred,truth):
	pred.sort(key = lambda x: x[1])
	truth.sort(key = lambda x: x[1])
	#itertools.chain(*x) flattens a list of lists
	predvs = np.array(list(itertools.chain(*[[x[1]]*x[0] for x in pred])))
	truthvs = np.array(list(itertools.chain(*[[x[1]]*x[0] for x in truth])))
	se = (truthvs - predvs)**2
	return sum(1-se)/float(len(truthvs))

# To cut down on some redundancy, this function accepts a filename, the role of the file (e.g. "prediction file") and
# a function pointer. The function pointer will be called on the file, alone with any additional arguments provided to
# the verify function.
def verify(filename,role,func,*args):
	try:
		f = open(filename)
		pred_data = f.read(10**6)
		f.close()
		pred = func(pred_data,*args)
	except IOError as e:
		print "Error opening " + role
		print e
		return None
	except ValidationError as e:
		print role + "does not validate"
		print e
		return None
	return pred;

def validatesubC1(specfile,pred1A,pred1B,pred1C):
	try:
		f = open(specfile)
		spec_data = f.read(10**6)
		f.close()
	except IOError:
		print "Error opening spec file. Exiting"
		return ["NA","NA","NA"]
	try:
		nssms = int(spec_data)
	except ValueError:
		print "Error reading spec file. Exiting"
		return ["NA","NA","NA"]

	pred1A = verify(pred1A, "input file for Challenge 1A", validate1A);

	if pred1A == None:
		return ["Invalid", "NA", "NA"]

	pred1B = verify(pred1B, "input file for Challenge 1B", validate1B);

	if pred1B == None:
		return ["Valid", "Invalid", "NA"]

	pred1C = verify(pred1C, "input file for Challenge 1C", validate1C, pred1B, nssms);

	if pred1C == None:
		return ["Valid", "Valid", "Invalid"]

	return ["Valid", "Valid", "Valid"]

def scoresubC1(specfile,pred1Afile,truth1Afile,pred1Bfile,truth1Bfile,pred1Cfile,truth1Cfile):
	try:
		f = open(specfile)
		spec_data = f.read(10**6)
		f.close()
	except IOError:
		print "Error opening spec file. Exiting"
		return ["NA","NA","NA"]
	try:
		nssms = int(spec_data)
	except ValueError:
		print "Error reading spec file. Exiting"
		return ["NA","NA","NA"]

	pred1A = verify(pred1Afile, "prediction file for Challenge 1A", validate1A);
	truth1A = verify(truth1Afile, "truth file for Challenge 1A", validate1A);
	[score1A, score1B, score1C] = ["NA","NA","NA"]

	if pred1A == None or truth1A == None:
		return [score1A, score1B, score1C]

	score1A = calculate1A(pred1A, truth1A)

	pred1B = verify(pred1Bfile, "prediction file for Challenge 1B", validate1B);
	truth1B = verify(truth1Bfile, "truth file for Challenge 1B", validate1B);

	if pred1B == None or truth1B == None:
		return [score1A, score1B, score1C]

	score1B = calculate1B(pred1B, truth1B)

	pred1C = verify(pred1Cfile, "prediction file for Challenge 1C", validate1C, pred1B, nssms);
	truth1C = verify(truth1Cfile, "truth file for Challenge 1C", validate1C, truth1B, nssms);

	if pred1C == None or truth1C == None:
		return [score1A, score1B, score1C]

	score1C = calculate1C(pred1C, truth1C)

	return [score1A, score1B, score1C]

def scoresubC1Alt(challenge,predfile,truthfile,validatefunc,calculatefunc, specfile=None):
	if specfile is not None:
		try:
			f = open(specfile)
			spec_data = f.read(10**6)
			f.close()
		except IOError:
			print "Error opening spec file. Exiting"
			return "NA"
		try:
			nssms = int(spec_data)
		except ValueError:
			print "Error reading spec file. Exiting"
			return "NA"

	pred = verify(predfile, "prediction file for Challenge %s" % (challenge), validatefunc)
	truth = verify(truthfile, "truth file for Challenge %s" % (challenge), validatefunc)

	if pred == None or truth == None:
		return "NA"

	score = calculatefunc(pred, truth)

	return score


evaluationMapping = {
	"1A" : (validate1A, calculate1A),
	"1B" : (validate1B, calculate1B),
	"1C" : (validate1C, calculate1C)
}

def main(scoring_args,output_file):
	scores = scoresubC1(*scoring_args)
	scores = dict(zip(['subC1A','subC1B','subC1C'],scores))
	try:
		f = open(output_file,'w')
		json.dump(scores,f)
		f.close()
	except IOError:
		print "Error opening/writing output file"


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("challenge")
	parser.add_argument("predfile")
	parser.add_argument("truthfile")
	parser.add_argument("outputfile")

	args = parser.parse_args()

	score = scoresubC1Alt(
		args.challenge, args.predfile, args.truthfile,
		evaluationMapping[args.challenge][0],
		evaluationMapping[args.challenge][1]
	)
	with open(args.outputfile, "w") as handle:
		jtxt = json.dumps( { args.challenge : score } )
		handle.write(jtxt)
