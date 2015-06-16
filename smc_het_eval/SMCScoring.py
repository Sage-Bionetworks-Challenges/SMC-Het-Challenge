import math
import numpy as np
import itertools
import sys
import json
import argparse
import StringIO

class ValidationError(Exception):
	def __init__(self, value):
		self.value = value
	def __str__(self):
		return repr(self.value)

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
	return (truth + 1 - min(truth+1,abs(pred-truth))) / float(truth+1)

def validate1C(data, nssms):
	data = data.split('\n')
	data = filter(None,data)
	data = [x.strip() for x in data]
	if len(data) < 1:
		raise ValidationError("Number of lines is less than 1")
	elif len(data) > 10:
		raise ValidationError("Number of lines is greater than 10")

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
	reported_nssms = sum([int(x[1]) for x in data2])
	if reported_nssms != nssms:
		raise ValidationError("Total number of reported mutations is %d.  Should be %d" % (reported_nssms,nssms))
	return zip([int(x[1]) for x in data2], [float(x[2]) for x in data2])

def calculate1C(pred,truth):
	pred.sort(key = lambda x: x[1])
	truth.sort(key = lambda x: x[1])
	#itertools.chain(*x) flattens a list of lists
	predvs = np.array(list(itertools.chain(*[[x[1]]*x[0] for x in pred])))
	truthvs = np.array(list(itertools.chain(*[[x[1]]*x[0] for x in truth])))
	se = abs(truthvs - predvs)
	return sum(1-se)/float(len(truthvs))

def validate2A(data,ssmlist,return_ccm=True):
	data = data.split('\n')
	data = filter(None,data)
	if len(data) != len(ssmlist):
		raise ValidationError("Input file contains a different number of lines than the specification file")
	data = [x.split('\t') for x in data]
	for i in range(len(data)):
		if len(data[i]) != 2:
			raise ValidationError("Number of tab separated columns in line %d is not 2" % (i+1))
	for i in range(len(data)):
		if data[i][0] != ssmlist[i]:
			raise ValidationError("Starting at line %d, mutation names are not consistent with specification file" % (i+1))
	cluster_entries = []
	for i in range(len(data)):
		try:
			cluster_n = int(data[i][1])
			cluster_entries.append(cluster_n)
			data[i][1] = cluster_n
		except ValueError:
			raise ValidationError("Cluster ID in line %d (ssm %s) can not be cast as an integer" % (i+1,data[i][0]))
	used_clusters = sorted(list(set(cluster_entries)))
	expected_clusters = range(1,len(set(cluster_entries)) +1)

	if used_clusters != expected_clusters:
		raise ValidationError("Cluster IDs used (%s) is not what is expected (%s)" % (str(used_clusters), str(expected_clusters)))

	c_m = np.zeros((len(data),len(set(cluster_entries))))
	for i in range(len(data)):
		c_m[i,data[i][1]-1] = 1
	if not return_ccm:
		return c_m
	else:
		ccm = np.dot(c_m,c_m.T)
		return ccm

def validate2Afor3A(data,ssmlist):
	return validate2A(data,ssmlist,False)

def calculate2(pred,truth):
	n = truth.shape[0]
	indices = np.triu_indices(n,k=1)
	count = (n**2 - n )/2.0
	res = np.sum(np.abs(pred[indices] - truth[indices])) 
	res = res / count
	return 1 - res

def validate2B(data,ssmlist):
	data = StringIO.StringIO(data)
	try:
		ccm = np.loadtxt(data,ndmin=2)
	except ValueError:
		raise ValidationError("Entry in co-clustering matrix could not be cast as a float")

	if ccm.shape != (len(ssmlist),len(ssmlist)):
		raise ValidationError("Shape of co-clustering matrix %s is wrong.  Should be %s" % (str(ccm.shape), str((len(ssmlist),len(ssmlist)))))
	if not np.allclose(ccm.diagonal(),np.ones((len(ssmlist)))):
		raise ValidationError("Diagonal entries of co-clustering matrix not 1")
	if np.any(np.isnan(ccm)):
		raise ValidationError("Co-clustering matrix contains NaNs")
	if np.any(np.isinf(ccm)):
		raise ValidationError("Co-clustering matrix contains non-finite entries")
	if np.any(ccm > 1):
		raise ValidationError("Co-clustering matrix contains entries greater than 1")
	if np.any(ccm < 0):
		raise ValidationError("Co-clustering matrix contains entries less than 0")
	if not np.allclose(ccm.T, ccm):
		raise ValidationError("Co-clustering matrix is not symmetric")
	return ccm

def validate3A(data, cas, ssmlist):
	predK = cas.shape[1]
	cluster_assignments = np.argmax(cas,1) + 1

	data = data.split('\n')
	data = filter(None,data)
	if len(data) != predK:
		raise ValidationError("Input file contains a different number of lines (%d) than expected (%d)")
	data = [x.split('\t') for x in data]
	for i in range(len(data)):
		if len(data[i]) != 2:
			raise ValidationError("Number of tab separated columns in line %d is not 2" % (i+1))
		try:
			data[i][0] = int(data[i][0])
			data[i][1] = int(data[i][1])
		except ValueError:
			raise ValidationError("Entry in line %d could not be cast as integer" % (i+1))
	
	if [x[0] for x in data] != range(1,predK+1):
		raise ValidationError("First column must have %d entries in acending order starting with 1" % predK)

	for i in range(len(data)):
		if data[i][1] not in set(range(predK+1)):
			raise ValidationError("Parent node label in line %d is not valid." % (i+1))

	# Form decendent of dict.  Each entry, keyed by cluster number, consists of a list of nodes that are decendents of the key. 
	decendent_of = dict()
	for i in range(predK+1):
		decendent_of[i] = []
	for child,parent in data:
		decendent_of[parent] += [child] + decendent_of[child]
		# gps (grandparents) are the list of nodes that are ancestors of the immediate parent
		gps = [x for x in decendent_of.keys() if parent in decendent_of[x]]
		for gp in gps:
			decendent_of[gp] += [child] + decendent_of[child]

	# Check that root has all nodes as decendants (equivalent to checking if the tree is connected)
	if set(decendent_of[0]) != set(range(1,predK+1)):
		print data
		print decendent_of
		raise ValidationError("Root of phylogeny not ancestor of all clusters / Tree is not connected")

	# Form AD matrix
	n = len(cluster_assignments)
	ad = np.zeros((n,n))
	for i in range(n):
		for j in range(n):
			if cluster_assignments[j] in decendent_of[cluster_assignments[i]]:
				ad[i,j] = 1
	return ad

def validate3B(data, ccm, ssmlist):
	data = StringIO.StringIO(data)
	k = ccm.shape[0]
	try:
		ad = np.loadtxt(data,ndmin=2)
	except ValueError:
		raise ValidationError("Entry in AD matrix could not be cast as a float")

	if ad.shape != ccm.shape:
		raise ValidationError("Shape of AD matrix %s is wrong.  Should be %s" % (str(ad.shape), str(ccm.shape)))
	if not np.allclose(ad.diagonal(),np.zeros(ccm.shape[0])):
		raise ValidationError("Diagonal entries of AD matrix not 0")
	if np.any(np.isnan(ad)):
		raise ValidationError("AD matrix contains NaNs")
	if np.any(np.isinf(ad)):
		raise ValidationError("AD matrix contains non-finite entries")
	if np.any(ad > 1):
		raise ValidationError("AD matrix contains entries greater than 1")
	if np.any(ad < 0):
		raise ValidationError("AD matrix contains entries less than 0")
	indices = np.triu_indices(k,1)
	if np.any(ad[indices] + ad.T[indices] + ccm[indices] > 1):
		raise ValidationError("For some i,j the sum of AD(i,j) + AD(j,i) + CCM(i,j) > 1.")

	return ad

def calculate3A(pred_ca, pred_ad, truth_ca, truth_ad):
	return calculate3(np.dot(pred_ca,pred_ca.T),pred_ad,np.dot(truth_ca,truth_ca.T),truth_ad)

def calculate3(pred_ccm, pred_ad, truth_ccm, truth_ad):
	n = truth_ccm.shape[0]
	indices = np.triu_indices(n,k=1)
	cc_res = np.sum(np.abs(pred_ccm[indices] - truth_ccm[indices])*truth_ccm[indices])
	ad_res = np.sum(np.abs(pred_ad.flatten() - truth_ad.flatten()) * truth_ad.flatten())
	truth_cous = 1 - truth_ccm[indices] - truth_ad[indices] - truth_ad.T[indices]
	pred_cous = 1 - pred_ccm[indices] - pred_ad[indices] - pred_ad.T[indices]
	cous_res = np.sum(np.abs(pred_cous - truth_cous) * truth_cous)
	count = (n**2 - n )*2.0
	res = (cc_res + ad_res + cous_res ) / count
	return 1 - res


def validateNssms(data):
	return int(data)

def validateSsmList(data):
	ssmlist = data.split('\n')
	nssms = int(ssmlist[0])
	ssmlist = ssmlist[1:]
	ssmlist = filter(None,ssmlist)
	if len(ssmlist) != nssms:
		raise ValidationError("Error in Spec file")
	return ssmlist

def verify(filename,role,func,*args):
	try:
		f = open(filename)
		pred_data = f.read(10**6)
		f.close()
		pred = func(pred_data,*args)
	except (IOError,TypeError) as e:
		print "Error opening " + role
		print e
		return None
	except (ValidationError,ValueError) as e:
		print role + " does not validate"
		print e
		return None
	return pred


challengeMapping = { 	'1A': {'spec': None, 'val_funcs':[validate1A],'score_func':calculate1A},
						'1B': {'spec': None, 'val_funcs':[validate1B],'score_func':calculate1B},
						'1C': {'spec': validateNssms, 'val_funcs':[validate1C],'score_func':calculate1C},
						'2A': {'spec': validateSsmList, 'val_funcs':[validate2A],'score_func':calculate2},
						'2B': {'spec': validateSsmList, 'val_funcs':[validate2B],'score_func':calculate2},
						'3A': {'spec': validateSsmList, 'val_funcs':[validate2Afor3A,validate3A],'score_func':calculate3A},
						'3B': {'spec': validateSsmList, 'val_funcs':[validate2B,validate3B],'score_func':calculate3},
					}

def verifyChallenge(challenge,predfiles,specfile):
	if challengeMapping[challenge]['spec']:
		spec = [verify(specfile,"specification file for Challenge %s" % (challenge), challengeMapping[challenge]['spec'])]
		if spec == [None]:
			print "Could not read specification file. Exiting"
			return "NA"
	else:
		spec = []
	if len(predfiles) != len(challengeMapping[challenge]['val_funcs']):
		print "Not enough input files for Challenge %s" % challenge
		return "Invalid"

	out = []
	for (predfile,valfunc) in zip(predfiles,challengeMapping[challenge]['val_funcs']):
		args = out + spec
		out.append(verify(predfile, "prediction file for Challenge %s" % (challenge),valfunc,*args))
		if out[-1] == None:
			return "Invalid"
	return "Valid"


def scoreChallenge(challenge,predfiles,truthfiles,specfile):
	if challengeMapping[challenge]['spec']:
		spec = [verify(specfile,"specification file for Challenge %s" % (challenge), challengeMapping[challenge]['spec'])]
		if spec == [None]:
			print "Could not read specification file. Exiting"
			return "NA"	
	else:
		spec=[]
	if len(predfiles) != len(challengeMapping[challenge]['val_funcs']) or len(truthfiles) != len(challengeMapping[challenge]['val_funcs']):
		print "Not enough input files for Challenge %s" % challenge
		return "NA"

	tout = []
	pout = []
	for predfile,truthfile,valfunc in zip(predfiles,truthfiles,challengeMapping[challenge]['val_funcs']):
		targs = tout + spec
		tout.append(verify(truthfile, "truth file for Challenge %s" % (challenge),valfunc,*targs))
		pargs = pout + spec
		pout.append(verify(predfile, "prediction file for Challenge %s" % (challenge),valfunc,*pargs))
		if tout[-1] == None or pout[-1] == None:
			return "NA"

	return challengeMapping[challenge]['score_func'](*(pout + tout))

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("challenge")
	parser.add_argument("--predfiles",nargs="+")
	parser.add_argument("--truthfiles",nargs="*")
	parser.add_argument("--specfile",nargs="?")
	parser.add_argument("outputfile")
	parser.add_argument('-v', action='store_true', default=False)

	args = parser.parse_args()

	if args.v:
		res = verifyChallenge(args.challenge,args.predfiles,args.specfile)
	else:
		res = scoreChallenge(args.challenge,args.predfiles,args.truthfiles,args.specfile)

	with open(args.outputfile, "w") as handle:
		jtxt = json.dumps( { args.challenge : res } )
		handle.write(jtxt)
