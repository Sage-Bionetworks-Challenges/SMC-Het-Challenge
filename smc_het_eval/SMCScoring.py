import math
import numpy as np
import itertools
import sys
import json
import argparse
import StringIO
import scipy.stats
import sklearn.metrics as mt
import copy as cp


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

def validate2A(data,nssms,return_ccm=True):
    data = data.split('\n')
    data = filter(None,data)
    if len(data) != nssms:
        raise ValidationError("Input file contains a different number of lines than the specification file")
    cluster_entries = []
    for i in range(len(data)):
        try:
            cluster_n = int(data[i])
            cluster_entries.append(cluster_n)
            data[i] = cluster_n
        except ValueError:
            raise ValidationError("Cluster ID in line %d (ssm %s) can not be cast as an integer" % (i+1,data[i][0]))
    used_clusters = sorted(list(set(cluster_entries)))
    expected_clusters = range(1,len(set(cluster_entries)) +1)

    if used_clusters != expected_clusters:
        raise ValidationError("Cluster IDs used (%s) is not what is expected (%s)" % (str(used_clusters), str(expected_clusters)))

    c_m = np.zeros((len(data),len(set(cluster_entries))))
    for i in range(len(data)):
        c_m[i,data[i]-1] = 1
    if not return_ccm:
        return c_m
    else:
        ccm = np.dot(c_m,c_m.T)
        return ccm

def validate2Afor3A(data,nssms):
    return validate2A(data,nssms,False)

def calculate2_quaid(pred,truth):
	n = truth.shape[0]
	indices = np.triu_indices(n,k=1)
	ones = np.sum(np.abs(pred[indices] - truth[indices]) * truth[indices])
	ones_count = np.count_nonzero(truth[indices])
	if ones_count > 0:
		ones_score = 1 - ones/float(ones_count)
	else:
		ones_score = -1

	zeros = np.sum(np.abs(pred[indices] - truth[indices]) * (1 - truth[indices]))
	zeros_count = len(truth[indices]) - ones_count
	if zeros_count > 0:
		zeros_score = 1 - zeros/float(zeros_count)
	else:
		zeros_score = -1

	if ones_score == -1:
		return zeros_score
	elif zeros_score == -1:
		return ones_score
	else:
		try:
			return 2.0/(1.0/ones_score + 1.0/zeros_score)
		except Warning:
			print ones_score, zeros_score
			return 0

def calculate2(pred,truth):
	return calculate2_orig(pred,truth)

def calculate2_orig(pred,truth):
	n = truth.shape[0]
	indices = np.triu_indices(n,k=1)
	count = (n**2 - n )/2.0
	res = np.sum(np.abs(pred[indices] - truth[indices]))
	res = res / count
	return 1 - res

def calculate2_sqrt(pred,truth):
	n = truth.shape[0]
	indices = np.triu_indices(n,k=1)
	count = (n**2 - n )/2.0
	res = np.sum(np.abs(pred[indices] - truth[indices]))
	res = res / count
	return np.sqrt(1 - res)

def calculate2_simpleKL(pred,truth,rnd=0.01):
    pred = np.abs(pred - rnd)
    n = truth.shape[0]
    indices = np.triu_indices(n,k=1)
    res = 0
    for i in range(len(indices[0])):

        if truth[indices[0][i],indices[1][i]]:
            res += np.log(pred[indices[0][i],indices[1][i]])
        else:
            res += np.log(1-pred[indices[0][i],indices[1][i]])
    return abs(res)

def calculate2_pseudoV(pred,truth,rnd=0.01):
    pred_cp = cp.deepcopy(pred)
    truth_cp = cp.deepcopy(truth)
    # make matrix upper triangular
    pred_cp = np.triu(pred_cp)
    truth_cp = np.triu(truth_cp)
    # Avoid dividing by zero
    # Note: it is ok to do this after making the matrix upper triangular
    # since the bottom triangle of the matrix will not affect the score
    pred_cp[pred_cp==0] = rnd
    truth_cp[truth_cp==0] = rnd

    # normalize data
    pred_cp = pred_cp / np.sum(pred_cp,axis=1)[:,np.newaxis]
    truth_cp = truth_cp / np.sum(truth_cp,axis=1)[:,np.newaxis]

    return np.sum(truth_cp * np.log(truth_cp/pred_cp))

def calculate2_sym_pseudoV(pred, truth, rnd=0.01):
    pred[pred==0] = rnd
    truth[truth==0] = rnd
    pred = pred / np.sum(pred,axis=1)[:,np.newaxis]
    truth = truth / np.sum(truth,axis=1)[:,np.newaxis]
    return np.sum(truth * np.log(truth/pred)) + np.sum(pred * np.log(pred/truth))

def calculate2_pearson(pred, truth):
    n = truth.shape[0]
    inds = np.triu_indices(n,k=1)
    return scipy.stats.pearsonr(pred[inds],truth[inds])[0]

def calculate2_aupr(pred,truth):
    n = truth.shape[0]
    inds = np.triu_indices(n,k=1)
    precision, recall, thresholds = mt.precision_recall_curve(truth[inds],pred[inds])
    aucpr = mt.auc(recall, precision)
    return aucpr

# Matthews Correlation Coefficient
# don't just use upper triangular matrix because you get na's with the AD matrix
def calculate2_mcc(pred,truth):
    n = truth.shape[0]

    tp = float(sum(pred[truth==1] == 1))
    tn = float(sum(pred[truth==0] == 0))
    fp = float(sum(pred[truth==0] == 1))
    fn = float(sum(pred[truth==1] == 0))
    print tp,tn,fp,fn
    return (tp*tn - fp*fn)/np.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))

def validate2B(data,nssms):
    data = StringIO.StringIO(data)
    try:
        ccm = np.loadtxt(data,ndmin=2)
    except ValueError:
        raise ValidationError("Entry in co-clustering matrix could not be cast as a float")

    if ccm.shape != (nssms,nssms):
        raise ValidationError("Shape of co-clustering matrix %s is wrong.  Should be %s" % (str(ccm.shape), str((nssms,nssms))))
    if not np.allclose(ccm.diagonal(),np.ones((nssms))):
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

#### SUBCHALLENGE 3 #########################################################################################

def validate3A(data, cas, nssms):
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

def validate3B(data, ccm, nssms):
    data = StringIO.StringIO(data)
    k = ccm.shape[0]
    try:
        ad = np.loadtxt(data,ndmin=2)
    except ValueError:
        raise ValidationError("Entry in AD matrix could not be cast as a float")

    if ad.shape != (nssms,nssms):
        raise ValidationError("Shape of AD matrix %s is wrong.  Should be %s" % (str(ad.shape), str(ccm.shape)))
    if not np.allclose(ad.diagonal(),np.zeros(ad.shape[0])):
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

def calculate3(pred_ccm, pred_ad, truth_ccm, truth_ad, method="orig_no_cc"):
    func = {"orig" : calculate3_orig,
               "orig_no_cc": calculate3_orig_no_cc,
               "pseudoV": calculate3_pseudoV,
    }.get(method, calculate3_other_no_cc)
    return func(pred_ccm, pred_ad, truth_ccm, truth_ad, verbose=True, method=method)

# Include a method field to make code more easily generalizable
def calculate3_orig(pred_ccm, pred_ad, truth_ccm, truth_ad, verbose=False, method="orig"):
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

# Include a method field to make code more easily generalizable
def calculate3_orig_no_cc(pred_ccm, pred_ad, truth_ccm, truth_ad, verbose=False, method="orig_no_cc"):
    n = truth_ccm.shape[0]
    indices = np.triu_indices(n,k=1)
    ad_res = np.sum(np.abs(pred_ad.flatten() - truth_ad.flatten()) * truth_ad.flatten())
    truth_cous = 1 - truth_ccm[indices] - truth_ad[indices] - truth_ad.T[indices]
    pred_cous = 1 - pred_ccm[indices] - pred_ad[indices] - pred_ad.T[indices]
    cous_res = np.sum(np.abs(pred_cous - truth_cous) * truth_cous)
    count = (n**2 - n )*2.0 - np.sum(truth_ccm)
    res = (ad_res + cous_res ) / count
    return 1 - res
    
# Use the pseudoV measure
# Include a method field to make code more easily generalizable
def calculate3_pseudoV(pred_ccm, pred_ad, truth_ccm, truth_ad, rnd=0.01, verbose=False, method="pseudoV"):
    # calculate the result without the co-clustering matrix
    res = calculate3_other_no_cc(pred_ccm, pred_ad, truth_ccm, truth_ad, rnd=rnd, verbose=verbose, method="pseudoV_no_cc",)
        
    cpred_ccm = cp.deepcopy(pred_ccm)
    ctruth_ccm = cp.deepcopy(truth_ccm)
    
    # Calculate the pseudo V measure for each
    cc_res = calculate2_pseudoV(cpred_ccm, ctruth_ccm, rnd)
    
    res =  (cc_res + res )
    if verbose:
        print("CC: %s" % str(cc_res))
    return res
    
# Use one of the SC2 metrics without using the co-clustering matrix
def calculate3_other_no_cc(pred_ccm, pred_ad, truth_ccm, truth_ad, rnd=0.01, verbose=False, method="pseudoV"):
    methods_SC2 = {"pseudoV_no_cc": calculate2_pseudoV,
               "simpleKL_no_cc": calculate2_simpleKL,
               "sqrt_no_cc": calculate2_sqrt,
               "sym_pseudoV_no_cc": calculate2_sym_pseudoV,
               "pearson_no_cc": calculate2_pearson,
               "aupr_no_cc":calculate2_aupr,
                "mcc_no_cc": calculate2_mcc,
    }
    # Get the cousin matrices
    truth_cous = 1 - truth_ccm - truth_ad - truth_ad.T
    pred_cous = 1 - pred_ccm - pred_ad - pred_ad.T
    if verbose:
        if(np.amax(truth_cous) > 1 or np.amin(truth_cous) < 0):
            print("Cousin Truth is wrong...")
        if(np.amax(pred_cous) > 1 or np.amin(pred_cous) < 0):
            print("Cousin Predicted is wrong...")

    # Calculate the metric measure for each matrix
    func = methods_SC2[method]
    if method in ("pseudoV_no_cc",
               "simpleKL_no_cc",
               "sym_pseudoV_no_cc"):
        ad_res = func(pred_ad, truth_ad, rnd)
        ad_res_t = func(np.transpose(pred_ad), np.transpose(truth_ad), rnd)
        cous_res = func(pred_cous, truth_cous, rnd)
    else:
        ad_res = func(pred_ad, truth_ad)
        ad_res_t = func(np.transpose(pred_ad), np.transpose(truth_ad))
        cous_res = func(pred_cous, truth_cous)

    res =  0
    for r in (ad_res, ad_res_t, cous_res):
        if not math.isnan(r):
            res += r
    if verbose:
        print("%s for Matrices\nAD: %s, AD Transpose: %s, Cousin: %s\nReuslt: %s" %
              (method, str(ad_res),str(ad_res_t),str(cous_res), str(res)))
    return res

def parseVCFSimple(data):
    data = data.split('\n')
    data = [x for x in data if x != '']
    data = [x for x in data if x[0] != '#']
    if len(data) == 0:
        raise ValidationError("Input VCF contains no SSMs")
    return [[len(data)],[len(data)]]

def parseVCFScoring(data):
    data = data.split('\n')
    data = [x for x in data if x != '']
    data = [x for x in data if x[0] != '#']
    if len(data) == 0:
        raise ValidationError("Input VCF contains no SSMs")
    total_ssms = len(data)
    tp_ssms = len([x for x in data if x[-4:] == "True"])
    mask = [x[-4:] == "True" for x in data]
    mask = [i for i,x in enumerate(mask) if x]
    return [[total_ssms],[tp_ssms],mask]

def filterFPs(matrix, mask):
    if matrix.shape[0] == matrix.shape[1]:
        return matrix[mask,:][:,mask]
    else:
        return matrix[mask,:]

def verify(filename,role,func,*args):
    try:
        f = open(filename)
        pred_data = f.read(10**6)
        f.close()
        pred = func(pred_data,*args)
    except (IOError,TypeError) as e:
        print "Error opening " + role
        print e
        print filename,func
        return None
    except (ValidationError,ValueError) as e:
        print role + " does not validate"
        print e
        return None
    return pred


challengeMapping = {     '1A': {'val_funcs':[validate1A],'score_func':calculate1A,'vcf_func':None, 'filter_func':None},
                        '1B': {'val_funcs':[validate1B],'score_func':calculate1B,'vcf_func':None, 'filter_func':None},
                        '1C': {'val_funcs':[validate1C],'score_func':calculate1C,'vcf_func':parseVCFSimple, 'filter_func':None},
                        '2A': {'val_funcs':[validate2A],'score_func':calculate2,'vcf_func':parseVCFScoring, 'filter_func':filterFPs},
                        '2B': {'val_funcs':[validate2B],'score_func':calculate2,'vcf_func':parseVCFScoring, 'filter_func':filterFPs},
                        '3A': {'val_funcs':[validate2Afor3A,validate3A],'score_func':calculate3A,'vcf_func':parseVCFScoring, 'filter_func':filterFPs},
                        '3B': {'val_funcs':[validate2B,validate3B],'score_func':calculate3,'vcf_func':parseVCFScoring, 'filter_func':filterFPs},
                    }

def verifyChallenge(challenge,predfiles,vcf):
    if challengeMapping[challenge]['vcf_func']:
        nssms = verify(vcf,"input VCF", parseVCFSimple)
        if nssms == None:
            print "Could not read input VCF. Exiting"
            return "NA"
    else:
        nssms = [[],[]]

    if len(predfiles) != len(challengeMapping[challenge]['val_funcs']):
        print "Not enough input files for Challenge %s" % challenge
        return "Invalid"

    out = []
    for (predfile,valfunc) in zip(predfiles,challengeMapping[challenge]['val_funcs']):
        args = out + nssms[0]
        out.append(verify(predfile, "prediction file for Challenge %s" % (challenge),valfunc,*args))
        if out[-1] == None:
            return "Invalid"
    return "Valid"


def scoreChallenge(challenge,predfiles,truthfiles,vcf):
    if challengeMapping[challenge]['vcf_func']:
        nssms = verify(vcf,"input VCF", challengeMapping[challenge]['vcf_func'])
        if nssms == None:
            print "Could not read input VCF. Exiting"
            return "NA"
    else:
        nssms = [[],[]]
    if len(predfiles) != len(challengeMapping[challenge]['val_funcs']) or len(truthfiles) != len(challengeMapping[challenge]['val_funcs']):
        print "Not enough input files for Challenge %s" % challenge
        return "NA"

    tout = []
    pout = []
    for predfile,truthfile,valfunc in zip(predfiles,truthfiles,challengeMapping[challenge]['val_funcs']):
        targs = tout + nssms[1]
        tout.append(verify(truthfile, "truth file for Challenge %s" % (challenge),valfunc,*targs))
        pargs = pout + nssms[0]
        pout.append(verify(predfile, "prediction file for Challenge %s" % (challenge),valfunc,*pargs))
        if tout[-1] == None or pout[-1] == None:
            return "NA"
    if challengeMapping[challenge]['filter_func']:
        pout = [challengeMapping[challenge]['filter_func'](x,nssms[2]) for x in pout]
    return challengeMapping[challenge]['score_func'](*(pout + tout))

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("challenge")
	parser.add_argument("--predfiles",nargs="+")
	parser.add_argument("--truthfiles",nargs="*")
	parser.add_argument("--vcf")
	parser.add_argument("outputfile")
	parser.add_argument('-v', action='store_true', default=False)

	args = parser.parse_args()

	if args.pred_config is not None and args.truth_config is not None:
		with open(args.pred_config) as handle:
			pred_config = {}
			for line in handle:
				v = json.loads(line)
				if isinstance(v,dict):
					pred_config = dict(pred_config **v)
		with open(args.truth_config) as handle:
			truth_config = {}
			for line in handle:
				v = json.loads(line)
				if isinstance(v,dict):
					truth_config = dict(truth_config **v)
		out = {}
		for challenge in pred_config:
			if challenge in truth_config:
				predfile = pred_config[challenge]
				vcf = truth_config[challenge]['vcf']
				truthfiles = truth_config[challenge]['truth']
				if args.v:
					res = verifyChallenge(challenge,predfiles,vcf)
				else:
					res = scoreChallenge(challenge,predfiles,truthfiles,vcf)
				out[challenge] = res
		with open(args.outputfile, "w") as handle:
			jtxt = json.dumps( out )
			handle.write(jtxt)
	else:
		if args.v:
			res = verifyChallenge(args.challenge,args.predfiles,args.vcf)
		else:
			res = scoreChallenge(args.challenge,args.predfiles,args.truthfiles,args.vcf)

		with open(args.outputfile, "w") as handle:
			jtxt = json.dumps( { args.challenge : res } )
			handle.write(jtxt)
