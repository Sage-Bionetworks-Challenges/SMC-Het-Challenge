import math
import numpy as np
import itertools
import json
import argparse
import StringIO
import scipy.stats
import sys
import sklearn.metrics as mt
import metric_behavior as mb
from functools import reduce


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

def calculate1A(pred,truth, err='abs'):
    if err is 'abs':
        return 1 - abs(truth - pred)
    elif err is 'sqr':
        return 1 - ((truth - pred) ** 2)
    else:
        raise KeyError('Invalid error penalty for scoring SC 1A. Choose one of "abs" or "sqr".')

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

def calculate1B(pred,truth, method='normalized'):
    if method is 'normalized':
        return (truth + 1 - min(truth+1,abs(pred-truth))) / float(truth+1)
    elif method is 'orig':
        return abs(truth - pred) / float(truth)
    else:
        raise KeyError('Invalid method for scoring SC 1B. Choose one of "orig" or "normalized".')

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
        raise ValidationError("Total number of reported mutations is %d. Should be %d" % (reported_nssms,nssms))
    return zip([int(x[1]) for x in data2], [float(x[2]) for x in data2])

def calculate1C(pred,truth, err='abs'):
    pred.sort(key = lambda x: x[1])
    truth.sort(key = lambda x: x[1])
    #itertools.chain(*x) flattens a list of lists
    predvs = np.array(list(itertools.chain(*[[x[1]]*x[0] for x in pred])))
    truthvs = np.array(list(itertools.chain(*[[x[1]]*x[0] for x in truth])))

    # calculate the score using the given error penalty
    if err is 'abs':
        se = abs(truthvs - predvs)
    elif err is 'sqr':
        se = ((truthvs - predvs) ** 2)
    else:
        raise KeyError('Invalid error penalty for scoring SC 1C. Choose one of "abs" or "sqr".')

    return sum(1-se)/float(len(truthvs))

def validate2A(data,nssms,return_ccm=True):
    data = data.split('\n')
    data = filter(None,data)
    if len(data) != nssms:
        raise ValidationError("Input file contains a different number of lines than the specification file. Input: %s lines Specification: %s lines" % (len(data), nssms))
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

def calculate2(pred,truth, full_matrix=True, method='default', pseudo_counts=None):
    '''Calculate the score for SubChallenge 2

    :param pred: predicted co-clustering matrix
    :param truth: true co-clustering matrix
    :param full_matrix: logical for whether to use the full matrix or just the upper triangular matrices when calculating the score
    :param method: scoring metric used, default is average of Pseudo V,
    :param pseudo_counts: logical for how many psuedo counts to add to the matrices
    :return: subchallenge 2 score for the predicted co-clustering matrix
    '''
    larger_is_worse_methods = ['pseudoV', 'sym_pseudoV'] # methods where a larger score is worse

    pc_pred = add_pseudo_counts(pred, num=pseudo_counts) # add pseudo counts to the matrices
    pc_truth = add_pseudo_counts(truth, num=pseudo_counts) # use np.copy so that original values are not altered
    ncluster = add_pseudo_counts(mb.get_ccm('NCluster', truth), num=pseudo_counts) # predicted CCM when every mutations is in its own cluster
    onecluster = add_pseudo_counts(mb.get_ccm('OneCluster', truth), num=pseudo_counts) # predicted CCM when all mutations are in the same cluster

    func_dict = {"orig" : calculate2_orig,
            "sqrt" : calculate2_sqrt,
            "pseudoV": calculate2_pseudoV,
            "sym_pseudoV": calculate2_sym_pseudoV,
            "spearman": calculate2_spearman,
            "pearson": calculate2_pearson,
            "aupr": calculate2_aupr,
            "mcc": calculate2_mcc
    }

    func = func_dict.get(method, None)
    if func is None:
        scores = []

        for m in ['pseudoV', 'pearson', 'mcc']:
            score = func_dict[m](pc_pred, pc_truth, full_matrix=full_matrix)

            # normalize the scores to be between (worst of OneCluster and NCluster scores) and (Truth score)
            ncluster_score = func_dict[m](ncluster, pc_truth, full_matrix=full_matrix)
            onecluster_score = func_dict[m](onecluster, pc_truth, full_matrix=full_matrix)
            if method in larger_is_worse_methods: # normalize scores where a larger score is worse
                # normalize the scores to be between 0 and 1 where 1 is the true matrix
                # and zero is the worse score of the NCluster matrix and the OneCluster matrix
                worst_score = max(ncluster_score, onecluster_score)
                score = 1 - (score / worst_score)
            else: # normalize scores where a smaller score is worse
                worst_score = min(ncluster_score, onecluster_score)
                score = (score - worst_score) / (1 - worst_score)
            scores.append(score)

        return np.mean(scores)
    else:
        score = func(pc_pred, pc_truth, full_matrix=full_matrix)
        ncluster_score = func(ncluster, pc_truth, full_matrix=full_matrix)
        onecluster_score = func(onecluster, pc_truth, full_matrix=full_matrix)
        if method in larger_is_worse_methods: # normalize the scores to be between 0 and 1 where 1 is the true matrix
            worst_score = max(ncluster_score, onecluster_score) # and zero is the worse score of the NCluster matrix
            score = 1 - (score / worst_score)                   # and the OneCluster matrix - similar to above
        else:
            worst_score = min(ncluster_score, onecluster_score)
            score = (score - worst_score) / (1 - worst_score)
        return score

def calculate2_orig(pred,truth, full_matrix=True):
    n = truth.shape[0]
    if full_matrix:
        pred_cp = np.copy(pred)
        truth_cp = np.copy(truth)
        count = (n**2 - n )

    else: # make matrix upper triangular
        inds = np.triu_indices(n,k=1)
        pred_cp = pred[inds]
        truth_cp = truth[inds]
        count = (n**2 - n )/2.0
    res = np.sum(np.abs(pred_cp - truth_cp))
    res = res / count
    return 1 - res


def calculate2_sqrt(pred,truth, full_matrix=True):
    n = truth.shape[0]
    if full_matrix:
        pred_cp = np.copy(pred)
        truth_cp = np.copy(truth)
        count = (n**2 - n)
    else: # make matrix upper triangular
        pred_cp = np.triu(pred)
        truth_cp = np.triu(truth)
        count = (n**2 - n )/2.0
    res = np.sum(np.abs(pred_cp - truth_cp))
    res = res / count
    return np.sqrt(1 - res)

def calculate2_simpleKL_norm(pred,truth,rnd=0.01):
    """Normalized version of the pseudo V measure where the return values are between 0 and 1
    with 0 being the worst score and 1 being the best

    :param pred:
    :param truth:
    :param rnd: small value to replace 0 entries in both matrices with. Used to avoid dividing by zero
    :return:
    """
    return 1 - calculate2_simpleKL(pred,truth,rnd=rnd) / 4000

# Out of date!!
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


def calculate2_pseudoV_norm(pred,truth,rnd=0.01, max_val=4000, full_matrix=True):
    """Normalized version of the pseudo V measure where the return values are between 0 and 1
    with 0 being the worst score and 1 being the best

    :param pred:
    :param truth:
    :param rnd: small value to replace 0 entries in both matrices with. Used to avoid dividing by zero
    :param max_val: maximum pseudoV value for this scenario - any prediction that has a pseudoV score >= max_val
        will be given a score of 0
    :return:
    """
    pv_val = calculate2_pseudoV(pred,truth,rnd=rnd, full_matrix=full_matrix)
    return max(1 -  pv_val/ max_val, 0)


def calculate2_pseudoV(pred,truth,rnd=0.01, full_matrix=True, sym=False):
    if full_matrix:
        pred_cp = np.copy(pred)
        truth_cp = np.copy(truth)
    else: # make matrix upper triangular
        pred_cp = np.triu(pred)
        truth_cp = np.triu(truth)

    # Avoid dividing by zero by rounding everything less than rnd up to rnd
    # Note: it is ok to do this after making the matrix upper triangular
    # since the bottom triangle of the matrix will not affect the score
    pred_cp = (1 - rnd) * pred_cp + rnd
    truth_cp = (1 - rnd) * truth_cp + rnd

    # normalize data
    pred_cp = pred_cp / np.sum(pred_cp,axis=1)[:,np.newaxis]
    truth_cp = truth_cp / np.sum(truth_cp,axis=1)[:,np.newaxis]

    if sym:
        return np.sum(truth_cp * np.log(truth_cp/pred_cp)) + np.sum(pred_cp * np.log(pred_cp/truth_cp))
    else:
        return np.sum(truth_cp * np.log(truth_cp/pred_cp))

def calculate2_sym_pseudoV_norm(pred,truth,rnd=0.01, max_val=8000, full_matrix=True):
    """Normalized version of the symmetric pseudo V measure where the return values are between 0 and 1
    with 0 being the worst score and 1 being the best

    :param pred:
    :param truth:
    :param rnd: small value to replace 0 entries in both matrices with. Used to avoid dividing by zero
    :param max_val: maximum pseudoV value for this scenario - any prediction that has a pseudoV score >= max_val
        will be given a score of 0
    :return:
    """
    spv_val = calculate2_sym_pseudoV(pred,truth,rnd=rnd, full_matrix=full_matrix)
    return max(1 - spv_val / max_val, 0)

def calculate2_sym_pseudoV(pred, truth, rnd=0.01, full_matrix=True):
    return calculate2_pseudoV(pred, truth, rnd=rnd, full_matrix=full_matrix, sym=True)

def calculate2_spearman(pred, truth, full_matrix=True):
    # use only the upper triangular matrix of the truth and
    # prediction matrices
    n = truth.shape[0]
    if full_matrix:
        pred_cp = pred.flatten()
        truth_cp = truth.flatten()
    else:
        inds = np.triu_indices(n,k=1)
        pred_cp = pred[inds]
        truth_cp = truth[inds]

    # implement spearman coefficient since scipy implementation
    # uses the covariance of the ranks, which could be zero
    # find the rank order of both sets of data
    predr = scipy.stats.rankdata(pred_cp)
    truthr = scipy.stats.rankdata(truth_cp)
    d = truthr - predr
    n = len(d)

    d = np.divide(d, np.sqrt(n)) # avoid overflow warnings
    row = 1 - (6 * sum(np.square(d))) / ((np.square(n) - 1))

    return row


def calculate2_pearson(pred, truth, full_matrix=True):
    n = truth.shape[0]
    if full_matrix:
        pred_cp = pred.flatten()
        truth_cp = truth.flatten()
    else:
        inds = np.triu_indices(n,k=1)
        pred_cp = pred[inds]
        truth_cp = truth[inds]
    return scipy.stats.pearsonr(pred_cp,truth_cp)[0]

def calculate2_aupr(pred,truth, full_matrix=True):
    n = truth.shape[0]
    if full_matrix:
        pred_cp = pred.flatten()
        truth_cp = truth.flatten()
    else:
        inds = np.triu_indices(n,k=1)
        pred_cp = pred[inds]
        truth_cp = truth[inds]
    precision, recall, thresholds = mt.precision_recall_curve(truth_cp,pred_cp)
    aucpr = mt.auc(recall, precision)
    return aucpr

# Matthews Correlation Coefficient
# don't just use upper triangular matrix because you get na's with the AD matrix
def calculate2_mcc(pred,truth, full_matrix=True):
    n = truth.shape[0]
    if full_matrix:
        pred_cp = pred
        truth_cp = truth
    else:
        inds = np.triu_indices(n,k=1)
        pred_cp = pred[inds]
        truth_cp = truth[inds]

    tp = float(sum(pred_cp[truth_cp==1] >= 0.5)) # Use 0.5 as the threshold for turning a probabalistic matrix into a binary matrix
    tn = float(sum(pred_cp[truth_cp==0] < 0.5))
    fp = float(sum(pred_cp[truth_cp==0] >= 0.5))
    fn = float(sum(pred_cp[truth_cp==1] < 0.5))

    # To avoid divide-by-zero cases
    denom_terms = [(tp+fp),(tp+fn),(tn+fp),(tn+fn)]
    for index, term in enumerate(denom_terms):
        if term == 0:
            denom_terms[index] = 1
    denom = np.sqrt(reduce(np.multiply, denom_terms, 1))

    if tp == 0 and fn == 0:
        num = (tn - fp)
    elif tn == 0 and fp == 0:
        num = (tp - fn)
    else:
        num = (tp*tn - fp*fn)

    return num / float(denom)

def validate2B(filename,nssms):
    try:
        if filename.endswith('.gz'):
            ccm = np.loadtxt(str(filename),ndmin=2)
        else:
            data = StringIO.StringIO(filename)
            ccm = np.loadtxt(data, ndmin=2)
    except ValueError as e:
        raise ValidationError("Entry in co-clustering matrix could not be cast as a float. Error message: %s" % e.message)

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

    # Form descendant of dict.  Each entry, keyed by cluster number, consists of a list of nodes that are decendents of the key.
    descendant_of = dict()
    for i in range(predK+1):
        descendant_of[i] = []
    for child,parent in data:
        descendant_of[parent] += [child] + descendant_of[child]
        # gps (grandparents) are the list of nodes that are ancestors of the immediate parent
        gps = [x for x in descendant_of.keys() if parent in descendant_of[x]]
        for gp in gps:
            descendant_of[gp] += [child] + descendant_of[child]

    # Check that root has all nodes as decendants (equivalent to checking if the tree is connected)
    if set(descendant_of[0]) != set(range(1,predK+1)):
        raise ValidationError("Root of phylogeny not ancestor of all clusters / Tree is not connected. " +
                              "Phelogeny matrix: %s, Descendant_of Dictionary %s" %
                              (data, descendant_of))

    # Form AD matrix
    n = len(cluster_assignments)
    ad = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            if cluster_assignments[j] in descendant_of[cluster_assignments[i]]:
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

def calculate3(pred_ccm, pred_ad, truth_ccm, truth_ad, method="sym_pseudoV_nc", weights=None, verbose=False, pseudo_counts=True, full_matrix=True, in_mat=2):
    """Calculate the score for subchallenge 3 using the given metric or a weighted average of the
    given metrics, if more than one are specified.

    :param pred_ccm: predicted co-clustering matrix
    :param pred_ad: predicted ancestor-descendant matrix
    :param truth_ccm: true co-clustering matrix
    :param truth_ad: trus ancestor-descendant matrix
    :param method: method to use when evaluating the submission or list of methods to use
    :param weights: weights to use in averaging the scores of different methods.
    :param verbose: boolean for whether to display information about the score calculations
    Only used if 'method' is a list - in this case must be a list of numbers of the same length as 'method'.
    :param full_matrix: boolean for whether to use the full CCM/AD matrix when calculating the score
    :param in_mat: number representing which matrices to use in calculating the SC3 scoring metric
        Options:
            1 - use all input matrics i.e. CCM, ADM, ADM^T and CM
            2 - use all except co-clustering matrix (CCM)
            3 - use all except ancestor descendant matrix (ADM)
            4 - use all except ADM^T
            5 - use all except cousin matrix (CM)
    :return: score for the given submission to subchallenge 3 using the given metric
    """
    larger_is_worse_methods = ['sym_pseudoV_nc', 'sym_pseudoV', 'pseudoV_nc', 'pseudoV', "simpleKL_nc", 'simpleKL'] # methods where a larger score is worse

    if pseudo_counts:
        if isinstance(pseudo_counts, int):
            pc_pred_ccm, pc_pred_ad = add_pseudo_counts(np.copy(pred_ccm), np.copy(pred_ad), num=pseudo_counts) # add pseudo counts to the matrices
            pc_truth_ccm, pc_truth_ad = add_pseudo_counts(np.copy(truth_ccm), np.copy(truth_ad), num=pseudo_counts) # use np.copy so that original values are not altered
            ncluster_ccm, ncluster_ad = add_pseudo_counts(mb.get_ccm('NClusterOneLineage', truth_ccm), mb.get_ad('NClusterOneLineage', truth_ad), num=pseudo_counts) # predicted matrices for each mutation being in their own cluster
            onecluster_ccm, onecluster_ad = add_pseudo_counts(mb.get_ccm('OneCluster', truth_ccm), mb.get_ad('OneCluster', truth_ad), num=pseudo_counts) # predicted matrices for all mutations being in the same cluster
        else:
            pc_pred_ccm, pc_pred_ad = add_pseudo_counts(np.copy(pred_ccm), np.copy(pred_ad))
            pc_truth_ccm, pc_truth_ad = add_pseudo_counts(np.copy(truth_ccm), np.copy(truth_ad))
            ncluster_ccm, ncluster_ad = add_pseudo_counts(mb.get_ccm('NClusterOneLineage', truth_ccm), mb.get_ad('NClusterOneLineage', truth_ad))
            onecluster_ccm, onecluster_ad = add_pseudo_counts(mb.get_ccm('OneCluster', truth_ccm), mb.get_ad('OneCluster', truth_ad))
    else:
        pc_pred_ccm, pc_pred_ad, pc_truth_ccm, pc_truth_ad = pred_ccm, pred_ad, truth_ccm, truth_ad
        ncluster_ccm, ncluster_ad = mb.get_ccm('NClusterOneLineage', truth_ccm), mb.get_ad('NClusterOneLineage', truth_ad)
        onecluster_ccm, onecluster_ad = mb.get_ccm('OneCluster', truth_ccm), mb.get_ad('OneCluster', truth_ad)

    if isinstance(method, list):
        res = [calculate3_onemetric(pc_pred_ccm, pc_pred_ad, pc_truth_ccm, pc_truth_ad,
                                    method=m, verbose=verbose, in_mat=in_mat) for m in method] # calculate the score for each method

        # normalize the scores to be between (worst of NCluster score and OneCluster score) and (Truth score)
        ncluster_score = [calculate3_onemetric(ncluster_ccm, ncluster_ad, pc_truth_ccm, pc_truth_ad,
                                               method=m, verbose=verbose, full_matrix=full_matrix, in_mat=in_mat) for m in method]
        onecluster_score = [calculate3_onemetric(onecluster_ccm, onecluster_ad, pc_truth_ccm, pc_truth_ad,
                                                 method=m, verbose=verbose, full_matrix=full_matrix, in_mat=in_mat) for m in method]
        for i in range(len(method)):
            if method[i] in larger_is_worse_methods: # normalization for methods where a larger score is worse
                worst_score = max(ncluster_score[i], onecluster_score[i]) # worst of NCluster and OneCluster scores
                res[i] = 1 - (res[i] / worst_score) # normalize the score
            else: # normalization for methods where a smaller score is worse
                worst_score = min(ncluster_score[i], onecluster_score[i])
                res[i] = (res[i] - worst_score) / (1 - worst_score)


        if weights is None: # if weights are not specified or if they cannot be normalized then default to equal weights
            weights = [1] * len(method)
        elif sum(weights) == 0:
            Warning('Weights sum to zero so they are invalid, defaulting to equal weights')
            weights = [1] * len(method)

        weights = np.array(weights) / float(sum(weights)) # normalize the weights
        score = sum(np.multiply(res, weights))
    else:
        score =  calculate3_onemetric(pc_pred_ccm, pc_pred_ad, pc_truth_ccm, pc_truth_ad,
                                      method=method, verbose=verbose, full_matrix=full_matrix, in_mat=in_mat)

        # normalize the score to be between (worst of NCluster score and OneCluster score) and (Truth score) - similar to above
        ncluster_score = calculate3_onemetric(ncluster_ccm, ncluster_ad, pc_truth_ccm, pc_truth_ad,
                                              method=method, verbose=verbose, full_matrix=full_matrix, in_mat=in_mat)
        onecluster_score = calculate3_onemetric(onecluster_ccm, onecluster_ad, pc_truth_ccm, pc_truth_ad,
                                                method=method, verbose=verbose, full_matrix=full_matrix, in_mat=in_mat)
        if method in larger_is_worse_methods:
            worst_score = max(ncluster_score, onecluster_score)
            score = 1 - (score / worst_score)
        else:
            worst_score = min(ncluster_score, onecluster_score)
            score = (score - worst_score) / (1 - worst_score)
    return score

# dictionary of method names and their corresponding metric functions
method_funcs = {"pseudoV": calculate2_pseudoV,
               "simpleKL": calculate2_simpleKL,
               "sqrt": calculate2_sqrt,
               "sym_pseudoV": calculate2_sym_pseudoV,
               "pearson": calculate2_pearson,
                "spearman": calculate2_spearman,
               "aupr":calculate2_aupr,
                "mcc": calculate2_pearson,
                "orig": calculate2_orig
    }

def calculate3_onemetric(pred_ccm, pred_ad, truth_ccm, truth_ad, rnd=0.01, method="orig_nc", verbose=False, full_matrix=True, in_mat=2):
    """Calculate the score for subchallenge 3 using the given metric

    :param pred_ccm: predicted co-clustering matrix
    :param pred_ad: predicted ancestor-descendant matrix
    :param truth_ccm: true co-clustering matrix
    :param truth_ad: trus ancestor-descendant matrix
    :param method: method to use when evaluating the submission
    :param verbose: boolean for whether to display information about the score calculations
    :param full_matrix: boolean for whether to use the full CCM/AD matrix when calculating the score
    :param in_mat: number representing which matrices to use in calculating the SC3 scoring metric
        Options:
            1 - use all input matrics i.e. CCM, ADM, ADM^T and CM
            2 - use all except co-clustering matrix (CCM)
            3 - use all except ancestor descendant matrix (ADM)
            4 - use all except ADM^T
            5 - use all except cousin matrix (CM)
    :return: score for the given submission to subchallenge 3 using the given metric
    """
    # Get the cousin matrices
    truth_cous = 1 - truth_ccm - truth_ad - truth_ad.T
    pred_cous = 1 - pred_ccm - pred_ad - pred_ad.T
    if verbose:
        if(np.amax(truth_cous) > 1 or np.amin(truth_cous) < 0):
            Warning("Cousin Truth is wrong. Maximum matrix entry is greater than 1 or minimum matrix entry is less than 0")
        if(np.amax(pred_cous) > 1 or np.amin(pred_cous) < 0):
            Warning("Cousin Predicted is wrong. Maximum matrix entry is greater than 1 or minimum matrix entry is less than 0")

    # Calculate the metric measure for each specified matrix
    func = method_funcs[method]
    results = []
    ccm_res, ad_res, ad_res_t, cous_res = [float('nan')] * 4
    if method in ("pseudoV",
               "simpleKL",
               "sym_pseudoV"):
        if in_mat != 2:
            ccm_res = func(pred_ccm, truth_ccm, rnd, full_matrix=full_matrix)
            results.append(ccm_res)
        if in_mat != 3:
            ad_res = func(pred_ad, truth_ad, rnd, full_matrix=full_matrix)
            results.append(ad_res)
        if in_mat != 4:
            ad_res_t = func(np.transpose(pred_ad), np.transpose(truth_ad), rnd, full_matrix=full_matrix)
            results.append(ad_res_t)
        if in_mat != 5:
            cous_res = func(pred_cous, truth_cous, rnd, full_matrix=full_matrix)
            results.append(cous_res)
    else:
        if in_mat != 2:
            ccm_res = func(pred_ccm, truth_ccm, full_matrix=full_matrix)
            results.append(ccm_res)
        if in_mat != 3:
            ad_res = func(pred_ad, truth_ad, full_matrix=full_matrix)
            results.append(ad_res)
        if in_mat != 4 or method in ('mcc',
                                     'pearson',
                                     'spearman'):
            ad_res_t = func(np.transpose(pred_ad), np.transpose(truth_ad), full_matrix=full_matrix)
            results.append(ad_res_t)
        if in_mat != 5:
            cous_res = func(pred_cous, truth_cous, full_matrix=full_matrix)
            results.append(cous_res)

    res =  0
    n = 0
    for r in results: # TODO: fix the NA's
        if not math.isnan(r):
            n += 1
            res += r
    if n > 0:
        res = res / float(n)

    if verbose:
        print("%s for Matrices\nCC: %s, AD: %s, AD Transpose: %s, Cousin: %s\nResult: %s" %
              (method, str(ccm_res), str(ad_res),str(ad_res_t),str(cous_res), str(res)))
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
    mask = [x[-4:] == "True" for x in data]
    mask = [i for i,x in enumerate(mask) if x]
    tp_ssms = len(mask)
    return [[total_ssms],[tp_ssms],mask]

def filterFPs(matrix, mask):
    if matrix.shape[0] == matrix.shape[1]:
        return matrix[np.ix_(mask, mask)]
    else:
        return matrix[mask,:]

def add_pseudo_counts(ccm,ad=None,num=None):
    """Add a small number of fake mutations or 'pseudo counts' to the co-clustering and ancestor-descendant matrices for
    subchallenges 2 and 3, each in their own, new cluster. This ensures that there are not cases where
    either of these matrices has a variance of zero. These fake mutations must be added to both the predicted
    and ruth matrices.

    :param ccm: co-clustering matrix
    :param ad: ancestor-descendant matrix (optional, to be compatible with subchallenge 2)
    :param num: number of pseudo counts to add
    :return: modified ccm and ad matrices
    """
    size = np.array(ccm.shape)[1]

    if num is None:
        num = np.sqrt(size)
    elif num == 0:
        return ccm, ad
    print(sys.getsizeof(ccm))

    ccm=np.vstack((ccm, np.zeros(shape=[num, size])))
    ccm=np.hstack((ccm, np.zeros(shape=[size+num, num])))

    print('Success')

    new_ccm = np.identity(size + num)
    new_ccm[:size, :size] = np.copy(ccm)
    ccm = new_ccm

    if ad is not None:
        new_ad = np.zeros([size + num]*2)
        new_ad[:size, :size] = np.copy(ad)
        new_ad[(size+num/2):(size+3*num/4),:(size)] = 1 # one quarter of the pseudo counts are ancestors of (almost) every other cluster
        new_ad[:(size),(size+3*num/4):(size+num)] = 1 # one quarter of the pseudo counts are descendants of (almost) every other cluster
        ad = new_ad                                         # half of the pseudo counts are cousins to all other clusters
        return ccm, ad

    return ccm

def verify(filename,role,func,*args):
    global err_msgs
    try:
        if filename.endswith('.gz'): #pass compressed files directly to 2B or 3B validate functions
            pred = func(filename,*args)
        else:
            f = open(filename)
            pred_data = f.read()
            f.close()
            pred = func(pred_data,*args)
    except (IOError,TypeError) as e:
        err_msgs.append("Error opening %s, from function %s using file %s in : %s" %  (role, func, filename, e.strerror))
        return None
    except (ValidationError,ValueError) as e:
        err_msgs.append("%s does not validate: %s" % (role, e.value))
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
    global err_msgs
    if challengeMapping[challenge]['vcf_func']:
        nssms = verify(vcf,"input VCF", parseVCFSimple)
        if nssms == None:
            err_msgs.append("Could not read input VCF. Exiting")
            return "NA"
    else:
        nssms = [[],[]]

    if len(predfiles) != len(challengeMapping[challenge]['val_funcs']):
        err_msgs.append("Not enough input files for Challenge %s" % challenge)
        return "Invalid"

    out = []
    for (predfile,valfunc) in zip(predfiles,challengeMapping[challenge]['val_funcs']):
        args = out + nssms[0]
        out.append(verify(predfile, "prediction file for Challenge %s" % (challenge),valfunc,*args))
        if out[-1] == None:
            return "Invalid"
    return "Valid"


def scoreChallenge(challenge,predfiles,truthfiles,vcf):
    global err_msgs
    if challengeMapping[challenge]['vcf_func']:
        nssms = verify(vcf,"input VCF", challengeMapping[challenge]['vcf_func'])
        if nssms == None:
            err_msgs.append("Could not read input VCF. Exiting")
            return "NA"

    else:
        nssms = [[],[]]
    if len(predfiles) != len(challengeMapping[challenge]['val_funcs']) or len(truthfiles) != len(challengeMapping[challenge]['val_funcs']):
        err_msgs.append("Not enough input files for Challenge %s" % challenge)
        return "NA"

    tout = []
    pout = []
    for predfile,truthfile,valfunc in zip(predfiles,truthfiles,challengeMapping[challenge]['val_funcs']):
        if truthfile.endswith('.gz') and challenge not in ['2B', '3B']:
            err_msgs.append('Incorrect format, must input a text file for challenge %s' % challenge)
            return "NA"
        targs = tout + nssms[1]
        tout.append(verify(truthfile, "truth file for Challenge %s" % (challenge),valfunc,*targs))

        if predfile.endswith('.gz') and challenge not in ['2B', '3B']:
            err_msgs.append('Incorrect format, must input a text file for challenge %s' % challenge)
            return "NA"
        pargs = pout + nssms[0]
        pout.append(verify(predfile, "prediction file for Challenge %s" % (challenge),valfunc,*pargs))
        if tout[-1] == None or pout[-1] == None:
            return "NA"
    if challengeMapping[challenge]['filter_func']:
        print('Filtering Challenge %s' % challenge)
        pout = [challengeMapping[challenge]['filter_func'](x,nssms[2]) for x in pout]
    return challengeMapping[challenge]['score_func'](*(pout + tout))


if __name__ == '__main__':
    global err_msgs
    err_msgs = []

    parser = argparse.ArgumentParser()
    parser.add_argument("--pred-config", default=None)
    parser.add_argument("--truth-config", default=None)
    parser.add_argument("-c", "--challenge", default=None)
    parser.add_argument("--predfiles",nargs="+")
    parser.add_argument("--truthfiles",nargs="*")
    parser.add_argument("--vcf")
    parser.add_argument("-o", "--outputfile")
    parser.add_argument('-v', action='store_true', default=False)

    args = parser.parse_args()

    if args.pred_config is not None and args.truth_config is not None:
        with open(args.pred_config) as handle:
            pred_config = {}
            for line in handle:
                try:
                    v = json.loads(line)
                    if isinstance(v,dict):
                        pred_config = dict(pred_config, **v)
                except ValueError as e:
                    pass
        with open(args.truth_config) as handle:
            truth_config = {}
            for line in handle:
                try:
                    v = json.loads(line)
                    if isinstance(v,dict):
                        truth_config = dict(truth_config, **v)
                except ValueError as e:
                    pass
        out = {}
        print "pred", pred_config
        print "truth", truth_config
        for challenge in pred_config:
            if challenge in truth_config:
                predfile = pred_config[challenge]
                vcf = truth_config[challenge]['vcf']
                truthfiles = truth_config[challenge]['truth']
                if args.v:
                    res = verifyChallenge(challenge,predfile,vcf)
                else:
                    res = scoreChallenge(challenge,predfile,truthfiles,vcf)
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

    if len(err_msgs) > 0:
        for msg in err_msgs:
            print msg
        raise ValidationError("Errors encountered. If running in Galaxy see stdout for more info. The results of any successful evaluations are in the Job data.")
