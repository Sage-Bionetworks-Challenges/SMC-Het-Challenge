import numpy as np
import math
import itertools
from scipy.special import comb, gammaln

def calcSame(num_of_ones, num_of_mutations, rnd=1e-50):
    num_of_zeros = num_of_mutations - num_of_ones
    sum_of_row = num_of_zeros * rnd + num_of_ones

    return num_of_zeros * (rnd/sum_of_row) * np.log(rnd/sum_of_row) + num_of_ones * (1/sum_of_row) * np.log(1/sum_of_row)


def calcDifferent(num_of_ones_in_pred, num_of_ones_in_truth, num_of_mutations, tp, rnd=1e-50):
    fp = num_of_ones_in_pred - tp
    fn = num_of_ones_in_truth - tp
    tn = num_of_mutations - tp - fp - fn

    num_of_zeros_in_pred = num_of_mutations - num_of_ones_in_pred
    num_of_zeros_in_truth = num_of_mutations - num_of_ones_in_truth

    sum_of_pred_row = num_of_ones_in_pred + num_of_zeros_in_pred*rnd
    sum_of_truth_row = num_of_ones_in_truth + num_of_zeros_in_truth*rnd

    return (((1/sum_of_pred_row)*np.log(1/sum_of_truth_row) + (1/sum_of_truth_row)*np.log(1/sum_of_pred_row))*tp +
        ((rnd/sum_of_pred_row)*np.log(1/sum_of_truth_row) + (1/sum_of_truth_row)*np.log(rnd/sum_of_pred_row))*fn +
        ((1/sum_of_pred_row)*np.log(rnd/sum_of_truth_row) + (rnd/sum_of_truth_row)*np.log(1/sum_of_pred_row))*fp +
        ((rnd/sum_of_pred_row)*np.log(rnd/sum_of_truth_row) + (rnd/sum_of_truth_row)*np.log(rnd/sum_of_pred_row))*tn)

def ccm_permute_N_cluster(ad_true, rnd=1e-50):
    num_of_mutations = ad_true.shape[0]
    num_of_descendants_in_cluster_list = []
    

    for i in range(num_of_mutations):
        descendants = np.count_nonzero(ad_true[i])
        if not (descendants in num_of_descendants_in_cluster_list):
            num_of_descendants_in_cluster_list.append(descendants)
    
    num_of_clusters = len(num_of_descendants_in_cluster_list)
    num_of_mutations_in_cluster = np.zeros((num_of_clusters, 1))
    num_of_descendants_in_cluster = np.zeros((num_of_clusters, 1))

    for i in range(num_of_mutations):
        descendants = np.count_nonzero(ad_true[i])
        ind = num_of_descendants_in_cluster_list.index(descendants)
        num_of_descendants_in_cluster[ind, 0] = descendants
        num_of_mutations_in_cluster[ind, 0] = num_of_mutations_in_cluster[ind, 0]+1

    p = 0;
    for i in range(num_of_mutations):
        p += calcSame(i, num_of_mutations, rnd=rnd)
    #p /= num_of_mutations
    #print p
    # q*log(q)
    q = 0
    for i in range(num_of_clusters):
        q += num_of_mutations_in_cluster[i, 0] * calcSame(num_of_descendants_in_cluster[i, 0], num_of_mutations, rnd=rnd)

    #print q
    r = 0
    for i in range(num_of_clusters):
        for j in range(num_of_mutations):
            #print i, " ", j
            minTP = max(0, num_of_descendants_in_cluster[i, 0]-num_of_mutations+j+1)
            maxTP = min(num_of_descendants_in_cluster[i, 0], j)
            #print "minTP", minTP
            #print "maxTP", maxTP
            for TP in range(int(minTP), int(maxTP)+1):
                r += (num_of_mutations_in_cluster[i, 0]*
                    np.exp( gammaln(num_of_descendants_in_cluster[i, 0] + 1) -
                        gammaln(TP+1) - 
                        gammaln(num_of_descendants_in_cluster[i, 0] - TP + 1) + 
                        gammaln(num_of_mutations-1-num_of_descendants_in_cluster[i, 0]+1) -
                        gammaln(j-TP+1) - 
                        gammaln(num_of_mutations-1-num_of_descendants_in_cluster[i, 0]-j+TP+1) -
                        gammaln(num_of_mutations+1) +
                        gammaln(j+1) + 
                        gammaln(num_of_mutations-j) ) *
                    calcDifferent(j, num_of_descendants_in_cluster[i, 0], num_of_mutations, TP, rnd=rnd))
    print "DEBUG p, q, r: ", p, q, r
    return p+q-r


def om_permute_N_cluster(om, num_of_descendants_in_cluster, rnd=1e-50):
    num_of_mutations = 0
    for row in range(om.shape[0]):
        for column in range(om.shape[1]):
            num_of_mutations += om[row, column]
    #print num_of_mutations
    num_of_clusters = om.shape[0]
    #print num_of_clusters
    num_of_mutations_in_cluster = np.zeros((num_of_clusters, 1))
    # p*log(p)
    for row in range(om.shape[0]):
        num_of_mutations_in_cluster[row, 0] = np.sum(om[row])

    p = 0;
    for i in range(num_of_mutations):
        p += calcSame(i, num_of_mutations, rnd=rnd)
    #p /= num_of_mutations
    #print p
    # q*log(q)
    q = 0
    for i in range(num_of_clusters):
        q += num_of_mutations_in_cluster[i, 0] * calcSame(num_of_descendants_in_cluster[i, 0], num_of_mutations, rnd=rnd)

    #print q
    r = 0
    for i in range(num_of_clusters):
        for j in range(num_of_mutations):
            #print i, " ", j
            minTP = max(0, num_of_descendants_in_cluster[i, 0]-num_of_mutations+j+1)
            maxTP = min(num_of_descendants_in_cluster[i, 0], j)
            #print "minTP", minTP
            #print "maxTP", maxTP
            for TP in range(int(minTP), int(maxTP)+1):
                r += (num_of_mutations_in_cluster[i, 0]*
                    np.exp( gammaln(num_of_descendants_in_cluster[i, 0] + 1) -
                        gammaln(TP+1) - 
                        gammaln(num_of_descendants_in_cluster[i, 0] - TP + 1) + 
                        gammaln(num_of_mutations-1-num_of_descendants_in_cluster[i, 0]+1) -
                        gammaln(j-TP+1) - 
                        gammaln(num_of_mutations-1-num_of_descendants_in_cluster[i, 0]-j+TP+1) -
                        gammaln(num_of_mutations+1) +
                        gammaln(j+1) + 
                        gammaln(num_of_mutations-j) ) *
                    calcDifferent(j, num_of_descendants_in_cluster[i, 0], num_of_mutations, TP, rnd=rnd))
    print p, q, r
    return p+q-r


def calculate2_pseudoV(pred, truth, rnd=1e-50, full_matrix=True, sym=False):
    if full_matrix:
        pred_cp = pred
        truth_cp = truth
    else: # make matrix upper triangular
        pred_cp = np.triu(pred)
        truth_cp = np.triu(truth)

    # Avoid dividing by zero by rounding everything less than rnd up to rnd
    # Note: it is ok to do this after making the matrix upper triangular
    # since the bottom triangle of the matrix will not affect the score

    size = np.array(pred_cp.shape)[1]
    res = 0 # result to be returned

    # do one row at a time to reduce memory usage
    for x in xrange(size):
        # (1 - rnd) will cast the pred_cp/truth_cp matrices automatically if they are int8
        pred_row = (1 - rnd) * pred_cp[x, ] + rnd
        truth_row = (1 - rnd) * truth_cp[x, ] + rnd

        pred_row /= np.sum(pred_row)
        truth_row /= np.sum(truth_row)
        if sym:
            res += np.sum(truth_row * np.log(truth_row/pred_row)) + np.sum(pred_row * np.log(pred_row/truth_row))
        else:
            res += np.sum(truth_row * np.log(truth_row/pred_row))
    return res

def calculate2_sym_pseudoV(pred, truth, rnd=1e-50, full_matrix=True):
    return calculate2_pseudoV(pred, truth, rnd=rnd, full_matrix=full_matrix, sym=True)

def calculate2_sym_pseudoV_average(arr, ad_truth):
    perms = itertools.permutations(arr)
    n = len(arr)
    total_score = 0

    pred = np.triu(np.ones((n, n)), k=1)
    for perm in perms:
        ad = np.zeros((len(arr), len(arr)))
        for i in range(len(perm)):
            for j in range(len(perm)):
                if ad_truth[perm[i]-1, perm[j]-1] == 1:
                    ad[i, j] = 1
        #print ad
        total_score += calculate2_sym_pseudoV(pred, ad)



    return total_score/math.factorial(n)




if __name__ == "__main__":
    ad_truth = np.matrix([[0, 1, 1], [0, 0, 1], [0, 0, 0]])
    arr = [1, 2, 2, 3, 3, 3, 3, 3]
    print calculate2_sym_pseudoV_average(arr, ad_truth)

    om = np.matrix(([1, 0, 0], [1, 1, 0], [0, 1, 4]))
    t = np.matrix(([7], [5], [0]))
    print om_permute_N_cluster(om, t)


    a1 = np.copy(np.matrix(([0, 1, 1, 1, 1, 1, 1, 1],
                [0, 0, 0, 1, 1, 1, 1, 1],
                [0, 0, 0, 1, 1, 1, 1, 1],
                [0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0])))
    b = np.triu(np.ones((8, 8)), k=1)
    #print calculate2_sym_pseudoV(a1, b)
    print ccm_permute_N_cluster(a1)

    ad_truth = np.matrix([[0, 1, 1], [0, 0, 0], [0, 0, 0]])
    arr = [1, 1, 2, 2, 2, 3, 3]
    print calculate2_sym_pseudoV_average(arr, ad_truth)

    om = np.matrix(([2, 0, 0], [0, 2, 1], [0, 0, 2]))
    t = np.matrix(([5], [0], [0]))
    
    print om_permute_N_cluster(om, t)

