import numpy as np
from permutations import*

import gc

def om_validate2A (pred_data, truth_data, nssms_x, nssms_y, filter_mut=None, mask=None, subchallenge="2A"):
    '''
    Creates overlapping matrix for SubChallenge 2 and 3
    :param pred_data: inputed data from prediction file
    :param truth_data: inputed data from truth file
    :param nssms_x: number of mutations prediction file (specified by vcf)
    :param filter_mut: list of mutations to filter in prediction file
    :param mask: mask applied
    :subchallenge: subchallenge scored
    :return: overlapping matrix and (for subchallenge 3) a list which specifies the cluster of each mutation
    '''
    pred_data = pred_data.split('\n')
    pred_data = filter(None, pred_data)
    pred_data = [x for i, x in enumerate(pred_data) if i in mask] if mask else pred_data
     
    if len(pred_data) != nssms_x:
        raise ValidationError("Prediction file contains a different number of lines than the specification file. Input: %s lines. Specification: %s lines" % (len(pred_data), nssms_x))
    pred_cluster_entries = set()

    for i in xrange(len(pred_data)):
        try:
            pred_data[i] = int(pred_data[i])
            pred_cluster_entries.add(pred_data[i])
        except ValueError:
            raise ValidationError("Cluster ID in line %d (ssm %s) can not be cast to an int", (i+1, pred_data[i][0]))

    num_pred_clusters = max(pred_cluster_entries)

    truth_data = truth_data.split('\n')
    truth_data = filter(None, truth_data)
    truth_data = [x for i, x in enumerate(truth_data) if i in mask] if mask else truth_data

    if len(truth_data) != nssms_y:
        raise ValidationError("Truth file contains a different number of lines than the specification file. Input: %s lines. Specification: %s lines" % (len(truth_data), nssms_y))

    truth_cluster_entries = set()
    for i in xrange(len(truth_data)):
        try:
            truth_data[i] = int(truth_data[i])
            truth_cluster_entries.add(truth_data[i])
        except ValueError:
            raise ValidationError("Cluster ID in line %d (ssm %s) can not be cast to an int", (i+1, truth_data[i][0]))

    num_truth_clusters = max(truth_cluster_entries)

    om = np.zeros((num_truth_clusters, num_pred_clusters), dtype=int)

    # print len(filter_mut)
    new_pred_data = []
    if filter_mut is not None:
        for i in range(len(pred_data)):     
            if i in filter_mut:
                new_pred_data.append(pred_data[i])
    else:
        new_pred_data = pred_data

    for i in range(len(new_pred_data)):
        # print(new_pred_data[i]), 
        om[truth_data[i]-1, new_pred_data[i]-1] += 1

    if subchallenge is "3A":
        return om, truth_data

    return om

def om_calculate2A(om, full_matrix=True, method='default', add_pseudo=True, pseudo_counts=None, rnd=1e-50):
    '''
    Calculate the score for SubChallenge 2
    :param om: overlapping matrix
    :param full_matrix: logical for whether to use the full matrix or just the upper triangular matrices when calculating the score
    :param method: scoring metric used, default is average of Pseudo V,
    :param pseudo_counts: logical for how many psuedo counts to add to the matrices
    :return: subchallenge 2 score for the predicted co-clustering matrix
    '''

    larger_is_worse_methods = ['pseudoV', 'sym_pseudoV', 'js_divergence'] # methods where a larger score is worse
    import gc
    func_dict = {
        "orig"           : om_calculate2_orig,
        "sqrt"           : om_calculate2_sqrt,
        "spearman"       : om_calculate2_spearman,
        "aupr"           : om_calculate2_aupr,
        "pseudoV"        : om_calculate2_pseudoV,
        "sym_pseudoV"    : om_calculate2_sym_pseudoV,
        "js_divergence"  : om_calculate2_js_divergence,
        "mcc"            : om_calculate2_mcc
    }
    func = func_dict.get(method, None)
    tp, fp, tn, fn = calculate_overlap_matrix(om)
    if add_pseudo:
        tp, fp, tn, fn = add_pseudo_counts_om_eff(tp, fp, tn, fn)
    if func is None:
        scores = []
        worst_scores = []

        # pearson and mcc give the same score for 2A
        functions = ['js_divergence','mcc','mcc']

        for m in functions:
            gc.collect()
            if m is 'pseudoV' or m is 'sym_pseudoV' or m is 'js_divergence':
                scores.append(func_dict[m](om, full_matrix=full_matrix, modify=add_pseudo, pseudo_counts=pseudo_counts, rnd=rnd))
            else:
                scores.append(func_dict[m](tp, fp, tn, fn, full_matrix=full_matrix, rnd=rnd))

            # normalize the scores to be between (worst of OneCluster and NCluster scores) and (Truth score)
        for m in functions:
            gc.collect()
            worst_scores.append(get_worst_score_om(om, func_dict[m], larger_is_worse=(m in larger_is_worse_methods), rnd=rnd))
#        print "DEBUG raw scores", scores
#        print "DEBUG worst scores", worst_scores
        for i, m in enumerate(functions):
            if m in larger_is_worse_methods:
                scores[i] = set_to_zero(1 - (scores[i] / worst_scores[i]))
            else:
                scores[i] = set_to_zero((scores[i] - worst_scores[i]) / (1 - worst_scores[i]))
        return [scores, np.mean(scores)]

    else:
        # if it is pseudo_count, immediately modify true
        if func is func_dict['pseudoV'] or func is func_dict['sym_pseudoV'] or func is func_dict['js_divergence']:
            if add_pseudo:
                score = func(om, full_matrix=full_matrix, modify=True, pseudo_counts=pseudo_counts, rnd=rnd)
            else:
                score = func(om, full_matrix=full_matrix, modify=False, pseudo_counts=pseudo_counts, rnd=rnd)
        else:
            score = func(tp, fp, tn, fn, full_matrix=full_matrix, rnd=rnd)

        if method in larger_is_worse_methods: # normalize the scores to be between 0 and 1 where 1 is the true matrix
            worst_score = get_worst_score_om(om, func, larger_is_worse=True, rnd=rnd) # and zero is the worse score of the NCluster matrix
            score = set_to_zero(1 - (score / worst_score))                   # and the OneCluster matrix - similar to above
        else:
            worst_score = get_worst_score_om(om, func, larger_is_worse=False, rnd=rnd)
            score = set_to_zero((score - worst_score) / (1 - worst_score))
        return score

# nssms should be the length of the pred file as well as the length of the truth file
def calculate_overlap_matrix(om):
    '''
    Calculates number of true postives, false postives, true negatives and false negatives from om
    :param om: overlapping matrix
    :return: number of true postives, false postives, true negatives and false negatives
    '''
    # tp is the number of true postives, p is the number of ones in the truth matrix and t is the total number of entries (size of the file)
    tp = 0
    p = 0
    t = 0

    # calculate number of true positives and total positives
    for row in range(om.shape[0]):
        p += np.sum(om[row])**2
        for column in range(om.shape[1]):
            tp += om[row, column]**2
            t += om[row, column]

    # fn is the number of false negatives       
    fn = p - tp

    # n is the number of zeros in the truth matrix
    n = t**2 - p
    
    # pp is the number of predicted postives
    pp = 0
    for column in range(om.shape[1]):
        pp += np.sum(om[:,column])**2

    fp = pp - tp

    tn = n - fp

    return tp, fp, tn, fn

#### METRICS DEVELOPED FOR OPTIMIZATION ####################################################################

# this function is the same as calculate2_orig but is customized for overlapping matrix
def om_calculate2_orig(tp, fp, tn, fn, full_matrix=True):
    if full_matrix:
        tp -= int(np.around(np.sqrt(tp+fn+fp+tn)))
    else:
        tp = int((tp - np.around(np.sqrt(tp+fn+fp+tn)))/2.0)
        fn /= 2 
        fp /= 2
        tn /= 2
    res = fp + fn
    count = tp + fn + fp + tn

    return 1-float(res)/count

# this function is the same as calculate2_sqrt but is customized for overlapping matrix
def om_calculate2_sqrt(tp, fp, tn, fn, full_matrix=True):
    if full_matrix:
        tp -= int(np.around(np.sqrt(tp+fn+fp+tn)))
    else:
        tp = int((tp - np.around(np.sqrt(tp+fn+fp+tn)))/2)
        fn /= 2 
        fp /= 2
        tn /= 2

    res = fp + fn
    count = tp + fn + fp + tn

    return np.sqrt(1-float(res)/count)

def om_calculate2_js_divergence_helper(om, t, row, column, modify=False, pseudo_counts=None, rnd=1e-50, P=None, T=None, challenge='2A'):

    """Helper function for JS-divergence score
    :param om: shared overlap matrix
    :param t: total number of entries
    :param row: row to use
    :param column: column to use
    :param modify: used to determine if pseudo_counts should be added
    :param pseudo_counts: number of pseudo_counts that will be added to overlapping matrix
    :return: one-sided score for subchallenge 2
    """
    # get smallest possible float -- for 0 values
    small = np.nextafter(0.,1.)

    if challenge is '2A':
        total_truth = np.sum(om[row,:]) * 1.0
        total_pred = np.sum(om[:,column]) * 1.0
    elif challenge is '3A':
        total_truth = T[row, 0] * 1.0
        total_pred = P[column, 0] * 1.0
        
        if total_truth == 0 or total_pred == 0:
            return 0

    tp = om[row,column] * 1.0
    fn = total_truth - tp
    fp = total_pred - tp
    tn = t + tp - total_truth - total_pred

    if modify:
        tn += pseudo_counts

    if tp < 0 or fn < 0 or fp < 0 or tn < 0:
        raise ValidationError("True positive or false negative should not be negative values")

    total_truth_pseudo = tp + fn + (fp + tn)*rnd
    total_pred_pseudo = tp + fp + (fn + tn)*rnd

    mid_tp_cases = 0.5 * (1/total_truth_pseudo + 1/total_pred_pseudo)
    mid_fn_cases = 0.5 * (1/total_truth_pseudo + rnd/total_pred_pseudo)
    mid_fp_cases = 0.5 * (rnd/total_truth_pseudo + 1/total_pred_pseudo)
    mid_tn_cases = 0.5 * (rnd/total_truth_pseudo + rnd/total_pred_pseudo)
    # mid_fp_cases is not 0, but is multiplied by a 0 term since truth file doesn't have it; so can ignore
    # mid_tn_cases is 0

    res = (
        total_truth_pseudo * np.log(1/total_truth_pseudo) -
        tp * np.log(mid_tp_cases) -
        fn * np.log(mid_fn_cases) -
        fp * (rnd+small) * np.log(mid_fp_cases+small) -
        tn * (rnd+small) * np.log(mid_tn_cases+small)
        )

    res /= total_truth
    return res
 
def om_calculate2_js_divergence(om, rnd=1e-50, full_matrix=True, sym=True, modify=False, pseudo_counts=None):

    """Calculates the Jensen-Shannon divergence score for subchallenge 2
    :param om: shared overlap matrix
    :param rnd: small value to replace 0 entries in both matrices with. Used to avoid dividing by zero
    :param full_matrix: used to determine if full matrix should be used
    :param sym: score will be the same if prediction and truth file were swapped
    :param modify: used to determine if pseudo_counts should be added
    :param pseudo_counts: number of pseudo_counts that will be added to overlapping matrix
    :return: score for subchallenge 2
    """
    res = 0
    t = 0

    for row in range(om.shape[0]):
        t += np.sum(om[row])

    if modify:
        if pseudo_counts is None:
            pseudo_counts = int(np.floor(np.sqrt(t)))

    for row in range(om.shape[0]):
        for column in range(om.shape[1]):
            for count in range(om[row, column]):
                sym1 = om_calculate2_js_divergence_helper(om, t, row, column, modify, pseudo_counts, rnd=rnd, challenge='2A')
                res += sym1

                if sym:
                    sym2 = om_calculate2_js_divergence_helper(np.transpose(om), t, column, row, modify, pseudo_counts, rnd=rnd, challenge='2A')
                    res += sym2

    return res

def om_calculate2_pseudoV(om, rnd=1e-50, full_matrix=True, sym=False, modify=False, pseudo_counts=None):

    """Calculates the pseudoV score for subchallenge 2
    :param srm: shared relative matrix (this could be shared ancestor matrix, shared cousin matrix, or shared descendent matrix)
    :param rnd: small value to replace 0 entries in both matrices with. Used to avoid dividing by zero
    :param full_matrix: used to determine if full matrix should be used
    :param sym: score will be the same if prediction and truth file were swapped
    :param modify: used to determine if pseudo_counts should be added
    :param pseudo_counts: number of pseudo_counts that will be added to overlapping matrix
    :return: score for subchallenge 2
    """
    res = 0
    t = 0

    for row in range(om.shape[0]):
        t += np.sum(om[row])

    if modify:
        if pseudo_counts is None:
            pseudo_counts = int(np.floor(np.sqrt(t)))

    for row in range(om.shape[0]):
        for column in range(om.shape[1]):
            for count in range(om[row, column]):
                tp = om[row, column]
                fn = np.sum(om[row]) - tp
                fp = np.sum(om[:,column]) - tp
                tn = t + tp - np.sum(om[row]) - np.sum(om[:,column])

                if modify:
                    tn += pseudo_counts

                if tp < 0 or fn < 0 or fp < 0 or tn < 0:
                    raise ValidationError("True positive, false negative, false postive and true negative should not be negative values")

                sum_of_truth_row = tp + fn + (fp + tn)*rnd
                sum_of_pred_row = tp + fp + (fn + tn)*rnd

                # get smallest possible float -- for 0 values
                small = np.nextafter(0.,1.)
                if tp != 0 or fp != 0 or tn != 0 or fn != 0:
                    sym2 = (
                        np.log(sum_of_pred_row) * (tp + rnd*fp + fn + rnd*tn) -
                        np.log(sum_of_truth_row) * (tp + rnd*fp + fn + rnd*tn) +
                        np.log(rnd+small) * ((rnd+small)*fp - fn)
                        )

                    sym2 /= sum_of_truth_row
                    res += sym2

                if sym:
                    if tp != 0 or fp != 0 or tn != 0 or fn != 0:
                        sym1 = (
                            np.log(sum_of_truth_row) * (tp + rnd*fn + fp + rnd*tn) -
                            np.log(sum_of_pred_row) * (tp + rnd*fn + fp + rnd*tn) +
                            np.log(rnd+small) * ((rnd+small)*fn - fp)
                            )

                        sym1 /= sum_of_pred_row
                        res += sym1

        # no need to consider the other pseudo_count rows, since sym1 and sym2 evaluate to zero for these rows

    return res

def calculate3_js_divergence(srm, om, P, T, rnd=1e-50, sym=True):
    """Calculates the JS divergence score for subchallenge 3
    :param srm: shared relative matrix (this could be shared ancestor matrix, shared cousin matrix, or shared descendent matrix)
    :param om: overlapping matrix
    :param P: matrix that specifies number of ancestors/descendents each cluster has in the prediction file
    :param T: matrix that specifies number of ancestors/descendents each cluster has in the truth file
    :param rnd: small value to replace 0 entries in both matrices with. Used to avoid dividing by zero
    :param sym: score will be the same if prediction and truth file were swapped
    :return: score for subchallenge 3
    """
    res = 0
    t = 0
    for row in range(om.shape[0]):
        t += np.sum(om[row])

    for row in range(om.shape[0]):
        for column in range(om.shape[1]):
            if om[row, column] != 0:
                for count in range(om[row, column]):
                    sym1 = om_calculate2_js_divergence_helper(srm, t, row, column, modify=False, rnd=rnd, P=P, T=T, challenge='3A')
                    res += sym1

                    if sym:
                        sym2 = om_calculate2_js_divergence_helper(np.transpose(srm), t, column, row, modify=False, rnd=rnd, P=T, T=P, challenge='3A')
                        res += sym2

    return res


def calculate3_pseudoV(srm, om, P, T, rnd=1e-50, sym=False):
    """Calculates the pseudoV score for subchallenge 3
    :param srm: shared relative matrix (this could be shared ancestor matrix, shared cousin matrix, or shared descendent matrix)
    :param om: overlapping matrix
    :param P: matrix that specifies number of ancestors/descendents each cluster has in the prediction file
    :param T: matrix that specifies number of ancestors/descendents each cluster has in the truth file
    :param rnd: small value to replace 0 entries in both matrices with. Used to avoid dividing by zero
    :param sym: score will be the same if prediction and truth file were swapped
    :return: score for subchallenge 3
    """
    res = 0
    t = 0
    for row in range(om.shape[0]):
        t += np.sum(om[row])

    for row in range(om.shape[0]):
        for column in range(om.shape[1]):
            if om[row, column] != 0:
                for i in range(om[row, column]):
                    tp = srm[row, column]
                    fn = T[row, 0]-tp
                    fp = P[column, 0]-tp
                    tn = t-tp-fn-fp

                    if tp < 0 or fn < 0 or fp < 0 or tn < 0:
                        raise ValidationError("True positive, false negative, false postive and true negative should not be negative values")

                    sum_of_truth_row = tp + fn + (fp + tn)*rnd
                    sum_of_pred_row = tp + fp + (fn + tn)*rnd

                    # get smallest possible float -- for 0 values
                    small = np.nextafter(0.,1.)
                    sum_of_truth_row += small
                    sum_of_pred_row += small

                    if tp != 0 or fp != 0 or tn != 0 or fn != 0:
                        sym2 = (
                            np.log(sum_of_pred_row) * (tp + rnd*fp + fn + rnd*tn) -
                            np.log(sum_of_truth_row) * (tp + rnd*fp + fn + rnd*tn) +
                            np.log(rnd+small) * ((rnd+small)*fp - fn)
                            )

                        sym2 /= sum_of_truth_row
                        res += sym2

                    if sym:
                        if tp != 0 or fp != 0 or tn != 0 or fn != 0:
                            sym1 = (
                                np.log(sum_of_truth_row) * (tp + rnd*fn + fp + rnd*tn) -
                                np.log(sum_of_pred_row) * (tp + rnd*fn + fp + rnd*tn) +
                                np.log(rnd+small) * ((rnd+small)*fn - fp)
                                )
    
                            sym1 /= sum_of_pred_row
                            res += sym1

    return res

def om_calculate2_pseudoV_norm(om, rnd=1e-50, max_val=4000, full_matrix=True, modify=False, pseudo_counts=None):
    pv_val = calculate2_pseudoV(om, rnd=rnd, full_matrix=full_matrix, modify=False, pseudo_counts=pseudo_counts)
    return max(1 -  pv_val/ max_val, 0) 

def om_calculate2_sym_pseudoV_norm(om, rnd=1e-50, max_val=8000, full_matrix=True, modify=False, pseudo_counts=None):
    spv_val = om_calculate2_sym_pseudoV(om, rnd=rnd, full_matrix=full_matrix, modify=False, pseudo_counts=pseudo_counts)
    return max(1 - spv_val / max_val, 0)

def om_calculate2_sym_pseudoV(om, rnd=1e-50, full_matrix=True, modify=False, pseudo_counts=None):
    return om_calculate2_pseudoV(om, rnd=rnd, full_matrix=full_matrix, sym=True, modify=modify, pseudo_counts=pseudo_counts)

# outputs the same result as calculate2_spearman, but customized for overlapping matrix
def om_calculate2_spearman(tp, fp, tn, fn, full_matrix = True):
    if (not full_matrix):
        tp = int((tp - np.around(np.sqrt(tp+fn+fp+tn)))/2)
        fn /= 2 
        fp /= 2
        tn /= 2

    # number of ones in pred file
    pos_pred = tp + fp
    # number of ones in truth file
    pos_truth = tp + fn
    # number of zeros in pred file
    neg_pred = tn + fn
    # number of zeros in truth file
    neg_truth = tn + fp

    rank_zero_pred = (neg_pred+1)/2.0
    rank_one_pred = (pos_pred+1)/2.0 + neg_pred

    rank_zero_truth = (neg_truth+1)/2.0
    rank_one_truth = (pos_truth+1)/2.0 + neg_truth

    n = np.sqrt(tp+fn+fp+tn)

    tp_score = ((rank_one_truth - rank_one_pred)/n)**2
    tn_score = ((rank_zero_truth - rank_zero_pred)/n)**2
    fn_score = ((rank_one_truth - rank_zero_pred)/n)**2
    fp_score = ((rank_zero_truth - rank_one_pred)/n)**2

    sum_of_scores = tp_score*tp + tn_score*tn + fn_score*fn + fp_score*fp
 
    row = 1 - (6 * sum_of_scores)/((tp+fn+fp+tn)**2-1)
    return row

# does the same thing as calculate2_aupr, but for overlapping matrix
def om_calculate2_aupr(tp, fp, tn, fn, full_matrix = True):
    if (not full_matrix):
        tp = int((tp - np.around(np.sqrt(tp+fn+fp+tn)))/2)
        fn /= 2 
        fp /= 2
        tn /= 2

    precision = []
    recall = []
    precision.append((tp+fn)/float(tp+fp+tn+fn))
    precision.append(tp / float(tp + fp))
    precision.append(1)
    recall.append(1)
    recall.append(tp / float(tp + fn))
    recall.append(0)
    import sklearn.metrics as mt
    aucpr = mt.auc(np.asarray(recall), np.asarray(precision))
    return aucpr

# does the same thing as calculate2_mcc, but customized for overlapping matrix
def om_calculate2_mcc(tp, fp, tn, fn, full_matrix=True, rnd=1e-50):
    if (not full_matrix):
        tp = int((tp - np.around(np.sqrt(tp+fn+fp+tn)))/2)
        fn /= 2 
        fp /= 2
        tn /= 2

    denom_terms = [(tp+fp), (tp+fn), (tn+fp), (tn+fn)]
    # print tp, fp, tn, fn

    for index, term in enumerate(denom_terms):
        if term == 0:
            denom_terms[index] = 1
    denom = np.sqrt(reduce(np.multiply, denom_terms, 1))
    # print "denom: ", denom

    if tp == 0 and fn == 0:
        num = (tn - fp)
    elif tn == 0 and fp == 0:
        num = (tp - fn)
    else:
        num = (tp*tn - fp*fn)
    # print num / float(denom) 
    return num / float(denom)

def om_validate3A(data_3A, predK, mask=None):
    """Constructs a matrix that describes the relationship between the clusters
    :param data_3A: inputted data for subchallenge 3A
    :param predK: number of clusters
    :param mask: mask used
    :return:
    """
    # read in the data
    data_3A = data_3A.split('\n')
    data_3A = filter(None, data_3A)
    if len(data_3A) != predK:
        raise ValidationError("Input file contains a different number of lines (%d) than expected (%d)")
    data_3A = [x.split('\t') for x in data_3A]
    for i in range(len(data_3A)):
        if len(data_3A[i]) != 2:
            raise ValidationError("Number of tab separated columns in line %d is not 2" % (i+1))
        try:
            data_3A[i][0] = int(data_3A[i][0])
            data_3A[i][1] = int(data_3A[i][1])
        except ValueError:
            raise ValidationError("Entry in line %d could not be cast as integer" % (i+1))

    if [x[0] for x in data_3A] != range(1, predK+1):
        raise ValidationError("First column must have %d entries in acending order starting with 1" % predK)

    for i in range(len(data_3A)):
        if data_3A[i][1] not in set(range(predK+1)):
            raise ValidationError("Parent node label in line %d is not valid." % (i+1))

    # Since cluster zero is not included in file
    ad_cluster = np.zeros((len(data_3A)+1, len(data_3A)+1), dtype=int)
    # file starts at one
    for i in range(len(data_3A)):
        ad_cluster[data_3A[i][1]][data_3A[i][0]] = 1
    # fill in a matrix which tells you whether or not one cluster is a descendant of another
    for i in range(len(data_3A)+1):
        for j in range(len(data_3A)+1):
            if (ad_cluster[j][i] == 1):
                for k in range(len(data_3A)+1):
                    if(ad_cluster[k][j] == 1):
                        ad_cluster[k][i] = 1
                    if(ad_cluster[i][k] == 1):
                        ad_cluster[j][k] = 1

    # check if all nodes are connected. If there are not, we could possibly run the above code again 
    if (not np.array_equal(np.nonzero(ad_cluster[0])[0], map(lambda x: x+1, range(len(data_3A))))):
        for i in range(len(data_3A)+1):
            for j in range(len(data_3A)+1):
                if (ad_cluster[j][i] == 1):
                    for k in range(len(data_3A)+1):
                        if(ad_cluster[k][j] == 1):
                            ad_cluster[k][i] = 1
                        if(ad_cluster[i][k] == 1):
                            ad_cluster[j][k] = 1

    if (not np.array_equal(np.nonzero(ad_cluster[0])[0], map(lambda x: x+1, range(len(data_3A))))):
        raise ValidationError("Root of phylogeny not ancestor of all clusters / Tree is not connected.")

    # print ad_cluster
    ad_cluster = np.delete(ad_cluster, 0, 0)
    ad_cluster = np.delete(ad_cluster, 0, 1)

    return ad_cluster

def construct_relative_matrix(om, pred, truth, mode="ancestor"):
    """Constructs the shared relative matrix
    :param om: overlapping matrix
    :param pred: matrix that specifies the relationship of the clusters of the predicted file
    :param truth: matrix that specifies the relationship of the clusters of the truth file
    :return: a matrix describing the number of shared relatives between each cluster
    """
    if (truth.shape[0] != om.shape[0]):
        raise ValidationError("Row size of overlapping matrix does not match size of matrix of truth file. Size of ad: %d. Size of om: %d." % (truth.shape[0], om.shape[0]))

    if (pred.shape[0] != om.shape[1]):
        raise ValidationError("Column size of overlapping matrix does not match size of matrix of predicition file.  Size of ad: %d. Size of om: %d." % (pred.shape[0], om.shape[1]))    

    srm = np.zeros((om.shape[0], om.shape[1]), dtype=int)
    for i in range(truth.shape[1]):
        for j in range(pred.shape[1]):
            true_relative_list = []
            pred_relative_list = []
            # np.nonzero returns an array of arrays that contain a list of indices that are nonzero
            if mode is "ancestor" or mode is "cousin":
                for row in range(truth.shape[1]):
                    if truth[row, i] == 1:
                        true_relative_list.append(row)
                for row in range(pred.shape[1]):
                    if pred[row, j] == 1:
                        pred_relative_list.append(row)
            else:
                for column in range(truth.shape[1]):
                    if truth[i, column] == 1:
                        true_relative_list.append(column)
                for column in range(pred.shape[1]):
                    if pred[j, column] == 1:
                        pred_relative_list.append(column)
            # print i, true_relative_list
            # print j, pred_relative_list
            for k in range(len(true_relative_list)):
                for l in range(len(pred_relative_list)): 
                    # add sam[i][j] by the shared mutations in their ancestor list
                    srm[i, j] += om[true_relative_list[k], pred_relative_list[l]]
    return srm  

def construct_relative_matrix_opt(om, truth, modification=None):
    """Constructs the shared relative matrix but optimizes for NClusterOneLineage
    :param om: overlapping matrix
    :param truth: matrix that specifies the relationship of the clusters of the truth file
    :return: a matrix describing the number of shared relatives between each cluster
    """
    srm = np.zeros((om.shape[0], om.shape[1]), dtype=int)
    for i in range(truth.shape[1]):
        for j in range(om.shape[1]):
            true_relative_list = []
            for column in range(truth.shape[1]):
                if truth[i, column] == 1:
                    true_relative_list.append(column)
            if modification is "transpose":
                pred_relative_list = np.arange(j)
                # print "j: ", pred_relative_list
            else:
                pred_relative_list = np.arange(j+1, om.shape[1])
            for k in range(len(true_relative_list)):
                for l in range(len(pred_relative_list)): 
                    # add sam[i][j] by the shared mutations in their ancestor list
                    srm[i, j] += om[true_relative_list[k], pred_relative_list[l]]
    return srm

def construct_c_cluster(ad_cluster):
    """Constructs the cousin matrix which describes the relationship between different clusters
    :param ad_cluster: ancestor-descendent matrix which describes the relationship between different clusters
    :return: 
    """
    c_cluster = np.ones((ad_cluster.shape[0], ad_cluster.shape[1]))
    old_cluster = np.copy(ad_cluster)
    for row in range(old_cluster.shape[0]):
        for column in range(old_cluster.shape[1]):
            if old_cluster[row, column] == 1:
                old_cluster[column, row] = 1

    # The idea here is that the ancestor-descendent matrix, the descendent-ancestor matrix, 
    # the cousin matrix, and the co-clustering matrix should sum to a matrix of ones
    c_cluster = c_cluster - old_cluster - np.identity(old_cluster.shape[0])
    return c_cluster

def construct_related_mutations_matrix(om, pred, truth, mode="ancestor"):
    """Constructs the related mutations matrices which describes the how many relatives each cluster has
    :param om: overlap matrix
    :param pred: matrix which describes the relationship between different clusters in prediction file
    :param truth: matrix which describes the relationship between different clusters in truth file
    :param mode: describes the relation that is described with the matrices (i.e. ancestor, descendent, cousin)
    :return: 
    """
    # T is a matrix with length equal to the number of clusters in the truth file and width 1; each entry is equal 
    # to the number of ancestor mutations each cluster have
    T = np.zeros((truth.shape[0], 1))
    # P is a matrix with length equal to the number of clusters in the pred file and width 1; each entry is equal 
    # to the number of ancestor mutations each cluster have
    P = np.zeros((pred.shape[0], 1))
    for i in range(truth.shape[0]):
        for j in range(truth.shape[0]):
            # if om[j, i]
            if mode is "ancestor" or mode is "cousin":
                if (truth[j, i] == 1):
                    T[i, 0] += np.sum(om[j])
            else:
                if (truth[i, j] == 1):
                    T[i, 0] += np.sum(om[j])

    for i in range(pred.shape[0]):
        for j in range(pred.shape[0]):
            # if om[j, i]
            if mode is "ancestor" or mode is "cousin":
                if (pred[j, i] == 1):
                    P[i, 0] += np.sum(om[:, j])
            else:
                if (pred[i, j] == 1):
                    P[i, 0] += np.sum(om[:, j])

    return P, T

# Note that this function is not used right now.... But would be needed if we decided to use another metric
# in place of the pseudoV metric
def calculate_descendant_matrix(srm, om, P, T):
    """Calculates the true postive, true negatives, false postives and false negatives for subchallenge 3
    :param srm: shared relative matrix
    :param om: overlap matrix
    :param P: matrix which describes the number of relatives each cluster has in the prediction file
    :param T: matrix which describes the number of relatives each cluster has in the truth file
    :return: 
    """
    tp = 0
    fp = 0
    fn = 0
    tn = 0
    t = 0
    for row in range(om.shape[0]):
        for column in range(om.shape[1]):
            tp += om[row, column]*srm[row, column]

    for column in range(om.shape[1]):
        fp += np.sum(om[:,column]) * P[column, 0]
    fp -= tp

    for row in range(om.shape[0]):
        fn += np.sum(om[row]) * T[row, 0]
    fn -= tp

    for row in range(om.shape[0]):
        for column in range(om.shape[1]):
            t+=om[row, column]

    tn = t**2-tp-fp-fn

    return tp, fp, tn, fn

# This method is a wrapper class for calculate3_pseudoV; ad_pred and ad_truth are adjusted according to the modification
def calculate3A_final(om, ad_pred, ad_truth, modification=None, rnd=1e-50):
    """Calculates the pseudoV score given an overlap matrix, and matrices which describes the relationship of clusters in both
    pred and truth files
    :param om: overlap matrix
    :param ad_pred: matrix which describes the relationshop between clusters has in the prediction file
    :param ad_truth: matrix which describes the relationship between clusters in the truth file
    :param modification: none, transpose or cousin
    :return: 
    """
    if modification is None:
        pred = ad_pred
        truth = ad_truth
    elif modification is "transpose":
        pred = ad_pred.T
        truth = ad_truth.T
    elif modification is "cousin":
        pred = construct_c_cluster(ad_pred)
        truth = construct_c_cluster(ad_truth)
    else:
        raise ValidationError('Incorrect modificiation')
    P, T = construct_related_mutations_matrix(om, pred, truth, mode="descendant")
    srm = construct_relative_matrix(om, pred, truth, mode="descendant")

    return calculate3_js_divergence(srm, om, P, T, rnd=rnd)
#    return calculate3_pseudoV(srm, om, P, T, rnd=rnd)

def calculate3A_worst(om, ad_pred, ad_truth, scenario="OneCluster", modification=None, truth_data=None, rnd=1e-50):
    """Calculates the worst score given an overlap matrix and matrices which describes the relationship of clusters in both
    pred and truth files
    :param om: overlap matrix
    :param ad_pred: matrix which describes the relationshop between clusters has in the prediction file
    :param ad_truth: matrix which describes the relationship between clusters in the truth file
    :param modification: none, transpose or cousin
    :truth_data: the order in which the mutations were inputted; this is necessary to calculate the score for N Cluster N Lineages
    :return: 
    """
    if modification is None:
        pred = ad_pred
        truth = ad_truth
        P, T = construct_related_mutations_matrix(om, pred, truth, mode="descendant")
        worst_srm, worst_om, worst_P = get_bad_ad_and_om(om, truth, scenario=scenario, truth_data=truth_data)
    elif modification is "transpose":
        pred = ad_pred.T
        truth = ad_truth.T
        P, T = construct_related_mutations_matrix(om, pred, truth, mode="descendant")
        worst_srm, worst_om, worst_P = get_bad_ad_and_om(om, truth, modification="transpose", scenario=scenario, truth_data=truth_data)
    elif modification is "cousin":
        pred = construct_c_cluster(ad_pred)
        truth = construct_c_cluster(ad_truth)
        P, T = construct_related_mutations_matrix(om, pred, truth, mode="descendant")
        worst_srm, worst_om, worst_P = get_bad_c_and_om(om, scenario=scenario)
    else:
        raise ValidationError('Incorrect modficiation')

    return calculate3_js_divergence(worst_srm, worst_om, worst_P, T, rnd=rnd)
#    return calculate3_pseudoV(worst_srm, worst_om, worst_P, T, rnd=rnd)

# This method is equivalent to the calculate3Final in the orignal script
def calculate3A(om, truth_data, ad_pred, ad_truth, rnd=1e-50):
    """Calculates the score given an overlap matrix and matrices which describes the relationship of clusters in both
    pred and truth files
    :param om: overlap matrix
    :param ad_pred: matrix which describes the relationshop between clusters has in the prediction file
    :param ad_truth: matrix which describes the relationship between clusters in the truth file
    :truth_data: the order in which the mutations were inputted; this is necessary to calculate the score for N Cluster N Lineages
    :return: the score for subchallenge 3
    """
    n_scores = []
    n_scores.append(set_to_zero(calculate3A_worst(om, ad_pred, ad_truth, scenario="NCluster", truth_data=truth_data, rnd=rnd)))
    n_scores.append(set_to_zero(calculate3A_worst(om, ad_pred, ad_truth, scenario="NCluster", modification="transpose", truth_data=truth_data, rnd=rnd)))
    n_scores.append(set_to_zero(calculate3A_worst(om, ad_pred, ad_truth, scenario="NCluster", modification="cousin", rnd=rnd)))
    del truth_data
    gc.collect()
    n_scores_permuted = []
    
    P, T = construct_related_mutations_matrix(om, ad_pred, ad_truth, mode="descendant")
    n_scores_permuted.append(set_to_zero(om_permute_N_cluster(om, T, rnd=rnd)))
    P, T = construct_related_mutations_matrix(om, ad_pred.T, ad_truth.T, mode="descendant")
    n_scores_permuted.append(set_to_zero(om_permute_N_cluster(om, T, rnd=rnd)))
    n_scores_permuted.append(n_scores[2])

    scores = []
    scores.append(set_to_zero(calculate3A_final(om, ad_pred, ad_truth, rnd=rnd)))
    scores.append(set_to_zero(calculate3A_final(om, ad_pred, ad_truth, modification="transpose", rnd=rnd)))
    scores.append(set_to_zero(calculate3A_final(om, ad_pred, ad_truth, modification="cousin", rnd=rnd)))

    one_scores = []
    one_scores.append(set_to_zero(calculate3A_worst(om, ad_pred, ad_truth, scenario="OneCluster", rnd=rnd)))
    one_scores.append(set_to_zero(calculate3A_worst(om, ad_pred, ad_truth, scenario="OneCluster", modification="transpose", rnd=rnd)))
    one_scores.append(set_to_zero(calculate3A_worst(om, ad_pred, ad_truth, scenario="OneCluster", modification="cousin", rnd=rnd)))

    score = sum(scores) / 3.0
    one_score = sum(one_scores) / 3.0
    n_score = sum(n_scores) / 3.0

    n_score_permuted = sum(n_scores_permuted)/3.0

#    return [set_to_zero((1 - (score / max(one_score, n_score)))), set_to_zero((1 - (score / max(one_score, n_score_permuted))))]  
    return  set_to_zero((1 - (score / max(one_score, n_score_permuted))))

# adds num pseudo counts to the om; default is square root of the number of mutations
def add_pseudo_counts_om(om, num=None):
    """Adds pseudo counts to the sample
    :param om: overlap matrix
    :param num: number of pseudo counts to add
    :return: an overlap matrix with pseudo counts added
    """

    if num is None:
        N = 0
        for row in range(om.shape[0]):
            for column in range(om.shape[1]):
                N += om[row, column]
        num = np.floor(np.sqrt(N))

    old_width = om.shape[0]
    om_width = int(old_width + num)

    old_length = om.shape[1]
    om_length = int(old_length + num)

    new_om = np.zeros((om_width, om_length), dtype = int)

    for row in range(old_width):
        for column in range(old_length):
            new_om[row, column] = om[row, column]

    for i in range(int(np.floor(np.sqrt(N)))):
        new_om[old_width+i, old_length+i] = 1

    return new_om

# adds num pseudo counts to the om; default is square root of the number of mutations
def add_pseudo_counts_om_eff(tp, fp, tn, fn, num=None):
    """Adds pseudo counts to the sample
    :param tp, fp, tn, fn: true postives, false postives, true negatives, false negatives
    :param num: number of pseudo counts to add
    :return: new true positives, false postives, true negatives, false negatives after pseudo counts are added
    """


    N = np.floor(np.sqrt(tp+fp+tn+fn))
    K = np.floor(np.sqrt(N))
    if num is not None:
        K = num
    tp += K
    tn += K**2 + 2*N*K - K
    return tp, fp, tn, fn 

# The equivalent to get_worst_score, but for overlap matrices
def get_worst_score_om(om, scoring_func, larger_is_worse=True, rnd=1e-50):
    """calculates worst score for subchallenge 2A given a function
    :param overlap matrix
    :param scoring_func: the method to be used
    :param larger_is_worse: whether or not larger means a worse score
    :return: the worst score for subchallenge 2A
    """
    if larger_is_worse:
        return max(get_bad_score_om(om, scoring_func, 'OneCluster', rnd=rnd),
                    get_bad_score_om(om, scoring_func, 'NCluster', rnd=rnd))
    else:
        return min(get_bad_score_om(om, scoring_func, 'OneCluster', rnd=rnd),
                    get_bad_score_om(om, scoring_func, 'NCluster', rnd=rnd))


# Equivalent to get_bad_score in original function
def get_bad_score_om(om, score_func, scenario='OneCluster', pseudo_counts=None, rnd=1e-50):
    if score_func is om_calculate2_pseudoV or score_func is om_calculate2_pseudoV_norm or score_func is om_calculate2_sym_pseudoV or score_func is om_calculate2_js_divergence:
        bad_om = get_bad_om(om, scenario)
        return score_func(bad_om, modify=True, pseudo_counts=pseudo_counts, rnd=rnd)
    else:
        tp, fp, tn, fn = add_pseudo_counts_om_eff(*calculate_overlap_matrix(get_bad_om(om, scenario)))
        return score_func(tp, fp, tn, fn, rnd=rnd)

# Eqiuivalent to get_bad_ccm 
def get_bad_om(om, scenario='OneCluster'):
    """constructs the worst om
    :param om: overlap matrix
    :param scenario: the scenario that will be used (OneCluster or NCluster)
    :return: the worst overlap matrix fir a given scenario
    """
    num_mutations = 0
    for row in range(om.shape[0]):
        num_mutations += np.sum(om[row])

    if scenario is 'NCluster':
        worst_matrix = np.zeros([om.shape[0], num_mutations], dtype=int)
        start = 0
        for row in range(om.shape[0]):
            cluster_length = np.sum(om[row])
            if cluster_length > num_mutations:
                raise ValidationError('Number of mutations in cluster %i exceeds total number of mutations' % (row+1))
            for column in range(start, start+cluster_length):
                worst_matrix[row][column] = 1
            start += cluster_length
        return worst_matrix 
    elif scenario is 'OneCluster':
        worst_matrix = np.zeros([om.shape[0], 1], dtype=int)
        for row in range(om.shape[0]):
            worst_matrix[row] = np.sum(om[row])
        return worst_matrix
    else:
        raise ValueError('Scenario must be one of OneCluster or NCluster')

# returns a series a matrices that is needed to calculate the worst score for subchallenge 3A
def get_bad_c_and_om(om, scenario='OneCluster'):
    if scenario is 'OneCluster':
        worst_om = get_bad_om(om, scenario=scenario)
        worst_c = np.zeros((om.shape[0], 1))
        P = np.zeros((1, 1))
    elif scenario is 'NCluster':
        worst_om = get_bad_om(om, scenario=scenario)
        worst_c = np.zeros((worst_om.shape[0], worst_om.shape[1]))
        P = np.zeros((worst_om.shape[1], 1))
    else:
        raise ValueError('Scenario must be one of OneCluster or NCluster')

    return worst_c, worst_om, P 

# returns a series a matrices that is needed to calculate the worst score for subchallenge 3A
def get_bad_ad_and_om(om, ad_truth, modification=None, scenario='OneCluster', truth_data=None):
    if scenario is 'OneCluster':
        worst_om = get_bad_om(om, scenario=scenario)
        worst_ad = np.zeros((om.shape[0], 1))
        P = np.zeros((1, 1))
    elif scenario is 'NCluster':
        if truth_data is None:
            raise ValidationError("Constructing the overlap matrix for NCluster requires the cluster information of every single mutation!")
        worst_om = get_NCluster_om_3A(om, truth_data)
        P = np.matrix(np.arange(worst_om.shape[1]))
        P.resize(worst_om.shape[1], 1)
        if modification is not "transpose":
            P = P[::-1]
        worst_ad = construct_relative_matrix_opt(worst_om, ad_truth, modification=modification)
    else:
        raise ValueError('Scenario must be one of OneCluster or NCluster')
    return worst_ad, worst_om, P

# returns the overlap matrix of subchallenge 3A, N Cluster N lineages
def get_NCluster_om_3A(om, truth_data):
    worst_om = np.zeros((om.shape[0], len(truth_data)), dtype=int)
    for i in range(len(truth_data)):
        worst_om[truth_data[i]-1, i] = 1

    return worst_om

# creates an overlapping matrix from predicted and truth co-clustering matrices; useful for debugging
def create_om(ccm_pred, ccm_truth):
    mutations_pred = np.zeros((ccm_pred.shape[0], 1), dtype=int)
    mutations_truth = np.zeros((ccm_truth.shape[0], 1), dtype=int)
    cluster_count = 1
    for row in range(ccm_pred.shape[0]):
        if mutations_pred[row, 0] == 0:
            for i in range(ccm_pred.shape[1]):
                if ccm_pred[row, i] == 1:
                    mutations_pred[i, 0] = cluster_count
            cluster_count+=1

    cluster_count = 1
    for row in range(ccm_truth.shape[0]):
        if mutations_truth[row, 0] == 0:
            for i in range(ccm_truth.shape[1]):
                if ccm_truth[row, i] == 1:
                    mutations_truth[i, 0] = cluster_count
            cluster_count+=1

    num_truth_clusters = mutations_truth.max()
    num_pred_clusters = mutations_pred.max()

    om = np.zeros((num_truth_clusters, num_pred_clusters), dtype=int)

    for i in range(mutations_pred.shape[0]):
        om[mutations_truth[i,0]-1][mutations_pred[i,0]-1] += 1

    return om

def pack(*args):
    mylist = []
    for element in args:
        mylist.append(element)
    return mylist


def adj_final(res):
    if ((res-1) < 0.00001 and res > 1):
        res = 1
    return res

def set_to_zero(res):
    if(res < 0):
        res = 0
    return res
