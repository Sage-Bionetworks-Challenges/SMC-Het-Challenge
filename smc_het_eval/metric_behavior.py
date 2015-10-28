from SMCScoring import *
import numpy as np
import csv

tsv_dir = './scoring_metric_data/text_files/' # directory to save tsv's to

def scoring1A_behavior():
    guesses = np.array(range(0,101))/100.0
    truths = np.array(range(10,110,10))/100.0
    res = []
    for i in range(len(truths)):
        for j in range(len(guesses)):
            res.append( [truths[i],guesses[j],calculate1A(guesses[j],truths[i])] )
    res = [map(str,x) for x in res]
    res = ['\t'.join(x) for x in res]
    f = open('scoring1A_behavior.tsv', 'w')
    f.write('\n'.join(res))
    f.close()


def scoring1B_behavior():
    guesses = np.array(range(1,11))
    truths = np.array(range(1,6))
    res = [] 
    for i in range(len(truths)):
        for j in range(len(guesses)):
            res.append( [truths[i],guesses[j],calculate1B(guesses[j],truths[i])] )
    res = [map(str,x) for x in res]
    res = ['\t'.join(x) for x in res]
    f = open('scoring1B_behavior.tsv', 'w')
    f.write('\n'.join(res))
    f.close()

def scoring1C_behavior():
    # Baseline Truth
    t_phis = np.array([.85,.5,.3])
    t_nssms = np.array([200,200,200])
    t_entry = zip(t_nssms,t_phis)

    # Zero-mean noise in phi
    res = []
    concentration = [2, 5, 10, 100, 1000]
    for c in concentration:
        for i in range(100):
            phis = []
            for p in t_phis:
                phis.append(np.random.beta(p*c,(1-p)*c))
            res.append([c,phis,calculate1C(t_entry,zip(t_nssms,phis))])
    res = [map(str,x) for x in res]
    res = ['\t'.join(x) for x in res]
    f = open('scoring1C_phi_ZM_behavior.tsv', 'w')
    f.write('\n'.join(res))
    f.close()

    # Systematic over/underestimation of phi
    sys_errors = np.array([.01,.03,.05,0.075,.10])
    sys_errors = np.concatenate((sys_errors,-sys_errors))
    res = []
    for sys_error in sys_errors:
        res.append([sys_error, calculate1C(t_entry,zip(t_nssms,t_phis+sys_error))])
    res = [map(str,x) for x in res]
    res = ['\t'.join(x) for x in res]
    f = open('scoring1C_phi_sys_behavior.tsv', 'w')
    f.write('\n'.join(res))
    f.close()

    res = []
    concentration = [2, 5, 10, 100, 1000]
    for c in concentration:
        for i in range(100):
            rd = np.random.dirichlet([c,c,c]) * sum(t_nssms)
            rd = map(round,rd)
            remainder = sum(t_nssms) - sum(rd)
            #Randomly assign remainder
            rd[np.random.randint(0,len(rd))] += remainder
            rd = map(int,rd)
            res.append([c,rd,calculate1C(t_entry,zip(rd,t_phis))])
    res = [map(str,x) for x in res]
    res = ['\t'.join(x) for x in res]
    f = open('scoring1C_nssm_behavior.tsv', 'w')
    f.write('\n'.join(res))
    f.close()            


    res = []
    # Collapse first two clusters
    phis = [(.85+.5)/2.0, .3]
    nssms = [400,200]
    entry = zip(nssms,phis)
    res.append(["Collapse12", calculate1C(t_entry,entry)])
    
    # Collapse last two clusters
    phis = [.85, (.5+.3)/2.0]
    nssms = [200,400]
    entry = zip(nssms,phis)
    res.append(["Collapse23", calculate1C(t_entry,entry)])

    # Collapse all clusters
    phis = [.55]
    nssms=[600]
    entry = zip(nssms,phis)
    res.append(["Collapse123", calculate1C(t_entry,entry)])

    # Assume all SSMs are clonal
    phis = [.85]
    nssms=[600]
    entry = zip(nssms,phis)
    res.append(["All_Clonal", calculate1C(t_entry,entry)])

    # For splits, phis were calculated as +/- 0.05 from center.  
    # 0.05 was obtained empirically by calculating the mean of the top and bottom half of a binomial sample with depth = 50. e.g.:
    # s = np.random.binomial(50,.85,(1000))
    # np.mean(s[s<np.median(s)]) / 50.0
    # 0.90
    # np.mean(s[s>np.median(s)]) / 50.0
    # 0.80

    # Split cluster 1
    phis = [.9,.8, .5, .3]
    nssms = [100,100,200,200]
    entry = zip(nssms,phis)
    res.append(["Split1", calculate1C(t_entry,entry)])

    # Split cluster 2
    phis = [.85, .55,.45, .3]
    nssms = [200,100,100,200]
    entry = zip(nssms,phis)
    res.append(["Split2", calculate1C(t_entry,entry)])
    # Split cluster 3
    phis = [.85, .5,.35,.25]
    nssms = [200,200,100,100]
    entry = zip(nssms,phis)
    res.append(["Split3", calculate1C(t_entry,entry)])
    
    res = [map(str,x) for x in res]
    res = ['\t'.join(x) for x in res]
    f = open('scoring1C_cases.tsv', 'w')
    f.write('\n'.join(res))
    f.close()        

def scoring2A_behavior():
    # True CCM:
    t_clusters = np.zeros((600,3))
    t_clusters[0:200,0] = 1
    t_clusters[200:400,1] = 1
    t_clusters[400:,2] = 1
    t_ccm = np.dot(t_clusters,t_clusters.T)

    # Cases:
    res = []
    # Split a cluster
    clusters = np.zeros((600,4))
    clusters[0:100,0] = 1
    clusters[100:200,3] = 1
    clusters[200:400,1] = 1
    clusters[400:,2] = 1
    ccm = np.dot(clusters,clusters.T)
    res.append(["SplitCluster",calculate2(ccm,t_ccm)])

    # Merge 2 clusters
    clusters = np.zeros((600,2))
    clusters[0:400,0] = 1
    clusters[400:,1] = 1
    ccm = np.dot(clusters,clusters.T)
    res.append(["MergeCluster",calculate2(ccm,t_ccm)])
    
    #All one cluster
    res.append(["OneCluster",calculate2(np.ones(t_ccm.shape),t_ccm)])
    # Each ssm own cluster
    res.append(["NClusters",calculate2(np.identity(t_ccm.shape[0]),t_ccm)])
    
    # Small extra cluster with some mutations from each cluster
    clusters = np.zeros((600,4))
    clusters[:,:-1] = np.copy(t_clusters)
    clusters[100,:] = np.array([0,0,0,1])
    clusters[300,:] = np.array([0,0,0,1])
    clusters[500,:] = np.array([0,0,0,1])
    ccm = np.dot(clusters,clusters.T)
    res.append(["SmallExtra", calculate2(ccm,t_ccm)])
    
    # Big extra cluster with some mutations from each cluster
    clusters = np.zeros((600,4))
    clusters[:,:-1] = np.copy(t_clusters)
    clusters[100:133,:] = np.array([0,0,0,1])
    clusters[300:333,:] = np.array([0,0,0,1])
    clusters[500:533,:] = np.array([0,0,0,1])
    ccm = np.dot(clusters,clusters.T)
    res.append(["BigExtra", calculate2(ccm,t_ccm)])

    res = [map(str,x) for x in res]
    res = ['\t'.join(x) for x in res]
    f = open('scoring_metric_data/scoring2A_cases.tsv', 'w')
    f.write('\n'.join(res))
    f.close()
    
    # Same thing but with more groups (all of the same size)
    # True CCM:
    t_clusters = np.zeros((2000,10))
    t_clusters[0:200,0] = 1
    t_clusters[200:400,1] = 1
    t_clusters[400:600,2] = 1
    t_clusters[600:800,3] = 1
    t_clusters[800:1000,4] = 1
    t_clusters[1000:1200,5] = 1
    t_clusters[1200:1400,6] = 1
    t_clusters[1400:1600,7] = 1
    t_clusters[1600:1800,8] = 1
    t_clusters[1800:2000,9] = 1
    
    t_ccm = np.dot(t_clusters,t_clusters.T)

    # Cases:
    res_more_cl = []
    # Split a cluster
    clusters = np.zeros((2000,11))
    clusters[:,:-1] = np.copy(t_clusters)
    clusters[1900:2000,9] = 0
    clusters[1900:2000,10] = 1
    ccm = np.dot(clusters,clusters.T)
    res_more_cl.append(["SplitCluster",calculate2(ccm,t_ccm)])

    # Merge 2 clusters
    clusters = np.copy(t_clusters[:,:-1])
    clusters[1800:2000,8] = 1
    ccm = np.dot(clusters,clusters.T)
    res_more_cl.append(["MergeCluster",calculate2(ccm,t_ccm)])
    
    #All one cluster
    res_more_cl.append(["OneCluster",calculate2(np.ones(t_ccm.shape),t_ccm)])
    # Each ssm own cluster
    res_more_cl.append(["NClusters",calculate2(np.identity(t_ccm.shape[0]),t_ccm)])
    
    # Small extra cluster with some mutations from each cluster
    clusters = np.zeros((2000,11))
    clusters[:,:-1] = np.copy(t_clusters)
    clusters[100,:] = np.array([0,0,0,0,0,0,0,0,0,0,1])
    clusters[300,:] = np.array([0,0,0,0,0,0,0,0,0,0,1])
    clusters[500,:] = np.array([0,0,0,0,0,0,0,0,0,0,1])
    clusters[700,:] = np.array([0,0,0,0,0,0,0,0,0,0,1])
    clusters[900,:] = np.array([0,0,0,0,0,0,0,0,0,0,1])
    clusters[1100,:] = np.array([0,0,0,0,0,0,0,0,0,0,1])
    clusters[1300,:] = np.array([0,0,0,0,0,0,0,0,0,0,1])
    clusters[1500,:] = np.array([0,0,0,0,0,0,0,0,0,0,1])
    clusters[1700,:] = np.array([0,0,0,0,0,0,0,0,0,0,1])
    clusters[1900,:] = np.array([0,0,0,0,0,0,0,0,0,0,1])
    ccm = np.dot(clusters,clusters.T)
    res_more_cl.append(["SmallExtra", calculate2(ccm,t_ccm)])
    
    # Big extra cluster with some mutations from each cluster
    clusters = np.zeros((2000,11))
    clusters[:,:-1] = np.copy(t_clusters)
    clusters[100:110,:] = np.array([0,0,0,0,0,0,0,0,0,0,1])
    clusters[300:310,:] = np.array([0,0,0,0,0,0,0,0,0,0,1])
    clusters[500:510,:] = np.array([0,0,0,0,0,0,0,0,0,0,1])
    clusters[700:710,:] = np.array([0,0,0,0,0,0,0,0,0,0,1])
    clusters[900:910,:] = np.array([0,0,0,0,0,0,0,0,0,0,1])
    clusters[1100:1110,:] = np.array([0,0,0,0,0,0,0,0,0,0,1])
    clusters[1300:1310,:] = np.array([0,0,0,0,0,0,0,0,0,0,1])
    clusters[1500:1510,:] = np.array([0,0,0,0,0,0,0,0,0,0,1])
    clusters[1700:1710,:] = np.array([0,0,0,0,0,0,0,0,0,0,1])
    clusters[1900:1910,:] = np.array([0,0,0,0,0,0,0,0,0,0,1])
    ccm = np.dot(clusters,clusters.T)
    res_more_cl.append(["BigExtra", calculate2(ccm,t_ccm)])    
    
    res_more_cl = [map(str,x) for x in res_more_cl]
    res_more_cl = ['\t'.join(x) for x in res_more_cl]
    f = open('scoring_metric_data/scoring2A_big_cases.tsv', 'w')
    f.write('\n'.join(res_more_cl))
    f.close()
    
    # Random re-assignment to arbitrary cluster with p=0.01,.03,.05,.1,.15,.25,.5 x 100
    res = []
    p_errors = [0.01,0.03,0.05,.1,.15,.25,.5]
    for p_err in p_errors:
        for i in range(100):
            clusters = np.copy(t_clusters)
            for j in range(t_ccm.shape[0]):
                if np.random.random() < p_err:
                    cluster = np.argmax(clusters[j,:])
                    if np.random.random() < 0.5:
                        cluster -= 1
                    else:
                        cluster += 1
                    cluster = cluster % 3
                    clusters[j,:] = 0
                    clusters[j,cluster] = 1
            ccm = np.dot(clusters,clusters.T)
            res.append([p_err,calculate2(ccm,t_ccm)])
    res = [map(str,x) for x in res]
    res = ['\t'.join(x) for x in res]
    f = open('scoring2A_random_reassignment.tsv', 'w')
    f.write('\n'.join(res))
    f.close()
    
    # Random re-assignment to closest cluster with p=0.01,.03,.05,.1,.15,.25,.5 x 100
    res = []
    p_errors = [0.01,0.03,0.05,.1,.15,.25,.5]
    for p_err in p_errors:
        for i in range(100):
            clusters = np.copy(t_clusters)
            for j in range(t_ccm.shape[0]):
                if np.random.random() < p_err:
                    cluster = np.argmax(clusters[j,:])
                    if cluster == 2 or cluster == 0:
                        cluster = 1
                    else:
                        if np.random.random() < 0.5:
                            cluster = 2
                        else:
                            cluster = 0
                    clusters[j,:] = 0
                    clusters[j,cluster] = 1
            ccm = np.dot(clusters,clusters.T)
            res.append([p_err,calculate2(ccm,t_ccm)])
    res = [map(str,x) for x in res]
    res = ['\t'.join(x) for x in res]
    f = open('scoring2A_closest_reassignment.tsv', 'w')
    f.write('\n'.join(res))
    f.close()
    
def scoring2B_behavior():
    t_clusters = np.zeros((600,3))
    t_clusters[0:200,0] = 1
    t_clusters[200:400,1] = 1
    t_clusters[400:,2] = 1
    t_ccm = np.dot(t_clusters,t_clusters.T)

    n_uniq = len(np.triu_indices(t_ccm.shape[0],k=1)[0])
    res = []
    concentrations = [1000,100,50,25,10,5,3,1]
    for c in concentrations:
        for i in range(50):
            ccm = np.copy(t_ccm)
            ccm[np.triu_indices(t_ccm.shape[0],k=1)] -= np.random.beta(1,c,n_uniq)
            #ccm[np.tril_indices(t_ccm.shape[0],k=-1)] = ccm[np.triu_indices(t_ccm.shape[0],k=1)]
            ccm[np.tril_indices(t_ccm.shape[0],k=-1)] = 0
            ccm = ccm + ccm.T
            np.fill_diagonal(ccm,1)
            ccm = np.abs(ccm)
            res.append([c,calculate2(ccm,t_ccm)])
    res = [map(str,x) for x in res]
    res = ['\t'.join(x) for x in res]
    f = open('scoring_metric_data/scoring2B_beta.tsv', 'w')
    f.write('\n'.join(res))
    f.close()  
    
def scoring3A_behavior(method="orig_nc", verbose=False, weights=None, save=True, pc_amount='more', full_matrix=True):
    '''Scoring behaviour of subchallenge 3 metrics

    Attributes:
    :param method: method or list of methods to use when evaluating each subchallenge scenario
    :param verbose: boolean for whether to output details of the scoring metrics
    :param weights: weights to pass into the scoring function for calculating the weighted average of multiple metrics
    :param save: boolean for whether or not to save the results of the scoring behaviour to a tsv file
    :param more_pc: booelane for whether to include more pseudo counts for each matrix
    '''
    # True values for each attribute
    t_ccm, t_clusters = get_ccm("Truth")
    t_ad = get_ad("Truth")

    res = list() # results of using the given method to score each scenario

    if pc_amount is 'more':
        n_pc = np.sqrt(t_ccm.shape[0]) # number of pseudo counts to include
        pc_ext = "_more_pc" # file name extension for the amount of pseudo counts to include
    elif pc_amount is 'less':
        n_pc = np.log(t_ccm.shape[0])
        pc_ext = ""
    elif pc_amount is 'none':
        n_pc = 0
        pc_ext = "_no_pc"

    for scenario in scenarios:
        if verbose:
            print '\nScenario %s' % scenario

        if scenario == 'Truth':
            res.append(['Truth',
                        calculate3(t_ccm,t_ad,t_ccm,t_ad,
                                   method=method, verbose=verbose, weights=weights, pseudo_counts=n_pc, full_matrix=full_matrix)])
        else:
            ccm = get_ccm(scenario, t_ccm=t_ccm)
            ad = get_ad(scenario, t_ad=t_ad)
            res.append([scenario,
                        calculate3(ccm,ad,t_ccm,t_ad,
                                   method=method, verbose=verbose, weights=weights, pseudo_counts=n_pc, full_matrix=full_matrix)])


    if save:
        if full_matrix:
            tri_ext = "_full"
        else:
            tri_ext = "_triu"
        if isinstance(method, list):
            f = open(tsv_dir + 'scoring3A_all_cases_' + '_'.join(method) + pc_ext + tri_ext + '.tsv', 'w')
        else:
            f = open(tsv_dir + 'scoring3A_all_cases_' + method + pc_ext + tri_ext + '.tsv', 'w')
        out_res = [map(str,x) for x in res]
        out_res = ['\t'.join(x) for x in out_res]
        f.write('\n'.join(out_res))
        f.close()

    return res

def scoring3A_behavior_all(verbose=True):
    for method in ['pseudoV_nc',
                    'pseudoV',
                    'orig_nc',
                    'orig',
                    'mcc_nc',
                    'pearson_nc',
                    'spearman_nc',
                    'aupr_nc',
                    'sqrt_nc',
                    'sym_pseudoV_nc',
                   ['pseudoV_nc', 'mcc_nc', 'pearson_nc'],
                   ['pseudoV_nc', 'pearson_nc', 'sym_pseudoV_nc'],
                   ['aupr_nc', 'sqrt_nc', 'sym_pseudoV_nc'],
                   ['aupr_nc', 'sqrt_nc', 'sym_pseudoV_nc', 'pearson_nc']]:
        for fm in [True, False]:
            for pc in ['none', 'less', 'more']:
                scoring3A_behavior(method=method, verbose=verbose,pc_amount=pc, full_matrix=fm)
                print 'Done %s - More PC: %s - Full Matrix: %s' % (method,pc,fm)

def scoring3A_weight_behavior(methods=["pseudoV_nc", "pearson_nc", "sym_pseudoV_nc"], verbose=False, res=None):
    '''Create the data on how the weights used in subchallenge 3 affect the score using the given scoring methods

    Attributes:
    :param methods: list of methods to use when evaluating each subchallenge scenario
    :param verbose: boolean for whether to output details of the scoring metrics
    '''
    # True values for each attribute
    if res is None:
        res = np.transpose(np.asarray([[row[1] for row in scoring3A_behavior(method, verbose=verbose)] for method in methods]))
    print res

    wght_res = {'Case':scenarios}
    n_method = len(methods)
    weights = [0,0.5,1]
    weight_list = get_weights(n_method, weights)
    print weight_list
    for wght in weight_list:
        norm_wght = wght / float(sum(wght))
        scores = np.sum(norm_wght * res, 1)
        if set(wght) == {0,weights[-1]} or set(wght) == {weights[-1]}:
            key = '+'.join([methods[i] for i in range(n_method) if wght[i] == weights[-1]])
            wght_res[key] = scores
        else:
            wght_res[str(wght)] = scores

    with open(tsv_dir + 'weights3A_all_cases_' + '_'.join(methods) + '.tsv', 'wb') as f:
        fields = sorted(wght_res.keys())
        print fields
        w = csv.DictWriter(f, delimiter='\t', fieldnames=fields)
        w.writeheader()
        for i in range(len(scenarios)):
            w.writerow({field:wght_res[field][i] for field in fields})

    return wght_res, res

def get_weights(n_objects, weights, cur_weight=None, list_weights= list()):
    """Generate a list of all possible unique combinations of weight assignments from the given list
    of weights. Does not make assumptions about the contents of the list, so some weight combinations
    may be duplicates if the list of weights contains elements that are multiples of each other.

    :param n_objects: number of objects that need to be assigned weights
    :param weights: list of possible weights
    :param cur_weight: current iteration of weight assignments
    :param list_weights: list of all weight assignments generated so far
    :return: list of unique weight assignments
    """
    if n_objects is 0:
        if len(set(cur_weight)) > 1 or cur_weight[0] is weights[-1]: # only include the case where all the weights are the same once
            if not (len(set(cur_weight)) <= 2 and 0 in set(cur_weight)) or weights[-1] in set(cur_weight): # only include the cases where one or more of the weights
                                                                                                        # are zero and all the other weights are the same once
                wght = np.copy(cur_weight)
                list_weights.append(wght)
    else:
        if cur_weight is None:
            cur_weight = [0] * n_objects
        for i in weights:
            cur_weight[n_objects-1] = i
            get_weights(n_objects-1, weights, cur_weight, list_weights)

    return np.asarray(list_weights)


def overall_score(p_cell, t_cell, p_ncluster, t_ncluster, p_1c, t_1c,  p_ccm, t_ccm, p_ad, t_ad, weights = [2000, 2000, 3000, 6000, 7000], verbose=False):
    '''Get the overall score for a prediction

    Attributes:
    :params p_cell, t_cell: predicted and true cellularity
    :params p_ncluster, t_ncluster: predicted and true number of clusters
    :params p_1c, t_1c: zippped variable with the size of each clusters and the cellullar frequency of each cluster
    :params p_ccm, t_ccm: co-clustering matrix for the mutations
    :params p_ad, t_ad: ancestry-descendant matrix for the mutations
    :param weights: weights for scores from each challenge when combining to give an overall score

    '''
    s1A = calculate1A(p_cell, t_cell)
    s1B = calculate1B(p_ncluster, t_ncluster)
    s1C = calculate1C(p_1c, t_1c)
    s2 = calculate2(p_ccm, t_ccm)
    s3 = calculate3(p_ccm, p_ad, t_ccm, t_ad)
    scores = [s1A,s1B,s1C,s2,s3]

    soverall = np.sum(np.multiply(scores, weights))

    scores.append(soverall)
    if verbose:
        print('Scores:\nSC1 - Part A: %s, Part B: %s, Part C: %s\nSC2 - %s\nSC3 - %s\nOverall: %s' % tuple(scores))

    return soverall

def score_all(l_p_cell, t_cell,
              l_p_ncluster, t_ncluster,
              l_p_1c, t_1c,
              l_p_ccm, t_ccm,
              l_p_ad, t_ad,
              weights = [2000, 2000, 3000, 6000, 7000],
              verbose = False):
    """Calculate the overall score across all sub-challenges for multiple scenarios.

    :param l_p_cell: list of predicted cellularity for each scenario
    :param t_cell: true cellularity
    :param l_p_ncluster: list of predicted number of mutation clusters for each scenario
    :param t_ncluster:
    :param l_p_1c: list of zip objects containing the predicted number of mutations in each cluster and the predicted
        cellular frequency of each cluster. one zip object for each scenario.
    :param t_1c: zip object containing the true number of mutations in each cluster and the true cellular frequency
        of each cluster
    :param l_p_ccm: list of predicted co-clustering matrices for each scenario
    :param t_ccm: true co-clustering matrix
    :param l_p_ad: list of predicted ancestor-descendant matrices for each scenario
    :param t_ad: true ancestor-descendant matrix
    :param weights: weights for each part of each sub-challenge tp be used when calculating the overall score
    :param verbose: boolean for whether to print information about the function execution
    :return: list of overall scores for each of the given scenarios
    """
    scores = dict()
    for i in range(len(l_p_cell)):
        if verbose:
            print 'Scenario: %s' % scenarios[i]
        score = overall_score(
            l_p_cell[i], t_cell,
              l_p_ncluster[i], t_ncluster,
              l_p_1c[i], t_1c,
              l_p_ccm[i], t_ccm,
              l_p_ad[i], t_ad,
              weights,
              verbose)
        scores[scenarios[i]] = score

    if verbose:
        print ("Overall Scores")
        print scores

    return scores

def rank(scores):
    '''Calculate the rank order of a list of scores (either from a single sub-challenge or overall scores)

    :param scores: list of scores
    :return: list of relative ranking of each score
    '''
    scores = np.asarray(scores)
    order = scores.argsort()
    rank = order.argsort()
    return rank + 1

def get_ccm(scenario, t_ccm=None, t_clusters=None, size_clusters=100):
    '''Find the co-clustering matrix for the given scenario

    Attributes:
    :param scenario: string representing the clustering scenario being evaluated
    :params t_ccm, t_clusters: - optional value for the true co-clustering matrix and the true cluster assignments,
            to avoid computing it multiple times
    :param return: (if scenario is 'Truth') the true co-clustering matrix and the clusters used to generate this matrix
                    (otherwise) the co-clustering matrix fro the given scenario
    '''
    # TODO: make this more generalizable by calculating CCM from size and number of clusters

    if t_clusters is None:
        t_clusters = np.zeros((6*size_clusters,6))
        t_clusters[0:size_clusters,0] = 1 #cluster 1
        t_clusters[size_clusters:2*size_clusters,1] = 1 #cluster 2
        t_clusters[2*size_clusters:3*size_clusters,2] = 1 #...
        t_clusters[3*size_clusters:4*size_clusters,3] = 1
        t_clusters[4*size_clusters:5*size_clusters,4] = 1
        t_clusters[5*size_clusters:6*size_clusters,5] = 1
    if t_ccm is None:
        t_ccm = np.dot(t_clusters,t_clusters.T)

    if scenario in "Truth":
        return t_ccm, t_clusters
    elif "ParentIs" in scenario:
        return t_ccm
    elif scenario is "OneCluster":
        return np.ones(t_ccm.shape)
    elif "NCluster" in scenario:
        return np.identity(t_ccm.shape[0])
    elif scenario is "SplitClusterBotSame":
        clusters = np.zeros((6*size_clusters,7))
        clusters[:,:-1] = np.copy(t_clusters)
        clusters[5.5 * size_clusters:6*size_clusters,5] = 0
        clusters[5.5*size_clusters:6*size_clusters,6] = 1
        return np.dot(clusters, clusters.T)
    elif scenario is "SplitClusterBotDiff":
        clusters = np.zeros((6*size_clusters,7))
        clusters[:,:-1] = np.copy(t_clusters)
        clusters[5.5*size_clusters:6*size_clusters,5] = 0
        clusters[5.5*size_clusters:6*size_clusters,6] = 1
        return np.dot(clusters, clusters.T)
    elif scenario is "SplitClusterMidOneChild":
        clusters = np.zeros((6*size_clusters,7))
        clusters[:,:-1] = np.copy(t_clusters)
        clusters[2.5*size_clusters:3*size_clusters,2] = 0
        clusters[2.5*size_clusters:3*size_clusters,6] = 1
        return np.dot(clusters, clusters.T)
    elif scenario is "SplitClusterMidMultiChild":
        clusters = np.zeros((6*size_clusters,7))
        clusters[:,:-1] = np.copy(t_clusters)
        clusters[1.5*size_clusters:2*size_clusters,1] = 0
        clusters[1.5*size_clusters:2*size_clusters,6] = 1
        return np.dot(clusters, clusters.T)

    elif scenario is "MergeClusterBot":
        clusters = np.copy(t_clusters[:,:-1])
        clusters[5*size_clusters:6*size_clusters,4] = 1 #fix cluster 5 (originally cluster 6)
        clusters[4*size_clusters:5*size_clusters,4] = 0 #merge clusters 4 and 5 (from true phylogeny)
        clusters[3*size_clusters:4*size_clusters,3] = 1
        return np.dot(clusters, clusters.T)
    elif scenario == "MergeClusterMid&BotOneChild":
        clusters = np.copy(t_clusters[:,:-1])
        clusters[5*size_clusters:6*size_clusters, 2] = 1 #merge clusters 3 and 6
        return np.dot(clusters, clusters.T)
    elif scenario == "MergeClusterMid&BotMultiChild":
        clusters = np.copy(t_clusters[:,:-1])
        clusters[5*size_clusters:6*size_clusters, 4] = 1 #fix cluster 5 (originally cluster 6)
        clusters[4*size_clusters:5*size_clusters, 1] = 1 #merge clusters 2 and 5 (from true phylogeny)
        clusters[4*size_clusters:5*size_clusters, 4] = 0
        return np.dot(clusters, clusters.T)
    elif scenario == "MergeClusterTop&Mid":
        clusters = np.zeros((6*size_clusters,5))
        clusters[0:2*size_clusters, 0] = 1 #merged cluster from clusters 1 and 2 in true phylogeny
        clusters[2*size_clusters:3*size_clusters, 1] = 1
        clusters[3*size_clusters:4*size_clusters, 2] = 1
        clusters[4*size_clusters:5*size_clusters, 3] = 1
        clusters[5*size_clusters:6*size_clusters, 4] = 1
        return np.dot(clusters, clusters.T)
    elif "SmallExtra" in scenario:
        clusters = np.zeros((6*size_clusters,7))
        clusters[:,:-1] = np.copy(t_clusters)
        for i in range(0,6):
            clusters[size_clusters*i,i] = 0
            clusters[size_clusters*i,6] = 1
        return np.dot(clusters,clusters.T)
    elif "BigExtra" in scenario:
        clusters = np.zeros((6*size_clusters,7))
        clusters[:,:-1] = np.copy(t_clusters)
        for i in range(0,6):
            clusters[size_clusters*i:size_clusters*i+15,i] = 0
            clusters[size_clusters*i:size_clusters*i+15,6] = 1
        return np.dot(clusters,clusters.T)
    else:
        raise LookupError("Invalid scenario")

def get_ad(scenario, t_ad=None, size_clusters=100):
    '''Find the ancestry-descendant matrix for the given scenario

    Attributes:
    :param scenario: string representing the clustering scenario being evaluated
    :param t_ad: optional value for the true AD matrix, to avoid comupting it multiple times
    '''
    if t_ad is None:
        t_ad = np.zeros((6*size_clusters,6*size_clusters))
        t_ad[0:size_clusters, size_clusters:] = 1
        t_ad[size_clusters:2*size_clusters, 3*size_clusters:5*size_clusters] = 1
        t_ad[2*size_clusters:3*size_clusters, 5*size_clusters:6*size_clusters] = 1

    if scenario in ["Truth", "SplitClusterBotSame", "MergeClusterBot"]:
        return t_ad
    elif scenario is "SplitClusterBotDiff":
        ad = np.copy(t_ad)
        ad[5*size_clusters:5.5*size_clusters,5.5*size_clusters:6*size_clusters] = 1.
        return ad
    elif scenario is "SplitClusterMidOneChild":
        ad = np.copy(t_ad)
        ad[2.5*size_clusters:3*size_clusters,5*size_clusters:6*size_clusters] = 0
        return ad
    elif scenario is "SplitClusterMidMultiChild":
        ad = np.copy(t_ad)
        ad[1.5*size_clusters:2*size_clusters,3*size_clusters:5*size_clusters] = 0
        return ad
    elif scenario == "MergeClusterMid&BotOneChild":
        ad = np.copy(t_ad)
        ad[2*size_clusters:3*size_clusters, 5*size_clusters:6*size_clusters] = 0
        return ad
    elif scenario == "MergeClusterMid&BotMultiChild":
        ad = np.copy(t_ad)
        ad[size_clusters:2*size_clusters, 4*size_clusters:5*size_clusters] = 0
        ad[4*size_clusters:5*size_clusters, 3*size_clusters:4*size_clusters] = 1
        return ad
    elif scenario == "MergeClusterTop&Mid":
        ad = np.zeros((6*size_clusters,6*size_clusters))
        ad[0:2*size_clusters, 2*size_clusters:] = 1
        ad[2*size_clusters:3*size_clusters, 5*size_clusters:] = 1
        return ad
    elif scenario is "ParentIsSibling":
        ad = np.copy(t_ad)
        ad[3*size_clusters:4*size_clusters,4*size_clusters:5*size_clusters] = 1
        return ad
    elif scenario is "ParentIsGrandparent":
        ad = np.copy(t_ad)
        ad[size_clusters:2*size_clusters,4*size_clusters:5*size_clusters] = 0
        return ad
    elif scenario is "ParentIsAunt":
        ad = np.copy(t_ad)
        ad[size_clusters:2*size_clusters,4*size_clusters:5*size_clusters] = 0
        ad[2*size_clusters:3*size_clusters,4*size_clusters:5*size_clusters] = 1
        return ad
    elif scenario is "ParentIsCousin":
        ad = np.copy(t_ad)
        ad[size_clusters:2*size_clusters,4*size_clusters:5*size_clusters] = 0
        ad[2*size_clusters:3*size_clusters,4*size_clusters:5*size_clusters] = 1
        ad[5*size_clusters:6*size_clusters,4*size_clusters:5*size_clusters] = 1
        return ad
    elif scenario is "ParentIsSiblingWithChildren":
        ad = np.copy(t_ad)
        ad[2*size_clusters:3*size_clusters, range(size_clusters,2*size_clusters)+range(3*size_clusters,5*size_clusters)] = 1 #adjust cluster 3's ancestry
        return ad
    elif scenario is "ParentIsNieceWithChildren":
        ad = np.copy(t_ad)
        ad[2*size_clusters:3*size_clusters, range(size_clusters,2*size_clusters)+range(3*size_clusters,5*size_clusters)] = 1 #adjust cluster 3's ancestry
        ad[5*size_clusters:6*size_clusters, range(size_clusters,2*size_clusters)+range(3*size_clusters,5*size_clusters)] = 1 #adjust cluster 6's ancestry
        return ad
    elif scenario is "OneCluster":
        return np.zeros(t_ad.shape)
    elif scenario is "NClusterOneLineage":
        return np.triu(np.ones(t_ad.shape), k=1)
    elif scenario is "NClusterTwoLineages":
        ad = np.triu(np.ones(t_ad.shape), k=1)
        ad[2:3*size_clusters+2,3*size_clusters+2:] = 0
        return ad
    elif scenario is "NClusterCorrectLineage":
        ad = np.triu(np.ones(t_ad.shape), k=1)
        ad[size_clusters:2*size_clusters,range(2*size_clusters,3*size_clusters)+range(5*size_clusters,6*size_clusters)] = 0 # equivalent of cluster 2 from true AD matrix
        ad[2*size_clusters:3*size_clusters,3*size_clusters:5*size_clusters] = 0 # cluster 3 from true AD matrix
        ad[3*size_clusters:4*size_clusters,4*size_clusters:] = 0 # cluster 4 from true AD matrix
        ad[4*size_clusters:5*size_clusters,5*size_clusters:6*size_clusters] = 0 # cluster 5 from true AD matrix
        return ad
    elif scenario is "SmallExtraNewBot":
        ad = np.copy(t_ad)
        for i in range(0,6):
            ad[:,size_clusters*i] = 0
            ad[range(1,size_clusters)+range(size_clusters+1,2*size_clusters)+range(3*size_clusters+1,4*size_clusters),size_clusters*i] = 1
            ad[size_clusters*i,:] = 0
        return ad
    elif scenario is "SmallExtraCurBot":
        ad = np.copy(t_ad)
        for i in range(0,6):
            ad[:,size_clusters*i] = 0
            ad[range(1,size_clusters)+range(2*size_clusters+1,3*size_clusters),size_clusters*i] = 1
            ad[size_clusters*i,:] = 0
        return ad
    elif scenario is "SmallExtraMid":
        ad = np.copy(t_ad)
        for i in range(0,6):
            ad[:,size_clusters*i] = 0
            ad[range(1,size_clusters),size_clusters*i] = 1
            ad[size_clusters*i,:] = 0
        return ad
    elif scenario is "SmallExtraTop":
        ad = np.copy(t_ad)
        for i in range(0,6):
            ad[:,size_clusters*i] = 0
            ad[size_clusters*i,
               range(1,size_clusters)+
               range(size_clusters+1,2*size_clusters)+
               range(2*size_clusters+1,3*size_clusters)+
               range(3*size_clusters+1,4*size_clusters)+range(4*size_clusters+1,5*size_clusters)+range(5*size_clusters+1,6*size_clusters)] = 1
        return ad
    elif scenario is "BigExtraNewBot":
        ad = np.copy(t_ad)
        for i in range(0,6):
            ad[:,size_clusters*i:size_clusters*i+15] = 0
            ad[range(15,size_clusters)+
               range(size_clusters+15,2*size_clusters)+
               range(3*size_clusters+15,4*size_clusters),
            size_clusters*i:size_clusters*i+15] = 1
            ad[size_clusters*i:size_clusters*i+15,:] = 0
        return ad
    elif scenario is "BigExtraCurBot":
        ad = np.copy(t_ad)
        for i in range(0,6):
            ad[:,size_clusters*i:size_clusters*i+15] = 0
            ad[range(15,size_clusters)+range(2*size_clusters+15,3*size_clusters),size_clusters*i:size_clusters*i+15] = 1
            ad[size_clusters*i:size_clusters*i+15,:] = 0
        return ad
    elif scenario is "BigExtraMid":
        ad = np.copy(t_ad)
        for i in range(0,6):
            ad[:,size_clusters*i:size_clusters*i+15] = 0
            ad[range(15,size_clusters),(size_clusters*i):(size_clusters*i+15)] = 1
            ad[size_clusters*i:size_clusters*i+15,:] = 0
        return ad
    elif scenario is "BigExtraTop":
        ad = np.copy(t_ad)
        for i in range(0,6):
            ad[:,size_clusters*i:size_clusters*i+15] = 0
            ad[size_clusters*i:size_clusters*i+15,
            range(15,size_clusters)+
            range(1*size_clusters+15,2*size_clusters)+
            range(2*size_clusters+15,3*size_clusters)+range(3*size_clusters+15,4*size_clusters)+
            range(4*size_clusters+15,5*size_clusters)+range(5*size_clusters+15,6*size_clusters)] = 1
        return ad
    else:
        raise LookupError("Invalid scenario")

def get_cluster_size(scenario, t_size=100):
    '''Find the size of each cluster for the given scenario

    Attributes:
    :param scenario: string representing the clustering scenario being evaluated
    :param t_size: True size of each cluster. It is assumed the the true clusters are all the same size.
    '''
    t_n = 6 # true number of clusters
    if scenario is "Truth" or "ParentIs" in scenario:
        return np.repeat(t_size, t_n)
    elif scenario is "OneCluster":
        return np.array([t_n * t_size])
    elif "NCluster" in scenario:
        return np.repeat(1,t_n * t_size)
    elif "Split" in scenario:
        return np.concatenate((np.repeat(t_size,t_n-1),np.array([t_size / 2, t_size / 2])))
    elif "Merge" in scenario:
        return np.concatenate((np.array([2 * t_size]),np.repeat(t_size,t_n - 2)))
    elif "SmallExtra" in scenario:
        return np.concatenate((np.repeat(t_size - 1,t_n), np.array([t_n])))
    elif "BigExtra" in scenario:
        return np.concatenate((np.repeat(t_size - 15,t_n), np.array([15*t_n])))
    else:
        raise LookupError("Invalid scenario")

def get_cf(scenario, size_clusters):
    '''Find the cellular frequency for each cluster for the given scenario

    Attributes:
    :param scenario: string representing the clustering scenario being evaluated
    :param size_cluster: the size of eah cluster in the given scenario
    '''
    return size_clusters / float(sum(size_clusters)) # default is frequency relative to size of cluster


def get_cellularity(scenario):
    '''Find the cellularity for the given scenario

    Attributes:
    :param scenario: string representing the clustering scenario being evaluated
    '''
    t_cell = 0.7 # true cellularity
    if scenario is "Truth":
        return t_cell
    else:
        return t_cell + np.random.normal(scale=0.05)

# list of scenarios for testing the metric behaviour for sub-challenge 2 and 3 scores and the overall score
scenarios = ["Truth", "OneCluster", "NClusterOneLineage", "NClusterTwoLineages", "NClusterCorrectLineage",
                "ParentIsNieceWithChildren", "ParentIsSiblingWithChildren", "ParentIsCousin","ParentIsAunt", "ParentIsGrandparent", "ParentIsSibling",
                "BigExtraTop", "BigExtraMid", "BigExtraCurBot", "BigExtraNewBot",
                "SmallExtraTop", "SmallExtraMid", "SmallExtraNewBot", 'SmallExtraCurBot',
                "SplitClusterMidMultiChild", "SplitClusterMidOneChild", "SplitClusterBotSame", "SplitClusterBotDiff",
                "MergeClusterTop&Mid", "MergeClusterMid&BotMultiChild", "MergeClusterMid&BotOneChild","MergeClusterBot"]

def scoringtotal_behavior(verbose=False):
    '''Scoring behaviour of all subchallenge metrics together

    Attributes:
    :param verbose: boolean for whether to output details of the scoring metrics
    '''
    # True values for each attribute
    t_cell = get_cellularity("Truth")
    t_cluster_size = get_cluster_size("Truth")
    t_ncluster = len(t_cluster_size)
    t_cf = get_cf("Truth", t_cluster_size)
    t_1C = zip(t_cluster_size, t_cf)

    t_ccm = get_ccm("Truth")
    t_ad = get_ad("Truth")

    # list of predicted value for each scenario
    l_p_cell = list()
    l_p_ncluster = list()
    l_p_1C = list()
    l_p_ccm = list()
    l_p_ad = list()

    for scenario in scenarios:
        l_p_cell.append(get_cellularity(scenario))
        p_cluster_size = get_cluster_size(scenario)
        p_cf = get_cf(scenario, p_cluster_size)
        l_p_ncluster.append(len(p_cluster_size))
        l_p_1C.append(zip(p_cluster_size, p_cf))

        l_p_ccm.append(get_ccm(scenario))
        l_p_ad.append(get_ad(scenario))

    scores =  score_all(l_p_cell, t_cell,
              l_p_ncluster, t_ncluster,
              l_p_1C, t_1C,
              l_p_ccm, t_ccm,
              l_p_ad, t_ad, verbose=verbose)

    f = open(tsv_dir + "all_scores.csv", 'wb')
    wr = csv.writer(f)
    wr.writerow(scores.keys())
    wr.writerow(scores.values())
    f.close()

    return scores



if __name__ == '__main__':
    #scoring1C_behavior()

    methods = ("orig", "orig_nc", "pseudoV", "pseudoV_nc", "simpleKL_nc",
                 "sqrt_nc", "sym_pseudoV_nc", "pearson_nc", "aupr_nc", "mcc_nc",
                ["pseudoV_nc", "mcc_nc", "spearman_nc"], ["pseudoV_nc", "mcc_nc", "pearson_nc"])
    for m in methods:
        print('Calculating metric %s...' % m)
        scoring3A_behavior(method=m, verbose=True)
    '''
    scoring1A_behavior()
    scoring1B_behavior()
    scoring1C_behavior()
    scoring2A_behavior()
    scoring2B_behavior()
'''
