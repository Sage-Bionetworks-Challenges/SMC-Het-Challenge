from SMCScoring import *
import numpy as np

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
    
def scoring3A_behavior():
    # Case 2: six clusters, one at top level, two at 2nd level,
    # three at 3rd level (one lineage with one cluster, the other with two)
    # Phylogeny Tree:
    #       [1]
    #      |   |
    #     [2] [3]
    #    |  |   |
    #   [4][5] [6]

    # True co-cluster and ancestry matrices
    t_clusters = np.zeros((600,6))
    t_clusters[0:100,0] = 1 #cluster 1
    t_clusters[100:200,1] = 1 #cluster 2
    t_clusters[200:300,2] = 1 #...
    t_clusters[300:400,3] = 1
    t_clusters[400:500,4] = 1
    t_clusters[500:600,5] = 1
    t_ccm = np.dot(t_clusters,t_clusters.T)
    
    t_ad = np.zeros((600,600))
    t_ad[0:100, 100:] = 1
    t_ad[100:200, 300:500] = 1
    t_ad[200:300, 500:600] = 1
    
    # Cases:
    res = []  
    
    # Splits:
    res_split = []
    
    # Split third level (new cluster is on the third level)
    # Incorrect Phylogeny Tree:
    #        [1]
    #      |     |
    #     [2]   [3]
    #    |  |  |   |
    #   [4][5][6] [7]
    clusters = np.zeros((600,7))
    clusters[:,:-1] = np.copy(t_clusters)
    clusters[550:600,5] = 0
    clusters[550:600,6] = 1
    ccm = np.dot(clusters, clusters.T)
    res_split.append(["SplitClusterBotSame",calculate3(ccm,t_ad,t_ccm,t_ad)])
    
    # Split third level (new cluster is on new bottom level)
    # Incorrect Phylogeny Tree:
    #        [1]
    #      |    |
    #     [2]  [3]
    #    |  |   |   
    #   [4][5] [6]
    #           |    
    #          [7]
    clusters = np.zeros((600,7))
    clusters[:,:-1] = np.copy(t_clusters)
    clusters[550:600,5] = 0
    clusters[550:600,6] = 1
    ccm = np.dot(clusters, clusters.T)
    
    ad = np.copy(t_ad)
    ad[500:550,550:600] = 1
    res_split.append(["SplitClusterBotDiff",calculate3(ccm,ad,t_ccm,t_ad)])
    
    #Split middle level with one child (new cluster is on same level)
    # Incorrect Phylogeny Tree:
    #          [1]
    #      |    |    |
    #     [2]  [3]  [7]
    #    |  |   |   
    #   [4][5] [6]
    clusters = np.zeros((600,7))
    clusters[:,:-1] = np.copy(t_clusters)
    clusters[250:300,2] = 0
    clusters[250:300,6] = 1
    ccm = np.dot(clusters, clusters.T)
    
    ad = np.copy(t_ad)
    ad[250:300,500:600] = 0
    res_split.append(["SplitClusterMidOneChild",calculate3(ccm,ad,t_ccm,t_ad)])
    
    #Split middle level with two children (new cluster is on same level)
    # Incorrect Phylogeny Tree:
    #          [1]
    #      |    |    |
    #     [7]  [2]  [3]
    #         |  |   |   
    #        [4][5] [6]
    clusters = np.zeros((600,7))
    clusters[:,:-1] = np.copy(t_clusters)
    clusters[150:200,2] = 0
    clusters[150:200,6] = 1
    ccm = np.dot(clusters, clusters.T)
    
    ad = np.copy(t_ad)
    ad[150:200,300:500] = 0
    res_split.append(["SplitClusterMidMultiChild",calculate3(ccm,ad,t_ccm,t_ad)])
    
    res_split = [map(str,x) for x in res_split]
    res_split = ['\t'.join(x) for x in res_split]
    f = open('scoring_metric_data/scoring3A_split_cases.tsv', 'w')
    f.write('\n'.join(res_split))
    f.close()
    
    # Merges:
    res_merge = []
    
    # Merge third level (new cluster is on the third level)
    # Incorrect Phylogeny Tree (* = merged cluster):
    #        [1]
    #      |     |
    #     [2]   [3]
    #      |     |
    #     [4*]   [5] 
    clusters = np.copy(t_clusters[:,:-1])
    clusters[500:600,4] = 1 #fix cluster 5 (originally cluster 6)
    clusters[400:500,4] = 0 #merge clusters 4 and 5 (from true phylogeny)
    clusters[300:400,3] = 1
    ccm = np.dot(clusters, clusters.T)
    res_merge.append(["MergeClusterBot",calculate3(ccm,t_ad,t_ccm,t_ad)])
    
    # Merge second and third level with one child (new cluster is on the second level)
    # Incorrect Phylogeny Tree (* = merged cluster):
    #        [1]
    #      |     |
    #     [2]   [3*]
    #    |   |
    #   [4] [5] 
    clusters = np.copy(t_clusters[:,:-1])
    clusters[500:600, 2] = 1 #merge clusters 3 and 6
    ccm = np.dot(clusters, clusters.T)
    
    ad = np.copy(t_ad)
    ad[200:300, 500:600] = 0
    res_merge.append(["MergeClusterMid&BotOneChild",calculate3(ccm,ad,t_ccm,t_ad)])
    
    # Merge second and third level with two children (new cluster is on the second level)
    # Incorrect Phylogeny Tree (* = merged cluster):
    #        [1]
    #      |     |
    #     [2*]  [3]
    #      |     |
    #     [4]   [5] 
    clusters = np.copy(t_clusters[:,:-1])
    clusters[500:600, 4] = 1 #fix cluster 5 (originally cluster 6)
    clusters[400:500, 1] = 1 #merge clusters 2 and 5 (from true phylogeny)
    clusters[400:500, 4] = 0
    ccm = np.dot(clusters, clusters.T)
    
    ad = np.copy(t_ad)
    ad[100:200, 400:500] = 0
    ad[400:500, 300:400] = 1
    res_merge.append(["MergeClusterMid&BotMultiChild",calculate3(ccm,ad,t_ccm,t_ad)])
    
    # Merge first and second level with two children (new cluster is on the first level)
    # Incorrect Phylogeny Tree (* = merged cluster):
    #         [1*]
    #      |   |    |
    #     [3] [4]  [2]
    #               |
    #              [5] 
    clusters = np.zeros((600,5))
    clusters[0:200, 0] = 1 #merged cluster from clusters 1 and 2 in true phylogeny
    clusters[200:300, 1] = 1
    clusters[300:400, 2] = 1
    clusters[400:500, 3] = 1
    clusters[500:600, 4] = 1
    ccm = np.dot(clusters, clusters.T)
    
    ad = np.zeros((600,600))
    ad[0:200, 200:] = 1
    ad[200:300, 500:] = 1
    res_merge.append(["MergeClusterTop&Mid",calculate3(ccm,ad,t_ccm,t_ad)])
    
    res_merge = [map(str,x) for x in res_merge]
    res_merge = ['\t'.join(x) for x in res_merge]
    f = open('scoring_metric_data/scoring3A_merge_cases.tsv', 'w')
    f.write('\n'.join(res_merge))
    f.close()
    
    # Incorrect Parent:
    res_parent = []
    
    # Cluster's parent is really a sibling
    # Incorrect Phylogeny Tree (* = incorrect cluster):
    #        [1]
    #      |     |
    #     [2]   [3]
    #      |     |
    #     [4]   [6] 
    #      |    
    #     [5*]    
    ad = np.copy(t_ad)
    ad[300:400,400:500] = 1
    res_parent.append(["ParentIsSibling",calculate3(t_ccm,ad,t_ccm,t_ad)])
    
    # Cluster's parent is really a grandparent
    # Incorrect Phylogeny Tree (* = incorrect cluster):
    #           [1]
    #      |     |    |
    #     [2]   [3]  [5]
    #      |     |
    #     [4]   [6]   
    ad = np.copy(t_ad)
    ad[100:200,400:500] = 0
    res_parent.append(["ParentIsGrandparent",calculate3(t_ccm,ad,t_ccm,t_ad)])
    
    # Cluster's parent is really an aunt/uncle
    # Incorrect Phylogeny Tree (* = incorrect cluster):
    #        [1]
    #      |      |
    #     [2]    [3]
    #      |   |    |
    #     [4] [5*] [6]    
    ad[200:300,400:500] = 1
    res_parent.append(["ParentIsAunt",calculate3(t_ccm,ad,t_ccm,t_ad)])
    
    # Cluster's parent is really a cousin
    # Incorrect Phylogeny Tree (* = incorrect cluster):
    #        [1]
    #      |     |
    #     [2]   [3]
    #      |     |
    #     [4]   [6] 
    #            |    
    #           [5*]     
    ad[500:600,400:500] = 1
    res_parent.append(["ParentIsCousin",calculate3(t_ccm,ad,t_ccm,t_ad)])
    
    # NOTE: Keep ordering of cases the same in order for the ancestry matrix to be correct    
    
    # Cluster's parent is really a sibling, and incorrect cluster has children
    # Incorrect Phylogeny Tree (* = incorrect cluster):
    #        [1]
    #         |
    #        [3]
    #      |      |
    #     [2*]   [6] 
    #    |    |    
    #   [4]  [5]    
    ad = np.copy(t_ad)
    ad[200:300, range(100,200)+range(300,600)] = 1 #adjust cluster 3's ancestry
    res_parent.append(["ParentIsSiblingWithChildren",calculate3(t_ccm,ad,t_ccm,t_ad)])
    
    # Cluster's parent is really a niece, and incorrect cluster has children
    # Incorrect Phylogeny Tree (* = incorrect cluster):
    #        [1]
    #         |
    #        [3]
    #         |
    #        [6] 
    #         |
    #        [2*]
    #       |    |    
    #      [4]  [5]    
    ad[200:300, range(100,200)+range(300,600)] = 1 #adjust cluster 3's ancestry
    res_parent.append(["ParentIsNieceWithChildren",calculate3(t_ccm,ad,t_ccm,t_ad)])
    
    res_parent = [map(str,x) for x in res_parent]
    res_parent = ['\t'.join(x) for x in res_parent]
    f = open('scoring_metric_data/scoring3A_parent_cases.tsv', 'w')
    f.write('\n'.join(res_parent))
    f.close()
    
    # Other Issues
    res_other = []
    
    # OneCluster: everything is in one big cluster
    # Incorrect Phylogeny Tree:
    #    [1]
    res_other.append(["OneCluster", calculate3(np.ones(t_ccm.shape),np.zeros(t_ad.shape),t_ccm,t_ad)])
    
    # NClustersOneLineage: everything in it's own cluster, all in one long lineage
    # Incorrect Phylogeny Tree:
    #  [1]
    #   |
    #  [2]
    #  ...
    #   |
    # [600]
    res_other.append(["NClusterOneLineage", calculate3(np.identity(t_ccm.shape[0]),np.triu(np.ones(t_ad.shape), k=1),t_ccm,t_ad)])
    
    # NClustersTwoLineages: everything in it's own cluster, split into two germlines
    # Incorrect Phylogeny Tree:
    #      [1]
    #   |      |
    #  [2]   [302]
    #  ...    ...
    #   |      |
    # [300]  [600]
    #   |
    # [301]
    ad = np.triu(np.ones(t_ad.shape), k=1)
    ad[2:302,302:] = 0
    res_other.append(["NClusterTwoLineages", calculate3(np.identity(t_ccm.shape[0]),ad,t_ccm,t_ad)])
    
    # NClustersCorrectLineage: everything in it's own cluster, but with
    # the same structure as the true matrix
    # Incorrect Phylogeny Tree:
    #         [1]
    #         ...
    #          |
    #        [100] 
    #     |         |
    #   [101]     [201]
    #    ...       ...
    #     |         |
    #   [200]     [300]
    #  |      |     |
    #[301]  [401] [501]
    # ...    ...   ...
    #  |      |     |
    #[400]  [500] [600]
    ad = np.triu(np.ones(t_ad.shape), k=1)
    ad[100:200,range(200,300)+range(500,600)] = 0 # equivalent of cluster 2 from true AD matrix
    ad[200:300,300:500] = 0 # cluster 3 from true AD matrix
    ad[300:400,400:] = 0 # cluster 4 from true AD matrix
    ad[400:500,500:600] = 0 # cluster 5 from true AD matrix
    res_other.append(["NClusterCorrectLineage", calculate3(np.identity(t_ccm.shape[0]),ad,t_ccm,t_ad)])
    
    # SmallExtraNewBot: small extra cluster with some mutations from each cluster
    # new cluster is in a new bottom level of the lineage
    # Incorrect Phylogeny Tree (* = incorrect cluster):
    #        [1]
    #      |      |
    #     [2]    [3]
    #    |   |    |
    #   [4] [5]  [6]
    #    |
    #  [7*]
    clusters = np.zeros((600,7))
    clusters[:,:-1] = np.copy(t_clusters)
    for i in range(0,6):
        clusters[100*i,i] = 0
        clusters[100*i,6] = 1
    ccm = np.dot(clusters,clusters.T)
    
    ad = np.copy(t_ad)
    for i in range(0,6):
        ad[:,100*i] = 0
        ad[range(1,100)+range(101,200)+range(301,400),100*i] = 1
        ad[100*i,:] = 0
    res_other.append(["SmallExtraNewBot", calculate3(ccm,ad,t_ccm,t_ad)])
    
    # SmallExtraCurBot: small extra cluster with some mutations from each cluster
    # new cluster is in the current bottom level of the lineage
    # Incorrect Phylogeny Tree (* = incorrect cluster):
    #        [1]
    #      |       |
    #     [2]     [3]
    #    |   |   |   |
    #   [4] [5] [6] [7*]

    # NOTE: co-clustering matrix does not change for any of the different small
    # extra cluster instances, just the AD matrix    
    ad = np.copy(t_ad)
    for i in range(0,6):
        ad[:,100*i] = 0
        ad[range(1,100)+range(201,300),100*i] = 1
        ad[100*i,:] = 0
    res_other.append(["SmallExtraCurBot", calculate3(ccm,ad,t_ccm,t_ad)])
    
    # SmallExtraMid: small extra cluster with some mutations from each cluster
    # new cluster is in the middle level of the lineage
    # Incorrect Phylogeny Tree (* = incorrect cluster):
    #           [1]
    #      |     |   |
    #     [2]   [3] [7*]
    #    |   |   |
    #   [4] [5] [6]   
    ad = np.copy(t_ad)
    for i in range(0,6):
        ad[:,100*i] = 0
        ad[range(1,100),100*i] = 1
        ad[100*i,:] = 0
    res_other.append(["SmallExtraMid", calculate3(ccm,ad,t_ccm,t_ad)])
    
    # SmallExtraTop: small extra cluster with some mutations from each cluster
    # new cluster is in a new top level of the lineage
    # Incorrect Phylogeny Tree (* = incorrect cluster):
    #        [7*]
    #         |
    #        [1]
    #      |     |
    #     [2]   [3] 
    #    |   |   |
    #   [4] [5] [6]   
    ad = np.copy(t_ad)
    for i in range(0,6):
        ad[:,100*i] = 0
        ad[100*i,range(1,100)+range(101,200)+range(201,300)+range(301,400)+range(401,500)+range(501,600)] = 1
    res_other.append(["SmallExtraTop", calculate3(ccm,ad,t_ccm,t_ad)])
    
    # BigExtraNewBot: big extra cluster with some mutations from each cluster
    # new cluster is in a new bottom level of the lineage
    # Incorrect Phylogeny Tree (* = incorrect cluster):
    #        [1]
    #      |      |
    #     [2]    [3]
    #    |   |    |
    #   [4] [5]  [6]
    #    |
    #  [7*]
    clusters = np.zeros((600,7))
    clusters[:,:-1] = np.copy(t_clusters)
    for i in range(0,6):
        clusters[100*i:100*i+15,i] = 0
        clusters[100*i:100*i+15,6] = 1
    ccm = np.dot(clusters,clusters.T)
    
    ad = np.copy(t_ad)
    for i in range(0,6):
        ad[:,100*i:100*i+15] = 0
        ad[range(15,100)+range(115,200)+range(315,400),100*i:100*i+15] = 1
        ad[100*i:100*i+15,:] = 0
    res_other.append(["BigExtraNewBot", calculate3(ccm,ad,t_ccm,t_ad)])
    
    # BigExtraCurBot: small extra cluster with some mutations from each cluster
    # new cluster is in the current bottom level of the lineage
    # Incorrect Phylogeny Tree (* = incorrect cluster):
    #        [1]
    #      |       |
    #     [2]     [3]
    #    |   |   |   |
    #   [4] [5] [6] [7*]

    # NOTE: co-clustering matrix does not change for any of the different small
    # extra cluster instances, just the AD matrix    
    ad = np.copy(t_ad)
    for i in range(0,6):
        ad[:,100*i:100*i+15] = 0
        ad[range(15,100)+range(215,300),100*i:100*i+15] = 1
        ad[100*i:100*i+15,:] = 0
    res_other.append(["BigExtraCurBot", calculate3(ccm,ad,t_ccm,t_ad)])
    
    # BigExtraMid: small extra cluster with some mutations from each cluster
    # new cluster is in the middle level of the lineage
    # Incorrect Phylogeny Tree (* = incorrect cluster):
    #           [1]
    #      |     |   |
    #     [2]   [3] [7*]
    #    |   |   |
    #   [4] [5] [6]   
    ad = np.copy(t_ad)
    for i in range(0,6):
        ad[:,100*i:100*i+15] = 0
        ad[range(15,100),100*i:100*i+15] = 1
        ad[100*i:100*i+15,:] = 0
    res_other.append(["BigExtraMid", calculate3(ccm,ad,t_ccm,t_ad)])
    
    # BigExtraTop: small extra cluster with some mutations from each cluster
    # new cluster is in a new top level of the lineage
    # Incorrect Phylogeny Tree (* = incorrect cluster):
    #        [7*]
    #         |
    #        [1]
    #      |     |
    #     [2]   [3] 
    #    |   |   |
    #   [4] [5] [6]   
    ad = np.copy(t_ad)
    for i in range(0,6):
        ad[:,100*i:100*i+15] = 0
        ad[100*i:100*i+15,range(15,100)+range(115,200)+range(215,300)+range(315,400)+range(415,500)+range(515,600)] = 1
    res_other.append(["BigExtraTop", calculate3(ccm,ad,t_ccm,t_ad)])
    
    res_other = [map(str,x) for x in res_other]
    res_other = ['\t'.join(x) for x in res_other]
    f = open('scoring_metric_data/scoring3A_other_cases.tsv', 'w')
    f.write('\n'.join(res_other))
    f.close()
    
    res = res_split + res_merge + res_parent + res_other
    f = open('scoring_metric_data/scoring3A_all_cases.tsv', 'w')
    f.write('\n'.join(res))
    f.close()
    
    

if __name__ == '__main__':
    scoring1A_behavior()
    scoring1B_behavior()
    scoring1C_behavior()
    scoring2A_behavior()
    scoring2B_behavior()
