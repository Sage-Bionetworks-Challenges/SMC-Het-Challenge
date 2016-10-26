from scoring_harness_optimized import *
from SMCScoring import*# calculate1C
import numpy as np
import itertools
import os
import json
import sys

def test_metrics():
    print "Testing metrics"
    entry1 = np.matrix([[2, 1, 0], [0, 0, 2]])
    a = [True, False]
    b = ["Full", "Half"]
    ans1 = calculate_overlap_matrix(entry1)
    assert ans1 == (9, 0, 12, 4)
    print "     1. Information derived from overlap matrix is correct"

    orig1 = {"Full" : 0.80, "Half" : 0.80}
    sqrt1 = {"Full" : 0.89 , "Half" : 0.89}
    mcc1 = {"Full" : 0.72, "Half" : 0.61}
    spearman1 = {"Full" : 0.80, "Half" : 0.76}
    aupr1 = {"Full" : 0.93, "Half" : 0.85}
    pseudoV1 = {"Full" : 4.24, "Half" : None}


    for p, q in zip(a, b):
        assert round(om_calculate2_orig(*ans1, full_matrix=p), 2) == orig1[q]
        assert round(om_calculate2_sqrt(*ans1, full_matrix=p), 2) == sqrt1[q]
        assert round(om_calculate2_mcc(*ans1, full_matrix=p), 2) == mcc1[q]
        assert round(om_calculate2_spearman(*ans1, full_matrix=p), 2) == spearman1[q]
        assert round(om_calculate2_aupr(*ans1, full_matrix=p), 2) == aupr1[q]
    
    assert round(om_calculate2_pseudoV(entry1), 2) == pseudoV1["Full"]
    print "     2. Metrics working correctly"
    print "Finished testing metrics"

def test_scaling():
    print "Testing scaling function"
    
    scenario = ["OneCluster", "NCluster"]
    larger_is_worse = [True, False]

    entry1 = np.matrix([[2, 1, 0], [0, 0, 2]])
    orig1 = {"OneCluster" : 0.71, "NCluster" : 0.81}
    sqrt1 = {"OneCluster" : 0.85 , "NCluster" : 0.90}
    mcc1 = {"OneCluster" : 0.60, "NCluster" : 0.61}
    spearman1 = {"OneCluster" : 0.72, "NCluster" : 0.80}
    aupr1 = {"OneCluster" : 0.78, "NCluster" : 0.81}
    pseudoV1 = {"OneCluster" : 4.24, "NCluster" : None}

    for p in scenario:
        assert round(get_bad_score_om(entry1, om_calculate2_orig, scenario=p), 2) == orig1[p]
        assert round(get_bad_score_om(entry1, om_calculate2_sqrt, scenario=p), 2) == sqrt1[p]
        assert round(get_bad_score_om(entry1, om_calculate2_mcc, scenario=p), 2) == mcc1[p]
        assert round(get_bad_score_om(entry1, om_calculate2_spearman, scenario=p), 2) == spearman1[p]
        assert round(get_bad_score_om(entry1, om_calculate2_aupr, scenario=p), 2) == aupr1[p]
    print "     1. Score outputed by each scenario for each function are correct"

    assert round(get_worst_score_om(entry1, om_calculate2_orig, larger_is_worse=False), 2) == orig1["OneCluster"]
    assert round(get_worst_score_om(entry1, om_calculate2_sqrt, larger_is_worse=False), 2) == sqrt1["OneCluster"]
    assert round(get_worst_score_om(entry1, om_calculate2_mcc, larger_is_worse=False), 2) == mcc1["OneCluster"]
    assert round(get_worst_score_om(entry1, om_calculate2_spearman, larger_is_worse=False), 2) == spearman1["OneCluster"]
    assert round(get_worst_score_om(entry1, om_calculate2_aupr, larger_is_worse=False), 2) == aupr1["OneCluster"]
    print "     2. Scaling working correctly"
    print "Finished testing scaling"

def xstr(num):
    if num is None:
        return "None"
    else:
        return str(num)

def test_add_pseudo():
    test_cases = [None, 1, 2, 3, 4, 5]

    entry1 = (9, 0, 12, 4)
    entry2 = (4, 5, 9, 7)
    entry3 = (2, 8, 7, 8)

    ans1 = {"None": (11, 0, 34, 4), 
            "1": (10, 0, 22, 4), 
            "2": (11, 0, 34, 4), 
            "3": (12, 0, 48, 4),
            "4": (13, 0, 64, 4),
            "5": (14, 0, 82, 4)
            }

    ans2 = {"None": (6, 5, 31, 7), 
            "1": (5, 5, 19, 7), 
            "2": (6, 5, 31, 7), 
            "3": (7, 5, 45, 7),
            "4": (8, 5, 61, 7),
            "5": (9, 5, 79, 7)
            }

    ans3 = {"None": (4, 8, 29, 8), 
            "1": (3, 8, 17, 8), 
            "2": (4, 8, 29, 8), 
            "3": (5, 8, 43, 8),
            "4": (6, 8, 59, 8),
            "5": (7, 8, 77, 8)
            }

    for element in test_cases:
        assert add_pseudo_counts_om_eff(entry1[0], entry1[1], entry1[2], entry1[3], element) == ans1[xstr(element)]
        assert add_pseudo_counts_om_eff(entry2[0], entry2[1], entry2[2], entry2[3], element) == ans2[xstr(element)]
        assert add_pseudo_counts_om_eff(entry3[0], entry3[1], entry3[2], entry3[3], element) == ans3[xstr(element)]

    print "Add pseudo count functionality seems to be working correctly"

# cases are taken from metric_behavior.py. The answer should be correct....hopefully.... (They seem right)
def test_calculate1C():
    scores = {"abs": [0.88, 0.93, 0.83, 0.82, 0.98, 0.98, 0.98], 
              "sqr": [0.98, 0.99, 0.95, 0.94, 1.00, 1.00, 1.00]}

    method = ["abs", "sqr"]
    t_phis = np.array([.85,.5,.3])
    t_nssms = np.array([200,200,200])
    t_entry = zip(t_nssms,t_phis)

    for m in method:
        phis = [(.85+.5)/2.0, .3]
        nssms = [400,200]
        entry = zip(nssms,phis)
        assert round(calculate1C(entry, t_entry, m), 2) == scores[m][0] 
        
        phis = [.85, (.5+.3)/2.0]
        nssms = [200,400]
        entry = zip(nssms,phis)
        assert round(calculate1C(entry, t_entry, m), 2) == scores[m][1] 

        phis = [.55]
        nssms=[600]
        entry = zip(nssms,phis)
        assert round(calculate1C(entry, t_entry, m), 2) == scores[m][2]

        phis = [.85]
        nssms=[600]
        entry = zip(nssms,phis)
        assert round(calculate1C(entry, t_entry, m), 2) == scores[m][3]

        phis = [.9,.8, .5, .3]
        nssms = [100,100,200,200]
        entry = zip(nssms,phis)
        assert round(calculate1C(entry, t_entry, m), 2) == scores[m][4]

        phis = [.85, .55,.45, .3]
        nssms = [200,100,100,200]
        entry = zip(nssms,phis)
        assert round(calculate1C(entry, t_entry, m), 2) == scores[m][5]

        phis = [.85, .5,.35,.25]
        nssms = [200,200,100,100]
        entry = zip(nssms,phis)
        assert round(calculate1C(entry, t_entry, m), 2) == scores[m][6]

    print "Subchallenge 1C scoring functionalities seem to be working correctly"

def test_calculate2A():
    entry1 = np.matrix([[2, 1, 0], [0, 0, 2]])
    assert round(om_calculate2A(entry1), 2) == 0.53

    entry2 = np.matrix([[2, 0, 0],[0, 2, 1],[0, 0, 2]])
    assert round(om_calculate2A(entry2), 2) == 0.52

    entry3 = np.matrix([[1, 0, 0],[1, 1, 0],[0, 1, 4]])
    assert round(om_calculate2A(entry3), 2) == 0.31

    entry4 = np.matrix([[5, 3, 0], [0, 0, 0], [0, 0, 0]])
    assert round(om_calculate2A(entry4), 2) == 0.35

    entry5 = np.matrix([[1, 1, 1], [1, 1, 1], [1, 1, 1]])
    assert round(om_calculate2A(entry5), 2) == 0

    print "Subchallenge 2A scoring functionalities seem to be working correctly"

def test_calculate3A():
    entry1_om = np.matrix([[2, 1, 0], [0, 0, 2]])
    entry1_ad_truth = np.matrix([[0, 1], [0, 0]])
    entry1_ad_pred = np.matrix([[0, 1, 1], [0, 0, 0], [0, 0, 0]])
    entry1_truth_data = [1, 1, 1, 2, 2]

    entry2_om = np.matrix([[2, 0, 0],[0, 2, 1],[0, 0, 2]])
    entry2_ad_pred = np.matrix([[0, 1, 1], [0, 0, 1], [0, 0, 0]])
    entry2_ad_truth = np.matrix([[0, 1, 1],[0, 0, 0],[0, 0, 0]])
    entry2_truth_data = [1, 1, 2, 2, 2, 3, 3]

    entry3_om = np.matrix([[1, 1, 0], [1, 1, 0], [0, 0, 4]])
    entry3_ad_pred = np.matrix([[0, 1, 1], [0, 0, 0], [0, 0, 0]])
    entry3_ad_truth = np.matrix([[0, 1, 1], [0, 0, 0], [0, 0, 0]])
    entry3_truth_data = [2, 1, 2, 1, 3, 3, 3, 3]

    assert calculate3A(entry1_om, entry1_truth_data, entry1_ad_pred, entry1_ad_truth) == 0
    assert round(calculate3A(entry2_om, entry2_truth_data, entry2_ad_pred, entry2_ad_truth), 2) == 0.33
    assert round(calculate3A(entry3_om, entry3_truth_data, entry3_ad_pred, entry3_ad_truth), 2) == 0.06

    print "Subchallenge 3A scoring functionality seems to be working correctly"

def test_validate2A():
    predfile2 = ['complex/pred_2A_3A_complex.txt', 
                'simple/pred_2A_3A_simple.txt',
                'simple/truth_2A_3A_simple.txt',
                'different/pred_2A_different_cluster.txt',
                'no_fp/predfile_2A_no_fp.txt',
                'simple/predfile_2A_simple.txt'
                ]
    truthfile2 = ['complex/truth_2A_3A_complex.txt',
                'simple/truth_2A_3A_simple.txt',
                'simple/truth_2A_3A_simple.txt',
                'different/truth_2A_different_cluster.txt',
                'no_fp/truthfile_2A_no_fp.txt',
                'simple/truthfile_2A_simple.txt'
                ]
    nssms = [6, 5, 5, 5, 5, 7]            
    complex_om = np.matrix([[1, 0, 1, 0], [0, 2, 0, 0], [0, 0, 0, 1], [0, 0, 0, 1]])
    simple_om_1 = np.matrix([[1, 0, 0], [0, 2, 1], [1, 0, 0]])
    simple_om_2 = np.matrix([[1, 0, 0], [0, 3, 0], [0, 0, 1]])
    different_om = np.matrix([[1, 0, 1, 0], [0, 2, 0, 0], [0, 0, 0, 1]])
    no_fp_om = np.matrix([[2, 1, 0], [0, 0, 2]])
    simple_om_3 = np.matrix([[2, 1, 0, 0],[0, 0, 2, 0],[0, 0, 0, 1], [0, 0, 0, 1]])
    om_matrics = [complex_om, simple_om_1, simple_om_2, different_om, no_fp_om, simple_om_3]

    for element in range(5):
        np.testing.assert_array_equal(
            om_validate2A(readFile(predfile2[element]), readFile(truthfile2[element]), nssms[element], nssms[element]), 
            om_matrics[element]
            )
    print "Subchallenge 2A validation seems to be working correctly"

def test_validate3A():
    file3 = ['complex/pred_3A_complex.txt', 'simple/pred_3A_simple.txt', 'complex/truth_3A_complex.txt', 'simple/truth_3A_simple.txt']
    simple_truth = np.matrix(([[0, 1, 1], [0, 0, 0], [0, 0, 0]]))
    simple_pred = np.matrix([[0, 1, 1], [0, 0, 1], [0, 0, 0]])
    complex_truth = np.matrix([[0, 1, 1, 1], [0, 0, 0, 1], [0, 0, 0, 0], [0, 0, 0, 0]])
    complex_pred = np.matrix([[0, 1, 1, 1], [0, 0, 1, 1], [0, 0, 0, 0], [0, 0, 0, 0]])

    ad_matrices = [complex_pred, simple_pred, complex_truth, simple_truth]

    for element in range(4):
        np.testing.assert_array_equal(om_validate3A(readFile(file3[element]), ad_matrices[element].shape[0]), ad_matrices[element])

    print "Subchallenge 3A validation seems to be working correctly"

def test_construct_srm():
    simple_truth = np.matrix([[0, 1, 1], [0, 0, 0], [0, 0, 0]])
    simple_pred = np.matrix([[0, 1, 1], [0, 0, 1], [0, 0, 0]])
    simple_om = np.matrix([[1, 0, 0], [0, 2, 1], [1, 0, 0]])
    relative_matrix_1 = np.matrix([[0, 0, 0], [0, 1, 1], [0, 1, 1]])

    np.testing.assert_array_equal(construct_relative_matrix(simple_om, simple_pred, simple_truth), relative_matrix_1)

    complex_pred = np.matrix([[0, 1, 1, 1], [0, 0, 1, 1], [0, 0, 0, 0], [0, 0, 0, 0]])
    complex_truth = np.matrix([[0, 1, 1, 1], [0, 0, 0, 1], [0, 0, 0, 0], [0, 0, 0, 0]])    
    complex_om = np.matrix([[1, 0, 1, 0], [0, 2, 0, 0], [0, 0, 0, 1], [0, 0, 0, 1]])
    relative_matrix_2 = np.matrix([[0, 0, 0, 0],[0, 1, 1, 1],[0, 1, 1, 1], [0, 1, 3, 3]])

    np.testing.assert_array_equal(construct_relative_matrix(complex_om, complex_pred, complex_truth), relative_matrix_2)

    cousin_om = np.matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 1, 1, 0], [0, 0, 0, 3]])
    cousin_pred = np.matrix([[0, 1, 1, 1], [0, 0, 1, 1], [0, 0, 0, 0], [0, 0, 0, 0]])
    cousin_truth = np.matrix([[0, 1, 1, 1], [0, 0, 0, 1], [0, 0, 0, 0], [0, 0, 0, 0]])
    relative_matrix_3 = np.matrix([[0, 0, 0, 0], [0, 1, 1, 1], [0, 1, 1, 1],[0, 1, 2, 2]])

    np.testing.assert_array_equal(construct_relative_matrix(cousin_om, cousin_pred, cousin_truth), relative_matrix_3)

    print "Shared relative matrix seems to be constructed properly"

def test_calculate3_pseudoV():
    om1 = np.matrix([[2, 0, 0], [0, 2, 1], [0, 0, 2]])
    srm1 = np.matrix([[0, 0, 0], [0, 2, 2], [0, 2, 2]])
    ad_pred1 = np.matrix([[0, 1, 1], [0, 0, 1], [0, 0, 0]])
    ad_truth1 = np.matrix([[0, 1, 1], [0, 0, 0], [0, 0, 0]])
    pdm1 = np.matrix([[0], [2], [4]])
    tdm1 = np.matrix([[0], [2], [2]])

    np.testing.assert_array_equal(construct_related_mutations_matrix(om1, ad_pred1, ad_truth1)[0], pdm1)
    np.testing.assert_array_equal(construct_related_mutations_matrix(om1, ad_pred1, ad_truth1)[1], tdm1)

    assert round(calculate3_pseudoV(srm1, om1, pdm1, tdm1), 2) == 6.72
    assert round(calculate3A_pseudoV_final(om1, ad_pred1, ad_truth1), 2) == 5.14
    assert round(calculate3A_pseudoV_final(om1, ad_pred1, ad_truth1, modification="transpose"), 2) == 6.72
    assert round(calculate3A_pseudoV_final(om1, ad_pred1, ad_truth1, modification="cousin"), 2) == 14.67

    om2 = np.matrix([[2, 0, 0], [0, 2, 1], [0, 0, 2]])
    ad_pred2 = np.zeros((3, 3))
    ad_truth2 = np.matrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]])
    pdm2 = np.zeros((3,1))
    tdm2 = np.matrix([[0], [2], [3]])

    np.testing.assert_array_equal(construct_related_mutations_matrix(om2, ad_pred2, ad_truth2)[0], pdm2)
    np.testing.assert_array_equal(construct_related_mutations_matrix(om2, ad_pred2, ad_truth2)[1], tdm2)

    print "The method to calculate the sym-pseudo-V-measure and its helper functions seem to be working correctly"

def readFile(txt):
    f= open(txt)
    data = f.read()
    f.close()
    return data

if __name__ == "__main__":
    os.chdir('./test_data')
    test_metrics()
    test_scaling()
    test_add_pseudo()
    test_construct_srm()
    test_calculate3_pseudoV()
    test_validate2A()
    test_validate3A()
    #test_calculate1C()
    test_calculate2A()   
    test_calculate3A()
