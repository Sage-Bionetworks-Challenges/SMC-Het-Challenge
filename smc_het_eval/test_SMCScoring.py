import pytest
from SMCScoring import *
import numpy as np
import itertools
import os
import json
import sys

def test_calculate1A():
	# No difference gives perfect score
	assert np.testing.assert_allclose(calculate1A(.8,.8),1.0) == None
	# Metric is symmetric
	assert np.testing.assert_allclose(calculate1A(.5,.8),0.7) == None
	assert np.testing.assert_allclose(calculate1A(.8,.5),0.7) == None
	# Worst possible score
	assert np.testing.assert_allclose(calculate1A(0.0,1.0),0.0) == None

def test_validate1A():
	# Empty files
	with pytest.raises(ValidationError) as e:
		validate1A("")
	assert 'contains zero' in str(e.value)

	with pytest.raises(ValidationError) as e:
		validate1A("\n")
	assert 'contains zero' in str(e.value)

	# Multiple values
	with pytest.raises(ValidationError) as e:
		validate1A("0.8\n0.4\n")
	assert 'more than' in str(e.value)

	with pytest.raises(ValidationError) as e:
		validate1A("0.8\n0.4")
	assert 'more than' in str(e.value)	

	# Can't be cast as a float
	with pytest.raises(ValidationError) as e:
		validate1A("a0.8")
	assert 'could not' in str(e.value)	

	# NaN / Infinite value
	with pytest.raises(ValidationError) as e:
		validate1A("NaN")
	assert 'NaN' in str(e.value)	

	with pytest.raises(ValidationError) as e:
		validate1A("-Inf")
	assert 'finite' in str(e.value)	

	with pytest.raises(ValidationError) as e:
		validate1A("Inf")
	assert 'finite' in str(e.value)	

	# Negative cellularities
	with pytest.raises(ValidationError) as e:
		validate1A("-1.2")
	assert 'was <' in str(e.value)	

	# Cellularity exceeding 100%
	with pytest.raises(ValidationError) as e:
		validate1A("1.2")
	assert 'was >' in str(e.value)


	assert validate1A("0.8") == 0.8
	assert validate1A("0.8\n") == 0.8
	# Test DOS line endings
	assert validate1A("0.8\r\n") == 0.8
	assert validate1A(".3") == 0.3
	# Test extremes of valid range
	assert validate1A("0") == 0.0
	assert validate1A("1") == 1.0

def test_calculate1B():
	# No difference gives zero error
	assert np.testing.assert_allclose(calculate1B(3,3),1.0) == None
	# Arbitrary examples
	assert np.testing.assert_allclose(calculate1B(2,3),0.75) == None
	assert np.testing.assert_allclose(calculate1B(3,2),2.0/3.0) == None
	assert np.testing.assert_allclose(calculate1B(5,4),0.8) == None

def test_validate1B():
	# Empty files
	with pytest.raises(ValidationError) as e:
		validate1B("")
	assert 'contains zero' in str(e.value)

	with pytest.raises(ValidationError) as e:
		validate1B("\n")
	assert 'contains zero' in str(e.value)

	#Multiple entries
	with pytest.raises(ValidationError) as e:
		validate1B("4\n2\n")
	assert 'more than' in str(e.value)

	with pytest.raises(ValidationError) as e:
		validate1B("4\n2")
	assert 'more than' in str(e.value)	

	#Things that can't be cast as integer
	with pytest.raises(ValidationError) as e:
		validate1B("a5")
	assert 'could not' in str(e.value)	
	# Floating point values
	with pytest.raises(ValidationError) as e:
		validate1B("3.0")
	assert 'could not' in str(e.value)	

	with pytest.raises(ValidationError) as e:
		validate1B("NaN")
	assert 'could not' in str(e.value)	

	with pytest.raises(ValidationError) as e:
		validate1B("-Inf")
	assert 'could not' in str(e.value)	

	with pytest.raises(ValidationError) as e:
		validate1B("Inf")
	assert 'could not' in str(e.value)	

	# Invalid number of clusters
	with pytest.raises(ValidationError) as e:
		validate1B("0")
	assert 'less than' in str(e.value)	

	with pytest.raises(ValidationError) as e:
		validate1B("-1")
	assert 'less than' in str(e.value)	

	with pytest.raises(ValidationError) as e:
		validate1B("10000")
	assert 'greater than' in str(e.value)	


	assert validate1B("4") == 4
	assert validate1B("2\n") == 2
	# Test DOS line endings
	assert validate1B("2\r\n") == 2

def test_calculate1C():
	entry1 = [(100,.9),(200,.5),(100,.1)]
	entry2 = [(200,.5),(100,.9),(100,.1)]
	entry3 = [(200,.45),(100,.95),(100,.09)]
	entry4 = [(210,.5),(95,.9),(95,.1)]
	entry5 = [(400,.5)]
	entry6 = [(100,.1),(200,.5),(100,.9)]
	#Identical
	assert np.testing.assert_allclose(calculate1C(entry1,entry1),1.0) == None
	assert np.testing.assert_allclose(calculate1C(entry5,entry5),1.0) == None
	#Identical but different order
	assert np.testing.assert_allclose(calculate1C(entry2,entry1),1.0) == None
	assert np.testing.assert_allclose(calculate1C(entry6,entry1),1.0) == None
	#Slightly different phis
	assert np.testing.assert_allclose(calculate1C(entry3,entry1),0.96) == None
	#Slightly different nssms
	assert np.testing.assert_allclose(calculate1C(entry4,entry1),0.99) == None
	#Totally different
	assert np.testing.assert_allclose(calculate1C(entry5,entry1),0.8) == None

def test_validate1C():
	# Empty files
	with pytest.raises(ValidationError) as e:
		validate1C("",1)
	assert 'Number of lines' in str(e.value)

	with pytest.raises(ValidationError) as e:
		validate1C("\n",1)
	assert 'Number of lines' in str(e.value)

	# No SC1B - Valid
	assert validate1C("1\t3\t.9\n2\t2\t.5",5) == [(3,.9),(2,.5)]
	# No SC1B - too many lines:
	with pytest.raises(ValidationError) as e:
		validate1C("1\t1\t0.9\n2\t1\t.5\n3\t1\t.5\n4\t1\t.5\n5\t1\t.5\n6\t1\t.5\n7\t1\t.5\n8\t1\t.5\n9\t1\t.5\n10\t1\t.5\n11\t1\t.5\n",11)
	assert 'Number of lines' in str(e.value)	
	# No SC1B - Empty
	with pytest.raises(ValidationError) as e:
		validate1C("",1)
	assert 'Number of lines' in str(e.value)

	# too many ssms
	with pytest.raises(ValidationError) as e:
		validate1C("1\t3\t0.9\n2\t2\t.5",4)
	assert 'Total number' in str(e.value)
	# too few ssms
	with pytest.raises(ValidationError) as e:
		validate1C("1\t3\t0.9\n2\t2\t.5",6)
	assert 'Total number' in str(e.value)
	# wrong number of columns
	with pytest.raises(ValidationError) as e:
		validate1C("1\t30.9\n2\t2\t.5",5)
	assert 'tab separated columns in line 1' in str(e.value)

	with pytest.raises(ValidationError) as e:
		validate1C("3\t0.9\n2\t.5",5)
	assert 'tab separated columns in line 1' in str(e.value)
	# No cluster ids
	with pytest.raises(ValidationError) as e:
		validate1C("1.0\t30.9\n2\t2\t.5",5)
	assert 'tab separated columns' in str(e.value)
	# bad cluster ids
	with pytest.raises(ValidationError) as e:
		validate1C("1.0\t3\t0.9\n2\t2\t.5",5)
	assert 'Cluster ID in line 1' in str(e.value)

	with pytest.raises(ValidationError) as e:
		validate1C("a\t3\t0.9\nb\t2\t.5",5)
	assert 'Cluster ID in line 1' in str(e.value)

	#zero-indexed clusters
	with pytest.raises(ValidationError) as e:
		validate1C("0\t3\t0.9\n1\t2\t.5",5)
	assert 'Cluster ID in line 1 is not' in str(e.value)
	
	# negative ssms
	with pytest.raises(ValidationError) as e:
		validate1C("1\t3\t0.9\n2\t-2\t.5",5)
	assert 'Number of mutations in line 2' in str(e.value)

	#zero ssms
	with pytest.raises(ValidationError) as e:
		validate1C("1\t0\t0.9\n2\t2\t.5",6)
	assert 'Number of mutations in line 1' in str(e.value)
	#string nssms
	with pytest.raises(ValidationError) as e:
		validate1C("1\tone\t0.9\n2\t2\t.5",6)
	assert 'Number of mutations in line 1' in str(e.value)
	#float nssms
	with pytest.raises(ValidationError) as e:
		validate1C("1\t3\t0.9\n2\t2.0\t.5",6)
	assert 'Number of mutations in line 2' in str(e.value)
	# negative CF
	with pytest.raises(ValidationError) as e:
		validate1C("1\t3\t0.9\n2\t2\t-.5",6)
	assert 'Cellular Frequency for cluster 2' in str(e.value)
	#infinite CF
	with pytest.raises(ValidationError) as e:
		validate1C("1\t3\tInf\n2\t2\t.5",6)
	assert 'Cellular Frequency for cluster 1' in str(e.value)
	#NaN CF
	with pytest.raises(ValidationError) as e:
		validate1C("1\t3\t.9\n2\t2\tNaN",6)
	assert 'Cellular Frequency for cluster 2' in str(e.value)
	# >1 CCF
	with pytest.raises(ValidationError) as e:
		validate1C("1\t3\t.9\n2\t2\t1.5",6)
	assert 'Cellular Frequency for cluster 2' in str(e.value)
	#not a float CF
	with pytest.raises(ValidationError) as e:
		validate1C("1\t3\t.9\n2\t2\tabc",6)
	assert 'Cellular Frequency for cluster 2 can not' in str(e.value)

	#Valid entries
	assert validate1C("1\t3\t.9\n2\t2\t.5",5) == [(3,.9),(2,.5)]
	assert validate1C("1\t3\t.9\n",3) == [(3,.9)]

	assert validate1C("1\t3\t.9\n2\t2\t.5\n",5) == [(3,.9),(2,.5)]
	#DOS line endings
	assert validate1C("1\t3\t.9\r\n2\t2\t.5",5) == [(3,.9),(2,.5)]

def test_score1A():
	try:
		os.remove('nosuchfile.txt')
	except OSError:
		pass
	s = ['valid.VCF','invalid.VCF','nosuchfile.txt']
	a = [['valid1A.txt'],['invalid1A.txt'],['nosuchfile.txt']]
	for p in itertools.product(a,a,s):
		params = list(p)
		res = scoreChallenge('1A',*params)
		if params[0] == a[0] and params[1] == a[0]:
			assert res == 1.0
		else:
			assert res == "NA"


def test_score1B():
	try:
		os.remove('nosuchfile.txt')
	except OSError:
		pass
	s = ['valid.VCF','invalid.VCF','nosuchfile.txt']
	b = [['valid1B.txt'],['invalid1B.txt'],['nosuchfile.txt']]
	for p in itertools.product(b,b,s):
		params = list(p)
		res = scoreChallenge('1B',*params)
		if params[0] == b[0] and params[1] == b[0]:
			assert res == 1.0
		else:
			assert res == "NA"		

def test_score1C():
	try:
		os.remove('nosuchfile.txt')
	except OSError:
		pass
	s = ['valid.VCF','invalid.VCF','nosuchfile.txt']
	c = [['valid1C.txt'],['invalid1C.txt'],['nosuchfile.txt']]
	for p in itertools.product(c,c,s):
		params = list(p)
		res = scoreChallenge('1C',*params)
		if params[0] == c[0] and params[1] == c[0] and params[2] == s[0]:
			assert res == 1.0
		else:
			assert res == "NA"	

def test_verify1A():
	try:
		os.remove('nosuchfile.txt')
	except OSError:
		pass

	s = ['valid.VCF','invalid.VCF','nosuchfile.txt']
	a = [['valid1A.txt'],['invalid1A.txt'],['nosuchfile.txt']]
	for p in itertools.product(a,s):
		params = list(p)
		res = verifyChallenge('1A',*params)
		if params[0] == a[0]:
			assert res == "Valid"
		else:
			assert res == "Invalid"

def test_verify1B():
	try:
		os.remove('nosuchfile.txt')
	except OSError:
		pass
	s = ['valid.VCF','invalid.VCF','nosuchfile.txt']
	b = [['valid1B.txt'],['invalid1B.txt'],['nosuchfile.txt']]
	for p in itertools.product(b,s):
		params = list(p)
		res = verifyChallenge('1B',*params)
		if params[0] == b[0]:
			assert res == "Valid"
		else:
			assert res == "Invalid"

def test_verify1C():
	try:
		os.remove('nosuchfile.txt')
	except OSError:
		pass
	s = ['valid.VCF','invalid.VCF','nosuchfile.txt']
	c = [['valid1C.txt'],['invalid1C.txt'],['nosuchfile.txt']]
	for p in itertools.product(c,s):
		params = list(p)
		res = verifyChallenge('1C',*params)
		if params[1] != s[0]:
			assert res == "NA"
		elif params[0] == c[0]:
			assert res == "Valid"
		else:
			assert res == "Invalid"

def test_calculate2():
	# Two mutations in one cluster, two in second
	c = np.zeros((4,2))
	c[0:2,0] = 1
	c[2:4,1] = 1
	c = np.dot(c,c.T)
	
	# Identical
	assert calculate2(c,c) == 1.0
	
	#Inverted
	c2 = np.abs(c-1)
	c2[np.diag_indices(4)] = 1
	assert calculate2(c,c2) == 0.0

	# Differences first 3 SSMs in first cluster, 4th ssm in second cluster
	c3 = np.zeros((4,2))
	c3[0:3,0] = 1
	c3[3:4,1] = 1
	c3 = np.dot(c3,c3.T)
	assert calculate2(c,c3) == 0.5

	# Metric doesn't count the diagnonal
	c4 = c+0
	c4[np.diag_indices(4)] = 0
	assert calculate2(c,c4) == 1.0

	c2[np.diag_indices(4)] = 0
	assert calculate2(c,c2) == 0.0
	

def test_validate2A():
	entry = "1\n2\n1\n2\n"
	# Empty file
	with pytest.raises(ValidationError) as e:
		validate2A("",4)
	assert 'contains a different number' in str(e.value)
	# Wrong number of ssms
	entry1 = "1\n2\n2\n"
	with pytest.raises(ValidationError) as e:
		validate2A(entry1,4)
	assert 'contains a different number' in str(e.value)	

	entry2 = "0\n1\n1\n0\n"
	# Zero indexing
	with pytest.raises(ValidationError) as e:
		validate2A(entry2,4)
	assert 'Cluster IDs used' in str(e.value)	

	# Cluster names are not integer
	entry5 = "1.0\n2.0\n2.0\n1.0\n"
	with pytest.raises(ValidationError) as e:
		validate2A(entry5,4)
	assert 'Cluster ID in line 1' in str(e.value)

	# Correct entry
	c = np.zeros((4,2))
	c[[0,2],0] = 1
	c[[1,3],1] = 1
	ccm = np.dot(c,c.T)
	assert np.testing.assert_allclose(ccm,validate2A(entry,4)) == None

	# All one cluster
	one_cluster_entry = "1\n1\n1\n1\n"
	res = validate2A(one_cluster_entry,4)
	assert np.testing.assert_allclose(res,np.ones((4,4))) == None

	# Each SSM it's own cluster
	four_cluster_entry = "1\n2\n3\n4\n"
	res = validate2A(four_cluster_entry,4)
	assert np.testing.assert_allclose(res,np.identity(4)) == None

	# Each SSM it's own cluster reverse order
	four_cluster_entry2 = "4\n3\n2\n1\n"
	res = validate2A(four_cluster_entry2,4)
	assert np.testing.assert_allclose(res,np.identity(4)) == None

	#Ensure cluster assignment version works:
	res = validate2A(entry,4,False)
	assert np.testing.assert_allclose(res,c) == None

def test_validate2B():
	ssmlist = 4
	# Empty file
	with pytest.raises(ValidationError) as e:
		validate2B("",ssmlist)
	assert 'Shape of' in str(e.value)	

	# Wrong size: Vector
	vector_entry = '\t'.join(map(str,[0.5]*16))
	with pytest.raises(ValidationError) as e:
		validate2B(vector_entry,ssmlist)
	assert 'Shape of' in str(e.value)	

	# Wrong size: 3x3
	three_by_three = "1\t.5\t.5\n.5\t1\t.5\n.5\t.5\t1\n"
	with pytest.raises(ValidationError) as e:
		validate2B(three_by_three,ssmlist)
	assert 'Shape of' in str(e.value)

	# Not symmetric
	non_sym = "1\t1\t1\t1\n0\t1\t1\t1\n1\t1\t1\t1\n1\t1\t1\t1\n"
	with pytest.raises(ValidationError) as e:
		validate2B(non_sym,ssmlist)
	assert 'is not symmetric' in str(e.value)

	# One diagonal entries not 1
	diag_ne_1 = "1\t1\t1\t1\n1\t0\t1\t1\n1\t1\t1\t1\n1\t1\t1\t1\n"
	with pytest.raises(ValidationError) as e:
		validate2B(diag_ne_1,ssmlist)
	assert 'Diagonal entries' in str(e.value)

	# No diagonal entries 1
	diag_a_ne_1 = "0\t1\t1\t1\n1\t0\t1\t1\n1\t1\t0\t1\n1\t1\t1\t0\n"
	with pytest.raises(ValidationError) as e:
		validate2B(diag_a_ne_1,ssmlist)
	assert 'Diagonal entries' in str(e.value)

	# Non-floats
	non_float = "1\t1\t1\t1\n1\t1\t1\t1\n1\ta\t1\t1\n1\t1\t1\t1\n"
	with pytest.raises(ValidationError) as e:
		validate2B(non_float,ssmlist)
	assert 'could not be cast' in str(e.value)
	# NaNs
	nans = "1\t1\t1\t1\n1\t1\t1\t1\n1\t1\t1\t1\n1\t1\tNaN\t1\n"
	with pytest.raises(ValidationError) as e:
		validate2B(nans,ssmlist)
	assert 'contains NaNs' in str(e.value)
	# Infinite entries
	inf = "1\t1\t1\t1\n-Inf\t1\t1\t1\n1\t1\t1\t1\n1\t1\t1\t1\n"
	with pytest.raises(ValidationError) as e:
		validate2B(inf,ssmlist)
	assert 'contains non-finite' in str(e.value)
	# > 1
	ge_1 = "1\t1\t1.2\t1\n1\t1\t1\t1\n1\t1\t1\t1\n1\t1\t1\t1\n"
	with pytest.raises(ValidationError) as e:
		validate2B(ge_1,ssmlist)
	assert 'greater than' in str(e.value)
	# < 0
	le_0 = "1\t1\t1\t1\n1\t1\t-1\t1\n1\t1\t1\t1\n1\t1\t1\t1\n"
	with pytest.raises(ValidationError) as e:
		validate2B(le_0,ssmlist)
	assert 'less than' in str(e.value)
	# Correct entry
	correct = "1\t1\t1\t1\n1\t1\t1\t1\n1\t1\t1\t1\n1\t1\t1\t1\n"
	ccm = np.ones((4,4))
	res = validate2B(correct,ssmlist)
	assert np.testing.assert_allclose(ccm,res) == None


def test_score2A():
	try:
		os.remove('nosuchfile.txt')
	except OSError:
		pass
	s = ['valid.VCF','invalid.VCF','nosuchfile.txt']
	a = [['valid2A.txt'],['invalid2A.txt'],['nosuchfile.txt']]
	for p in itertools.product(a,a,s):
		params = list(p)
		res = scoreChallenge('2A',*params)
		if params[0] == a[0] and params[1] == a[0] and params[2] == s[0]:
			assert res == 1.0
		else:
			assert res == "NA"

def test_score2B():
	try:
		os.remove('nosuchfile.txt')
	except OSError:
		pass
	s = ['valid.VCF','invalid.VCF','nosuchfile.txt']
	b = [['valid2B.txt'],['invalid2B.txt'],['nosuchfile.txt']]
	for p in itertools.product(b,b,s):
		params = list(p)
		res = scoreChallenge('2B',*params)
		if params[0] == b[0] and params[1] == b[0] and params[2] == s[0]:
			assert res == 1.0
		else:
			assert res == "NA"		

def test_verify2A():
	try:
		os.remove('nosuchfile.txt')
	except OSError:
		pass
	s = ['valid.VCF','invalid.VCF','nosuchfile.txt']
	a = [['valid2A.txt'],['invalid2A.txt'],['nosuchfile.txt']]
	for p in itertools.product(a,s):
		params = list(p)
		res = verifyChallenge('2A',*params)
		if params[1] != s[0]:
			assert res == "NA"
		elif params[0] == a[0]:
			assert res == "Valid"
		else:
			assert res == "Invalid"

def test_verify2B():
	try:
		os.remove('nosuchfile.txt')
	except OSError:
		pass
	s = ['valid.VCF','invalid.VCF','nosuchfile.txt']
	b = [['valid2B.txt'],['invalid2B.txt'],['nosuchfile.txt']]
	for p in itertools.product(b,s):
		params = list(p)
		res = verifyChallenge('2B',*params)
		print params
		if params[1] != s[0]:
			assert res == "NA"
		if params[0] == b[0] and params[1] == s[0]:
			assert res == "Valid"
		if params[0] != b[0] and params[1] == s[0]:
			assert res == "Invalid"

def test_validate3A():
	CAs4 = np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
	valid_data = "1\t0\n2\t1\n3\t2\n4\t3\n"
	# Empty file
	with pytest.raises(ValidationError) as e:
		validate3A("",CAs4,"")
	assert "Input file contains a different" in str(e.value)		

	# Too many lines
	CAs3 = np.array([[1,0,0],[0,1,0],[0,0,1],[0,0,1]])
	with pytest.raises(ValidationError) as e:
		validate3A(valid_data,CAs3,"")
	assert "Input file contains a different" in str(e.value)	

	CAs5 = np.array([[1,0,0,0,0],[0,1,0,0,0],[0,0,1,0,0],[0,0,0,1,0]])
	# Too few lines
	with pytest.raises(ValidationError) as e:
		validate3A(valid_data,CAs5,"")
	assert "Input file contains a different" in str(e.value)	

	# Too many columns
	three_col_data = "1\t0\t.5\n2\t1\t.5\n3\t2\t.5\n4\t3\t.5\n"
	with pytest.raises(ValidationError) as e:
		validate3A(three_col_data,CAs4,"")
	assert "Number of tab separated" in str(e.value)	

	# Too few columns
	one_col_data = "43\n32\n21\n10\n"
	with pytest.raises(ValidationError) as e:
		validate3A(one_col_data,CAs4,"")
	assert "Number of tab separated" in str(e.value)	

	# Mis ordered cluster numbers
	mis_data = "4\t3\n3\t2\n2\t1\n1\t0\n"
	with pytest.raises(ValidationError) as e:
		validate3A(mis_data,CAs4,"")
	assert "First column must have" in str(e.value)	

	# Negative cluster numbers
	child_neg_data = "1\t0\n-2\t1\n3\t2\n4\t3\n"
	with pytest.raises(ValidationError) as e:
		validate3A(child_neg_data,CAs4,"")
	assert "First column must have" in str(e.value)	

	# Negative parent cluster
	parent_neg_data = "1\t0\n2\t1\n3\t-2\n4\t3\n"
	with pytest.raises(ValidationError) as e:
		validate3A(parent_neg_data,CAs4,"")
	assert "Parent node label in line 3" in str(e.value)

	# String cluster id
	child_string_data = '1\t0\n"2"\t1\n3\t2\n4\t3\n'
	with pytest.raises(ValidationError) as e:
		validate3A(child_string_data,CAs4,"")
	assert "Entry in line 2 could not be" in str(e.value)	

	# String parent id
	parent_string_data = '1\t"0"\n2\t"1"\n3\t"2"\n4\t"3"\n'
	with pytest.raises(ValidationError) as e:
		validate3A(parent_string_data,CAs4,"")
	assert "Entry in line 1 could not be" in str(e.value)		

	# Reverse parent / child columns
	rev_data = "0\t1\n1\t2\n2\t3\n3\t4\n"
	with pytest.raises(ValidationError) as e:
		validate3A(rev_data,CAs4,"")
	assert 'First column must have' in str(e.value)	
	
	# Disconnected tree
	discon_data = "1\t2\n2\t1\n3\t2\n4\t3\n"
	with pytest.raises(ValidationError) as e:
		validate3A(discon_data,CAs4,"")
	assert 'Root of phylogeny' in str(e.value)
	
	# Valid case
	expected_ADM = np.zeros((4,4))
	expected_ADM[np.triu_indices(4,1)] = 1
	assert np.testing.assert_allclose(validate3A(valid_data, CAs4,""),expected_ADM) == None

	#Valid but parent nodes re-ordered
	rev_valid_data = "1\t2\n2\t3\n3\t4\n4\t0\n"
	expected_rev_ADM = np.zeros((4,4))
	expected_rev_ADM[np.tril_indices(4,-1)] = 1
	res = validate3A(rev_valid_data, CAs4,"")
	assert np.testing.assert_allclose(res,expected_rev_ADM) == None

	# Valid test files
	f = open("valid3A.txt")
	test_data = f.read()
	f.close()
	test_CAs = np.array([[1,0],[0,1],[1,0],[0,1]])
	f = open("valid3B.txt")
	d = f.read()
	f.close()
	test_ADM = np.loadtxt(StringIO.StringIO(d))
	res = validate3A(test_data,test_CAs,"")
	assert np.testing.assert_allclose(res,test_ADM) == None 


def test_validate3B():

	#Valid example
	ssmlist = 4
	valid = "0\t1\t1\t1\n0\t0\t1\t1\n0\t0\t0\t1\n0\t0\t0\t0\n"
	expected_ADM = np.zeros((4,4))
	expected_ADM[np.triu_indices(4,1)] = 1
	res = validate3B(valid,np.identity(4),ssmlist)
	assert np.testing.assert_allclose(res,expected_ADM) == None
	
	# Empty file
	with pytest.raises(ValidationError) as e:
		validate3B("",np.identity(4),ssmlist)
	assert 'Shape of' in str(e.value)	

	# Wrong size: Vector
	vector_entry = "0\t1\t1\t1\t0\t0\t1\t1\t0\t0\t0\t1\t0\t0\t0\t1"
	with pytest.raises(ValidationError) as e:
		validate3B(vector_entry,np.identity(4),ssmlist)
	assert 'Shape of' in str(e.value)	

	# Wrong size: 3x3
	three_by_three = "0\t1\t1\n0\t0\t1\n0\t0\t0\n"
	with pytest.raises(ValidationError) as e:
		validate3B(three_by_three,np.identity(4),ssmlist)
	assert 'Shape of' in str(e.value)
	
	# Ones in diagonal
	one_in_diag = "0\t1\t1\t1\n0\t0\t1\t1\n0\t0\t0\t1\n0\t0\t0\t1\n"
	with pytest.raises(ValidationError) as e:
		validate3B(one_in_diag,np.identity(4),ssmlist)
	assert 'Diagonal entries' in str(e.value)

	# Non-zero in diagonal
	half_in_diag = "0.5\t1\t1\t1\n0\t0\t1\t1\n0\t0\t0\t1\n0\t0\t0\t0\n"
	with pytest.raises(ValidationError) as e:
		validate3B(half_in_diag,np.identity(4),ssmlist)
	assert 'Diagonal entries' in str(e.value)

	#NaN entry
	nan = "0\t1\t1\t1\n0\t0\t1\tNaN\n0\t0\t0\t1\n0\t0\t0\t0\n"
	with pytest.raises(ValidationError) as e:
		validate3B(nan,np.identity(4),ssmlist)
	assert 'contains NaNs' in str(e.value)

	#Infinite entry
	inf = "0\t1\t1\tInf\n0\t0\t1\t1\n0\t0\t0\t1\n0\t0\t0\t0\n"
	with pytest.raises(ValidationError) as e:
		validate3B(inf,np.identity(4),ssmlist)
	assert 'contains non-finite' in str(e.value)
	
	# > 1 entry
	gt1 = "0\t1.01\t1\t1\n0\t0\t1\t1\n0\t0\t0\t1\n0\t0\t0\t0\n"
	with pytest.raises(ValidationError) as e:
		validate3B(gt1,np.identity(4),ssmlist)
	assert 'greater than 1' in str(e.value)

	# < 0 entry
	le0 = "0\t1\t-0.005\t1\n0\t0\t1\t1\n0\t0\t0\t1\n0\t0\t0\t0\n"
	with pytest.raises(ValidationError) as e:
		validate3B(le0,np.identity(4),ssmlist)
	assert 'less than 0' in str(e.value)

	#string entry
	string = "0\t1\t1\t1\n0\t0\t1\tNaN\n0\t0\t0\ta\n0\t0\t0\t0\n"
	with pytest.raises(ValidationError) as e:
		validate3B(string,np.identity(4),ssmlist)
	assert 'could not be cast' in str(e.value)


	# ADi,j + ADj,i > 1
	sumgt1 = "0\t.5\t1\t1\n.6\t0\t1\t1\n0\t0\t0\t1\n0\t0\t0\t0\n"
	with pytest.raises(ValidationError) as e:
		validate3B(sumgt1,np.identity(4),ssmlist)
	assert 'the sum of' in str(e.value)

	# ADi,j + ADj,i + CCi,j > 1
	sum_ccm_gt1 = "0\t.5\t1\t1\n.3\t0\t1\t1\n0\t0\t0\t1\n0\t0\t0\t0\n"
	ccm = np.identity(4)
	ccm[0,1] = 1
	ccm[1,0] = 1

	with pytest.raises(ValidationError) as e:
		validate3B(sum_ccm_gt1,ccm,ssmlist)
	assert 'the sum of' in str(e.value)

def test_score3A():
	try:
		os.remove('nosuchfile.txt')
	except OSError:
		pass
	s = ['valid.VCF','invalid.VCF','nosuchfile.txt']
	two = [['valid2A.txt'],['invalid2A.txt'],['nosuchfile.txt'],[]]
	three = [['valid3A.txt'],['invalid3A.txt'],['nosuchfile.txt'],[]]
	for p in itertools.product(two,three,two,three,s):
		params = [p[0]+p[1],p[2]+p[3],p[4]]
		res = scoreChallenge('3A',*params)
		if params[0:2] == [two[0] + three[0]]*2 and params[-1]==s[0]:
			assert res == 1.0
		else:
			assert res == "NA"

def test_score3B():
	try:
		os.remove('nosuchfile.txt')
	except OSError:
		pass
	s = ['valid.VCF','invalid.VCF','nosuchfile.txt']
	two = [['valid2B.txt'],['invalid2B.txt'],['nosuchfile.txt'],[]]
	three = [['valid3B.txt'],['invalid3B.txt'],['nosuchfile.txt'],[]]
	for p in itertools.product(two,three,two,three,s):
		params = [p[0]+p[1],p[2]+p[3],p[4]]
		res = scoreChallenge('3B',*params)
		if params[0:2] == [two[0] + three[0]]*2 and params[-1]==s[0]:
			assert res == 1.0
		else:
			assert res == "NA"	

def test_verify3A():
	try:
		os.remove('nosuchfile.txt')
	except OSError:
		pass
	s = ['valid.VCF','invalid.VCF','nosuchfile.txt']
	two = [['valid2A.txt'],['invalid2A.txt'],['nosuchfile.txt'],[]]
	three = [['valid3A.txt'],['invalid3A.txt'],['nosuchfile.txt'],[]]
	for p in itertools.product(two,three,s):
		params = [p[0]+p[1],p[2]]
		res = verifyChallenge('3A',*params)
		if params[1] != s[0]:
			assert res == "NA"
		if params[0] == two[0] + three[0] and params[1] == s[0]:
			assert res == "Valid"
		if params[0] != two[0] + three[0] and params[1] == s[0]:
			assert res == "Invalid"

def test_verify3B():
	try:
		os.remove('nosuchfile.txt')
	except OSError:
		pass
	s = ['valid.VCF','invalid.VCF','nosuchfile.txt']
	two = [['valid2B.txt'],['invalid2B.txt'],['nosuchfile.txt'],[]]
	three = [['valid3B.txt'],['invalid3B.txt'],['nosuchfile.txt'],[]]
	for p in itertools.product(two,three,s):
		params = [p[0]+p[1],p[2]]
		res = verifyChallenge('3B',*params)
		if params[1] != s[0]:
			assert res == "NA"
		if params[0] == two[0] + three[0] and params[1] == s[0]:
			assert res == "Valid"
		if params[0] != two[0] + three[0] and params[1] == s[0]:
			assert res == "Invalid"

def test_integration():
	score_mapping = {	'1A': 'python SMCScoring.py 1A sc.json --predfiles valid1A.txt --truthfiles valid1A.txt --vcf valid.VCF',
						'1B': 'python SMCScoring.py 1B sc.json --predfiles valid1B.txt --truthfiles valid1B.txt --vcf valid.VCF',
						'1C': 'python SMCScoring.py 1C sc.json --predfiles valid1C.txt --truthfiles valid1C.txt --vcf valid.VCF',
						'2A': 'python SMCScoring.py 2A sc.json --predfiles valid2A.txt --truthfiles valid2A.txt --vcf valid.VCF',
						'2B': 'python SMCScoring.py 2B sc.json --predfiles valid2B.txt --truthfiles valid2B.txt --vcf valid.VCF',
						'3A': 'python SMCScoring.py 3A sc.json --predfiles valid2A.txt valid3A.txt --truthfiles valid2A.txt valid3A.txt --vcf valid.VCF',
						'3B': 'python SMCScoring.py 3B sc.json --predfiles valid2B.txt valid3B.txt --truthfiles valid2B.txt valid3B.txt  --vcf valid.VCF'
					}

	verify_mapping =  {	'1A': 'python SMCScoring.py 1A sc.json -v --predfiles valid1A.txt --vcf valid.VCF',
						'1B': 'python SMCScoring.py 1B sc.json -v --predfiles valid1B.txt --vcf valid.VCF',
						'1C': 'python SMCScoring.py 1C sc.json -v --predfiles valid1C.txt --vcf valid.VCF',
						'2A': 'python SMCScoring.py 2A sc.json -v --predfiles valid2A.txt --vcf valid.VCF',
						'2B': 'python SMCScoring.py 2B sc.json -v --predfiles valid2B.txt --vcf valid.VCF',
						'3A': 'python SMCScoring.py 3A sc.json -v --predfiles valid2A.txt valid3A.txt --vcf valid.VCF',
						'3B': 'python SMCScoring.py 3B sc.json -v --predfiles valid2B.txt valid3B.txt --vcf valid.VCF'
					}
	for challenge,cmd_string in score_mapping.iteritems():
		try:
			os.remove('sc.json')
		except OSError:
			pass
		print cmd_string
		os.system(cmd_string)
		f = open('sc.json')
		out = json.load(f)
		f.close()
		assert out == {challenge: 1.0}

	for challenge,cmd_string in verify_mapping.iteritems():
		try:
			os.remove('sc.json')
		except OSError:
			pass
		os.system(cmd_string)
		f = open('sc.json')
		out = json.load(f)
		f.close()
		assert out == {challenge: "Valid"}


