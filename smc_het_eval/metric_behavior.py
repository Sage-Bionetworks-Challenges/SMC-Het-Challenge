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

	clusters = np.zeros((600,4))
	clusters[:,:-1] = np.copy(t_clusters)
	clusters[100,:] = np.array([0,0,0,1])
	clusters[300,:] = np.array([0,0,0,1])
	clusters[500,:] = np.array([0,0,0,1])
	ccm = np.dot(clusters,clusters.T)
	res.append(["SmallExtra", calculate2(ccm,t_ccm)])

	res = [map(str,x) for x in res]
	res = ['\t'.join(x) for x in res]
	f = open('scoring2A_cases.tsv', 'w')
	f.write('\n'.join(res))
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
	f = open('scoring2B_beta.tsv', 'w')
	f.write('\n'.join(res))
	f.close()	

if __name__ == '__main__':
	scoring1A_behavior()
	scoring1B_behavior()
	scoring1C_behavior()
	scoring2A_behavior()
	scoring2B_behavior()
