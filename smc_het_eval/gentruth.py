#!/usr/bin/env python

import os
import yaml
import itertools
import SMCScoring
import numpy as np
import argparse

def getTruth1A(yml_ob,truth_vcf_prefix,vcf_location,filename=''):
	node_muts = get_node_muts(yml_ob, truth_vcf_prefix, vcf_location)
	cellularity = sum([x['phi'] for x in node_muts if x['parent'] == '' and len(x['mutations']) > 0])
	if filename:
		with open(filename,'w') as f:
			f.write('%f\n' % cellularity)
	return cellularity

def getTruth1B(yml_ob,truth_vcf_prefix,vcf_location,filename=''):
	node_muts = get_node_muts(yml_ob, truth_vcf_prefix, vcf_location)
	cancerous = len([x for x in node_muts if len(x['mutations']) > 0 ])
	if filename:
		with open(filename,'w') as f:
			f.write('%d\n' % cancerous)
	return cancerous

def getTruth2A(yml_ob,truth_vcf_prefix,vcf_location,filename=''):
	node_muts = get_node_muts(yml_ob, truth_vcf_prefix, vcf_location)
	called_muts = filter_true_positives(node_muts,vcf_location)
	real_nodes = [x for x in node_muts if len(x['mutations']) > 0]
	mut_map = {}
	for i,n in enumerate(real_nodes):
		for m in n['mutations']:
			mut_map[m] = i+1
	out = [mut_map[x] for x in called_muts]
	out = '\n'.join([str(x) for x in out])
	if filename:
		with open(filename,'w') as f:
			f.write(out)
	return out

def filter_true_positives(node_muts,vcf_location):
	all_true_muts = list(itertools.chain(*[x['mutations'] for x in node_muts]))
	called_muts = get_ids(vcf_location)
	return [x for x in called_muts if x in all_true_muts]

def getScoringVCF(yml_ob,truth_vcf_prefix,vcf_location,filename=''):
	node_muts = get_node_muts(yml_ob, truth_vcf_prefix, vcf_location)
	all_true_muts = list(itertools.chain(*[x['mutations'] for x in node_muts]))
	f = open(vcf_location)
	d = f.readlines()
	f.close()
	head = [x for x in d if x[0] == "#"]
	d = [x.split('\t') for x in d if x[0] != "#"]
	d = [[y.strip('\n') for y in x ] for x in d]
	d = [x + (['True'] if (x[0]+'_'+x[1]) in all_true_muts else ['False']) for x in d]
	out = ''.join(head) + '\n'.join(['\t'.join(x) for x in d])
	if filename:
		with open(filename,'w') as f:
			f.write(''.join(out))
	return out

def getTruth2B(yml_ob,truth_vcf_prefix,vcf_location,filename=''):
	data = getTruth2A(yml_ob,truth_vcf_prefix,vcf_location)
	ccm = SMCScoring.validate2A(data,data.count('\n')+1)
	if filename:
		np.savetxt(filename,ccm)
	return ccm

def getTruth3A(yml_ob,truth_vcf_prefix,vcf_location,filename=''):
	node_muts = get_node_muts(yml_ob, truth_vcf_prefix, vcf_location)
	real_nodes = [x for x in node_muts if len(x['mutations']) > 0]
	node_map = {'':0}
	for i,n in enumerate(real_nodes):
		node_map[n['prefix']] = i+1
	out = [(node_map[x['prefix']], node_map[x['parent']]) for x in real_nodes]
	out = '\n'.join(['%d\t%d' % x for x in out])
	if filename:
		with open(filename,'w') as f:
			f.write(out)
	return out


def getTruth3B(yml_ob,truth_vcf_prefix,vcf_location,filename=''):
	data_2 = getTruth2A(yml_ob,truth_vcf_prefix,vcf_location)
	cluster_assignments = SMCScoring.validate2Afor3A(data_2,data_2.count('\n')+1)
	data_3 = getTruth3A(yml_ob,truth_vcf_prefix,vcf_location)
	adm =  SMCScoring.validate3A(data_3,cluster_assignments,data_2.count('\n')+1)
	if filename:
		np.savetxt(filename,adm)
	return adm

def get_ids(filename):
	if not os.path.exists(filename):
		return []
	f = open(filename)
	d = f.readlines()
	f.close()
	d = [x for x in d if x[0] != "#"]
	d = [x.split('\t') for x in d]
	d = [x[0]+'_'+x[1] for x in d]
	return d

def construct_loc_node_map(yml_ob,vcf_prefix):
	chrs = [str(x) for x in range(1,23)] + ['X','Y']
	mut_types = ['snv','sv','indel']
	def decend(node, par_muts, par_prefix, ind):
		if not node.has_key('children') and (node['mut_type'] or par_muts):
			muts = [get_ids(vcf_prefix+'chr%s%s_c%d_%s.vcf' % (c,par_prefix,ind,s)) for c,s in  itertools.product(chrs,mut_types)]
			muts = itertools.chain(*muts)
			return [{'prefix':'%s_c%d' % (par_prefix,ind),'parent': par_prefix, 'mutations':set(muts)}]
		elif node.has_key('children'):
			ret = []
			lin_muts = None
			for i,ch in enumerate(node['children']):
				if ind == -1:
					pref = ''
					if not ch['mut_type']:
						continue
				else:
					pref = par_prefix + '_c%d' % ind
				muts = decend(ch,node['mut_type'],pref, i)
				if not ch['mut_type'] and not ch.has_key('children'):
					lin_muts = muts[0]
					lin_muts['prefix'] = lin_muts['parent']
					lin_muts['parent'] = '_'.join(lin_muts['parent'].split('_')[:-1])
				else:
					ret = ret + muts
			if lin_muts:
				for i,mut in enumerate(ret):
					if mut['mutations'] != lin_muts['mutations']:
						ret[i]['mutations'] = ret[i]['mutations'] - lin_muts['mutations']
				return [lin_muts] + ret
			return ret
		else:
			return [[]]
	node_muts = decend(yml_ob[1]['root'],None,'',-1)
	phis = get_phis(yml_ob)
	for i,n in enumerate(node_muts):
		node_muts[i]['phi'] = phis[i]
	return node_muts

def get_phis(yml_ob):
	def decend(node,par_perc):
		if not node.has_key('children'):
			return [node['percent'] * par_perc]
		else:
			ret = [node['percent'] * par_perc]
			for ch in node['children']:
				if ch['mut_type']:
					ret = ret + decend(ch,node['percent']*par_perc)
		return ret

	return decend(yml_ob[1]['root'],1.0)[1:]

def get_node_muts(yml_ob, truth_vcf_prefix, vcf_location):
	node_muts = construct_loc_node_map(yml_ob,truth_vcf_prefix)
	muts = get_ids(vcf_location)
	for node in node_muts:
		node['mutations'] = [x for x in node['mutations'] if x in muts]
	return node_muts

def getTruth1C(yml_ob, truth_vcf_prefix, vcf_location, filename=''):
	node_muts = get_node_muts(yml_ob, truth_vcf_prefix, vcf_location)
	real_nodes = [x for x in node_muts if len(x['mutations']) > 0]
	all_true_muts = list(itertools.chain(*[x['mutations'] for x in real_nodes]))
	called_muts = get_ids(vcf_location)
	false_pos_muts = [x for x in called_muts if x not in all_true_muts]
	real_nodes.append({'mutations':false_pos_muts,'phi':0.0})
	out = ['%d\t%d\t%f\n' % (i+1,len(x['mutations']),x['phi']) for i,x in enumerate(real_nodes)]
	if filename:
		with open(filename,'w') as f:
			f.write(''.join(out))
	return out

func_map = 	{	'1A':{'func': getTruth1A,'extension':'txt'},
				'1B':{'func': getTruth1B,'extension':'txt'},
				'1C':{'func': getTruth1C,'extension':'txt'},
				'2A':{'func': getTruth2A,'extension':'txt'},
				'2B':{'func': getTruth2B,'extension':'gz'},
				'3A':{'func': getTruth3A,'extension':'txt'},
				'3B':{'func': getTruth3B,'extension':'gz'},
				'scoring_vcf':{'func': getScoringVCF,'extension':'vcf'}
			}

def generate_truth(yaml_file,truth_vcf_prefix,vcf_location,output_prefix):
	with open(yaml_file, 'r') as stream:
		yaml_ob = yaml.load(stream)
	for challenge in func_map:
		filename = '%s.truth.%s.%s' % (output_prefix, challenge, func_map[challenge]['extension'])
		func_map[challenge]['func'](yaml_ob,truth_vcf_prefix,vcf_location,filename)


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("truth_yaml")
	parser.add_argument("truth_vcf_prefix")
	parser.add_argument("mutect_vcf")
	parser.add_argument("output_prefix")

	args = parser.parse_args()

	generate_truth(args.truth_yaml,args.truth_vcf_prefix,args.mutect_vcf,args.output_prefix)
