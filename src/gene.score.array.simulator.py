import numpy as np 
import pandas as pd 
import sys
import string
import time
import subprocess
from collections import Counter
import string
import random

def random_pheno_generator(size=6,chars=string.ascii_lowercase):
	return ''.join(random.choice(chars) for _ in range(size))

#First argument is the gene score distribution that you want to draw from, the second is the type of clusters to generate
#If 'large' only clusters with a large number of shared genes will be simulated
#If 'mixed' one cluster with only a few shared genes will be simulated
subprocess.call('mkdir NewSims_nothreshenforced',shell = True)

if len(sys.argv) < 3:
	sys.exit("Enter the ICD10 code of interest as the first argument, and either 'mixed' or 'large' as the second argument depending on desired number of significant genes in a cluster.")

class simulator():

	def __init__(self,type_of_clusters,num_draws,sim_status,percentage,sim_label):
		self.example_dist = pd.read_csv('merged.pegasus.results.'+ sys.argv[1] + '.txt', delimiter = '\t').set_index('Gene')
		self.genes = np.array(self.example_dist.index.tolist())
		self.num_clusters = int(np.random.uniform(2,int(num_draws) * 0.15))
		# self.num_clusters = int(np.random.uniform(2,3))
		self.phenos = np.array([random_pheno_generator() for i in range(num_draws)])
		self.clusters = {}
		self.unique_sig_genes = {}
		self.cluster_type = type_of_clusters
		self.percentage = float(percentage)/100
		self.draw_status = sim_status
		self.sim_label = sim_label
		if self.draw_status == 'limited':
			self.num_draws = num_draws
		else:
			self.num_draws == 100000000000
	def _gen_clusters_(self):
		self.possible_genes = list(self.genes)
		self.possible_phenos = list(self.phenos)
		total_genes = 175
		self.ref_count = {}
		for i in range(self.num_clusters):
			#Set size of clusters, both number of phenos and sig genes
			num_sig_shared_genes = int(total_genes*self.percentage)
			genes,phenos = self.cluster_sharing(num_sig_shared_genes,np.random.randint(2,8),self.possible_genes,self.possible_phenos)
			#Update sets of genes and phenos so that there is not overlap between the clusters (first run)
			self.possible_phenos = list(set(self.possible_phenos).difference(phenos))
			self.possible_genes = list(set(self.possible_genes).difference(genes))
			self.clusters['cluster' + str(i)] = {'Gene':list(genes),'Phenos':list(phenos)}
			for j in phenos:
				self.ref_count[str(j)] = len(genes)

		for i in self.phenos:
			if i not in self.ref_count.keys():
				self.ref_count[i] = 0


		self.unique_genes(self.phenos)

	def cluster_sharing(self,num_unique_genes,num_unique_phenos,possible_genes,possible_phenos):
		genes = set()
		while len(genes) < num_unique_genes:
			genes.add(np.random.choice(possible_genes))
		phenos = set()
		while len(phenos) < num_unique_phenos:
			phenos.add(np.random.choice(possible_phenos))
		return genes,phenos

	def draw_counter(self,gene_dict,selected_genes):
		if self.draw_status == 'limited':
			for i in selected_genes:
				gene_dict[i] +=1
			for x,y in gene_dict.items():
				if y >= self.num_draws:
					del gene_dict[x]
			return gene_dict
		else:
			return gene_dict

	#Generates a list of genes that are also significant for each phenotype, whether or not they have been assigned to a cluster
	def unique_genes(self,phenos):
		self.counter_dict = {}
		for i in self.possible_genes:
			self.counter_dict[i] = 0
		for pheno in phenos:
			self.number_siggenes = 175
			pheno_only_genes = np.random.choice(self.possible_genes, size = int(self.number_siggenes - self.ref_count[pheno]),replace = False)
			self.counter_dict = self.draw_counter(self.counter_dict,pheno_only_genes)
			self.unique_sig_genes[pheno] = list(set(pheno_only_genes))
			self.possible_genes = list(self.counter_dict.keys())

	def generate_matrix(self):
		all_scores = np.array(self.example_dist).flatten()
		small_scores = all_scores[all_scores <= 0.001]
		non_sig_scores = all_scores[all_scores > 0.001]
		data = np.zeros((len(self.phenos),len(self.genes)))
		for j in range(len(self.phenos)):
			data[j] = np.negative(np.log(np.array(np.random.choice(non_sig_scores,len(self.genes)))))
		scorematrix = pd.DataFrame(data.T,index = self.genes,columns = self.phenos)
		for key,value in self.clusters.items():
			for phenotype in value['Phenos']:
				for gene in value['Gene']:
					self.unique_sig_genes[phenotype].append(gene)
		#Fill in significant gene scores that are unique to each phenotype
		for key,value in self.unique_sig_genes.items():
			for x in value:
				scorematrix.loc[x,key] = np.negative(np.log(np.random.choice(small_scores)))

		return scorematrix

	def write(self,dataframe):

		if self.draw_status == 100000000000:
			y = str(self.percentage) + '_' + str(self.sim_label)
			subprocess.call('mkdir NewSims_nothreshenforced/Simulations' + y+str(self.num_clusters),shell = True)
			dataframe = dataframe*-1
			dataframe = 10**dataframe.astype(float)
			dataframe.index.name = 'Gene'
			dataframe.to_csv('NewSims_nothreshenforced/Simulations'+y+str(self.num_clusters)+'/Simulated.scores.using.' + sys.argv[1] + '.gene.dist.' + y + '.csv', header = True, index = True)
			for key,value in self.clusters.items():
				newfile = open('NewSims_nothreshenforced/Simulations'+y+str(self.num_clusters)+ '/' + str(key) + 'gene.and.pheno.info.txt','w')
				newfile.write('Shared Significant Genes:\n')
				newfile.write(','.join(value["Gene"]))
				newfile.write('\nPhenos:\n')
				newfile.write(','.join(value['Phenos']))
		
		else:
			y = str(self.percentage) + '_' + str(self.sim_label)
			subprocess.call('mkdir NewSims_nothreshenforced/Simulations' + y + '_num_draws_' + str(self.num_draws),shell = True)
			dataframe = dataframe*-1
			dataframe = 10**dataframe.astype(float)
			dataframe.index.name = 'Gene'
			dataframe.to_csv('NewSims_nothreshenforced/Simulations'+y+ '_num_draws_' + str(self.num_draws) + '/Simulated.scores.using.' + sys.argv[1] + '.gene.dist.' + str(self.num_clusters) + '.clusters.' + str(self.num_draws)+'.pos.draws.csv', header = True, index = True)
			for key,value in self.clusters.items():
				newfile = open('NewSims_nothreshenforced/Simulations'+y+ '_num_draws_' + str(self.num_draws)+ '/' + str(key) + 'gene.and.pheno.info.txt','w')
				newfile.write('Shared Significant Genes:\n')
				newfile.write(','.join(value["Gene"]))
				newfile.write('\nPhenos:\n')
				newfile.write(','.join(value['Phenos']))

	def test(self):
		self._gen_clusters_()
		self.write(self.generate_matrix())

def main():
	#each item of z is the number of phenotypes in a simulation
	for z in [25,50,75,100]:
		#The amount of shared significant architecture to be imposed on a cluster
		shared_percentage = [1,10,25,50,75]
		for g in shared_percentage:
			#How many simulations for each set of parameters should be run
			for j in range(1,1001):
				print('Generated ' + str(g) + '% with unlimited random draws simulation,' +str(z) + ' phenotypes: ' + str(j)) 
				limiteddraw = simulator(sys.argv[2],z,'limited',g,j)
				limiteddraw.test()
main()


