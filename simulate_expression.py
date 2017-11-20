import numpy
import argparse
import datetime
import sys
import os
from math import isclose
from random import shuffle

class MyParser(argparse.ArgumentParser):
    # overriding argument parser to print help message when
    # arguments are not given correctly
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

parser = MyParser(description='Simulate cellcount, genotype and epression data with different MAF and beta settings for two celltypes')
parser.add_argument('-g','--number_of_genotypes', 
                    help='Number of genotypes to simulate', 
                    required=True,
                    type=int)
parser.add_argument('-s','--number_of_samples', 
                    help='Number of samples to simulate', 
                    required=True,
                    type = int)
parser.add_argument('-o','--out_directory', 
                    help='Diretory to write simulated counts, expression and genotypes to', 
                    required=True)
parser.add_argument('-c','--cellcount_mean', help='Cellcount mean to use', 
                    required=False, 
                    default=60,
                    type = float)
parser.add_argument('-d','--cellcount_deviation', help='Cellcount deviation to use', 
                    required=False, 
                    default=5,
                    type = float)
parser.add_argument('-m','--maf', help='comma separated list of maf values to use (will split equally over the groups, starting from lowest)', 
                    required=False, 
                    default="0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5")

def simulate_cellcounts(cellcount_mean, cellcount_deviation, number_of_samples):
    cellcounts_1 = numpy.random.normal(cellcount_mean, cellcount_deviation, number_of_samples)
    cellcounts_2 = [100-x for x in cellcounts_1]
    return(cellcounts_1, cellcounts_2)

def simulate_genotypes(number_of_samples, maf):
    hets = int(number_of_samples*maf)
    alts = int(number_of_samples*maf/2)
    genotypes = [0]*number_of_samples

    for i in range(0, hets):
        genotypes[i] = 1
    for i in range(hets, hets+alts):
        genotypes[i] = 2
    
    simulated_maf = sum(genotypes)/(len(genotypes)*2)
    if not isclose(simulated_maf, maf, abs_tol=0.001):
        raise RuntimeError('Set maf was '+str(maf)+' but simulated maf was '+str(simulated_maf)+
                           '.\nThis is outside tolerated range of +/- 0.001')
    # convert the int genotypes to str so that they can be joined for writing later on
    genotypes = [str(x) for x in genotypes]
    shuffle(genotypes)
    return(genotypes)

def simulate_expression(betas, genotypes):
    expression_per_sample = []
    for index, genotype in enumerate(genotypes):
        # we simulate to test the linear model y ~ cc1 + cc2 + cc1*GT + cc2 * GT
        # so expression = b1 * cc1 + b2 * cc2 + b3 * CC1 * GT + b4 * cc2 * GT
        expression = (beta1 * cellcounts_1[index] + 
                  beta2 * cellcounts_2[index] +
                  beta3 * cellcounts_1[index] * genotype +
                  beta4 * cellcounts_2[index] * genotype)
        expression_per_sample.append(expression)

if __name__ == "__main__":
    args = vars(parser.parse_args())
    # parameter settings
    number_of_genotypes = args['number_of_genotypes']
    number_of_samples = args['number_of_samples']
    cellcount_mean = args['cellcount_mean']
    cellcount_deviation = args['cellcount_deviation']
    maf_list = [float(x) for x in args['maf'].replace(' ','').split(',')]
    out_directory = args['out_directory']

    print('number of genotypes:', number_of_genotypes)
    print('number of samples:',number_of_samples)
    print('cellcount mean:', cellcount_mean)
    print('cellcount deviation:', cellcount_deviation)
    print('maf values to use:',maf_list)
    print('Directory to write output to :',out_directory)


    now = datetime.datetime.now().strftime('%Y%B%d')
    cellcount_file = out_directory+'/simulated_cellcounts_'+str(now)+'.txt'
    genotype_file = out_directory+'/simulated_genotypes_'+str(now)+'.txt'
    expression_file = out_directory+'/simulated_expression_'+str(now)+'.txt'
    
    if not os.path.exists(out_directory):
        os.makedirs(out_directory)

    cellcounts_1, cellcounts_2 = simulate_cellcounts(cellcount_mean, cellcount_deviation, number_of_samples)
    
    with open(cellcount_file,'w') as out:
        out.write('\tcc1\tcc2')
        for index, cc1 in enumerate(cellcounts_1):
            out.write('\nsample_'+str(index)+'\t'+str(cc1)+'\t'+str(cellcounts_2[index]))
    print('written to:',cellcount_file)
    
    with open(genotype_file,'w') as out:
        for sample in range(0, number_of_samples, 1):
            out.write('\tsample_'+str(sample))
        x = 0
        for i in range(0, number_of_genotypes, 1):
            if x >= len(maf_list):
                x = 0
            maf = maf_list[x]
            x+=1
            genotypes = simulate_genotypes(number_of_samples, maf)
        
            out.write('\ngenotype_'+str(i)+'_maf'+str(maf)+'\t'+'\t'.join(genotypes))
    print('written to:',genotype_file)
    #    beta1 = beta2 = beta3 = beta4 = 1
    #    expressionValues = simulate_expression(betas,genotypes)
        


