import numpy
import argparse
import datetime
import sys
import os
from math import isclose
from random import shuffle
import itertools

class MyParser(argparse.ArgumentParser):
    # overriding argument parser to print help message when
    # arguments are not given correctly
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

parser = MyParser(description='Simulate cellcount, genotype and epression data with different MAF and beta settings for two celltypes')
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
    error_list = []
    for index, genotype in enumerate(genotypes):
        # we simulate to test the linear model y ~ cc1 + cc2 + cc1*GT + cc2 * GT
        # so expression = b1 * cc1 + b2 * cc2 + b3 * CC1 * GT + b4 * cc2 * GT
        error = numpy.random.normal(20, 5, 1)
        expression = (float(betas[0]) * cellcounts_1[index] + 
                      float(betas[1]) * cellcounts_2[index] +
                      float(betas[2]) * cellcounts_1[index] * float(genotype) +
                      float(betas[3]) * cellcounts_2[index] * float(genotype) +
                      float(error[0]))
        # change to str so that it can be joined for writing to file later
        expression_per_sample.append(str(expression))
        error_list.append(str(error))
    return(expression_per_sample, error_list)


if __name__ == "__main__":
    args = vars(parser.parse_args())
    # parameter settings
    number_of_samples = args['number_of_samples']
    cellcount_mean = args['cellcount_mean']
    cellcount_deviation = args['cellcount_deviation']
    maf_list = [0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]

    out_directory = args['out_directory']

    print('number of samples:',number_of_samples)
    print('cellcount mean:', cellcount_mean)
    print('cellcount deviation:', cellcount_deviation)
    print('maf values to use:',maf_list)
    print('Directory to write output to :',out_directory)


    now = datetime.datetime.now().strftime('%Y%B%d')
    cellcount_file = out_directory+'/simulated_cellcounts_'+str(now)+'.txt'
    genotype_file = out_directory+'/simulated_genotypes_'+str(now)+'.txt'
    expression_file = out_directory+'/simulated_expression_'+str(now)+'.txt'
    snp_expression_file = out_directory+'/snp_expression_'+str(now)+'.txt'
    error_file = out_directory+'/error_'+str(now)+'.txt'
    
    if not os.path.exists(out_directory):
        os.makedirs(out_directory)

    cellcounts_1, cellcounts_2 = simulate_cellcounts(cellcount_mean, cellcount_deviation, number_of_samples)
    
    with open(cellcount_file,'w') as out:
        out.write('\tcc1\tcc2')
        for index, cc1 in enumerate(cellcounts_1):
            out.write('\nsample_'+str(index)+'\t'+str(cc1)+'\t'+str(cellcounts_2[index]))
    print('written to:',cellcount_file)
    
    # do genotypes and expression in one go because for the expression values I need those specific genotypes
    with open(genotype_file,'w') as out_genotypes, open(expression_file, 'w') as out_expression:
        with open(snp_expression_file,'w') as out_snp_expression, open(error_file, 'w') as out_error:
            out_snp_expression.write('gene\tsnp')
            for sample in range(0, number_of_samples, 1):
                out_genotypes.write('\tsample_'+str(sample))
                out_expression.write('\tsample_'+str(sample))
                out_error.write('\tsample_'+str(sample))
            # loop over all the mafs and all the betas, and calculate
            # the expression level with different mafs and beta settings
            x = 0
            #for i in range(0,10, 1):
            print('Simulating for maf:')
            for maf in maf_list:
                print(maf)
                beta_list = [0, 0.25, 0.5, 0.75, 1]
                betas_to_combine = [list(beta_list) for i in range(0,4, 1)]
                all_beta_combinations = list(itertools.product(*betas_to_combine))
                
                number_of_combi = len(all_beta_combinations)
                print(number_of_combi,'beta combinations to simulate for')
                for betas in all_beta_combinations:
                    # simulate genotypes
                    genotypes = simulate_genotypes(number_of_samples, maf)
                    snp_name = 'genotype_'+str(x)+'__maf_'+str(maf)
                    out_genotypes.write('\n'+snp_name+'\t'+'\t'.join(genotypes))
            
                    # simulate expression
                    expression, error = simulate_expression(betas, genotypes)
                    gene_name = 'expression_'+str(x)+'__betas_'+str(betas[0])+'_'+str(betas[1])+'_'+str(betas[2])+'_'+str(betas[3])
                    out_expression.write('\n'+gene_name+'\t'+'\t'.join(expression))
                    out_error.write('\n'+gene_name+'\t'+'\t'.join(error))
            
                    # write snp-gene coupling file
                    out_snp_expression.write('\n'+gene_name+'\t'+snp_name)
                    x += 1
    print('written to:',genotype_file)
    print('written to:',expression_file)
    print('written to:',snp_expression_file)
    print('written to:',error_file)
    
    
    

