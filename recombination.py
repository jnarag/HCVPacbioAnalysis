from Bio import SeqIO

from Bio import AlignIO
import os
import numpy as np
import scipy.stats as stats
from scipy.stats import ks_2samp
from pylab import *
from scipy import optimize
import re
import pandas as pd



''' Estimating Recombination Rate - Adapted from Neher and Leitner (2010) '''

cutoff = 3
minsamplesize = 0
no_bins = 10


L = 1680 #alignment length
G = 1  # generations in neher hiv - 15 per month (convert time into months)
sample_length = 1# neher hiv 7 months

file ="LD_versus_genetic_distance_cutoff_"+str(cutoff)+".csv"
print "sample_length", sample_length
output_handle = open('%s'%file,'w')
output_handle.write("key,distance,D_prime,r,AB,Ab,aB,ab\n")

def getBiallelicSites(fastain, aln_over_time):

    file = SeqIO.parse('%s'%fastain, 'fasta')

    aln = AlignIO.read('%s'%fastain, 'fasta')

    biallelic_aln = {}


    for i in xrange(len(aln[1,:])):

        biallelic_aln[i+1] = aln[:,i] #i+1 is site number, which is a key in the biallelic_aln dict.


    aln_over_time.append(biallelic_aln)
    #return biallelic_aln


def count_biallelic_haplotypes(no_of_haplotypes, fastain, aln_over_time):

    #keys = (biallelic_aln).keys()
    #keys.sort()

    #file = fastain.replace(".fasta", "_LD_versus_genetic_distance_cutoff_"+str(cutoff)+".csv")
    #file = "LD_versus_genetic_distance_cutoff_"+str(cutoff)+".csv"

    aln = AlignIO.read('%s'%fastain, 'fasta')

    no_sites = len(aln[1,:])

    biallelic_3_haplotypes = {}

    biallelic_aln = {}
    for i in range(no_sites):

        biallelic_aln[i+1] = aln[:,i]

        if str(aln[0:,i]).isupper():
            A_i = aln[:,i].count('A')
            C_i = aln[:,i].count('C')
            G_i = aln[:,i].count('G')
            T_i = aln[:,i].count('T')

        else:
            A_i = aln[:,i].count('a')
            C_i = aln[:,i].count('c')
            G_i = aln[:,i].count('g')
            T_i = aln[:,i].count('t')

        counts_i = [A_i, C_i, G_i, T_i]
        major_af_i = 0
        minor_af_i = 0

        #print counts_i

        if(len(aln[:,i]) == sum(counts_i) and counts_i.count(0)==2):

            #site i is biallelic and gapless




            for c in range(len(counts_i)):

                if counts_i[c] > 0:

                    af = float(counts_i[c])/float(sum(counts_i))

                    if af < 0.5: minor_af_i = af; major_af_i = 1.0-af
                    else: major_af_i = af; minor_af_i = 1.0-major_af_i


        if minor_af_i*sum(counts_i) > cutoff+0.1:

            for ii in range(i+1, no_sites):

                # print "i, ii", i, ii
                # comparing pair of biallelic sites

                haplotype_set = set()
                key_i = i
                key_ii = ii

                #print 'y', key_i, key_ii
                #print key_i, key_ii

                no_of_seqs = lno_sites = len(aln[:,1]) # this is the same as no of sites - simplify

                # A_i = aln[:,i].count('A')
                # C_i = aln[:,i].count('C')
                # G_i = aln[:,i].count('G')
                # T_i = aln[:,i].count('T')
                #
                # counts_i = [A_i, C_i, G_i, T_i]
                #
                # major_af_i = 0
                # minor_af_i = 0
                #
                #
                # if(len(aln[:,i]) == sum(counts_i) and counts_i.count(0)==2):
                #
                #     #site i is biallelic and gapless
                #
                #
                #
                #
                #     for c in range(len(counts_i)):
                #
                #         if counts_i[c] > 0:
                #
                #             af = float(counts_i[c])/float(sum(counts_i))
                #
                #             if af < 0.5: minor_af_i = af; major_af_i = 1.0-af
                #             else: major_af_i = af; minor_af_i = 1.0-major_af_i
                #


                if str(aln[0:,i]).isupper():

                    A_ii = aln[:,ii].count('A')
                    C_ii = aln[:,ii].count('C')
                    G_ii = aln[:,ii].count('G')
                    T_ii = aln[:,ii].count('T')

                else:
                    A_ii = aln[:,ii].count('a')
                    C_ii = aln[:,ii].count('c')
                    G_ii = aln[:,ii].count('g')
                    T_ii = aln[:,ii].count('t')

                counts_ii = [A_ii, C_ii, G_ii, T_ii]


                major_af_ii = 0
                minor_af_ii = 0


                if(len(aln[:,ii])==sum(counts_ii) and counts_ii.count(0)==2):

                    # then site ii is a biallelic site with no gaps

                    # print i+1, aln[:,i]
                    # print ii+1, aln[:,ii]
                    #
                    for c in range(len(counts_ii)):

                        if counts_ii[c] > 0:

                            af = float(counts_ii[c])/float(sum(counts_ii))

                            if af < 0.5: minor_af_ii = af; major_af_ii = 1.0-af
                            else: major_af_ii = af; minor_af_ii = 1.0-major_af_ii


                    if minor_af_i*sum(counts_i) > cutoff+0.1 and minor_af_ii*sum(counts_ii) > cutoff+0.1:


                        #print "aln", aln[:, i], aln[:,ii]

                        dd = []

                        # find the unique biallelic haplotypes
                        for j in range(no_of_seqs):

                            #haplotype = biallelic_aln[key_i+1][j] + biallelic_aln[key_ii+1][j]
                            haplotype = aln[j,i]+aln[j,ii]
                            haplotype_set.add(haplotype)


                        if no_of_haplotypes is None:
                            #if(len(haplotype_set)) > 3:

                            #print keys[i], keys[ii], len(haplotype_set)


                            #biallelic_3_haplotypes.append([keys[i], keys[ii]])

                            key = str(key_i+1)+"_"+str(key_ii+1)
                            #print i, ii, "key", key
                            biallelic_3_haplotypes[key] = haplotype_set
                        else:


                            if(len(haplotype_set)) == no_of_haplotypes:

                                key = str(i+1)+"_"+str(ii+1)
                                biallelic_3_haplotypes[key] = haplotype_set
                                for j in range(no_of_seqs):
                                    haplotype = aln[j, i] + aln[j, ii]
                                    dd.append(haplotype)

                                LD_results = calcLD(haplotype_set, dd, len(aln[:,i]))
                                print key, abs(key_i-key_ii), LD_results[0], LD_results[1]


                                AB = str(LD_results[2])
                                Ab = str(LD_results[3])
                                aB = str(LD_results[4])
                                ab = str(LD_results[5])

                                output_handle.write(key+","+str(abs(key_i-key_ii))+","+str(LD_results[0])+","+str(LD_results[1])+str(","+AB+","+Ab+","+aB+","+ab)+"\n")
                                # for each in haplotype_set:
                                #
                                #     print each, float(dd.count(each))/float(len(aln[:,i]))
                                #
                                # print

    aln_over_time.append(biallelic_aln)
    return biallelic_3_haplotypes

def count_biallelic_haplotypes_bootstrap(no_of_haplotypes, fastain, aln_over_time, bootstrap_sites):

    #keys = (biallelic_aln).keys()
    #keys.sort()

    #file = fastain.replace(".fasta", "_LD_versus_genetic_distance_cutoff_"+str(cutoff)+".csv")
    #file = "LD_versus_genetic_distance_cutoff_"+str(cutoff)+".csv"

    aln = AlignIO.read('%s'%fastain, 'fasta')

    no_sites = len(aln[1,:])

    #bootstrap_sites = np.random.randint(no_sites, size = no_sites)


    biallelic_3_haplotypes = {}

    biallelic_aln = {}
    for i in range(no_sites):

        site1 = bootstrap_sites[i]

        print site1, no_sites, fastain
        biallelic_aln[i+1] = aln[:,site1]

        if str(aln[0:,site1]).isupper():
            A_i = aln[:,site1].count('A')
            C_i = aln[:,site1].count('C')
            G_i = aln[:,site1].count('G')
            T_i = aln[:,site1].count('T')

        else:
            A_i = aln[:,site1].count('a')
            C_i = aln[:,site1].count('c')
            G_i = aln[:,site1].count('g')
            T_i = aln[:,site1].count('t')

        counts_i = [A_i, C_i, G_i, T_i]
        major_af_i = 0
        minor_af_i = 0

        #print counts_i

        if(len(aln[:,site1]) == sum(counts_i) and counts_i.count(0)==2):

            #site i is biallelic and gapless




            for c in range(len(counts_i)):

                if counts_i[c] > 0:

                    af = float(counts_i[c])/float(sum(counts_i))

                    if af < 0.5: minor_af_i = af; major_af_i = 1.0-af
                    else: major_af_i = af; minor_af_i = 1.0-major_af_i


        if minor_af_i*sum(counts_i) > cutoff+0.1:

            for ii in range(i+1, no_sites):

                # print "i, ii", i, ii
                # comparing pair of biallelic sites

                site2 = bootstrap_sites[ii]

                haplotype_set = set()
                key_i = i
                key_ii = ii

                #print 'y', key_i, key_ii
                #print key_i, key_ii

                no_of_seqs = lno_sites = len(aln[:,1]) # this is the same as no of sites - simplify

                # A_i = aln[:,i].count('A')
                # C_i = aln[:,i].count('C')
                # G_i = aln[:,i].count('G')
                # T_i = aln[:,i].count('T')
                #
                # counts_i = [A_i, C_i, G_i, T_i]
                #
                # major_af_i = 0
                # minor_af_i = 0
                #
                #
                # if(len(aln[:,i]) == sum(counts_i) and counts_i.count(0)==2):
                #
                #     #site i is biallelic and gapless
                #
                #
                #
                #
                #     for c in range(len(counts_i)):
                #
                #         if counts_i[c] > 0:
                #
                #             af = float(counts_i[c])/float(sum(counts_i))
                #
                #             if af < 0.5: minor_af_i = af; major_af_i = 1.0-af
                #             else: major_af_i = af; minor_af_i = 1.0-major_af_i
                #


                if str(aln[0:,site2]).isupper():

                    A_ii = aln[:,site2].count('A')
                    C_ii = aln[:,site2].count('C')
                    G_ii = aln[:,site2].count('G')
                    T_ii = aln[:,site2].count('T')

                else:
                    A_ii = aln[:,site2].count('a')
                    C_ii = aln[:,site2].count('c')
                    G_ii = aln[:,site2].count('g')
                    T_ii = aln[:,site2].count('t')

                counts_ii = [A_ii, C_ii, G_ii, T_ii]


                major_af_ii = 0
                minor_af_ii = 0


                if(len(aln[:,site2])==sum(counts_ii) and counts_ii.count(0)==2):

                    # then site ii is a biallelic site with no gaps

                    # print i+1, aln[:,i]
                    # print ii+1, aln[:,ii]
                    #
                    for c in range(len(counts_ii)):

                        if counts_ii[c] > 0:

                            af = float(counts_ii[c])/float(sum(counts_ii))

                            if af < 0.5: minor_af_ii = af; major_af_ii = 1.0-af
                            else: major_af_ii = af; minor_af_ii = 1.0-major_af_ii


                    if minor_af_i*sum(counts_i) > cutoff+0.1 and minor_af_ii*sum(counts_ii) > cutoff+0.1:


                        #print "aln", aln[:, i], aln[:,ii]

                        dd = []

                        # find the unique biallelic haplotypes
                        for j in range(no_of_seqs):

                            #haplotype = biallelic_aln[key_i+1][j] + biallelic_aln[key_ii+1][j]
                            haplotype = aln[j,site1]+aln[j,site2]
                            haplotype_set.add(haplotype)


                        if no_of_haplotypes is None:
                            #if(len(haplotype_set)) > 3:

                            #print keys[i], keys[ii], len(haplotype_set)


                            #biallelic_3_haplotypes.append([keys[i], keys[ii]])

                            key = str(key_i+1)+"_"+str(key_ii+1)
                            #print i, ii, "key", key
                            biallelic_3_haplotypes[key] = haplotype_set
                        else:


                            if(len(haplotype_set)) == no_of_haplotypes:

                                key = str(i+1)+"_"+str(ii+1)
                                biallelic_3_haplotypes[key] = haplotype_set
                                for j in range(no_of_seqs):
                                    haplotype = aln[j, site1] + aln[j, site2]
                                    dd.append(haplotype)

                                LD_results = calcLD(haplotype_set, dd, len(aln[:,site1]))
                                print key, abs(key_i-key_ii), LD_results[0], LD_results[1]


                                AB = str(LD_results[2])
                                Ab = str(LD_results[3])
                                aB = str(LD_results[4])
                                ab = str(LD_results[5])

                                output_handle.write(key+","+str(abs(key_i-key_ii))+","+str(LD_results[0])+","+str(LD_results[1])+str(","+AB+","+Ab+","+aB+","+ab)+"\n")
                                # for each in haplotype_set:
                                #
                                #     print each, float(dd.count(each))/float(len(aln[:,i]))
                                #
                                # print

    aln_over_time.append(biallelic_aln)
    return biallelic_3_haplotypes

def count_pair_of_biallelic_sites_from_aln(file, aln_over_time):

    #print 'timepoint', timepoint
    #biallelic_aln = getBiallelicSites(file)
    #aln_over_time.append(biallelic_aln)

    print file
    return count_biallelic_haplotypes(None, file, aln_over_time)

def count_pair_of_biallelic_sites_from_aln_bootstrap(file, aln_over_time, bootstrap_sites):

    #print 'timepoint', timepoint
    #biallelic_aln = getBiallelicSites(file)
    #aln_over_time.append(biallelic_aln)

    #print file
    return count_biallelic_haplotypes_bootstrap(None, file, aln_over_time, bootstrap_sites)



recomb_haplotype_freq = {}

patient_hist = {}

bins_start = 0
bins_length = L*G/no_bins
print bins_length
bins = []



total_hist = {}

for i in range(no_bins):

    bins.append([bins_start, bins_start+bins_length])
    key = int((bins_start+bins_start+bins_length)/2)
    bins_start = bins_start+bins_length
    print bins_start

    patient_hist[key] = []
    total_hist[key] = []

def get_recomb_haplotypes(filedir, bins, datepos):

    filenames = os.listdir(filedir)
    filenames.sort(key=natural_keys)

    #bins = []#[[0,8000], [8000,16000], [16000,24000], [24000,32000], [32000,40000], [40000,48000], [48000, 56000]]
    biallelic_sites_3_haplotypes_by_time = {}
    hist = {}


    Mp1p2 = []
    rec_missinghap_copies = np.zeros((bins))
    rec_diversification_dis = np.zeros((bins))
    diversificationcount = np.zeros((bins))
    nonrec_diversification_count = np.zeros((bins))

    missinghap_counter = 0

    aln_over_time = []


    max = 0
    for i in range(len(filenames)):

        filename = filenames[i]

        if '_' in filename and ".fasta" in filename:
            parts = filename.split('_')
            parts2 = parts[len(parts)-1].split('.fasta')
            print parts2, parts2[len(parts2)-1]
            #date = float(parts2[len(parts2)-2])
            #print "date", date
            date = float(parts2[len(parts2)-2])


            # if date > 2.5:
            #     break
            #date = float(parts[len(parts)-2]) # subsampled dataset
            '''count pairs of biallelic sites'''
            biallelic_sites_3_haplotypes_by_time[date] = count_pair_of_biallelic_sites_from_aln(filedir+filename, aln_over_time)

    #print filedir

    timepoints = biallelic_sites_3_haplotypes_by_time.keys()

    #print timepoints

    timepoints.sort()
    missing_haplotypes_set = set()



    for i in range(len(timepoints)-1):
        # for each time point get the three haplotype set

        pair_of_biallelic_sites_haplotypes = biallelic_sites_3_haplotypes_by_time[timepoints[i]]


        for key in pair_of_biallelic_sites_haplotypes.keys():

            biallelic_sites = key.split('_')

            site1 = int(str(biallelic_sites[0]).strip())
            site2 = int(str(biallelic_sites[1]).strip())

            haplotypes_current = list(pair_of_biallelic_sites_haplotypes[key])

            if len(haplotypes_current) == 3:

                if (haplotypes_current[0][0]+haplotypes_current[1][1])not in haplotypes_current: missing_haplotype=[haplotypes_current[0][0]+haplotypes_current[1][1]]
                elif (haplotypes_current[0][0]+haplotypes_current[2][1]) not in haplotypes_current: missing_haplotype=[haplotypes_current[0][0]+haplotypes_current[2][1]]
                elif (haplotypes_current[1][0]+haplotypes_current[0][1]) not in haplotypes_current: missing_haplotype=[haplotypes_current[1][0]+haplotypes_current[0][1]]
                elif (haplotypes_current[1][0]+haplotypes_current[2][1]) not in haplotypes_current: missing_haplotype=[haplotypes_current[1][0]+haplotypes_current[2][1]]
                elif (haplotypes_current[2][0]+haplotypes_current[0][1]) not in haplotypes_current: missing_haplotype=[haplotypes_current[2][0]+haplotypes_current[0][1]]
                elif (haplotypes_current[2][0]+haplotypes_current[1][1]) not in haplotypes_current: missing_haplotype=[haplotypes_current[2][0]+haplotypes_current[1][1]]
                #for j in range(i+1, len(timepoints)):

                j = i+1
                aln_j = aln_over_time[j]
                aln_i = aln_over_time[i]

                no_of_biallelic_sites = len(aln_j) # no of biallelic sites

                #print aln_j.keys()

                aln_j_keys = aln_j.keys()
                aln_j_keys.sort()

                haplotype_new = set()


                n_seqs = len(aln_j[site1])
                p1 = float(aln_j[site1].count(missing_haplotype[0][0]))/float(n_seqs)
                p2 = float(aln_j[site2].count(missing_haplotype[0][0]))/float(n_seqs)
                Mp1p2.append(float(n_seqs)*p1*p2)


                if site1 in aln_j_keys and site2 in aln_j_keys:

                    # n_seqs = len(aln_j[site1])
                    # p1 = float(aln_j[site1].count(missing_haplotype[0][0]))/float(n_seqs)
                    # p2 = float(aln_j[site2].count(missing_haplotype[0][0]))/float(n_seqs)
                    # Mp1p2.append(float(n_seqs)*p1*p2)

                    for n in range(0, n_seqs):
                        haplotype = aln_j[site1][n]+aln_j[site2][n]
                        #print haplotype
                        haplotype_new.add(haplotype)

                    #missing_haplotype1 = list(haplotype_new - set(haplotypes_current))

                    time_diff = (float(timepoints[j])-float(timepoints[i])) #multiply by 12 if you want in months

                    dist = abs(float(site2)-float(site1))

                    if(time_diff*dist > max):

                        max = time_diff*dist

                    ii=min(np.floor(dist*time_diff*no_bins/L/sample_length), no_bins-1)
                    if missing_haplotype[0] in haplotype_new:
                        #print missing_haplotype, missing_haplotype1
                        missing_haplotype_key = missing_haplotype[0]+'_'+str(site1)+"_"+str(site2)
                        missinghap_counter+=1
                        # p1 = float(aln_j[site1].count(missing_haplotype[0][0]))/float(n_seqs)
                        # p2 = float(aln_j[site2].count(missing_haplotype[0][0]))/float(n_seqs)
                        # Mp1p2.append(float(n_seqs)*p1*p2)


                        #if missing_haplotype_key not in missing_haplotypes_set:

                        #print timepoints[i], missing_haplotype_key
                        missing_haplotypes_set.add(missing_haplotype_key)

                        count = 0

                        #print "count", (aln_j[site1]).count(missing_haplotype[0][0])
                        for n in range(0, n_seqs):



                            if aln_j[site1][n]+aln_j[site2][n] == missing_haplotype[0]:

                                count +=1



                        frequency = float(count)/float(n_seqs)
                        #print missing_haplotype[0], count, frequency, n_seqs, site1, site2




                        #dist_times_gen = (time_diff*365/2.0)*(float(site2)-float(site1))
                        #print timepoints[j], timepoints[i], site1, site2, site2-site1,  dist_times_gen, frequency, missing_haplotype[0], time_diff, haplotype_new, haplotypes_current, ii, float(n_seqs)*p1*p2


                        rec_missinghap_copies[ii] = rec_missinghap_copies[ii]+count

                        rec_diversification_dis[ii] += 1


                        #Mp1p2.append(float(n_seqs)*p1*p2)


                    elif(len(haplotype_new)>3):
                        nonrec_diversification_count[ii]=nonrec_diversification_count[ii]+1
                    if(len(haplotype_new)>4):
                        nonrec_diversification_count[ii] += 1

                diversificationcount[ii]+=1

        #print 'xx', len(pair_of_biallelic_sites_haplotypes)
    print 'max', max
    return Mp1p2, rec_diversification_dis, rec_missinghap_copies, diversificationcount, nonrec_diversification_count


def get_recomb_haplotypes_bootstrap(filedir, bins, datepos, bootstrap_sites):

    filenames = os.listdir(filedir)
    filenames.sort(key=natural_keys)

    #bins = []#[[0,8000], [8000,16000], [16000,24000], [24000,32000], [32000,40000], [40000,48000], [48000, 56000]]
    biallelic_sites_3_haplotypes_by_time = {}
    hist = {}


    Mp1p2 = []
    rec_missinghap_copies = np.zeros((bins))
    rec_diversification_dis = np.zeros((bins))
    diversificationcount = np.zeros((bins))
    nonrec_diversification_count = np.zeros((bins))

    missinghap_counter = 0

    aln_over_time = []


    max = 0
    for i in range(len(filenames)):

        filename = filenames[i]

        if '_' in filename and ".fas" in filename:
            parts = filename.split('_')
            parts2 = parts[len(parts)-1].split('.fas')
            date = float(parts2[len(parts2)-2])
            print "date", date
            #date = float(parts[len(parts)-2])

            #print parts2, parts2[len(parts2)-1]

            # if date > 2.5:
            #     break
            #date = float(parts[len(parts)-2]) # subsampled dataset
            '''count pairs of biallelic sites'''
            biallelic_sites_3_haplotypes_by_time[date] = count_pair_of_biallelic_sites_from_aln_bootstrap(filedir+filename, aln_over_time, bootstrap_sites)

    #print filedir

    timepoints = biallelic_sites_3_haplotypes_by_time.keys()

    #print timepoints

    timepoints.sort()
    missing_haplotypes_set = set()



    for i in range(len(timepoints)-1):
        # for each time point get the three haplotype set

        pair_of_biallelic_sites_haplotypes = biallelic_sites_3_haplotypes_by_time[timepoints[i]]

        #print len(biallelic_sites_haplotypes), biallelic_sites_haplotypes

        for key in pair_of_biallelic_sites_haplotypes.keys():

            biallelic_sites = key.split('_')

            site1 = int(str(biallelic_sites[0]).strip())
            site2 = int(str(biallelic_sites[1]).strip())

            #print site1, site2
            haplotypes_current = list(pair_of_biallelic_sites_haplotypes[key])

            if len(haplotypes_current) == 3:

                if (haplotypes_current[0][0]+haplotypes_current[1][1])not in haplotypes_current: missing_haplotype=[haplotypes_current[0][0]+haplotypes_current[1][1]]
                elif (haplotypes_current[0][0]+haplotypes_current[2][1]) not in haplotypes_current: missing_haplotype=[haplotypes_current[0][0]+haplotypes_current[2][1]]
                elif (haplotypes_current[1][0]+haplotypes_current[0][1]) not in haplotypes_current: missing_haplotype=[haplotypes_current[1][0]+haplotypes_current[0][1]]
                elif (haplotypes_current[1][0]+haplotypes_current[2][1]) not in haplotypes_current: missing_haplotype=[haplotypes_current[1][0]+haplotypes_current[2][1]]
                elif (haplotypes_current[2][0]+haplotypes_current[0][1]) not in haplotypes_current: missing_haplotype=[haplotypes_current[2][0]+haplotypes_current[0][1]]
                elif (haplotypes_current[2][0]+haplotypes_current[1][1]) not in haplotypes_current: missing_haplotype=[haplotypes_current[2][0]+haplotypes_current[1][1]]
                #for j in range(i+1, len(timepoints)):

                j = i+1
                aln_j = aln_over_time[j]
                aln_i = aln_over_time[i]

                no_of_biallelic_sites = len(aln_j) # no of biallelic sites

                #print aln_j.keys()

                aln_j_keys = aln_j.keys()
                aln_j_keys.sort()

                haplotype_new = set()


                n_seqs = len(aln_j[site1])
                p1 = float(aln_j[site1].count(missing_haplotype[0][0]))/float(n_seqs)
                p2 = float(aln_j[site2].count(missing_haplotype[0][0]))/float(n_seqs)
                Mp1p2.append(float(n_seqs)*p1*p2)


                if site1 in aln_j_keys and site2 in aln_j_keys:

                    # n_seqs = len(aln_j[site1])
                    # p1 = float(aln_j[site1].count(missing_haplotype[0][0]))/float(n_seqs)
                    # p2 = float(aln_j[site2].count(missing_haplotype[0][0]))/float(n_seqs)
                    # Mp1p2.append(float(n_seqs)*p1*p2)

                    for n in range(0, n_seqs):
                        haplotype = aln_j[site1][n]+aln_j[site2][n]
                        #print haplotype
                        haplotype_new.add(haplotype)

                    #missing_haplotype1 = list(haplotype_new - set(haplotypes_current))

                    time_diff = (float(timepoints[j])-float(timepoints[i])) #multiply by 12 if you want in months

                    dist = abs(float(site2)-float(site1))

                    if(time_diff*dist > max):

                        max = time_diff*dist

                    ii=min(np.floor(dist*time_diff*no_bins/L/sample_length), no_bins-1)
                    if missing_haplotype[0] in haplotype_new:
                        #print missing_haplotype, missing_haplotype1
                        missing_haplotype_key = missing_haplotype[0]+'_'+str(site1)+"_"+str(site2)
                        missinghap_counter+=1
                        # p1 = float(aln_j[site1].count(missing_haplotype[0][0]))/float(n_seqs)
                        # p2 = float(aln_j[site2].count(missing_haplotype[0][0]))/float(n_seqs)
                        # Mp1p2.append(float(n_seqs)*p1*p2)


                        #if missing_haplotype_key not in missing_haplotypes_set:

                        #print timepoints[i], missing_haplotype_key
                        missing_haplotypes_set.add(missing_haplotype_key)

                        count = 0

                        #print "count", (aln_j[site1]).count(missing_haplotype[0][0])
                        for n in range(0, n_seqs):



                            if aln_j[site1][n]+aln_j[site2][n] == missing_haplotype[0]:

                                count +=1



                        frequency = float(count)/float(n_seqs)
                        #print missing_haplotype[0], count, frequency, n_seqs, site1, site2




                        #dist_times_gen = (time_diff*365/2.0)*(float(site2)-float(site1))
                        #print timepoints[j], timepoints[i], site1, site2, site2-site1,  dist_times_gen, frequency, missing_haplotype[0], time_diff, haplotype_new, haplotypes_current, ii, float(n_seqs)*p1*p2


                        rec_missinghap_copies[ii] = rec_missinghap_copies[ii]+count

                        rec_diversification_dis[ii] += 1


                        #Mp1p2.append(float(n_seqs)*p1*p2)


                    elif(len(haplotype_new)>3):
                        nonrec_diversification_count[ii]=nonrec_diversification_count[ii]+1
                    if(len(haplotype_new)>4):
                        nonrec_diversification_count[ii] += 1

                diversificationcount[ii]+=1



        #print 'xx', len(pair_of_biallelic_sites_haplotypes)
    print 'max', max
    return Mp1p2, rec_diversification_dis, rec_missinghap_copies, diversificationcount, nonrec_diversification_count

def haplotype_diversfication_timeavg(bins, patients, datepos):

    rec_diversification_dis = np.zeros((bins))
    rec_missinghap_copies = np.zeros((bins))
    Mp1p2 = []
    diversification_count = np.zeros((bins))
    nonrec_diversification_dis = np.zeros((bins))

    print len(rec_diversification_dis), len(rec_missinghap_copies), len(diversification_count)

    rec_err=np.zeros((bins))
    nonrec_err=np.zeros((bins))
    fs=16   #font size in plots
    ''' calculate center bin values, 764*7 are magic numbers, 15 is the number generations per month '''
    binpos=np.zeros((bins))

    #for HIV neher data
    #for ii in range(bins): binpos[ii]=764*15*0.5/bins+764*15*ii/bins

    #for HCV pacbio data
    for ii in range(bins): binpos[ii]=L*G*0.5/bins+L*G*ii/bins

    print 'binpos', binpos

    for i in range(len(patients)):

        print (i+1)
        print ()
        tempMp1p2, rec_d, rec_m, dc,  nrec_d, = get_recomb_haplotypes(patients[i], bins, datepos)

        rec_diversification_dis = rec_diversification_dis+rec_d
        rec_missinghap_copies = rec_missinghap_copies+rec_m
        diversification_count = diversification_count+dc
        nonrec_diversification_dis = nonrec_diversification_dis+nrec_d
        Mp1p2.extend(tempMp1p2)

    print 'missing_haps', rec_missinghap_copies
    print 'diversfication_count', diversification_count
    print 'nondiversification_count', nonrec_diversification_dis

    '''normalize the histograms'''
    for ii in range(bins):
        if (diversification_count[ii]>0):
            rec_diversification_dis[ii]=rec_diversification_dis[ii]/diversification_count[ii]
            nonrec_diversification_dis[ii]=nonrec_diversification_dis[ii]/diversification_count[ii]
            rec_missinghap_copies[ii]=rec_missinghap_copies[ii]/diversification_count[ii]
            '''error bars for the mean values assuming everything is independent'''
            rec_err[ii]=(rec_diversification_dis[ii]/diversification_count[ii])**0.5
            nonrec_err[ii]=(nonrec_diversification_dis[ii]/diversification_count[ii])**0.5

    print 'normalised counts'
    print 'missing_haps', rec_missinghap_copies
    print 'diversfication_count', diversification_count
    print 'nondiversification_count', nonrec_diversification_dis

    '''plot mean haplotype observation frequency as a function of distance'''
    c = 2
    errorbar(binpos, rec_diversification_dis, yerr=rec_err, label="recombinant haplotype, cutoff "+str(c+1),  color='red', linewidth=2)#, cutoff "+str(c+1))
    errorbar(binpos, nonrec_diversification_dis, yerr=rec_err,label="other haplotypes, cutoff "+str(c+1),  color='green', linewidth=2)#, cutoff "+str(c+1))


    baseline=0
    baselinecount=0
    for ii in range(bins):
        baseline+=diversification_count[ii]*nonrec_diversification_dis[ii]
        baselinecount+=diversification_count[ii]

    print "<Mp1p2>: ", mean(Mp1p2)
    print "1-exp(-<Mp1p2>): ", 1-exp(-mean(Mp1p2))
    print "1-<exp(-Mp1p2)>: ", 1-mean(exp(-array(Mp1p2)))

    # print Mp1p2
    # print len(Mp1p2)


    baseline/=baselinecount

    print 'baseline', baseline
    plot([0, L], [baseline, baseline], label="weighted mean",  color='black', linewidth=1)
    legend(loc=2)
    xlabel('Distance [bp]', fontsize=fs)
    ylabel('P(haplotype present)', fontsize=fs)
    ylim([0,0.5])
    xlim([0,L])
    text(-50,0.54,'A', fontsize=fs+4)
    ax=gca()
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fs)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fs)

    figname='diversification_distance_c_'+str(cutoff)
    #if minimum_sample_size>0: figname+='_minsize_'+str(minimum_sample_size)
    #if delta_t>1: figname+='_dt_'+str(delta_t)
    savefig(figname+'.pdf')
    close()


    ''' repeat with distance x time dependence '''
    ''' calculate center bin values, 764*7 are magic numbers, 15 is the number generations per month '''
    binpos=zeros((bins))
    for ii in range(bins): binpos[ii]=L*sample_length*G*0.5/bins+L*sample_length*G*ii/bins
    print binpos



    ''' fit an exponential to the data '''
    fitfunc = lambda p, x: p[0]*(1-exp(-x*p[1]))
    fitfunc2 = lambda p, x: p[2]+p[0]*(1-exp(-x*p[1]))
    errfunc= lambda p, x, y, e:(fitfunc2(p,x)-y)/e
    '''fit using optimization routine from scipy'''
    p1,cov,info,m,success=optimize.leastsq(errfunc,[0.3, 0.0005, baseline],args=(binpos, rec_diversification_dis, rec_err), full_output=1)
    print p1

    errorbar(binpos, rec_diversification_dis, yerr=rec_err, label="recombinant haplotypes", color='red', linewidth=2) #, cutoff "+str(c+1))
    lstring="%(o).2f+%(n).2f(1-exp(-%(r).2e d $\Delta$t))" % { 'o':p1[2],'n': p1[0], 'r': p1[1]}
    print lstring
    print 'recombination estimate:', p1[0]*p1[1]/(mean(Mp1p2)-p1[2])

    rec_estimate1=p1[0]*p1[1]/(mean(Mp1p2)-p1[2])
    plot(binpos, fitfunc2(p1, binpos), label=lstring, color='black', linewidth=2)
    plot([0,7000], [p1[2],p1[2]], color='black')

    print "bp", max(binpos)
    legend(loc=2)
    xlabel('Distance x time [bp x generation]', fontsize=fs)
    ylabel('P(haplotype present)', fontsize=fs)
    ylim([0,0.5])
    xlim([0,7000])
    text(-50,0.53,'B', fontsize=fs+4)
    ax=gca()
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fs)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fs)

    figname='diversification_distance_time_c_'+str(cutoff)
    #if minimum_sample_size>0: figname+='_minsize_'+str(minimum_sample_size)
   # if delta_t>1: figname+='_dt_'+str(delta_t)
    savefig(figname+'.pdf')
    close()

    p1,cov,info,m,success=optimize.leastsq(errfunc,[0.3, 0.0005, baseline],args=(binpos, rec_missinghap_copies, rec_err), full_output=1)
    print p1

    lstring="%(o).2f+%(n).2f(1-exp(-%(r).2e d $\Delta$t))" % { 'o':p1[2],'n': p1[0], 'r': p1[1]}
    print lstring
    print 'recombination estimate:', p1[1]
    rec_estimate2=p1[1]

    plot(binpos, rec_missinghap_copies)
    print "bp", max(binpos)

    #plot(binpos, fitfunc2(p1, binpos), label=lstring, color='black', linewidth=2)
    xlabel('Distance x time [bp x generation]', fontsize=fs)
    ylabel('<N_AB>')
    figname='hap_copies_distance_time_c_'+str(cutoff)
    #if minimum_sample_size>0: figname+='_minsize_'+str(minimum_sample_size)
    #if delta_t>1: figname+='_dt_'+str(delta_t)
    savefig(figname+'.pdf')
    close()

    return rec_estimate1, rec_estimate2 #rec_diversification_dis, diversification_count, Mp1p2, rec_missinghap_copies, nonrec_diversification_dis

def haplotype_diversfication_timeavg_bootstrap(bins, patients, datepos, n_bootstrap):

    recomb_est_A = []
    recomb_est_B = []

    for i in xrange(n_bootstrap):

        print "bootstrap ", i+1
        bootstrap_sites = np.random.randint(L, size = L)

        rec_diversification_dis = np.zeros((bins))
        rec_missinghap_copies = np.zeros((bins))
        Mp1p2 = []
        diversification_count = np.zeros((bins))
        nonrec_diversification_dis = np.zeros((bins))

        #print len(rec_diversification_dis), len(rec_missinghap_copies), len(diversification_count)

        rec_err=np.zeros((bins))
        nonrec_err=np.zeros((bins))
        fs=16   #font size in plots
        ''' calculate center bin values, 764*7 are magic numbers, 15 is the number generations per month '''
        binpos=np.zeros((bins))

        #for HIV neher data
        #for ii in range(bins): binpos[ii]=764*15*0.5/bins+764*15*ii/bins

        #for HCV pacbio data
        for ii in range(bins): binpos[ii]=L*G*0.5/bins+L*G*ii/bins

        print 'binpos', binpos

        for i in range(len(patients)):

            # print (i+1)
            # print ()
            tempMp1p2, rec_d, rec_m, dc,  nrec_d, = get_recomb_haplotypes_bootstrap(patients[i], bins, datepos, bootstrap_sites)

            rec_diversification_dis = rec_diversification_dis+rec_d
            rec_missinghap_copies = rec_missinghap_copies+rec_m
            diversification_count = diversification_count+dc
            nonrec_diversification_dis = nonrec_diversification_dis+nrec_d
            Mp1p2.extend(tempMp1p2)

        # print 'missing_haps', rec_missinghap_copies
        # print 'diversfication_count', diversification_count
        # print 'nondiversification_count', nonrec_diversification_dis

        '''normalize the histograms'''
        for ii in range(bins):
            if (diversification_count[ii]>0):
                rec_diversification_dis[ii]=rec_diversification_dis[ii]/diversification_count[ii]
                nonrec_diversification_dis[ii]=nonrec_diversification_dis[ii]/diversification_count[ii]
                rec_missinghap_copies[ii]=rec_missinghap_copies[ii]/diversification_count[ii]
                '''error bars for the mean values assuming everything is independent'''
                rec_err[ii]=(rec_diversification_dis[ii]/diversification_count[ii])**0.5
                nonrec_err[ii]=(nonrec_diversification_dis[ii]/diversification_count[ii])**0.5

        # print 'normalised counts'
        # print 'missing_haps', rec_missinghap_copies
        # print 'diversfication_count', diversification_count
        # print 'nondiversification_count', nonrec_diversification_dis

        '''plot mean haplotype observation frequency as a function of distance'''
        c = 2
        #errorbar(binpos, rec_diversification_dis, yerr=rec_err, label="recombinant haplotype, cutoff "+str(c+1),  color='red', linewidth=2)#, cutoff "+str(c+1))
        #errorbar(binpos, nonrec_diversification_dis, yerr=rec_err,label="other haplotypes, cutoff "+str(c+1),  color='green', linewidth=2)#, cutoff "+str(c+1))


        baseline=0
        baselinecount=0
        for ii in range(bins):
            baseline+=diversification_count[ii]*nonrec_diversification_dis[ii]
            baselinecount+=diversification_count[ii]

        # print "<Mp1p2>: ", mean(Mp1p2)
        # print "1-exp(-<Mp1p2>): ", 1-exp(-mean(Mp1p2))
        # print "1-<exp(-Mp1p2)>: ", 1-mean(exp(-array(Mp1p2)))

        # print Mp1p2
        # print len(Mp1p2)


        baseline/=baselinecount

        #print 'baseline', baseline
        # plot([0, L], [baseline, baseline], label="weighted mean",  color='black', linewidth=1)
        # legend(loc=2)
        # xlabel('Distance [bp]', fontsize=fs)
        # ylabel('P(haplotype present)', fontsize=fs)
        # ylim([0,0.5])
        # xlim([0,L])
        # text(-50,0.54,'A', fontsize=fs+4)
        # ax=gca()
        # for tick in ax.xaxis.get_major_ticks():
        #     tick.label1.set_fontsize(fs)
        # for tick in ax.yaxis.get_major_ticks():
        #     tick.label1.set_fontsize(fs)
        #
        # figname='diversification_distance_c_'+str(cutoff)
        # #if minimum_sample_size>0: figname+='_minsize_'+str(minimum_sample_size)
        # #if delta_t>1: figname+='_dt_'+str(delta_t)
        # savefig(figname+'.pdf')
        # close()


        ''' repeat with distance x time dependence '''
        ''' calculate center bin values, 764*7 are magic numbers, 15 is the number generations per month '''
        binpos=zeros((bins))
        for ii in range(bins): binpos[ii]=L*sample_length*G*0.5/bins+L*sample_length*G*ii/bins
        #print binpos



        # rec_diversification_dis=zeros((bins))
        # nonrec_diversification_dis=zeros((bins))
        # missingfixed=zeros((bins))
        # diversificationcount=zeros((bins))
        #''' loop over patients '''
        # for patient_id, p in self.patients.iteritems():
        #     ''' calc mean crossover freq as a func distance x time: type 0 (last argument) '''
        #     rec_d, nrec_d, mf, dc, tempMp1p2, rec_copies = p.haplotype_diversification_timeavg(delta_t=delta_t,bins=bins,cutoff=c, distance=False, minimum_sample_size=minimum_sample_size)
        #     rec_diversification_dis=rec_diversification_dis+rec_d;
        #     nonrec_diversification_dis=nonrec_diversification_dis+nrec_d;
        #     diversificationcount=diversificationcount+dc;
        #     missingfixed=missingfixed+mf
        #     rec_missinghap_copies=rec_missinghap_copies+rec_copies

        # err=zeros((bins))
        # for ii in range(bins):
        #     if (diversificationcount[ii]>0):
        #         rec_diversification_dis[ii]=rec_diversification_dis[ii]/diversificationcount[ii]
        #         nonrec_diversification_dis[ii]=nonrec_diversification_dis[ii]/diversificationcount[ii]
        #         missingfixed[ii]=missingfixed[ii]/diversificationcount[ii]
        #         rec_missinghap_copies[ii]=rec_missinghap_copies[ii]/diversificationcount[ii]
        #         '''error bars for the mean values assuming everything is independent'''
        #         err[ii]=((rec_diversification_dis[ii]+1)/diversificationcount[ii])**0.5
        #     else: err[ii]=10000

        ''' fit an exponential to the data '''
        fitfunc = lambda p, x: p[0]*(1-exp(-x*p[1]))
        fitfunc2 = lambda p, x: p[2]+p[0]*(1-exp(-x*p[1]))
        errfunc= lambda p, x, y, e:(fitfunc2(p,x)-y)/e
        '''fit using optimization routine from scipy'''
        p1,cov,info,m,success=optimize.leastsq(errfunc,[0.3, 0.0005, baseline],args=(binpos, rec_diversification_dis, rec_err), full_output=1)
        #print p1

        # errorbar(binpos, rec_diversification_dis, yerr=rec_err, label="recombinant haplotypes", color='red', linewidth=2) #, cutoff "+str(c+1))
        # lstring="%(o).2f+%(n).2f(1-exp(-%(r).2e d $\Delta$t))" % { 'o':p1[2],'n': p1[0], 'r': p1[1]}
        # print lstring
        # print 'recombination estimate:', p1[0]*p1[1]/(mean(Mp1p2)-p1[2])

        rec_estimate1=p1[0]*p1[1]/(mean(Mp1p2)-p1[2])
       #  plot(binpos, fitfunc2(p1, binpos), label=lstring, color='black', linewidth=2)
       #  plot([0,7000], [p1[2],p1[2]], color='black')
       #
       #  print "bp", max(binpos)
       #  legend(loc=2)
       #  xlabel('Distance x time [bp x generation]', fontsize=fs)
       #  ylabel('P(haplotype present)', fontsize=fs)
       #  ylim([0,0.5])
       #  xlim([0,7000])
       #  text(-50,0.53,'B', fontsize=fs+4)
       #  ax=gca()
       #  for tick in ax.xaxis.get_major_ticks():
       #      tick.label1.set_fontsize(fs)
       #  for tick in ax.yaxis.get_major_ticks():
       #      tick.label1.set_fontsize(fs)
       #
       #  figname='diversification_distance_time_c_'+str(cutoff)
       #  #if minimum_sample_size>0: figname+='_minsize_'+str(minimum_sample_size)
       # # if delta_t>1: figname+='_dt_'+str(delta_t)
       #  savefig(figname+'.pdf')
       #  close()

        p1,cov,info,m,success=optimize.leastsq(errfunc,[0.3, 0.0005, baseline],args=(binpos, rec_missinghap_copies, rec_err), full_output=1)
        # print p1
        #
        # lstring="%(o).2f+%(n).2f(1-exp(-%(r).2e d $\Delta$t))" % { 'o':p1[2],'n': p1[0], 'r': p1[1]}
        # print lstring
        # print 'recombination estimate:', p1[1]
        rec_estimate2=p1[1]

        # plot(binpos, rec_missinghap_copies)
        # print "bp", max(binpos)
        #
        # #plot(binpos, fitfunc2(p1, binpos), label=lstring, color='black', linewidth=2)
        # xlabel('Distance x time [bp x generation]', fontsize=fs)
        # ylabel('<N_AB>')
        # figname='hap_copies_distance_time_c_'+str(cutoff)
        # #if minimum_sample_size>0: figname+='_minsize_'+str(minimum_sample_size)
        # #if delta_t>1: figname+='_dt_'+str(delta_t)
        # savefig(figname+'.pdf')
        # close()

        recomb_est_A.append(rec_estimate1)
        recomb_est_B.append(rec_estimate2)
        #return rec_estimate1, rec_estimate2 #rec_diversification_dis, diversification_count, Mp1p2, rec_missinghap_copies, nonrec_diversification_dis

    print recomb_est_A, np.mean(recomb_est_A), np.var(recomb_est_A)
    print recomb_est_B, np.mean(recomb_est_B), np.var(recomb_est_B)


    df = pd.DataFrame.from_items([('RecombA',recomb_est_A), ('RecombB', recomb_est_B)])
    df.to_csv("/Users/jayna/Dropbox/BEAST_Katrina/recombination_bootstap_gag.csv")

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]

def calcLD(haplotype_set, dd, n):

    freq_set = {}
    for each in haplotype_set:

        freq = float(dd.count(each))/float(n);
        freq_set[each] = freq

    #find haplotype with highest freq

    max = 0.0

    AB = ''
    Ab = ''
    aB = ''
    ab = ''


    for each in haplotype_set:

        if freq_set[each] > max:

            max = freq_set[each]
            AB = each


    temp = haplotype_set.copy()
    temp.remove(AB)
    for each in temp:

        Ab_t = AB[0]+each[1]
        aB_t = each[0]+AB[1]

        if Ab_t in temp and aB_t in temp:
            Ab = Ab_t
            aB = aB_t
            ab = each
            break



    p_AB = freq_set[AB]
    p_Ab = freq_set[Ab]
    p_aB = freq_set[aB]
    p_ab = freq_set[ab]

    p_A = p_AB+p_Ab
    p_B = p_AB+p_aB
    p_a = p_aB+p_ab
    p_b = p_Ab+p_ab


    p1q2 = p_A*p_b
    p2q1 = p_a*p_B
    p1q1 = p_A*p_B
    p2q2 = p_a*p_b


    D = p_AB - float(p_A)*float(p_B)
    D_prime = 0.0

    if D > 0:

        Dmax = min(p1q2,p2q1)
        D_prime = D/Dmax

    else:


        Dmin = min(p1q1,p2q2)

        D_prime = abs(D/Dmin)



    r = D/sqrt((p_A*(1.0-p_A))*(p_B*(1.0-p_B)))

    return [D_prime, r, p_AB, p_Ab, p_aB, p_ab]




subjects = ['~/AMC_HCV_DATA/fasta_nt/p_37_sub/',
            '~//AMC_HCV_DATA/fasta_nt/p_53_sub/',
            '~//AMC_HCV_DATA/fasta_nt/p_4_sub/',
            '~//AMC_HCV_DATA/fasta_nt/p_61_sub/'
             ]



haplotype_diversfication_timeavg(10, subjects, 1)





