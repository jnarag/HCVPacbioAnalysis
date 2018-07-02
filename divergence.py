__author__ = 'jayna'

import numpy as np
import os

from matplotlib import pyplot as plt
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align import AlignInfo
from Bio.Alphabet import generic_dna

import pandas as pd
import statsmodels.formula.api as smf
from Bio.Seq import Seq

import matplotlib as mpl
import matplotlib.cm as cm

import random
import array
import re

''' Divergence over time - based on patristic distance from the founder sequence (consensus sequence
 from first time point '''

def divergence(fastain, translate, date_part, patient_id, sites):

    seqs_by_timepoint = split(fastain, date_part)

    mean_divergence = []
    median_divergence = []

    lower_divergence_25 = []
    upper_divergence_75 = []
    lower_divergence_5 = []
    upper_divergence_95 = []
    divergence_std = []
    mean_N_divergence = []
    median_N_divergence = []

    lower_N_divergence_25 = []
    upper_N_divergence_75 = []
    lower_N_divergence_5 = []
    upper_N_divergence_95 = []
    N_divergence_std = []
    mean_S_divergence = []
    median_S_divergence = []
    lower_S_divergence_25 = []
    upper_S_divergence_75 = []
    lower_S_divergence_5 = []
    upper_S_divergence_95 = []
    S_divergence_std = []
    dN = []
    dN_med = []
    dN_lower_25 = []
    dN_upper_75 = []
    dN_lower_5 = []
    dN_upper_95 = []
    dN_std = []
    dS = []
    dS_med = []
    dS_lower_25 = []
    dS_upper_75 = []
    dS_lower_5 = []
    dS_upper_95 = []
    dS_std = []
    patient = []

    # parts = str.split(fastain, "/")
    # parts2 = str.split(parts[len(parts)-1], "_")


    patient.append(patient_id)

    nonsyn_sites, syn_sites = number_of_N_and_S_sites(fastain, None)
    print nonsyn_sites, syn_sites

    sorted_timepoints = seqs_by_timepoint.keys()
    sorted_timepoints.sort(key=natural_keys)

    print sorted_timepoints
    first_timepoint = AlignIO.MultipleSeqAlignment(seqs_by_timepoint[sorted_timepoints[0]])

    consensus = AlignInfo.SummaryInfo(first_timepoint).dumb_consensus(threshold=0.01).upper()

    sampleTimes = []
    for t in sorted_timepoints:
        sampleTimes.append(float(t))


    #for f in filelist:
    for t in range(0,len(sorted_timepoints)):

        divergence = []
        divergence_N = []
        divergence_S = []
        divergence_dN = []
        divergence_dS = []
        # diff = 0


        seqs_at_t = seqs_by_timepoint[sorted_timepoints[t]]

        for each in seqs_at_t:

            parts = str.split(each.name, "_")
            freq = 1
            diff = 0
            diff_N = 0
            diff_S = 0

            seq = Seq(str(each.seq).upper().replace('-', 'N'))[sites[0]:sites[1]]

            codon_pos_start = 0
            codon_pos_end = 2

            A_i = str(seq).find('A')
            T_i = str(seq).find('T')
            G_i = str(seq).find('G')
            C_i = str(seq).find('C')

            start = [A_i, T_i, G_i, C_i]

            A_ii = str(seq).rfind('A')
            T_ii = str(seq).rfind('T')
            G_ii = str(seq).rfind('G')
            C_ii = str(seq).rfind('C')

            end = [A_ii, T_ii, G_ii, C_ii]

            start_i = min(start)
            end_i = max(end)

            if start_i > -1 and end_i > -1:
                # print start_i, end_i

                remainder_1 = start_i % 3
                remainder_2 = end_i % 3

                if remainder_1 != 0:
                    b = remainder_1 != codon_pos_start
                    # print start_i, start_i + (3-remainder_1)
                    start_i = start_i + (3 - remainder_1)

                if remainder_2 != 2:
                    # tprint end_i, end_i + (3-remainder_2)
                    end_i = end_i + (2 - remainder_2)

                seq = seq[start_i: end_i + 1]
                gaps = str(seq).count('N')

                seq_length = len(seq)
                aa_length = seq_length / 3



                conseq = Seq(str(consensus).replace('X', 'N'))[sites[0]:sites[1]]



                translated_seq = seq.translate()

                gaps_con = str(conseq).count('N')

                if gaps_con == seq_length:

                    print ("all gaps in conseq")
                    break

                else:
                    # if (seq_length >= length and (float(gaps) / float(seq_length)) < 0.05 and
                    #             (float(gaps_con) / float(seq_length)) < 0.05):

                        # print translated_seq, conseq.translate(),
                        count = 0
                        if (translate):
                            seq = each.seq.translate()

                        # count +=1

                        for a in range(seq_length):

                            i = a

                            if (str(conseq[i]) != "N"):

                                if (str(seq[i]) != "N"):

                                    count = count + 1

                                    if (conseq[i] != seq[i]):

                                        codon = []

                                        if (i % 3 == 0):
                                            cp = i
                                            cp_a = i + 1
                                            cp_b = i + 2

                                            codon = [cp, cp_a, cp_b]


                                        elif (i % 3 == 1):
                                            cp_a = i - 1
                                            cp = i
                                            cp_b = i + 1

                                            codon = [cp_a, cp, cp_b]

                                        else:

                                            cp_a = i - 2
                                            cp_b = i - 1
                                            cp = i

                                            codon = [cp_a, cp_b, cp]



                                        consensus_aa = conseq[codon[0]:(codon[2] + 1)].translate()
                                        current_aa = seq[codon[0]:(codon[2] + 1)].translate()

                                        #print(str(consensus_aa), str(current_aa))
                                        if 'X' in conseq[codon[0]:(codon[2] + 1)]:
                                            break

                                        if (str(consensus_aa) != str(current_aa)):

                                            diff_N += 1
                                        else:
                                            diff_S += 1

                                        diff += 1

                        # print diff/count, diff/seq_length


                        divergence.extend([float(diff) / float(count)]*freq)
                        divergence_N.extend([float(diff_N) / float(count)]*freq)
                        divergence_S.extend([float(diff_S) / float(count)]*freq)
                        divergence_dN.extend([float(diff_N) / float(nonsyn_sites)]*freq)
                        divergence_dS.extend([float(diff_S) / float(syn_sites)]*freq)


        if len(divergence) < 100:

            mean_divergence.append(float('nan'))
            median_divergence.append(float('nan'))
            lower_divergence_25.append(float('nan'))
            upper_divergence_75.append(float('nan'))
            lower_divergence_5.append(float('nan'))
            upper_divergence_95.append(float('nan'))
            divergence_std.append(float(1000))

            mean_N_divergence.append(float('nan'))
            median_N_divergence.append(float('nan'))
            lower_N_divergence_25.append(float('nan'))
            upper_N_divergence_75.append(float('nan'))
            lower_N_divergence_5.append(float('nan'))
            upper_N_divergence_95.append(float('nan'))
            N_divergence_std.append(float(1000))

            mean_S_divergence.append(float('nan'))
            median_S_divergence.append(float('nan'))
            lower_S_divergence_25.append(float('nan'))
            upper_S_divergence_75.append(float('nan'))
            lower_S_divergence_5.append(float('nan'))
            upper_S_divergence_95.append(float('nan'))
            S_divergence_std.append(float(1000))

            dN.append(float('nan'))
            dN_med.append(float('nan'))
            dN_lower_25.append(float('nan'))
            dN_upper_75.append(float('nan'))
            dN_lower_5.append(float('nan'))
            dN_upper_95.append(float('nan'))
            dN_std.append(float('nan'))

            dS.append(float('nan'))
            dS_med.append(float('nan'))
            dS_lower_25.append(float('nan'))
            dS_upper_75.append(float('nan'))
            dS_lower_5.append(float('nan'))
            dS_upper_95.append(float('nan'))
            dS_std.append(float(1000))


        else:

            #print divergence
            mean_divergence.append(np.mean(divergence))
            median_divergence.append(np.percentile(divergence, 50))
            lower_divergence_25.append(np.percentile(divergence, 25))
            upper_divergence_75.append(np.percentile(divergence, 75))
            lower_divergence_5.append(np.percentile(divergence, 5))
            upper_divergence_95.append(np.percentile(divergence, 95))
            divergence_std.append(np.std(divergence))

            mean_N_divergence.append(np.mean(divergence_N))
            median_N_divergence.append(np.percentile(divergence_N, 50))
            lower_N_divergence_25.append(np.percentile(divergence_N, 25))
            upper_N_divergence_75.append(np.percentile(divergence_N, 75))
            lower_N_divergence_5.append(np.percentile(divergence_N, 5))
            upper_N_divergence_95.append(np.percentile(divergence_N, 95))
            N_divergence_std.append(np.std(divergence_N))

            mean_S_divergence.append(np.mean(divergence_S))
            median_S_divergence.append(np.percentile(divergence_S, 50))
            lower_S_divergence_25.append(np.percentile(divergence_S, 25))
            upper_S_divergence_75.append(np.percentile(divergence_S, 75))
            lower_S_divergence_5.append(np.percentile(divergence_S, 5))
            upper_S_divergence_95.append(np.percentile(divergence_S, 95))
            S_divergence_std.append(np.std(divergence_S))

            dN.append(np.mean(divergence_dN))
            dN_med.append(np.percentile(divergence_dN, 50))
            dN_lower_25.append(np.percentile(divergence_dN, 25))
            dN_upper_75.append(np.percentile(divergence_dN, 75))
            dN_lower_5.append(np.percentile(divergence_dN, 5))
            dN_upper_95.append(np.percentile(divergence_dN, 95))
            dN_std.append(np.std(divergence_dN))

            dS.append(np.mean(divergence_dS))
            dS_med.append(np.percentile(divergence_dS, 50))
            dS_lower_25.append(np.percentile(divergence_dS, 25))
            dS_upper_75.append(np.percentile(divergence_dS, 75))
            dS_lower_5.append(np.percentile(divergence_dS, 5))
            dS_upper_95.append(np.percentile(divergence_dS, 95))
            dS_std.append(np.std(divergence_dS))


    df = pd.DataFrame.from_items(
        [('Times', sampleTimes), ('Divergence_median', median_divergence), ('Divergence_mean', mean_divergence),
         ('Divergence_N_med', median_N_divergence), ('Divergence_N_mean', mean_N_divergence),
         ('Divergence_S_med', median_S_divergence), ('Divergence_S_mean', mean_S_divergence),
         ('Divergence_lower_25', lower_divergence_25), ('Divergence_upper_75', upper_divergence_75),
         ('Divergence_lower_5', lower_divergence_5), ('Divergence_upper_95', upper_divergence_95),
         ("Divergence_std", divergence_std),
         ('Divergence_N_lower_25', lower_N_divergence_25), ('Divergence_N_upper_75', upper_N_divergence_75),
         ('Divergence_N_lower_5', lower_N_divergence_5), ('Divergence_N_upper_95', upper_N_divergence_95),
         ("Divergence_N_std", N_divergence_std),
         ('Divergence_S_lower_25', lower_S_divergence_25), ('Divergence_S_upper_75', upper_S_divergence_75),
         ('Divergence_S_lower_5', lower_S_divergence_5), ('Divergence_S_upper_95', upper_S_divergence_95),
         ("Divergence_S_std", S_divergence_std), ('Patients', patient*len(sampleTimes))])

    # csvfilename = filename.replace("filelist_", "divergence_results_by_birth_patristic_sites_"+str(sites_pos[0])+"_to_"+str(sites_pos[1])+"_")
    # csvfilename = csvfilename.replace(".txt",".csv")
    # df.to_csv(csvfilename)

    df = df.dropna()

    # print sampleTimes - (np.ones(len(sampleTimes))*sampleTimes[0])


    # print mean_divergence
    # print divergence_std

    if len(df) == 1:

        total_div = ['total', float('nan'), float('nan'), float('nan'), float('nan'), float('nan')]
        N_div = ['N', float('nan'), float('nan'), float('nan'), float('nan'), float('nan')]
        S_div = ['S', float('nan'), float('nan'), float('nan'), float('nan'), float('nan')]

    return df




