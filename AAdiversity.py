from Bio import AlignIO
import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO

import numpy as np



def plot_pairwise_diff(fastain, window_size, sliding_window):

    AA = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    aln = AlignIO.read('%s'%fastain, 'fasta')

    trans_aln = []

    for i in range(len(aln)):

        trans_aln.append(SeqIO.SeqRecord(Seq(str(aln[i].seq).replace('-','N')).translate()))

    trans_aln = AlignIO.MultipleSeqAlignment(trans_aln)

    seq_length = len(aln[1,:])

    n_windows = seq_length/window_size

    midpoint = []
    pwd = []
    raw_diff = []
    window_no = []

    count = 1
    end = window_size
    sliding_window_size = 100



    for each in range(n_windows):


        start = each * window_size
        end = (each + 1) * window_size


        print start, end, each+1
        sub_aln = trans_aln[:, each]

        aa_freq = []

        for a in AA:

            align_array = np.array(sub_aln, np.str)
            aa_freq.append(str(align_array).count(a))

        total_aa_freq = sum(aa_freq)

        diff = 0.0
        for i in range(0,len(aa_freq)-1):

            freq_i = float(aa_freq[i])/float(total_aa_freq)

            for j in range((i+1), len(aa_freq)):

                freq_j = float(aa_freq[j])/float(total_aa_freq-1)

                diff += freq_i*freq_j


        #print diff

        pwd.append(diff)
        window_no.append(each+1)


    return {"codon_no":window_no, "pwd":pwd}





#variable called window size controls how large the windows are. If window_size = 3, then you will estimate per codon pairwise diversity
window_size = 3

# sliding window is set to True
sliding_window = False

fastain1 = str('/Users/jayna/Dropbox/AMC_HCV_DATA/fasta_nt/All_patients_seq/FP7_patient_004_allseqs.fasta')
fastain2 = str('/Users/jayna/Dropbox/AMC_HCV_DATA/fasta_nt/All_patients_seq/FP7_patient_037_allseqs.fasta')
fastain3 = str('/Users/jayna/Dropbox/AMC_HCV_DATA/fasta_nt/All_patients_seq/FP7_patient_053.fasta')
fastain4 = str('/Users/jayna/Dropbox/AMC_HCV_DATA/fasta_nt/All_patients_seq/FP7_patient_061.fasta')

p4 = plot_pairwise_diff(fastain1, window_size, sliding_window)
p37 = plot_pairwise_diff(fastain2, window_size, sliding_window)
p53 = plot_pairwise_diff(fastain3, window_size, sliding_window)
p61 = plot_pairwise_diff(fastain4, window_size, sliding_window)

p4_d = pd.DataFrame(data=p4)
p37_d = pd.DataFrame(data=p37)
p53_d = pd.DataFrame(data=p53)
p61_d = pd.DataFrame(data=p61)

p4_d.to_csv("/Users/jayna/Dropbox/HCV_Pacbio_manuscript/figures/pairwise_diversity_codon_E1E2_p4_aa_new.csv")
p37_d.to_csv("/Users/jayna/Dropbox/HCV_Pacbio_manuscript/figures/pairwise_diversity_codon_E1E2_p37_aa_new.csv")
p53_d.to_csv("/Users/jayna/Dropbox/HCV_Pacbio_manuscript/figures/pairwise_diversity_codon_E1E2_p53_aa_new.csv")
p61_d.to_csv("/Users/jayna/Dropbox/HCV_Pacbio_manuscript/figures/pairwise_diversity_codon_E1E2_p61_aa_new.csv")

