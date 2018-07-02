from Bio import AlignIO
import pandas as pd
import random as r
from Bio.Seq import Seq



def plot_pairwise_diff(fastain, window_size, sliding_window):


    aln = AlignIO.read('%s'%fastain, 'fasta')

    indices = range(len(aln))

    r.shuffle(indices)

    random_array = []

    for i in range(500):

        random_array.append(aln[indices[i]])

    random_aln = AlignIO.MultipleSeqAlignment(random_array)


    seq_length = len(random_aln[1,:])

    n_windows = seq_length/window_size

    midpoint = []
    pwd = []
    raw_diff = []
    window_no = []

    count = 1
    end = window_size
    sliding_window_size = 100



    if(sliding_window):
        while end < seq_length:

            #print window_size
            start = (count-1)*sliding_window_size

            print count, start, end, end-start

            if(start >= end):
                break

            midpoint.append(float((start+end)/2))
            sub_aln = aln[:,start:end]
            diff = 0
            x = 0
            for s in range(len(sub_aln)-1):

                seq1 = sub_aln[s,:].seq
                seq2 = sub_aln[(s+1),:].seq

                for b in range(len(seq1)):

                    if(seq1[b] != '-'):

                        if(seq2[b] != '-'):

                            x+=1
                            if(str(seq1[b]).upper()!=str(seq2[b]).upper()):

                                diff+=1

            print "count", x, float(len(sub_aln)), float(x)/float(len(sub_aln)*window_size), float(diff)/float(x)

            pwd_i = float(diff)/float(x)#float(end-start)/float(len(sub_aln))
            pwd.append(pwd_i)
            raw_diff.append(float(diff) / float(len(sub_aln)))
            end += sliding_window_size
            count += 1



    else:


        for each in range(n_windows):

            diff = 0
            x = 0


            start = each * window_size
            end = (each + 1) * window_size


            print start, end, each+1
            sub_aln = random_aln[:, start:end]
            midpoint.append(float((start + end + 2) / 2))


            for j in range(len(sub_aln)):

                for k in range(len(sub_aln)-1):

                    seq1 = sub_aln[j, :].seq
                    seq2 = sub_aln[(k + 1), :].seq

                    seq1 = str(seq1).replace('-','N')
                    seq2 = str(seq2).replace('-','N')


                    seq1_aa = Seq(seq1).translate()
                    seq2_aa = Seq(seq2).translate()

                    #print seq1_aa, seq2_aa

                    for b in range(len(seq1_aa)):

                        if (seq1_aa[b] != 'X'):

                            if (seq2_aa[b] != 'X'):

                                x += 1
                                if (seq1_aa[b] != seq2_aa[b]):
                                    diff += 1

            #print x, diff

            #print (each+1), float(diff)/float(end-start)/float(len(sub_aln))

            if x>0:

                print each, float(diff)/x, x
                pwd_i = float(diff)/x
            else:
                pwd_i = float('nan')

            pwd.append(pwd_i)
            raw_diff.append(float(diff)/float(len(sub_aln)))
            window_no.append(each+1)



    print len(pwd), len(raw_diff)

    return {"codon_no":window_no, "mid": midpoint, "pwd":pwd}

#variable called window size, controls how large the windows are
window_size = 3

# sliding window is set to True
sliding_window = False

fastain1 = str('~/AMC_HCV_DATA/fasta_nt/All_patients_seq/FP7_patient_004_allseqs.fasta')
fastain2 = str('~/AMC_HCV_DATA/fasta_nt/All_patients_seq/FP7_patient_037_allseqs.fasta')
fastain3 = str('~/AMC_HCV_DATA/fasta_nt/All_patients_seq/FP7_patient_053.fasta')
fastain4 = str('~/AMC_HCV_DATA/fasta_nt/All_patients_seq/FP7_patient_061.fasta')

p4 = plot_pairwise_diff(fastain1, window_size, sliding_window)
p37 = plot_pairwise_diff(fastain2, window_size, sliding_window)
p53 = plot_pairwise_diff(fastain3, window_size, sliding_window)
p61 = plot_pairwise_diff(fastain4, window_size, sliding_window)

p4_d = pd.DataFrame(data=p4)
p37_d = pd.DataFrame(data=p37)
p53_d = pd.DataFrame(data=p53)
p61_d = pd.DataFrame(data=p61)

p4_d.to_csv("~/figures/pairwise_diversity_codon_E1E2_p4_aa.csv")
p37_d.to_csv("~/figures/pairwise_diversity_codon_E1E2_p37_aa.csv")
p53_d.to_csv("~/figures/pairwise_diversity_codon_E1E2_p53_aa.csv")
p61_d.to_csv("~/figures/pairwise_diversity_codon_E1E2_p61_aa.csv")

