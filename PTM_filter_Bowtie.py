# This script removes the reads generated from internal polyA priming events, thus false positive. 

import sys
import pysam


def get_name(ls_input):
    ls_name = []
    ls = open(ls_input, 'r')
    for line in ls:
        ls_name.append(line.rstrip('\n'))
    ls.close()
    return ls_name


def ptm_filter(infile, ref_genome_path, window_nbp):

    print "loading input ..."
    outfile = 'filtered_' + infile   # name of output
    # out_rvs = 'filtered_rvs_' + infile
    logfile = outfile[:-4] + '.log'   # name of log file
    ref_genome = pysam.FastaFile(ref_genome_path)   # load reference genome file

    # open files to read and write
    insam = open(infile, 'r')
    outsam = open(outfile, 'w')
    # rvssam = open(out_rvs, 'w')
    log = open(logfile, 'w')

    # init count variables
    unmapped_reads = 0
    filtered_reads = 0
    total_reads = 0
    remaining_reads = 0

    print "scanning input..."
    for line in insam:
        if line[0] == '@':
            outsam.write(line)
            #rvssam.write(line)
            continue  # skip and copy headers
        else:
            total_reads += 1
            ls_line = line.split()
            if not ls_line[2] == '*':  # check if the read is well-mapped
                interprimed = False
                count = 0   # init count variable for distinguish A-richness at downstream
                pos_p = int(ls_line[3]) + int(ls_line[5][:-1])  # POS in sam file specifies the 1-based left-most coordinates
                start_p = pos_p - 1                              # while ref_genome is 0-based
                end_r = int(ls_line[3]) - 1                      # so the end points for - strand differ from + strand

                chr_s = ls_line[2]
                strand_sig = ls_line[1]
                if int(strand_sig) == 0:  # if the read is mapped to the forward + strand
                    seq = ref_genome.fetch(chr_s, start_p, start_p+window_nbp)  # ref_genome.fetch only takes left end
                    if seq[:6] == 'AAAAAA':
                        interprimed = True
                    else:
                        for i in range(len(seq)):
                            if seq[i] == 'A':
                                count += 1
                elif int(strand_sig) == 16:  # if the read is mapped to the reverse - strand
                    if end_r-window_nbp >= 0:
                        seq = ref_genome.fetch(chr_s, end_r-window_nbp, end_r)
                        if seq[-6:] == 'TTTTTT':
                            interprimed = True
                        else:
                            for i in range(len(seq)):
                                if seq[i] == 'T':
                                    count += 1
                    #else:
                        #pass
                else:
                    print "there are other numbers except 0 or 16"

                if not interprimed and count < 7:
                    outsam.write(line)
                    remaining_reads += 1
                else:
                    filtered_reads += 1
            else:
                unmapped_reads += 1
                continue

    log.write("The process has gone through {0} reads in {1} \n".format(total_reads, infile))
    log.write("{0} unmapped reads have been removed \n".format(unmapped_reads))
    log.write("{0} reads have been filtered \n".format(filtered_reads))
    log.write("{0} reads remain in {1} \n".format(remaining_reads, outfile))
    log.write("----------------------------------------\n")

    insam.close()
    outsam.close()
    log.close()

mm9_genome_path = '/cbcl/szang/Chip-seq/mm9.RefGenome/NCBIM37_mm9.genome.fa'

ls_input = sys.argv[1]
ls_name = get_name(ls_input)
for name in ls_name:
    infile = name
    ptm_filter(infile, mm9_genome_path, 10)  # 10 is downstream window size
