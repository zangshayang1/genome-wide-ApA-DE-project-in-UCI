# This script increments the count of the reads that fall into a [-40, 40] range around a ApA site
# with a ApA-sites masterlist as reference. 

import sys


def get_name(ls_file):
    ls_name = []
    ls = open(ls_file, 'r')
    for line in ls:
        ls_name.append(line.rstrip('\n'))
    ls.close()
    return ls_name


def pas_incrememt(in_bed, masterlist):

    print "loading file..."
    argv3 = in_bed[:-4]+'_PAS_list.txt'  # name of output
    argv4 = in_bed[:-10]+'.log'  # name of log (the input name ends with "_NoDup.bed")

    infile = open(in_bed, 'r')
    reffile = open(masterlist, 'r')
    outfile = open(argv3, 'w')
    log = open(argv4, 'a')


    dict_coordinates = {}
    hits = 0

    print "collecting reads coordinates..."
    for line in infile:
        ls_line = line.split()
        if ls_line[5] == '+':
            coords = ls_line[0]+ls_line[5]+ls_line[2]    # the PAS on forward strand
        else:
            coords = ls_line[0]+ls_line[5]+ls_line[1]    # the PAS on reverse strand

        if dict_coordinates.get(coords) is None:
            dict_coordinates[coords] = 1
        else:
            dict_coordinates[coords] += 1


    print "reads info cached."

    print "scanning PAS-list..."
    for line in reffile:
        incremental = 0
        ls_line = line.split()
        pas = ls_line[4]

        negstr_idx = pas.find('-')
        posstr_idx = pas.find('+')
        strand_idx = max(negstr_idx, posstr_idx)
        strand_sig = pas[strand_idx]

        chr_s = pas[:strand_idx]
        pos = int(pas[strand_idx+1:])
        lb = pos - 40
        hb = pos + 40

        for p_point in range(lb, hb):
            p_key = chr_s+strand_sig+str(p_point)
            if dict_coordinates.get(p_key) is None:
                continue
            else:
                incremental += dict_coordinates[p_key]

        if incremental > 0:
        #    incremental += 1
            hits += 1
        if hits % 100 == 0:
            print "got {0} hit".format(hits)
        outfile.write(ls_line[0]+'\t'+ls_line[1]+'\t'+ls_line[2]+'\t'+ls_line[3]+'\t'+pas+'\t'+str(incremental)+'\n')

    print "storing PAS info..."
    print "completed."

    log.write('{0} hits were found.\n'.format(hits))
    log.write("----------------------------------------\n")

    infile.close()
    reffile.close()
    outfile.close()
    log.close()


ref_path = '/cbcl/szang/Chip-seq/3_PASseq/postBowtie_process/PAS_hits/tian_PAS_data.txt'

ls_input = sys.argv[1]
ls_name = get_name(ls_input)
for name in ls_name:
    pas_incrememt(name, ref_path)
