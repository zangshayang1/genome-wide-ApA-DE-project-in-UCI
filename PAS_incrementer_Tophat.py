# this script is the PAS_incrementer.py "Tophat" version.
# Tophat brings the output a little different from what came out of Bowtie.

import sys


def get_name(ls_file):
    ls_name = []
    ls = open(ls_file, 'r')
    for line in ls:
        ls_name.append(line.rstrip('\n'))
    ls.close()
    return ls_name


def masterlist_increment(masterlist, dict_coordinates, output_name):   # take a masterlist and a pre-cached dictionary as inputs
    ref = open(masterlist, 'r')
    outfile = open(output_name, 'w')
    print "scanning APA-masterlist..."

    hits = 0

    for line in ref:
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
            hits += 1
        if hits % 100 == 0:
            print "got {0} hit".format(hits)
        outfile.write(ls_line[0]+'\t'+ls_line[1]+'\t'+ls_line[2]+'\t'+ls_line[3]+'\t'+pas+'\t'+str(incremental)+'\n')
    print "completed."
    return hits


def count_bed_reads(filename):
    bed = open(filename, 'r')
    print "loading file..."

    dict_coordinates = {}
    print "collecting reads coordinates..."
    for line in bed:
        ls_line = line.split()
        if ls_line[5] == "+":
            coords = ls_line[0] + ls_line[5] + ls_line[2]    # the PAS on forward strand such as "chr+rightend"
        else:
            coords = ls_line[0] + ls_line[5] + ls_line[1]    # the PAS on reverse strand such as "chr-leftend"

        if dict_coordinates.get(coords) is None:
            dict_coordinates[coords] = 1
        else:
            dict_coordinates[coords] += 1
    print "reads info cached."

    bed.close()
    return dict_coordinates


def main(ls_input):
    ref_path = '/cbcl/szang/Chip-seq/3_PASseq/postBowtie_process/PAS_hits/tian_PAS_data.txt'
    ls_name = get_name(ls_input)
    for bed_name in ls_name:
        out_name = bed_name[:-4]+'_PAS_list.txt'  # name of output
        log_name = bed_name[:-4]+'.log'           # name of log (the input name ends with ".bed")
        dict_coords = count_bed_reads(bed_name)   # cache the coords from bed
        hits = masterlist_increment(ref_path, dict_coords, out_name)    # write output and return hits
        log = open(log_name, 'a')                 # write log
        log.write('{0} hits were found.\n'.format(hits))
        log.write("----------------------------------------\n")
        log.close()


###########
# EXECUTE
###########
main(sys.argv[1])
