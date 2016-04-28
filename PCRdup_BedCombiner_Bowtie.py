# this scripts removes the duplicated reads from the bed file.
# duplicates have the same coordinates and same quatitative barcode.

import sys


def get_name(ls_input):
    ls_name = []
    ls = open(ls_input, 'r')
    for line in ls:
        ls_name.append(line.rstrip('\n'))
    ls.close()
    return ls_name


def pcrdup_combine(infile):

    outfile = infile[:-4] + '_NoDup.bed'  # name of output
    outlog = infile[:-4] + '.log'    # name of log

    inbed = open(infile, 'r')
    outbed = open(outfile, 'w')
    log = open(outlog, 'a')   # let's keep log name unchanged

    templist = ['chr', 'start', 'end', 'barcode']  # init a list to store info of last line
    dup_count = 0
    total_reads = 0
    remaining_reads = 0

    for line in inbed:
        total_reads += 1
        ls_line = line.split()   # bed file is sorted so we just need to compare current line with the previous one
        if ls_line[0]==templist[0] and ls_line[1]==templist[1] and ls_line[2]==templist[2] and ls_line[3][-4:]==templist[3]:
            dup_count += 1
            continue
        else:
            templist[0] = ls_line[0]
            templist[1] = ls_line[1]
            templist[2] = ls_line[2]
            templist[3] = ls_line[3][-4:]
            outbed.write(line)
            remaining_reads += 1

    log.write("The process has gone through {0} reads in {1} \n".format(total_reads, infile))
    log.write("{0} PCR duplicates have been removed \n".format(dup_count))
    log.write("{0} reads remain in {1} \n".format(remaining_reads, outfile))
    log.write("----------------------------------------\n")

    inbed.close()
    outbed.close()
    log.close()


ls_input = sys.argv[1]
ls_name = get_name(ls_input)
for name in ls_name:
    infile = name
    pcrdup_combine(infile)
