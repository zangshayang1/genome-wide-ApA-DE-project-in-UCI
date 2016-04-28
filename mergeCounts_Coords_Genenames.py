# this script would combine the outputs from PAS_incrementer.py into one .txt file 
# that includes counting data from multiple samples for each row of ApA sites.


import sys


def get_inputs(ls_filename):
  ls_file = open(ls_filename, 'r')
  ls_inputs = []
  for line in ls_file:
    ls_inputs.append(line.rstrip('\n'))
  ls_file.close()
  return ls_inputs


def list_counts(filename):
  infile = open(filename, 'r')
  ls_counts = []
  for line in infile:
    ls_line = line.rstrip().split('\t')
    ls_counts.append(ls_line[-1])
  infile.close()
  return ls_counts

###########################

# you can specify pathway for reference and output
# you can manipulate sys.argvs here

###########################


def main():
  output_name = sys.argv[2]
  output_path = '/cbcl/szang/Chip-seq/3_PASseq/postBowtie_process/PAS_hits/'
  out = open(output_path + output_name, 'w')
  
  reference = '/cbcl/szang/Chip-seq/3_PASseq/postBowtie_process/PAS_hits/tian_PAS_data.txt'
  ref = open(reference, 'r')
  ls_coords = []
  ls_names = []
  for line in ref:
    ls_line = line.rstrip().split('\t')
    ls_names.append(ls_line[0])
    ls_coords.append(ls_line[4])
  ref.close()
 
  ls_inputs = get_inputs(sys.argv[1])
  ary_counts = []
  for filename in ls_inputs:
    ls_counts = list_counts(filename)  
    ary_counts.append(ls_counts)

  N = len(ls_names)
  M = len(ary_counts)
  for i in range(N):
    out.write(ls_names[i] + '\t')
    out.write(ls_coords[i])
    for j in range(M):
      if ary_counts[j][i].isdigit():
        out.write('\t' + ary_counts[j][i])
      else:
        raise ValueError('non-numeric element detected!')
        break
    out.write('\n')
  out.close()

  
#################################
# Execute
#################################
main()
