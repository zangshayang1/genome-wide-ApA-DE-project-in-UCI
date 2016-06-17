#!/usr/bin/env python
import argparse
from APAP.filter import run_filter
from APAP.counter import run_counter
from APAP.merger import run_merger
from APAP.cleaner import run_cleaner
from APAP.summer import run_summer
from APAP.rc2pcter import run_rc2pcter
from APAP.redcalculator import run_REDcalculator
from APAP.top2selecter import run_top2selecter
from APAP.hister import run_hister


parser = argparse.ArgumentParser(prog='apap')
subparsers = parser.add_subparsers()

#======================================
# add IMfilter subcommand
#======================================

parser_IMfilter = subparsers.add_parser('IMfilter',
                                        help='IMfilter removes internally mis-primed reads')

parser_IMfilter.add_argument('reference_genome',
                             help='FASTA file as reference genome')

parser_IMfilter.add_argument('infile_sam',
                             help='SAM file as input')

parser_IMfilter.add_argument('--window_size', default=10, type=int,
                             help='Please specify a positive number as the size of look-up window to \
                                  justify internally mis-primed event')

parser_IMfilter.add_argument('--deny_number', default=6, type=int,
                             help='Please specify a positive number less than look-up window size \
                                  as the threshold to deny internally mis-primed event \
                                  \t i.e. Default is 6. If downstream sequence includes 6 consecutive As or \
                                  more than 6 As are found within the look-up window, the read will be removed')

parser_IMfilter.set_defaults(func=run_filter)


'''
different versions of genome references can vary widely
'''



#======================================
# add PARcounter subcommand
#======================================

parser_PARcounter = subparsers.add_parser('PARcounter',
                                        help='PARcounter increment counts of reads at given APA positions.')

parser_PARcounter.add_argument('masterlist', type=str,
                             help='a master list of ApA sites')

parser_PARcounter.add_argument('infile_bed', type=str,
                             help='BED file as input')

parser_PARcounter.add_argument('sample_name', type=str,
                             help='Please specify a positive number as the size of the window to \
                                  include reads as it hits the corresponding ApA site.')

parser_PARcounter.add_argument('--window_size', default=40, type=int,
                             help='Please specify a positive number as the size of the window to \
                                  include reads as it hits the corresponding ApA site.')

parser_PARcounter.set_defaults(func=run_counter)



#======================================
# add another MergeAPAL
#======================================


parser_MergeAPAL = subparsers.add_parser('MergeAPAL',
                                        help='MergeAPAL merges a list of results from PARcounter into one read counts (RC) table.')

parser_MergeAPAL.add_argument('PARC_list',
                              help='Please provide a file containing the name of all the inputs, polyA read counts lists.')

parser_MergeAPAL.add_argument('--undefined', default=False,
                              help='The default setting is to merge PARcounter results from well-defined masterlist\
                                   If you want to merge un-identified ApA read counts from PARcounter, you need to\
                                   specify this parameter so it automatically turns to be True upon set.',
                              action='store_true')


parser_MergeAPAL.add_argument('--outputname', default='ApA_RC.txt',
                              help='Please specify a output name.')

parser_MergeAPAL.set_defaults(func=run_merger)


#======================================
# add another MLCleaner
#======================================


parser_MLCleaner = subparsers.add_parser('MLCleaner',
                                        help='MLCleaner removes the ApAs that:\
                                        1. are "upstreamAntisense";\
                                        2. are not associated with any geneNames;\
                                        3. are within a certain range of each other;\
                                        4. are repeated.')

parser_MLCleaner.add_argument('masterlist',
                              help='Please put in the masterlist that you want to clean.')

parser_MLCleaner.add_argument('--window_size', default=40, type=int,
                              help='You can specify a window size. So if there are two or more ApA sites found within this window, \
                                   only the one with the highest read counts across samples remain.')

parser_MLCleaner.add_argument('--hasHeader', default=True, type=bool,
                              help='You can specify if the masterlist has header or not.\
                                   If not, default is True.')

parser_MLCleaner.add_argument('--delrc', default=False,
                              help='Please specify this parameter if you want to clean off the read counts data.',
                              action='store_true')

parser_MLCleaner.set_defaults(func=run_cleaner)



#======================================
# add another PARCSummer
#======================================

parser_PARCSummer = subparsers.add_parser('PARCSummer',
                                        help='PARCSummer calculates the sums of RC values of the same gene.')

parser_PARCSummer.add_argument('PARC_Table',
                             help='Please provide a standard ApA read counts table.')

parser_PARCSummer.add_argument('--hasHeader', default=True, type=bool,
                              help='You can specify if the masterlist has header or not.\
                                   If not, default is True.')


parser_PARCSummer.set_defaults(func=run_summer)



#======================================
# add another rc2pct
#======================================

parser_rc2pcter = subparsers.add_parser('rc2pct',
                                        help='rc2pct converts read counts to percentages \
                                        using the sums of RC values of the same gene.')

parser_rc2pcter.add_argument('PARC_Table',
                             help='Please provide a standard ApA read counts table.')

parser_rc2pcter.add_argument('--decimal_point', default=4, type=int,
                           help='You can specify some decimal point that you would like percentages to be rounded\
                            about. Default is 4.')

parser_rc2pcter.add_argument('--hasHeader', default=True, type=bool,
                              help='You can specify if the masterlist has header or not.\
                                   If not, default is True.')


parser_rc2pcter.set_defaults(func=run_rc2pcter)


#======================================
# add another calcRED
#======================================

parser_calcRED = subparsers.add_parser('calcRED',
                                        help='calcRED defines the approximal and distal pAs from the top2 selected\
                                         read counts table and generates a RED score table.')

parser_calcRED.add_argument('PARC_Table',
                             help='Please provide a RC table in a standard format, that contains no greater\
                              than 2 ApA sites per gene.')

parser_calcRED.add_argument('--decimal_point', default=6, type=int,
                           help='You can specify some decimal point that you would like RED scores to be rounded about.')

parser_calcRED.add_argument('--hasHeader', default=True, type=bool,
                              help='You can specify if the RCtable has header or not.\
                                   Default is True.')

parser_calcRED.add_argument('--hasPCT', default=True, type=bool,
                              help='You can specify if the RCtable includes percentage info or not.\
                                   Default is True.')


parser_calcRED.set_defaults(func=run_REDcalculator)


#======================================
# add another Top2selecter
#======================================

parser_top2selecter = subparsers.add_parser('select2',
                                             help='This function selects the top2 enriched ApAs per gene from a standard\
                                             RCtable')

parser_top2selecter.add_argument('PARC_Table',
                                help='Please provide a RC table in a standard format, that contains no greater\
                                than 2 ApA sites per gene.')

parser_top2selecter.add_argument('--hasHeader', default=True, type=bool,
                                help='You can specify if the RCtable has header or not.\
                                Default is True.')

parser_top2selecter.add_argument('--hasPCT', default=True, type=bool,
                                 help='You can specify if the RCtable includes percentage info or not.\
                                 Default is True.')


parser_top2selecter.set_defaults(func=run_top2selecter)



#======================================
# add another hister
#======================================
parser_hister = subparsers.add_parser('hister',
                                      help='This function takes a labeled ApA logFC lists and\
                                            outputs three types of data:\
                                            1. sliced fasta sequences of a certain length centered at ApA sites\
                                                for the purpose of generating weblogos.\
                                            2. the counting data of a certain short sequence pattern appearing\
                                                at the each position of the slicing window, for further plotting purpose\
                                            3. the preliminary plots of the counting data.')
parser_hister.add_argument('in_list',
                           help='Please provide a filtered and labeled (distal/approx) ApA logFC list.\
                                if the provided list does not contain header, the first record will be ignored.')

parser_hister.add_argument('reference_genome',
                           help='Please provide a reference genome.')

parser_hister.add_argument('pattern',
                           help='Please provide the sequence pattern that you are interested.\
                                such as TATG or AA/TTAAA.')

parser_hister.add_argument('output_name',
                           help='Please provide a naming pattern that will be followed by all the outputs with\
                           appropriate suffixes according to their types.')

parser_hister.add_argument('--window_size', default = 200, type=int,
                           help='You can specify the size of the slicing window centered at ApAs\
                           (aka the length of the sequence segment in fasta outputs). The default is 200.')

parser_top2selecter.set_defaults(func=run_hister)

#======================================
# run
#======================================
if __name__ == '__main__':
    args = parser.parse_args()
    args.func(args)



