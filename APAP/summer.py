import myclasses
import time
import logger
import myio
import const
import numpy as np


class PARCSummerApp(myclasses.RCTableCacher):
    def __init__(self, rcTable, hasHeader):
        myclasses.RCTableCacher.__init__(self, rcTable, hasHeader)
        self._rcMatrix = None
        self._sumList_perGene = None
        self._sumList_perApA = None
        self._lgr = logger.set_logger(self._rcTable.name[:-4] + '.log')

    def _generate_GeneOnly_header(self):
        ls_header = self._header.rstrip().split('\t')
        # remove the first five colnames: GENES, IDs, TRANSCRIPTS, TYPES, COORDS
        ls_samples = ls_header[const.DEFINED_REFERENCE_COLN:]
        go_header = 'GENES'
        # stitch GENES to start, so I can accumulate go_header this way
        for colname in ls_samples:
            go_header += '\t' + colname
        go_header += '\n'
        return go_header

    def convert2rcMatrix(self):
        '''
        It functions on self._cached_dict and outputs a matrix in a shape of (num_apas, num_samples) for each gene in a dict
        '''
        rc_ary_per_gene = {}
        for gene in self._cached_dict:
            ary = []
            ls_apa = self._cached_dict[gene]
            for apa in ls_apa:
                ary.append(apa.readcounts())
            rc_ary_per_gene[gene] = ary
        return rc_ary_per_gene

    def calc_sum_per_gene(self):
        sum_list_per_gene = {}
        for gene in self._rcMatrix:
            rc_matrix = np.array(self._rcMatrix[gene])
            sum_matrix = rc_matrix.sum(0)
            sum_list_per_gene[gene] = sum_matrix
        return sum_list_per_gene

    def calc_sum_per_ApA(self):
        sum_list_per_apa = {}
        for gene in self._rcMatrix:
            rc_matrix = np.array(self._rcMatrix[gene])
            sum_matrix = rc_matrix.sum(1)
            sum_list_per_apa[gene] = sum_matrix
        return sum_list_per_apa

    def run(self):
        self._lgr.info("Start Summing Up Read Counts for Each Gene...")
        self._lgr.info("hasHeader? %s", self._hasHeader)

        if not self._hasHeader:
            self._cached_dict, self._header = self._cache_into_dict(linestart=0)
        else:
            self._cached_dict, self._header = self._cache_into_dict()
            # the default linestart = 1, hasPCT = False

        self._lgr.info("converting to rcMaxtrix...")
        self._rcMatrix = self.convert2rcMatrix()
        self._lgr.info("summing up for each gene...")
        self._sumList_perGene = self.calc_sum_per_gene()

        myresult = myio.RCSummerOutput(self._rcTable.name)

        self._header = self._generate_GeneOnly_header()
        myresult.add2content(self._header)
        self._lgr.info("the output header set to be: %s", self._header)

        self._lgr.info("writing output ...")
        for gene in self._sumList_perGene:
            sum_line = gene
            ls = self._sumList_perGene[gene]
            for s in ls:
                sum_line += '\t' + str(s)
            sum_line += '\n'
            myresult.add2content(sum_line)
        myresult.open2write(myresult.content)
        self._lgr.info("Done!")


def run_summer(args):
    '''
    this RUN func should create a APP object and deliver all the input and params from args
    :param args:
    :return:
    '''
    time_start = time.time()

    app = PARCSummerApp(args.PARC_Table,
                        args.hasHeader)
    app.run()

    time_end = time.time()
    print 'Run time: {0:.2f} seconds'.format(time_end - time_start)
    # sys.stdout.flush()

