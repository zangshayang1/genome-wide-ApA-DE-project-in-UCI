import time
import myio
import const
import summer


class RC2PCTerApp(summer.PARCSummerApp):
    def __init__(self, rcTable, dpoint, hasHeader):
        summer.PARCSummerApp.__init__(self, rcTable, hasHeader)
        self._dpoint = dpoint
        self._newheader = None

    def generate_pct_header(self):
        # self._header : str : "gene\tid\ttranscript\ttype\tcoords\tsample1\tsample2...."
        ls_samples = self._header.rstrip().split('\t')[const.DEFINED_REFERENCE_COLN:]
        ls_pct = []
        pct_header = ''
        for sample in ls_samples:
            ls_pct.append(sample + const.PCT_HEADER_SUFFIX)
        for pct in ls_pct:
            pct_header += '\t' + pct
        return pct_header

    def run(self):
        self._lgr.info("Start Converting Read Counts to Percentages ....")
        self._lgr.info("hasHeader? %s", self._hasHeader)

        if not self._hasHeader:
            self._cached_dict, self._header = self._cache_into_dict(linestart=0)
        else:
            self._cached_dict, self._header = self._cache_into_dict()

        self._lgr.info("converting to rcMaxtrix...")
        self._rcMatrix = self.convert2rcMatrix()
        self._lgr.info("summing up for each gene...")
        self._sumList_perGene = self.calc_sum_per_gene()
        self._newheader = self._header.rstrip() + self.generate_pct_header() + '\n'
        self._lgr.info("the output header set to be: %s", self._newheader)

        myresult = myio.RC2PCTerOutput(self._rcTable.name)
        myresult.add2content(self._newheader)

        self._lgr.info("calculating percentages...writing output ...")
        for gene in self._cached_dict:
            '''
            Behind each gene key, there is a apaList at axis=0 (v), a 2D-rcMatrix, and a sumList at axis=1 (h).
            Each time, a ApA got picked, rc data from rcMatrix will follow
            And pct will be calculated based on correspondent value in sumList
            This part is pretty neat.
            '''
            apaList = self._cached_dict[gene]
            rcMatrix = self._rcMatrix[gene]
            sumList = self._sumList_perGene[gene]
            sample_num = len(sumList)
            apa_num = len(apaList)
            for i in range(apa_num):
                myline = apaList[i].build_line().rstrip()
                for j in range(sample_num):
                    try:
                        pct = round(rcMatrix[i][j] / float(sumList[j]), self._dpoint)
                    except ZeroDivisionError:
                        pct = 'NA'
                    myline += '\t' + str(pct)
                myline += '\n'
                myresult.add2content(myline)

        myresult.open2write(myresult.content)
        self._lgr.info("Done!")


def run_rc2pcter(args):
    '''
    this RUN func should create a APP object and deliver all the input and params from args
    :param args:
    :return:
    '''
    time_start = time.time()

    app = RC2PCTerApp(args.PARC_Table,
                        args.decimal_point,
                        args.hasHeader)
    app.run()

    time_end = time.time()
    print 'Run time: {0:.2f} seconds'.format(time_end - time_start)
    # sys.stdout.flush()

