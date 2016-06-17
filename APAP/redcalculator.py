import top2selecter
import time
import utils
import myio
import const


class REDcalculatorApp(top2selecter.Top2SelectionApp):
    def __init__(self, rcTable, dpoint, hasHeader, hasPCT):
        top2selecter.Top2SelectionApp.__init__(self, rcTable, hasHeader, hasPCT)
        self._dpoint = dpoint
        self._red_dict = None
        self._redHeader = None

    def _define_and_calculate_RED(self):
        dict_red = {}
        for gene in self._top2_dict:
            ls_apa = self._top2_dict[gene]
            if len(ls_apa) > 2:
                raise Exception('The input has not gone through top2 filter.')
            elif len(ls_apa) <= 1:
                raise Exception('There is either single ApA or empty gene without any APAs cached.')
            else:
                apa1, apa2 = ls_apa[:]
                ls_red = []
                if utils.first_is_approximal(apa1.coords(), apa2.coords()) is None:
                    continue
                elif utils.first_is_approximal(apa1.coords(), apa2.coords()):
                    for approximal_pA, distal_pA in zip(apa1.readcounts(), apa2.readcounts()):
                        red = round(utils.calculate_red_score(approximal_pA, distal_pA), self._dpoint)
                        ls_red.append(red)
                else:
                    for distal_pA, approximal_pA in zip(apa1.readcounts(), apa2.readcounts()):
                        red = round(utils.calculate_red_score(approximal_pA, distal_pA), self._dpoint)
                        ls_red.append(red)
                dict_red[gene] = ls_red
        return dict_red

    def run(self):
        self._lgr.info("Start Calculating RED...")
        self._lgr.info("hasHeader? %s", self._hasHeader)
        self._lgr.info("hasPCT? %s", self._hasPCT)

        if self._hasHeader and self._hasPCT:
            self._cached_dict, self._header = self._cache_into_dict(hasPCT=True)
        elif self._hasHeader and not self._hasPCT:
            self._cached_dict, self._header = self._cache_into_dict() # default: linestart = 1, hasPCT = False
        elif not self._hasHeader and self._hasPCT:
            self._cached_dict, self._header = self._cache_into_dict(linestart=0, hasPCT=True)
        elif not self._hasHeader and not self._hasPCT:
            self._cached_dict, self._header = self._cache_into_dict(linestart=0)
        else:
            raise Exception('The read count table is in wrong format.')

        self._lgr.info("converting to rcMaxtrix...")
        self._rcMatrix = self.convert2rcMatrix()
        self._lgr.info("summing up for each ApA site...")
        self._sumList_perApA = self.calc_sum_per_ApA()
        self._lgr.info("selecting top2 ApAs...")
        self._top2_dict, self._singles = self._select_top2()
        self._lgr.info("calculating RED per gene...")
        self._red_dict = self._define_and_calculate_RED()
        self._header = const.RED_OUTPUT_HEADER
        self._lgr.info("the output header set to be: %s", self._header)

        myresult = myio.REDtableOutput(self._rcTable.name)
        myresult.add2content(self._header)

        self._lgr.info("writing output...")
        for gene in self._red_dict:
            ls_red = self._red_dict[gene]
            myline = str(gene)
            for red in ls_red:
                myline += '\t' + str(red)
            myline += '\n'
            myresult.add2content(myline)
        myresult.open2write(myresult.content)
        self._lgr.info("Done!")


def run_REDcalculator(args):
    '''
    this RUN func should create a APP object and deliver all the input and params from args
    :param args:
    :return:
    '''
    time_start = time.time()

    app = REDcalculatorApp(args.PARC_Table,
                            args.decimal_point,
                            args.hasHeader,
                            args.hasPCT)
    app.run()

    time_end = time.time()
    print 'Run time: {0:.2f} seconds'.format(time_end - time_start)
    # sys.stdout.flush()

