import time
import myio
import summer


class Top2SelectionApp(summer.PARCSummerApp):
    def __init__(self, rcTable, hasHeader, hasPCT):
        summer.PARCSummerApp.__init__(self, rcTable, hasHeader)
        self._hasPCT = hasPCT
        self._singles = None
        self._top2_dict = None

    def _select_top2(self):
        top2_dict = {}
        single_dict = {}
        for gene in self._sumList_perApA:
            npary_apa_sum = self._sumList_perApA[gene]  # it is a np.array
            ls_apa = self._cached_dict[gene]            # it is a list
            if len(ls_apa) < 2:                         # pick out the ApAs that are single within the gene
                single_dict[gene] = ls_apa
                continue
            elif len(ls_apa) != npary_apa_sum.shape[0]:
                raise Exception('During _select_top2, the sum_per_apa list and apa list has different lengths.')
            else:
                pass
            # what we have here is:
            # 1. convert np.ary to list, sort, and take the top 2 elements
            # 2. When top1 == top2, list.index() would give same the index for these two elements
            #     So, after indexing top1, we covered it with -1, and then index top2.
            #     This way, even we have sorted_top2 = [0, 0], ls_apa[idx1] and ls_apa[idx2] won't be the same apa.
            # 3. Take two different apas.

            ls_apa_sum = npary_apa_sum.tolist()
            sorted_top2 = sorted(ls_apa_sum)[-2:]
            top2, top1 = sorted_top2[:]

            idx1 = ls_apa_sum.index(top1)
            ls_apa_sum[idx1] = -1
            idx2 = ls_apa_sum.index(top2)

            new_list_apa = [ls_apa[idx1], ls_apa[idx2]]
            top2_dict[gene] = new_list_apa
        return top2_dict, single_dict

    def run(self):
        self._lgr.info("Start Selecting Top2 ApAs...")
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

        myresult = myio.Top2selectionOutput(self._rcTable.name)
        myresult.add2content(self._header)
        self._lgr.info("the output header set to be: %s", self._header)

        self._lgr.info("writing output...")
        for gene in self._top2_dict:
            ls_apa = self._top2_dict[gene]
            for apa in ls_apa:
                apaline = apa.build_line()
                myresult.add2content(apaline)
        myresult.open2write(myresult.content)
        self._lgr.info("Done!")


def run_top2selecter(args):
    '''
    this RUN func should create a APP object and deliver all the input and params from args
    :param args:
    :return:
    '''
    time_start = time.time()

    app = Top2SelectionApp(args.PARC_Table,
                            args.hasHeader,
                            args.hasPCT)
    app.run()

    time_end = time.time()
    print 'Run time: {0:.2f} seconds'.format(time_end - time_start)
    # sys.stdout.flush()

