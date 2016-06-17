import myclasses
import time
import logger
import myio
import const
import utils


class MLCleanerApp(myclasses.RCTableCacher):
    def __init__(self, rcTable, window_size, hasHeader, delrc):
        myclasses.RCTableCacher.__init__(self, rcTable, hasHeader)
        self._window_size = window_size
        self._sorted_dict = None
        self._nodup_dict = None
        self._merged_dict = None
        self._delrc = delrc
        self._lgr = logger.set_logger(self._rcTable.name[:-4] + '.log')

    def _cache_into_dict(self, linestart=1, hasPCT=False):
        # hasPCT still there only because it fully covers the "_cache_into_dict()" method in RCTableCacher class.

        infile = open(self._rcTable.__str__(), 'r')
        ls_content = infile.readlines()
        infile.close()

        if linestart == 1:
            if not self._delrc:
                header = ls_content[0]
            else:
                ls_header = ls_content[0].rstrip().split('\t')
                header = ''
                for i in range(const.DEFINED_REFERENCE_COLN):
                    header += ls_header[i] + '\t'
                header = header.rstrip() + '\n'
        else:
            header = None

        cache_dict = {}
        for i in range(linestart, len(ls_content)):
            line = ls_content[i]
            ls_line = line.rstrip().split('\t')

            # remove the lines without gene names or gene ids or its type is upstreamAntisense
            if '' in ls_line or ls_line[const.NA_COLN_upstreamAntisense] == 'upstreamAntisense':
                continue
            else:
                gene = ls_line[const.GENES_coli]
                id = ls_line[const.IDS_coli]
                transcript = ls_line[const.TRANSCRIPTS_coli]
                type = ls_line[const.TYPES_coli]
                coords = ls_line[const.COORDS_coli]
                if not self._delrc:
                    ls_rc = ls_line[const.DEFINED_REFERENCE_COLN:]
                else:
                    ls_rc = None
                apa = myclasses.APA(gene, id, transcript, type, coords, ls_rc)
                try:
                    cache_dict[gene].append(apa)
                except KeyError:
                    cache_dict[gene] = []
                    cache_dict[gene].append(apa)

        return cache_dict, header

    def _remove_APAdups(self):
        nodup_dict = {}
        # This func applies on sorted_cached_dict
        for gene in self._sorted_dict:
            ls_apa = self._sorted_dict[gene]
            new_coords = []
            newls_apa = []
            for apa in ls_apa:
                if apa.coords() in new_coords:
                    continue
                else:
                    new_coords.append(apa.coords())
                    newls_apa.append(apa)
            nodup_dict[gene] = newls_apa
        return nodup_dict

    def _sort_APAs_by_coords(self):
        sorted_dict = {}
        for gene in self._cached_dict:
            ls_coords = []
            newls_apa = []
            ls_apa = self._cached_dict[gene]
            for apa in ls_apa:
                ls_coords.append(apa.coords())
            sorted_coords = sorted(ls_coords)
            for i in range(len(sorted_coords)):
                for j in range(len(ls_apa)):
                    if ls_apa[j].coords() == sorted_coords[i]:
                        newls_apa.append(ls_apa[j])
            sorted_dict[gene] = newls_apa
        return sorted_dict

    def _merge_ApAs(self):
        '''
        This step removes the single ApAs
        '''
        merged_dict = {}
        for gene in self._nodup_dict:
            ls_apa = self._nodup_dict[gene]
            if len(ls_apa) >= 2:
                newls_apa = utils.close_merge(ls_apa, self._window_size)
                merged_dict[gene] = newls_apa
            else:
                continue
        return merged_dict

    def run(self):
        self._lgr.info("hasHeader? %s", self._hasHeader)

        if not self._hasHeader:
            self._cached_dict, self._header = self._cache_into_dict(linestart=0)
        else:
            self._cached_dict, self._header = self._cache_into_dict()  # the default linestart = 1

        self._lgr.info("resorting...")
        self._sorted_dict = self._sort_APAs_by_coords()
        self._lgr.info("removing dups...")
        self._nodup_dict = self._remove_APAdups()

        if not self._delrc:
            self._lgr.info("merging close ApAs...")
            self._merged_dict = self._merge_ApAs()
        else:
            self._merged_dict = None
            self._lgr.info("skip merging close ApAs \
            because --delrc is specified and no read counts can be taken into consideration.")

        myresult = myio.MLCleanerOutput(self._rcTable.name)
        myresult.add2content(self._header)
        self._lgr.info("output header is: %s", self._header)

        self._lgr.info("writing clean table...")
        if not self._delrc:
            for gene in self._merged_dict:
                ls_apa = self._merged_dict[gene]
                for apa in ls_apa:
                    apaline = apa.build_line()
                    myresult.add2content(apaline)
        else:
            for gene in self._nodup_dict:
                ls_apa = self._nodup_dict[gene]
                for apa in ls_apa:
                    apaline = apa.build_line()
                    myresult.add2content(apaline)

        myresult.open2write(myresult.content)
        self._lgr.info("Done!")


def run_cleaner(args):
    '''
    this RUN func should create a APP object and deliver all the input and params from args
    :param args:
    :return:
    '''
    time_start = time.time()

    app = MLCleanerApp(args.masterlist,
                        args.window_size,
                        args.hasHeader,
                        args.delrc)
    app.run()

    time_end = time.time()
    print 'Run time: {0:.2f} seconds'.format(time_end - time_start)
    # sys.stdout.flush()

