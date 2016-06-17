import myclasses
import time
import logger
import myio
import const
import utils


class MergeAPALApp:
    def __init__(self, ls_inputs, undefined, outputname):
        '''
        This app class takes a list of single-sample RC files and merge them into one.
        Default settings merge the list of files derived from defined masterlist.
        If you want to merge "Unidentified" UID_files, specify --undefined upon calling.
        '''
        self._ls_inputs = utils.list_inputs_from_file(ls_inputs)
        # self._ls_inputs contains a list of Path() objects that are through existing check
        self._undefined = undefined
        # default is True, when undefined is given in command line, it turns False
        self._outputname = outputname
        self._lgr = logger.set_logger('Merged_' + self._outputname[:-4] + '.log')
    '''
    The following two funcs are in one set, serve for one purpose
    '''

    def take_reference_from_inputs(self):
        '''
        If the inputs all came in as the standard 6-col single-sample ApA readcounts table,
            they should have the same first 5 columns, which goes down to the ref_header.
            Here I am referencing the first 5 columns from the Nth input.
        '''
        fp = self._ls_inputs[const.Nth_INPUT_AS_REFERENCE].__str__()
        self._lgr.info('Reference info taken from #%s file of input', str(const.Nth_INPUT_AS_REFERENCE+1))
        opened_fp = open(fp, 'r')
        header = opened_fp.readline().rstrip().split('\t')
        header = header[:const.DEFINED_REFERENCE_COLN]
        # Take the first 5 columns as reference header
        # It is a configurable const.
        ref_header = ''
        if len(header) == const.DEFINED_REFERENCE_COLN and \
            header[const.GENES_coli] == 'GENES' and\
            header[const.IDS_coli] == 'IDs' and\
            header[const.TRANSCRIPTS_coli] == 'TRANSCRIPTS' and \
            header[const.TYPES_coli] =='TYPES' and \
            header[const.COORDS_coli] == 'COORDS':
            # make new header
            for i in range(const.DEFINED_REFERENCE_COLN):
                ref_header += header[i] + '\t'
                # leave '\t' at the end of ref_header because more sample names will be added during running.
            ls_genes = []
            ls_geneids = []
            ls_transcripts = []
            ls_types = []
            ls_coords = []
            for line in opened_fp:
                ls_line = line.rstrip().split('\t')
                ls_genes.append(ls_line[const.GENES_coli])
                ls_geneids.append(ls_line[const.IDS_coli])
                ls_transcripts.append(ls_line[const.TRANSCRIPTS_coli])
                ls_types.append(ls_line[const.TYPES_coli])
                ls_coords.append(ls_line[const.COORDS_coli])

            ref_dict = {'genes': ls_genes,
                        'geneids':ls_geneids,
                        'transcripts':ls_transcripts,
                        'types':ls_types,
                        'coords': ls_coords
                        }
            opened_fp.close()
        else:
            raise Exception('The input should be standard single-sample ApA readcounts table.\
                            You can specify column indexes in const.py if they are switched.')

        return ref_dict, ref_header

    def cache_rc(self):
        '''
        This func takes the RC info from the list of inputs at the specific column number
        '''
        ary_rc = []
        for filepath in self._ls_inputs:
            # use a func that can extract a column of info from the input file
            # the colnum is set default as const.RC_coli = 5
            # if header is specified as True, then the colnames in the header will be included.
            ls_rc = utils.list_1column_from(filepath.__str__(), colnum=const.RC_coli, header=True)
            ary_rc.append(ls_rc)
        return ary_rc

    '''
    The following three funcs are in one set, serve for one purpose
    '''
    def collect_uidapa(self):
        '''
        If the inputs all came in as the standard 2-col single-sample un-id ApA readcounts table,
            they should have different content in the first column.
            They would be cached in a sequence.
        '''
        ls_counters = []
        ls_headers = []
        ls_coords = []
        ls_rc = []
        accum = 0       # this is an accumulative counter
        ls_inputnames = []
        for filepath in self._ls_inputs:
            opened_afile = open(filepath.__str__(), 'r')
            ls_inputnames.append(filepath.name)

            headerline = opened_afile.readline()
            if len(headerline.rstrip().split('\t')) == const.UID_STD_COLN:
                # take the sample name
                header = headerline.rstrip().split('\t')[const.UID_RC_coli]
                for line in opened_afile:
                    ls_line = line.rstrip().split('\t')
                    ls_coords.append(ls_line[const.UID_COORDS_coli])
                    ls_rc.append(ls_line[const.UID_RC_coli])     # the rc data here is str
                    accum += 1                # memorize how many lines of coords have gone through so far
                ls_headers.append(header)
                ls_counters.append(accum)
                opened_afile.close()
            else:
                raise Exception('The input should be standard single-sample UID-ApA readcounts table.')
        ref_dict = {'coords': ls_coords,
                    'rc': ls_rc,
                    'counters': ls_counters,
                    'headers': ls_headers
                    }
        fnumber = len(ls_inputnames)
        snumber = len(ls_headers)
        if fnumber == snumber:
            pass
        else:
            raise Exception()

        for i in range(fnumber):
            if i == 0:
                fname, sname, counts = ls_inputnames[i], ls_headers[i], ls_counters[i]
            else:
                fname, sname, counts = ls_inputnames[i], ls_headers[i], ls_counters[i] - ls_counters[i-1]
            self._lgr.info("%s uid coords were collected from %s (sample:%s)", str(counts), fname, sname)
        self._lgr.info("%s uid coords in total were collected.", str(len(ref_dict['coords'])))

        return ref_dict

    def repack_uidapa_from(self, dict_uidapa):
        # load info from the previously generated dictionary
        ls_headers = dict_uidapa['headers']
        ls_counters = dict_uidapa['counters']
        ls_coords = dict_uidapa['coords']
        ls_rc = dict_uidapa['rc']

        self._lgr.info("check - the length of ls_rc is: %s", str(len(ls_rc)))
        self._lgr.info("check - 'if length of ls_rc equals to the length of ls_coords: %s", str(len(ls_rc) == len(ls_coords)))

        # re-pack the info into a new dictionary, samplewise.
        new_dict = {}
        left_idx = 0

        for header, accum in zip(ls_headers, ls_counters):
            new_list = []

            ls_subcoords = ls_coords[left_idx:accum]
            ls_subrc = ls_rc[left_idx:accum]

            for coords, rc in zip(ls_subcoords, ls_subrc):
                _ls_rc = [rc]  # in UIDAPA class, the _ls_rc attribute stores a list of rc

                # This is fxcking horrible, I named this variable ls_rc the first time i did it
                # and spent the whole afternoon chasing after the "bug"...
                # because I name the entire rc list after that name
                uidapa = myclasses.UIDAPA(coords, _ls_rc)
                new_list.append(uidapa)
            # add the element
            new_dict[header] = new_list
            left_idx = accum

        return new_dict

    def screen_uidapa_from(self, dict_uidapa):

        qualified_coords = []

        sample_list = dict_uidapa.keys()
        std_list = dict_uidapa[sample_list[const.Nth_UIDAPA_SCREEN_INPUT_AS_STANDARD]]
        self._lgr.info("#%s sample was used as standard coords provider.", str(const.Nth_UIDAPA_SCREEN_INPUT_AS_STANDARD + 1))

        del sample_list[const.Nth_UIDAPA_SCREEN_INPUT_AS_STANDARD]
        # the std got deleted from being examined through following:

        # set the minimum number of times same coords need to show up before it is established as a new ApA site.
        # len(dict_uidapa) == the number of samples included.
        amin = len(sample_list) - const.MISS_ALLOWED_FOR_NEWAPA_ESTABLISHMENT
        self._lgr.info("the screen minimum is: %s", amin)

        hits = 0
        for uidapa in std_list:
            counter = 0
            for samplename in sample_list:
                hasit = False
                for eachapa in dict_uidapa[samplename]:
                    if uidapa.coords() == eachapa.coords():
                        hasit = True
                        break
                if hasit:
                    counter += 1
            if counter >= amin:
                qualified_coords.append(uidapa.coords())
                hits += 1
                if hits % const.FIVE_HUNDRED_HITS == 0:
                    print "{} hits now".format(hits)

        return qualified_coords

    def run(self):
        self._lgr.info("self._undefined is: %s", self._undefined)
        if not self._undefined:
            self._lgr.info("taking reference info...")
            ref_dict, ref_header = self.take_reference_from_inputs()
            myresult = myio.MergeAPALOutput(self._outputname)
            self._lgr.info("taking read counts info...")
            ary_rc = self.cache_rc()

            myheader = ref_header
            sample_num = len(ary_rc)
            sample_size = len(ary_rc[0]) - 1   # ary_rc[i][0] is sample name
            self._lgr.info("sample number is %s with size of %s each", str(sample_num), str(sample_size))
            for i in range(sample_num):
                myheader += ary_rc[i][0] + '\t'
            self._lgr.info("output header built as: %s", myheader)
            myresult.add2content(myheader.rstrip() + '\n')

            self._lgr.info("merging read counts info ...")
            for i in range(sample_size):
                myline = ref_dict['genes'][i] + '\t ' + ref_dict['geneids'][i] + '\t' + ref_dict['transcripts'][i] + '\t'\
                        + ref_dict['types'][i] + '\t' + ref_dict['coords'][i] + '\t'
                for j in range(sample_num):

                    myline += ary_rc[j][i+1] + '\t'
                myresult.add2content(myline.rstrip() + '\n')

            myresult.open2write(myresult.content)
            self._lgr.info("Done!")
        else:
            collected_uidapa = self.collect_uidapa()
            self._lgr.info("uidapa collected.")
            repacked_uidapa = self.repack_uidapa_from(collected_uidapa)
            self._lgr.info("uidapa repacked.")
            qualified_coords = self.screen_uidapa_from(repacked_uidapa)
            self._lgr.info("uidapa screened.")

            myresult = myio.MergeAPAL_UID_Output(self._outputname)
            '''
            Here I am creating an intermediate result in non-standard format.
            '''
            self._lgr.info("storing %s qualified new coords.", str(len(qualified_coords)))
            for coords in qualified_coords:
                myline = coords + '\n'
                myresult.add2content(myline)
            myresult.open2write(myresult.content)
            self._lgr.info("Done!")


def run_merger(args):
    '''
    this RUN func should create a APP object and deliver all the input and params from args
    :param args:
    :return:
    '''
    time_start = time.time()

    app = MergeAPALApp(args.PARC_list,
                        args.undefined,
                        args.outputname
                        )
    app.run()

    time_end = time.time()
    print 'Run time: {0:.2f} seconds'.format(time_end - time_start)
    # sys.stdout.flush()

