import myclasses
import time
import pathlib
import logger
import reference
import myio
import const
import utils


class PARcounterApp:
    '''
    This App class should take a master list as reference,
                                a bed file as a provider of reads info
                                window_size param
    '''
    def __init__(self, reference, infilename, sample_name, window_size):
        self._reference = reference
        self._infilename = pathlib.Path(infilename)
        self._sample_name = sample_name
        self._window_size = window_size
        self._cached_dict = None
        # the name of the input is supposed to be XXXXXX.bed
        self._lgr = logger.set_logger(self._infilename.name[:-4] + '.log')

    def cache_BED_input(self):
        opened_input = open(self._infilename.__str__(), 'r')
        # self._infilename points to a Path() object
        self._lgr.info('BED file opened.')
        dict_bedreads = {}
        for line in opened_input:
            # BED file doesn't contain header
            bedread = myclasses.LineUp(line).parse_line()
            try:
                dict_bedreads[bedread.form_apa_coords()].append(bedread)
            except KeyError:
                dict_bedreads[bedread.form_apa_coords()] = []
                dict_bedreads[bedread.form_apa_coords()].append(bedread)
        opened_input.close()
        return dict_bedreads

    def run(self):
        '''
        Whatever is the header of the reference, I will take it down as the first line of the output
        '''

        ref = reference.Reference(self._reference)
        masterlist = ref.show_masterlist()

        header = masterlist.readline().rstrip() + '\t' + self._sample_name + '\n'
        # so the header of the reference must be clean
        self._lgr.info("header of the output: %s", header)

        myresult = myio.PARcounterOutput(self._infilename.name)
        # self._infilename points to a Path() object

        myresult.add2content(header)

        # 1. cache the info from input BED file
        self._cached_dict = self.cache_BED_input()
        self._lgr.info("input BED cached.")

        # 2. increment on the reference
        hits = 0
        self._lgr.info("start matching with the given reference list...")
        for line in masterlist:
            ls_line = line.rstrip().split('\t')
            coords = ls_line[const.COORDS_coli]
            incrementals, hitted, self._cached_dict = utils.increment_reads_at(coords, self._window_size, self._cached_dict)
            # what is worth-noting
            # self._cached_dict is mutated within utils.increment_reads_at()
            # for the saking of write out previously unidentified ApA coords
            if hitted:
                hits += 1
                if hits % const.FIVE_HUNDRED_HITS == 0:
                    print "{0} hits".format(hits)

            new_line = ls_line[const.GENES_coli] + '\t'\
                       + ls_line[const.IDS_coli] + '\t'\
                       + ls_line[const.TRANSCRIPTS_coli] + '\t'\
                       + ls_line[const.TYPES_coli] + '\t'\
                       + ls_line[const.COORDS_coli] + '\t'\
                       + str(incrementals) + '\n'
            myresult.add2content(new_line)
        self._lgr.info("%s hits on the reference list were found!", str(hits))

        # write output
        myresult.open2write(myresult.content)

        ref.close_masterlist()

        self._lgr.info("PARcounter table generated.")
        '''
        The following outputs unidentified ApA coords and their read counts.
        '''
        if len(self._cached_dict) > 0:
            sideresult = myio.UnIdentifiedAPAsOutput(self._infilename.name)
            # self._infilename points to a Path() object

            sideheader = const.UID_HEADER + self._sample_name + '\n'
            sideresult.add2content(sideheader)
            counter_uidapa = 0
            for coords in self._cached_dict:
                counter_uidapa += 1
                # Behind each "coords" key in self._cached_dict, it is a list of BedRead objects
                # So the len(ls_bedreads) is the number of hits on that coords.
                ls_bedreads = self._cached_dict[coords]
                uid_apa_line = coords + '\t' + str(len(ls_bedreads)) + '\n'
                sideresult.add2content(uid_apa_line)
            self._lgr.info("start to output %s un-identified potential ApAs from %s ...", str(counter_uidapa), self._infilename.name)
            sideresult.open2write(sideresult.content)
            self._lgr.info("Done!")
        else:
            self._lgr.info("There is no un-identified potential ApAs left.")
            pass


def run_counter(args):
    '''
    this RUN func should create a APP object and deliver all the input and params from args
    :param args:
    :return:
    '''
    time_start = time.time()

    app = PARcounterApp(args.masterlist,
                      args.infile_bed,
                      args.sample_name,
                      args.window_size,
                      )
    app.run()

    time_end = time.time()
    print 'Run time: {0:.2f} seconds'.format(time_end - time_start)
    # sys.stdout.flush()

