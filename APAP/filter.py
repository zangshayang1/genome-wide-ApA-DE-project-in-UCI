import myclasses
import time
import pathlib
import logger
import reference
import myio


class IMPriming:
    '''
    need to build this class that takes a Reference object,
                                        a MyInput object that runs the following out:
                                            a strand from SamRead.strand
                                            a return from SamRead.pat_locator(),
                                            a chroms from SamRead.chroms
                                        a window_size: int,
                                        a deny_number: int,
                              and gives a logic return
    '''
    def __init__(self, refasta, strand, chroms, start, window_size, deny_number):
        self._refasta = refasta
        self._strand = strand
        self._chroms = chroms
        self._start = start
        self._window_size = window_size
        self._dn = deny_number
        self._imp = self._imp_check()    # internal method applied upon init

    def _imp_check(self):
        '''
        6 consecutive "A" justify an internal priming event.
        OR
        7/10 (deny_number / window_size) "A" justify one.

        The standard is not configurable outside this class.
        '''
        counter = 0

        if self._strand == '+':
            seq = self._refasta.fetch(self._chroms, self._start, self._start + self._window_size)
            flag_imp = True if seq[:self._dn] == 'AAAAAA' or seq[:self._dn] == 'aaaaaa' else False
                                                            # genome.fa may contain lowercase letters
                                                            # depending on which one you use
            for sgnd in seq:        # sgnd, short for single-nucleotide
                if sgnd == 'A' or sgnd == 'a':
                    counter += 1
        elif self._strand == '-':
            seq = self._refasta.fetch(self._chroms, self._start - self._window_size, self._start)
            flag_imp = True if seq[0-self._dn:] == 'TTTTTT' or seq[0-self._dn:] == 'tttttt' else False
            for sgnd in seq:
                if sgnd == 'T' or sgnd == 't':
                    counter += 1
        else:
            raise Exception("Other strand detected during 'imp_check'")

        if flag_imp:
            pass
        elif counter > self._dn:
            flag_imp = True
        else:
            pass

        return flag_imp


class IMFilterApp:
    '''
    This App class should take a referencefile, inputfile, window_size, deny_number params
    '''
    def __init__(self, reference, infilename, window_size, deny_number):
        # infilename could be a path/file, pack them in a Path() object will be better
        self.infilename = pathlib.Path(infilename)
        if not self.infilename.exists():
            raise Exception('The input SAM does not exist!')
        self.reference = pathlib.Path(reference)
        if not self.reference.exists():
            raise Exception('The input reference does not exist!')
        self.window_size = window_size
        self.deny_number = deny_number
        # the name of the infile should be XXXXXX.sam
        self._lgr = logger.set_logger(self.infilename.name[:-4] + '.log')

    def run(self):

        ref_path = self.reference.__str__()
        refasta = reference.Reference(ref_path).show_fasta()
        # Having reference an independent class leaves enough room for accommodating future expansion

        myresult = myio.IMFilterOutput(self.infilename.name)
        # another problem tho is args.infile_sam came along way with multiple "/"
        # so here you need to extract the basename

        opened_input = open(self.infilename.__str__(), 'r')

        self._lgr.info('referencing the genome fasta provided...')

        counter_all, counter_imp = 0, 0
        for line in opened_input:
            samline = myclasses.LineUp(line)
            if samline._identity == 'SAMheader':
                continue
                # The header is excluded from the output SAM
            elif samline._identity == 'SAMread':
                counter_all += 1
                samread = samline.parse_line()
                start = samread.pat_locator()
                p = IMPriming(refasta, samread._strand, samread._chroms, start, self.window_size, self.deny_number)
                if p._imp:
                    counter_imp += 1
                    continue
                else:
                    newline = samread.build_line()
                    myresult.add2content(newline)
            else:
                raise Exception('Unexpected format of lines appear in SAM input')
        self._lgr.info("creating filtered SAM file excluding headers ...")
        self._lgr.info("through {0} reads".format(counter_all))
        self._lgr.info("\t\t{0} ({1}) were removed due to internal priming.".format(counter_imp, float(counter_imp)/float(counter_all)))
        myresult.open2write(myresult.content)
        self._lgr.info("Done!")


def run_filter(args):
    '''
    this RUN func should create a APP object and deliver all the input and params from args
    :param args:
    :return:
    '''
    time_start = time.time()

    app = IMFilterApp(args.reference_genome,
                      args.infile_sam,
                      args.window_size,
                      args.deny_number
                      )
    app.run()

    time_end = time.time()
    print 'Run time: {0:.2f} seconds'.format(time_end - time_start)
    # sys.stdout.flush()

