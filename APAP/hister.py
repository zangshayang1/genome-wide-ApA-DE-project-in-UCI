import pathlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import time
import logger
import const
import utils
import reference
import myclasses


class ApALFC_HISTApp:
    def __init__(self, infile, refgenome, pattern, output_name_pattern, window_size=200):
        self._infile = pathlib.Path(infile)
        if not self._infile.exists():
            raise Exception('Can not find input file.')

        self._refgenome = pathlib.Path(refgenome)
        if not self._refgenome.exists():
            raise Exception('Can not find refgenome file')

        self._pattern = pattern
        if '/' in self._pattern:
            self._length_of_pattern = len(self._pattern) - 2
        else:
            self._length_of_pattern = len(self._pattern)

        self._window_size = window_size
        if not self._window_size % 2 == 0:
            raise Exception('Please specify an even number for window_size')

        self._output_name_pattern = output_name_pattern

        self._dict_apalfc = None
        self._opened_ref = None

        self._lgr = logger.set_logger(self._infile.name[:-3]+'log')

    def _open_reference(self):
        ref = reference.Reference(self._refgenome.__str__())
        opened_ref = ref.show_fasta()
        self._lgr.info("reference genome opened successfully.")
        # return an object that you can call .fetch() on
        return opened_ref

    def _cached_apalfclist(self):
        opened_infile = open(self._infile.__str__(), 'r')
        content = opened_infile.readlines()
        opened_infile.close()

        # skip the header
        ls_distal = []
        ls_approx = []
        for i in range(1, len(content)):
            ls_line = content[i].rstrip().split('\t')
            gene = ls_line[const.GENES_coli]
            id = ls_line[const.IDS_coli]
            transcript = ls_line[const.TRANSCRIPTS_coli]
            type = ls_line[const.TYPES_coli]
            coords = ls_line[const.COORDS_coli]
            relative_pos = ls_line[const.RELATIVE_POS_coli]

            apalfc = myclasses.ApALFC(gene, id, transcript, type, coords, relative_pos)
            if relative_pos == 'distal':
                ls_distal.append(apalfc)
            elif relative_pos == 'approximal':
                ls_approx.append(apalfc)

        dict_apalfc = {'distal': ls_distal,
                       'approx': ls_approx
                       }
        return dict_apalfc

    def get_data_on(self, which):
        if self._dict_apalfc is None:
            raise Exception('Input ApALFCList has not been cached ready.')
        if self._opened_ref is None:
            raise Exception('Reference has not been opened for reading.')
        ls = self._dict_apalfc[which]
        ls_fasta_info = []
        y = []
        for apalfc in ls:
            coords = apalfc.coords()
            chroms, strand, pos = utils.parse_apa_coords(coords)[:]
            half_size = self._window_size/2
            seq = self._opened_ref.fetch(chroms, pos - half_size, pos + half_size - 1 + self._length_of_pattern)
            # so the length of the seq input should be window_size + leng_of_pattern
            seq = seq.upper()
            # capitalize it before going to pattern matching
            if strand == '+':
                xlist = utils.count_pattern_on_xlist(seq, self._pattern)
                fasta_tup = (coords, seq)
            elif strand == '-':
                rseq = utils.reverse_compliment(seq)
                xlist = utils.count_pattern_on_xlist(rseq, self._pattern)
                fasta_tup = (coords, rseq)
            else:
                raise Exception("there are other type of strand info")
            y.append(xlist)
            ls_fasta_info.append(fasta_tup)
        return np.array(y).sum(0), ls_fasta_info

    def plot(self, x, y, distal_or_approx):
        fig = plt.figure()
        plt.plot(x, y)
        fig.savefig(self._output_name_pattern + '.' + distal_or_approx + '.pdf', dpi=100)

    def export_fasta(self, ls_fasta, distal_or_approx, lowercase=True):
        o = open(self._output_name_pattern + '.' + distal_or_approx + '.fasta', 'w')
        for tup in ls_fasta:
            o.write(">COORDS:{}\n".format(tup[0]))
            if lowercase:
                o.write("{}\n".format(tup[1]).lower())
            else:
                o.write("{}\n".format(tup[1]))
        o.close()
        return None

    def export_counting_result(self, x_idx, y_distal, y_approx):
        if not isinstance(x_idx, list):
            x_idx = x_idx.tolist()
        if not isinstance(y_distal, list):
            y_distal = y_distal.tolist()
        if not isinstance(y_approx, list):
            y_approx = y_approx.tolist()
        # write out data for plotting outside of this program
        o = open(self._output_name_pattern + '.hist_counts', 'w')
        header = 'idx\tdistal\tapprox\n'
        o.write(header)
        for i, c1, c2 in zip(x_idx, y_distal, y_approx):
            o.write(str(i) + '\t' + str(c1) + '\t' + str(c2) + '\n')
        o.close()

    def run(self):
        self._lgr.info("Start running...")

        self._dict_apalfc = self._cached_apalfclist()

        self._lgr.info("Input list info cached.")

        self._opened_ref = self._open_reference()

        self._lgr.info("start getting data ...")

        y_distal, fasta_distal = self.get_data_on('distal')

        self._lgr.info("the length of distal data: %s", len(y_distal))

        y_approx, fasta_approx = self.get_data_on('approx')

        self._lgr.info("the length of approx data: %s", len(y_approx))
        # both lengths should be 201

        # generate x-axis coordinates
        half_size = self._window_size/2
        x_idx = np.arange(-half_size, half_size)

        self._lgr.info("start exporting fasta...")

        self.export_fasta(fasta_distal, 'distal')

        self.export_fasta(fasta_approx, 'approx')

        self._lgr.info("Successful!")

        self._lgr.info("start exporting counting results...")

        self.export_counting_result(x_idx, y_distal, y_approx)

        self._lgr.info("Successful!")

        self._lgr.info("start plotting...")

        self.plot(x_idx, y_distal, 'distal')
        self.plot(x_idx, y_approx, 'approx')
        self._lgr.info("Successful!")


def run_hister(args):
    '''
    this RUN func should create a APP object and deliver all the input and params from args
    :param args:
    :return:
    '''
    time_start = time.time()

    app = ApALFC_HISTApp(args.in_list,
                         args.reference_genome,
                         args.pattern,
                         args.output_name,
                         args.window_size
                         )
    app.run()

    time_end = time.time()
    print 'Run time: {0:.2f} seconds'.format(time_end - time_start)
    # sys.stdout.flush()