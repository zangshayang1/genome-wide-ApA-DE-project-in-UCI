import utils
import pathlib
import const


class LineUp:
    def __init__(self, line):
        self._line = line
        self._identity = self._identify()
    # Use internal method to assign an identity to itself upon init
    # Enough room left to have the class expanded to nearly every type of input format

    def _identify(self):
        if self._line.startswith('@'):
            identity = 'SAMheader'
        elif self._line.startswith('chr'):
            identity = 'BEDread'
        else:
            identity = 'SAMread'
        return identity

    def parse_line(self):
        if self._identity == 'SAMread':
            ls_line = self._line.rstrip().split('\t')
            if len(ls_line) < const.SAM_ESSENTIAL_LENGTH:
                raise IndexError('Input SAM file is incomplete.')
            else:
                # The way every attribute got assigned is in line with the std format of SAM file
                # You need to make changes here if some SAM file comes in another format
                qname, flag, chroms, leftmost, mqual, cigar, rnext, pnext, tlen, seq, qual = ls_line[:const.SAM_ESSENTIAL_LENGTH]
                other = ls_line[const.SAM_ESSENTIAL_LENGTH:]
                samread = SamRead(qname, flag, chroms, leftmost, mqual, cigar, rnext, pnext, tlen, seq, qual, other)
                return samread
        elif self._identity == 'SAMheader':
            return self._line
        elif self._identity == 'BEDread':
            ls_line = self._line.rstrip().split('\t')
            if len(ls_line) < const.BED_ESSENTIAL_LENGTH:
                raise IndexError('Input BED file is incomplete.')
            else:
                # The way every attribute got assigned is in line with the std format of BED file
                # You need to make changes here if some BED file comes in another format
                chrom, leftmost, rightmost, name, mqual, strand = ls_line[:const.BED_ESSENTIAL_LENGTH]
                otherlist = ls_line[const.BED_ESSENTIAL_LENGTH:]
                bedread = BedRead(chrom, leftmost, rightmost, name, mqual, strand, otherlist)
            return bedread
        else:
            return None

    def get_id(self):
        return self._identity


class BedRead:
    def __init__(self, chrom, leftmost, rightmost, name, mqual, strand, otherlist):
        self._chrom = chrom
        self._leftmost = leftmost
        self._rightmost = rightmost
        self._name = name
        self._mqual = mqual
        self._strand = strand
        self._otherlist = otherlist

    def form_apa_coords(self):
        if self._strand == '+':
            return self._chrom + self._strand + self._rightmost    # the PAS on forward strand such as "chr+rightend"
        elif self._strand == '-':
            return self._chrom + self._strand + self._leftmost
        else:
            raise Exception('This BED line cannot be defined on strand')

    def build_line(self):
        bedline = self._chrom + '\t' + self._leftmost + '\t' + self._rightmost + '\t' + self._name + '\t' +\
            self._mqual + '\t' + self._strand
        for item in self._otherlist:
            bedline += '\t' + item
        bedline += '\n'
        return bedline


class SamRead:
    def __init__(self, qname, flag, chroms, leftmost, mqual, cigar, rnext, pnext, tlen, seq, qual, other):
        '''
        :param qname: query name
        :param flag: FLAG info
        :param chroms: chromosome name
        :param leftmost: 1-based leftmost mapping position
        :param mqual: mapping quality
        :param cigar: CIGAR info
        :param rnext: see samfile specifications
        :param pnext: see samfile specifications
        :param tlen: see samfile specifications
        :param seq: read sequence
        :param qual: read sequence quality
        '''
        self._qname = qname
        self._flag = int(flag)
        self._chroms = chroms
        self._leftmost = int(leftmost)
        self._mqual = mqual
        self._cigar = cigar
        self._rnext =rnext
        self._pnext = pnext
        self._tlen = tlen
        self._seq = seq
        self._qual = qual
        self._strand = self._get_strand() # internal method applied upon init
        self._other = other

    def _get_strand(self):
        if self._flag == const.SAMFLAG_FORWARD or self._flag == const.SAMFLAG_2ndALIGNMENT_FORWARD:
            strand = '+'
        elif self._flag == const.SAMFLAG_REVERSE or self._flag == const.SAMFLAG_2ndALIGNMENT_REVERSE:
            strand = '-'
        else:
            raise Exception("Failed to parse FLAG info")
        return strand

    def pat_locator(self):
        if self._strand == '+':
            pat_start = self._leftmost + utils.cigar_dist(self._cigar) - const.SAM_BASE  # SAM is 1-based
        elif self._strand == '-':
            pat_start = self._leftmost - const.SAM_BASE                           # but reference genome is 0-based
        else:
            raise Exception("Unrecognized strand info")
        return pat_start

    def build_line(self):
        samline = self._qname + '\t' + str(self._flag) + '\t' + self._chroms + '\t' + str(self._leftmost) + '\t' \
                  + self._mqual + '\t' + self._cigar + '\t' + self._rnext + '\t' + self._pnext + '\t' \
                  + self._tlen + '\t' + self._seq + '\t' + self._qual
        for item in self._other:
            samline = samline + '\t' + item
        samline += '\n'
        return samline


class APA:
    # here I leave *args for pct info if needed in the future
    def __init__(self, gene, id, transcript, type, coords, ls_rc, *args):
        self._gene = gene
        self._id = id
        self._transcript = transcript
        self._type = type
        self._coords = coords
        self._ls_rc = ls_rc
        self._ls_pct = None
        if not len(args) == 0:
            self._ls_pct = args[0]

    def gene(self):
        return self._gene

    def id(self):
        return self._id

    def transcript(self):
        return self._transcript

    def type(self):
        return self._type

    def coords(self):
        return self._coords

    def readcounts(self):
        numeric_rc = []
        for rc in self._ls_rc:
            numeric_rc.append(int(rc))
        return numeric_rc

    def percentages(self):
        numeric_pct = []
        for pct in self._ls_pct:
            numeric_pct.append(float(pct))
        return numeric_pct

    def position(self):
        tup = utils.parse_apa_coords(self._coords)      # this func returns a tup (chr, strand, postion)
        return tup[-1]

    def chromsome(self):
        tup = utils.parse_apa_coords(self._coords)      # this func returns a tup (chr, strand, postion)
        return tup[0]

    def strand(self):
        tup = utils.parse_apa_coords(self._coords)      # this func returns a tup (chr, strand, postion)
        return tup[1]

    def sum_of_rc(self):
        the_sum = 0
        try:
            the_sum = sum(self._ls_rc)
        except TypeError:
            for number in self._ls_rc:
                the_sum += int(number)
        return the_sum

    def build_line(self):
        apaline = self._gene + '\t' + self._id + '\t' + self._transcript + '\t' + self._type + '\t' + self._coords
        if self._ls_rc is None:
            apaline += '\n'
            return apaline
        else:
            for rc in self._ls_rc:
                apaline += '\t' + rc  # rc::str
            if self._ls_pct is None:
                pass
            else:
                for pct in self._ls_pct:
                    apaline += '\t' + pct
            apaline += '\n'
            return apaline


class UIDAPA:
    def __init__(self, coords, ls_rc):
        self._coords = coords
        self._ls_rc = ls_rc

    def position(self):
        tup = utils.parse_apa_coords(self._coords)   # this func returns a tup (chr, strand, postion)
        return tup[-1]

    def sum_of_rc(self):
        the_sum = 0
        try:
            the_sum = sum(self._ls_rc)
        except TypeError:
            for number in self._ls_rc:
                the_sum += int(number)
        return the_sum

    def coords(self):
        return self._coords

    def readcounts(self):
        return self._ls_rc


class RCTableCacher:
    def __init__(self, rcTable, hasHeader):
        self._rcTable = pathlib.Path(rcTable)
        if not self._rcTable.exists():
            raise Exception('Cannot find the input read counts table.')
        self._hasHeader = hasHeader

        self._header = None
        self._cached_dict = None

    def _cache_into_dict(self, linestart=1, hasPCT=False):
        infile = open(self._rcTable.__str__(), 'r')
        ls_content = infile.readlines()
        infile.close()

        if linestart == 1:
            header = ls_content[0]
        else:
            header = None

        cache_dict = {}
        for i in range(linestart, len(ls_content)):
            line = ls_content[i]
            ls_line = line.rstrip().split('\t')

            gene = ls_line[const.GENES_coli]
            id = ls_line[const.IDS_coli]
            transcript = ls_line[const.TRANSCRIPTS_coli]
            type = ls_line[const.TYPES_coli]
            coords = ls_line[const.COORDS_coli]

            if not hasPCT:
                ls_rc = ls_line[const.DEFINED_REFERENCE_COLN:]
                apa = APA(gene, id, transcript, type, coords, ls_rc)
            else:
                sep_idx = len(ls_line[const.DEFINED_REFERENCE_COLN:])/2 + const.DEFINED_REFERENCE_COLN
                ls_rc = ls_line[const.DEFINED_REFERENCE_COLN:sep_idx]
                ls_pct = ls_line[sep_idx:]
                apa = APA(gene, id, transcript, type, coords, ls_rc, ls_pct)
            try:
                cache_dict[gene].append(apa)
            except KeyError:
                cache_dict[gene] = []
                cache_dict[gene].append(apa)
        return cache_dict, header




class ApALFC:
    def __init__(self, gene, id, transcript, types, coords, relative_pos, *args):
        '''
        built this class to store ApA logFC data info
        *args were left for future development
        '''
        self._gene = gene
        self._id = id
        self._transcript = transcript
        self._types = types
        self._coords = coords
        self._relative_pos = relative_pos

    def relative_pos(self):
        return self._relative_pos

    def coords(self):
        return self._coords
