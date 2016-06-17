#import pysam as ps


class Reference:
    def __init__(self, fpath):
        self._fpath = fpath

    def show_fasta(self):
        self._fasta = ps.FastaFile(self._fpath)
        return self._fasta

    def show_masterlist(self):
        self._mlist = open(self._fpath, 'r')
        return self._mlist

    def close_masterlist(self):
        try:
            self._mlist.close()
        except NameError:
            print 'There is no masterlist opened.'
        return None