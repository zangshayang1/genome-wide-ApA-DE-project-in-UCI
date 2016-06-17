

class MyOutput:
    def __init__(self, name):
        self._name = name
        self.content = []

    def add2content(self, line):
        self.content.append(line)
        return None

    def open2write(self, content):
        o = open(self._name, 'w')
        for line in content:
            o.write(line)
        o.close()


class IMFilterOutput(MyOutput):
    def __init__(self, inputname):
        MyOutput.__init__(self, 'IMFiltered_' + inputname)


class PARcounterOutput(MyOutput):
    def __init__(self, inputname):
        MyOutput.__init__(self, 'PARcounted_' + inputname[:-3] + 'txt')


class UnIdentifiedAPAsOutput(MyOutput):
    def __init__(self, inputname):
        MyOutput.__init__(self, 'UIDAPA_from_' + inputname[:-3] + 'txt')


class MergeAPALOutput(MyOutput):
    def __init__(self, inputname):
        MyOutput.__init__(self, 'Merged_' + inputname)


class MergeAPAL_UID_Output(MyOutput):
    def __init__(self, outputname):
        MyOutput.__init__(self, 'New_ApAs_' + outputname)


class MLCleanerOutput(MyOutput):
    def __init__(self, inputname):
        MyOutput.__init__(self, 'Cleaned_' + inputname)


class RCSummerOutput(MyOutput):
    def __init__(self, inputname):
        MyOutput.__init__(self, 'summed_DGE_from_' + inputname)


class RC2PCTerOutput(MyOutput):
    def __init__(self, inputname):
        MyOutput.__init__(self, 'PCT_followed_' + inputname)


class Top2selectionOutput(MyOutput):
    def __init__(self, inputname):
        MyOutput.__init__(self, 'Top2_' + inputname)


class REDtableOutput(MyOutput):
    def __init__(self, inputname):
        MyOutput.__init__(self, 'RED_' + inputname)
