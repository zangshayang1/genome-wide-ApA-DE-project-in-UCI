import pathlib
import math


def increment_reads_at(coords, window_size, dict_bedreads):
    incrementals, hitted = 0, False
    tup_coords = parse_apa_coords(coords)
    chroms, strand, pos = tup_coords[:]      # this func in utilis returns chr, strand, pos(int) in a tuple

    lowerbound, higherbound = pos - window_size, pos + window_size
    for point in range(lowerbound, higherbound):
        point_coords = chroms + strand + str(point)
        if dict_bedreads.get(point_coords) is None:
            continue
        else:
            incrementals += len(dict_bedreads[point_coords])
            del dict_bedreads[point_coords]
            # dict_bedreads got mutated every time a hit found

    if incrementals > 0:
        hitted = True
    else:
        pass
    return incrementals, hitted, dict_bedreads


def parse_apa_coords(coords):
    neg_idx = coords.find('-')
    pos_idx = coords.find('+')
    strand_idx = max(neg_idx, pos_idx)
    strand = coords[strand_idx]
    chr = coords[:strand_idx]
    pos = int(coords[strand_idx+1:])
    tup = (chr, strand, pos)
    return tup


def cigar_dist(s):
    # This is brilliant
    left_idx = [0]
    right_idx = []
    for i in range(len(s)):
        # when i == "=" or capital letters
        if ord(s[i]) == 61 or 65 <= ord(s[i]) <= 90:
            right_idx.append(i)
            try:
                # i+1 indexes next numerics if i+1 is not out of bound
                left_idx.append(i+1)
            except IndexError:
                break
    ls_match = []
    for l, r in zip(left_idx, right_idx):
        ls_match.append(int(s[l:r]))
    return sum(ls_match)


class Check:
    def __init__(self, problem):
        self.problem = problem

    def pp(self):
        print "==================CHECK===================="
        print self.problem
        print "==================CHECK===================="


def list_inputs_from_file(ls_filename):
    f = open(ls_filename, 'r')
    ls_out = []
    for line in f:
        p = pathlib.Path(line.rstrip('\n'))
        if not p.exists():
            raise Exception('One or more input does not exist!')
        ls_out.append(p)
    f.close()
    return ls_out


def list_1column_from(filename, colnum, **kwargs):
    fp = pathlib.Path(filename)
    if fp.exists():
        f = open(filename, 'r')
        collist = []
        if 'header' in kwargs and kwargs['header'] == True:
            for line in f:
                ls_line = line.rstrip().split('\t')
                collist.append(ls_line[colnum])
        else:
            f.readline()
            for line in f:
                ls_line = line.rstrip().split('\t')
                collist.append(ls_line[colnum])
    else:
        raise Exception('Failed to find the input file')
    return collist


def close_merge(ls, r):
    '''
    This func deploys several APA objects' functions.

    This func is specifically designed for MLCleanerApp because it functions at one level higher than APA object
    but one level lower than MLCleaner object. I should've created a intermediate class to deal with ls_apa objects.

    :param ls: a list of sorted objects
    :param r: searching radius
    :return: a new list of remaining objects
    '''
    new_list = []

    i = 0
    while i < len(ls):

        need_filter = False
        start = ls[i].position()

        j = 1
        while i+j < len(ls):
            if start + r > ls[i+j].position():
                need_filter = True
                j += 1
            else:
                break

        if need_filter:
            big, idx = 0, 0
            for item in ls[i: i+j]:
                if item.sum_of_rc() >= big:
                    big = item.sum_of_rc()
                    idx = ls[i: i+j].index(item) + i
            new_list.append(ls[idx])
        else:
            new_list.append(ls[i])

        i = i + j

    return new_list


def calculate_red_score(approx, distal):
    if not isinstance(approx, int) or not isinstance(distal, int):
        raise Exception('The input pA1 or pA2 is not int.')
    if approx < 0 or distal < 0:
        raise Exception('The input pA1 or pA2 is negative.')
    approx = float(approx)     # so when making divisions, there won't be math.log(0,2) rendering math domain error.
    distal = float(distal)
    red = math.log((distal+1)/(approx+1), 2)
    return red


def first_is_approximal(coords1, coords2):
    tup1, tup2 = parse_apa_coords(coords1), parse_apa_coords(coords2)
    strand1, pos1 = tup1[1:]
    strand2, pos2 = tup2[1:]
    if (strand1 == '+' and strand2 == '+' and pos1 < pos2) or (strand1 == '-' and strand2 == '-' and pos1 > pos2):
        return True
    elif (strand1 == '+' and strand2 == '+' and pos1 > pos2) or (strand1 == '-' and strand2 == '-' and pos1 < pos2):
        return False
    elif (strand1 == '+' and strand2 == '-') or (strand1 == '-' and strand2 == '+'):
        print "The input two coords info are not in the same strand."
        print tup1, tup2
        print "====================================================="
        return None
    else:
        print "===================== ERROR ========================="
        print tup1, tup2
        print "===================== ERROR ========================="
        return None




'''
The following functions are built for hister specifically.
'''

def is_single_pattern(p):
    '''
    input a pattern
    return if it includes "/" which means an alternative sgnd.
    '''
    or_idx = p.find('/')
    if or_idx == -1:
        return True
    elif or_idx == 0:
        raise Exception('Invalid Pattern Found. / can not be the first position of the input pattern.')
    else:
        return False


def derive_two_patterns(p):
    '''
    input a pattern that includes alternatives
    :returns two patterns
    '''
    or_idx = p.find('/')
    if or_idx == -1:
        raise Exception('The pattern does not contain /.')
    interchangeable1 = p[or_idx-1]
    interchangeable2 = p[or_idx+1]
    p1 = p[:or_idx-1] + interchangeable1 + p[or_idx+2:]
    p2 = p[:or_idx-1] + interchangeable2 + p[or_idx+2:]
    return p1, p2


def count_pattern_on_xlist(seq, p):
    '''
    make a histogram-type incrementation once the given pattern appears
    utilize the "is_single_pattern func" at the root of the decision tree
    So if the input pattern has alternatives, both of them must be considered during incrementation.
    '''
    # the seq will be capitalized before becoming the input of this func.
    if not is_single_pattern(p):
        p1, p2 = derive_two_patterns(p)
        if not len(p1) == len(p2):
            raise Exception("the lengths of the two derived patterns are not the same.")
        f = len(seq)
        l = len(p1)
        xlist = [0] * (f - l + 1)
        for i in range(len(xlist)):
            if seq[i:i+l] == p1 or seq[i:i+l] == p2:
                xlist[i] += 1
    else:
        f = len(seq)
        l = len(p)
        xlist = [0] * (f - l + 1)
        for i in range(len(xlist)):
            if seq[i:i+l] == p:
                xlist[i] += 1
    return xlist


def reverse_compliment(seq):
    '''
    found out the reverse complimentary seq
    '''
    rseq = seq[::-1]
    new_seq = ''
    for i in range(len(rseq)):
        if rseq[i] == 'A':
            new_seq += 'T'
        elif rseq[i] == 'T':
            new_seq += 'A'
        elif rseq[i] == 'G':
            new_seq += 'C'
        elif rseq[i] == 'C':
            new_seq += 'G'
        elif rseq[i] == 'N':
            new_seq += 'N'
        else:
            raise Exception('There are other sgnd than just ATGC')
    return new_seq




