###############################################################################
# Date: 28/1/23
# Author: Sebastian Porras
# Aims: Adapted code from the Binfpy library to read and write FASTA seqs.
# Major difference is that Alphabets have been changed to ones compatible with
# GRASP and that Sequence.alphabet is now the alphabet name rather than list
# of all possible letters in that alphabet.
###############################################################################

import re
from sym import *
import math


class Sequence(object):
    """ A biological sequence. Stores the sequence itself (as a compact array), 
    the alphabet (i.e., type of sequence it is), and optionally a name and further 
    information. """

    def __init__(self, sequence, alphabet=None, name: str = '', info: str = '', gappy=False):
        """ Create a sequence with the sequence data. Specifying the alphabet,
        name and other information about the sequence are all optional.
        The sequence data is immutable (stored as a string).

        Parameters:
            sequence(list): The array of symbols that make up the sequence

            alphabet(Alphabet): The alphabet from which symbols come

            name(str): The name (identifier) of a sequence

            info(str): Other information (free text; e.g. annotations)

            length(int): The number of symbols that the sequence is composed of

            gappy(bool): True if the sequence has "gaps", i.e. positions that 

            represent deletions relative another sequence

        Example:
        >>> myseq = Sequence('MVSAKKVPAIAMSFGVSF')
        will create a sequence with no name, and assign one of the predefined
        alphabets on the basis of what symbols were used.
        >>> myseq.alphabet.symbols
        will output the standard protein alphabet:
        ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q',
        'R', 'S', 'T', 'V', 'W', 'Y'] 
        """

        self.sequence = sequence

        # Assign an alphabet
        # If no alphabet is provided, attempts to identify the alphabet from sequence
        self.alphabet = None
        if not alphabet is None:
            for sym in self.sequence:
                # error check: bail out
                if not sym in alphabet and (sym != '-' or not gappy):
                    raise RuntimeError(
                        'Invalid symbol: %c in sequence %s' % (sym, name))
            self.alphabet = alphabet
        else:
            for alphaName in preferredOrder:

                alpha = predefAlphabets[alphaName]
                valid = True

                for sym in self.sequence:
                    if not sym in alpha and (sym != '-' or not gappy):
                        valid = False
                        break
                if valid:
                    self.alphabet = alphaName
                    break

            if self.alphabet is None:
                raise RuntimeError(
                    'Could not identify alphabet from sequence: %s' % name)

        # Store other information
        self.name = name
        self.info = info
        self.length = len(self.sequence)
        self.gappy = gappy

    def __len__(self):
        """ Defines what the "len" operator returns for an instance of Sequence, e.g.
        >>> seq = Sequence('ACGGTAGGA', DNA_Alphabet)
        >>> print (len(seq))
        9
        """
        return len(self.sequence)

    def __str__(self):
        """ Defines what should be printed when the print statement is used on a Sequence instance """
        str = self.name + ': '
        for sym in self:
            str += sym
        return str

    def __iter__(self):
        """ Defines how a Sequence should be "iterated", i.e. what its elements are, e.g.
        >>> seq = Sequence('AGGAT', DNA_Alphabet)
        >>> for sym in seq:
                print (sym)
        will print A, G, G, A, T (each on a separate row)
        """
        tsyms = tuple(self.sequence)
        return tsyms.__iter__()

    def __contains__(self, item):
        """ Defines what is returned when the "in" operator is used on a Sequence, e.g.
        >>> seq = Sequence('ACGGTAGGA', DNA_Alphabet)
        >>> print ('T' in seq)
        True
            which is equivalent to 
        >>> print (seq.__contains__('T'))
        True
        >>> print ('X' in seq)
        False
        """
        for sym in self.sequence:
            if sym == item:
                return True
        return False

    def __getitem__(self, ndx):
        """ Retrieve a specified index (or a "slice" of indices) of the sequence data.
            Calling self.__getitem__(3) is equivalent to self[3] 
        """
        if type(ndx) is slice:
            return ''.join(self.sequence[ndx])
        else:
            return self.sequence[ndx]

    def writeFasta(self):
        """ Write one sequence in FASTA format to a string and return it. """
        if parseDefline(self.info)[0] == self.name:  # this sequence was previously "parsed" and info should hold the original header
            fasta = '>' + self.info + '\n'
        else:
            fasta = '>' + self.name + ' ' + self.info + '\n'
        data = ''.join(self.sequence)
        nlines = int(math.ceil((len(self.sequence) - 1) / 60 + 1))
        for i in range(nlines):
            lineofseq = ''.join(data[i*60: (i+1)*60]) + '\n'
            fasta += lineofseq
        return fasta

    def getDegapped(self):
        """ Create the sequence excluding gaps, and provide the corresponding indices for the gapped version, e.g.
        >>> gappy = Sequence('AC--TA-GA', DNA_Alphabet, name = 'myseq', gappy = True)
        >>> degapped, indices = gappy.getDegapped()
        >>> print(degapped)
            myseq: ACTAGA
        >>> print(indices)
            [0, 1, 4, 5, 7, 8]
        """
        idxs = []
        newseq = []
        for i in range(len(self.sequence)):
            if not self.sequence[i] == '-':
                newseq.append(self.sequence[i])
                idxs.append(i)
        return Sequence(newseq, self.alphabet, self.name, self.info, gappy=False), idxs

    def find(self, findme, gappy=False):
        """ Find the position of the specified symbol or sub-sequence """
        if gappy == False or self.gappy == False:
            return ''.join(self.sequence).find(findme)
        else:  # if the sequence is gappy AND the function is called with gappy = True THEN run the find on the de-gapped sequence
            degapped, idxs = self.getDegapped()
            idx = ''.join(degapped).find(findme)
            return idxs[idx] if idx >= 0 else -1


"""
Below are some useful methods for loading data from strings and files.
Recognize the FASTA format (nothing fancy).
"""


def readFasta(string, alphabet=None, ignore=False, gappy=False, parse_defline=True):
    """ Read the given string as FASTA formatted data and return the list of
        sequences contained within it.
        If alphabet is specified, use it, if None (default) then guess it.
        If ignore is False, errors cause the method to fail.
        If ignore is True, errors will disregard sequence.
        If gappy is False (default), sequence cannot contain gaps,
        if True gaps are accepted and included in the resulting sequences.
        If parse_defline is False, the name will be set to everything before the first space, else parsing will be attempted."""
    seqlist = []    # list of sequences contained in the string
    seqname = None  # name of *current* sequence
    seqinfo = None
    seqdata = []    # sequence data for *current* sequence
    for line in string.splitlines():    # read every line
        if len(line) == 0:              # ignore empty lines
            continue
        if line[0] == '>':  # start of new sequence
            if seqname:     # check if we've got one current
                try:
                    current = Sequence(seqdata, alphabet,
                                       seqname, seqinfo, gappy)
                    seqlist.append(current)
                except RuntimeError as errmsg:
                    if not ignore:
                        raise RuntimeError(errmsg)
            # now collect data about the new sequence
            seqinfo = line[1:].split()  # skip first char
            if len(seqinfo) > 0:
                try:
                    if parse_defline:
                        parsed = parseDefline(seqinfo[0])
                        seqname = parsed[0]
                        seqinfo = line[1:]
                    else:  # we are not parsing the sequence name so no need to duplicate it in the info
                        seqname = seqinfo[0]
                        if len(seqinfo) > 0:  # more than a name
                            edited_info = ''
                            for infopart in seqinfo[1:]:
                                edited_info += infopart + ' '
                            seqinfo = edited_info
                        else:
                            seqinfo = ''
                except IndexError as errmsg:
                    if not ignore:
                        raise RuntimeError(errmsg)
            else:
                seqname = ''
                seqinfo = ''
            seqdata = []
        else:               # we assume this is (more) data for current
            cleanline = line.split()
            for thisline in cleanline:
                seqdata.extend(tuple(thisline.strip('*')))
    # we're done reading the file, but the last sequence remains
    if seqname:
        try:
            lastseq = Sequence(seqdata, alphabet, seqname, seqinfo, gappy)
            seqlist.append(lastseq)
        except RuntimeError as errmsg:
            if not ignore:
                raise RuntimeError(errmsg)
    return seqlist


def parseDefline(string):
    """ Parse the FASTA defline (see http://en.wikipedia.org/wiki/FASTA_format)
        GenBank, EMBL, etc                gi|gi-number|gb|accession|locus
        SWISS-PROT, TrEMBL                sp|accession|name
        ...
        Return a tuple with
        [0] primary search key, e.g. UniProt accession, Genbank GI
        [1] secondary search key, e.g. UniProt name, Genbank accession
        [2] source, e.g. 'sp' (SwissProt/UniProt), 'tr' (TrEMBL), 'gb' (Genbank)
    """
    if len(string) == 0:
        return ('', '', '', '')
    s = string.split()[0]
    if re.match("^sp\|[A-Z][A-Z0-9]*\|\S+", s):
        arg = s.split('|')
        return (arg[1], arg[2], arg[0], '')
    elif re.match("^tr\|[A-Z][A-Z0-9]*\|\S+", s):
        arg = s.split('|')
        return (arg[1], arg[2], arg[0], '')
    elif re.match("^gi\|[0-9]*\|\S+\|\S+", s):
        arg = s.split('|')
        return (arg[1], arg[3], arg[0], arg[2])
    elif re.match("gb\|\S+\|\S+", s):
        arg = s.split('|')
        return (arg[1], arg[2], arg[0], '')
    elif re.match("emb\|\S+\|\S+", s):
        arg = s.split('|')
        return (arg[1], arg[2], arg[0], '')
    elif re.match("^refseq\|\S+\|\S+", s):
        arg = s.split('|')
        return (arg[1], arg[2], arg[0], '')
    elif re.match("[A-Z][A-Z0-9]*\|\S+", s):
        arg = s.split('|')
        return (arg[0], arg[1], 'UniProt', '')  # assume this is UniProt
    else:
        return (s, '', '', '')


def readFastaFile(filename, alphabet=None, ignore=False, gappy=False, parse_defline=True):
    """ Read the given FASTA formatted file and return the list of sequences
        contained within it. Note that if alphabet is NOT specified, it will take a
        separate guess for each sequence.
        If ignore is False, errors cause the method to fail.
        If ignore is True, errors will disregard sequence.
        If gappy is False (default), sequence cannot contain gaps,
        if True gaps are accepted and included in the resulting sequences.
        If parse_defline is False, the name will be set to everything before the first space, else parsing will be attempted."""
    fh = open(filename)
    seqlist = []
    batch = ''  # a batch of rows including one or more complete FASTA entries
    rowcnt = 0
    for row in fh:
        row = row.strip()
        if len(row) > 0:
            if row.startswith('>') and rowcnt > 0:
                more = readFasta(batch, alphabet, ignore, gappy, parse_defline)
                if len(more) > 0:
                    seqlist.extend(more)
                batch = ''
                rowcnt = 0
            batch += row + '\n'
            rowcnt += 1
    if len(batch) > 0:
        more = readFasta(batch, alphabet, ignore, gappy, parse_defline)
        if len(more) > 0:
            seqlist.extend(more)
    fh.close()
    return seqlist


def writeFastaFile(filename, seqs):
    """ Write the specified sequences to a FASTA file. """
    fh = open(filename, 'w')
    for seq in seqs:
        fh.write(seq.writeFasta())
    fh.close()
