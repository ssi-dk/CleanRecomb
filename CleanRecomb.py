#!/usr/bin/env python3

import argparse
import sys
import random


parser = argparse.ArgumentParser(
    description='Remove recombination from snp matrix')
parser.add_argument('alnfile', type=argparse.FileType(), help="alignment file")
parser.add_argument('-p', '--posfile', type=argparse.FileType(), default=None)
args = parser.parse_args()


class segment:
    ntest = 0
    totsnps = 0

    def __init__(self, h, start, end, snps):
        self.hash = h
        self.start = start
        self.end = end
        self.snps = snps
        segment.ntest += 1
        segment.totsnps += snps

    def test(self):
        try:
            statistic = 1.0*segment.ntest*(
                hashprob[self.hash]/segment.totsnps*1.0)**(self.snps-1)
        except KeyError:
            print(hashprob)
            sys.exit()
        if statistic < 0.05:
            print("Length: {}, stat: {}, ntest: {}, l: {}, hash: {}".format(
                len(self), statistic, segment.ntest, segment.totsnps,
                self.hash), file=sys.stderr)
        return statistic >= 0.05

    def __len__(self):
        return self.end-self.start

    def discard(self):
        hashprob[self.hash] -= self.snps
        segment.totsnps -= self.snps
        segment.ntest -= 1


def hash(pattern):
    """ Determines if patterns are the same.
    Only identical patterns should give the same hash """
    prev = ""
    n = 1
    hashstring = list()
    letters = list()
    for i in pattern:
        if i == prev:
            n += 1
        else:
            try:
                l = letters.index(prev)
            except ValueError:
                l = len(letters)
                letters.append(prev)
            hashstring.append((l, n))
            n = 1
        prev = i
    try:
        l = letters.index(prev)
    except ValueError:
        l = len(letters)
        letters.append(prev)
    hashstring.append((l, n))
    return tuple(hashstring[1:])


line = next(args.alnfile)
(norgs, alnlength) = line.split()
(norgs, alnlength) = (int(norgs), int(alnlength))
aln = list()
for i in range(alnlength):
    aln.append(list("N" * norgs))

forspalte = list()
org = 0
for line in args.alnfile:
    line = line.strip()
    forspalte.append(line[0:10])
    pos = 0
    for char in line[10:]:
        try:
            aln[pos][org] = char
        except IndexError:
            #print("norgs=%i\norg=%i\nalnlength=%i\ni=%i" % (
            #    norgs, org, alnlength, i), file=sys.stderr)
            pass
        pos += 1
    org += 1

""" Calculate frequencies of hash values """

hashes = [hash(profile) for profile in aln]
hashprob = {h: hashes.count(h) for h in set(hashes)}
l = len(aln)
last = hash("")

segments = list()
begin = None
snps = 1
for i in range(l):
    sh = hashes[i]
    if len(sh) == 1:  # Ignore sites without variation
        continue
    if sh == last:
        """ Current segment continues """
        snps += 1
    else:
        """ segment ends and is recorded """
        if begin:
            segments.append(segment(last, begin, i, snps))
        """ start a new segment """
        last = sh
        begin = i
        snps = 1

discardedPositions = set()
lastDiscarded = None
segments.sort(key=len, reverse=True)
for i in range(len(segments)):
    s = segments[i]
    if s.test() is False:
        s.discard()
        discardedPositions.update(range(s.start, s.end))
        print("Discarding {} SNPs over in positions [{},{}]".format(
            s.snps, s.start, s.end), file=sys.stderr)
        lastDiscarded = i
print("%i %i" % (norgs, l-len(discardedPositions)))

for i in range(norgs):
    output = forspalte[i]
    output += "".join(
        [aln[j][i] for j in range(l) if j not in discardedPositions])
    print(output)


# Read the snp positions, and output the ones that were not discarded.

if args.posfile:
    fo = open(args.posfile.name + ".clean", "w")
    i = 0
    for line in args.posfile:
        if i not in discardedPositions:
            print(line.strip(), file=fo)
        i += 1
    args.posfile.close()
    fo.close()
