#!/usr/bin/env python

from collections import namedtuple

Feature = namedtuple('Feature', ['id', 
                                 'chromosome',
                                 'start', 
                                 'end', 
                                 'strand'])

Seqinfo = namedtuple("Seqinfo", ["sequence",
                                 "compsequence",
                                 "id",
                                 "chromosome",
                                 "start",
                                 "end",
                                 "strand",
                                 "offset"])

Sequence = namedtuple("Sequence", [])


if __name__ == "__main__":
    pass