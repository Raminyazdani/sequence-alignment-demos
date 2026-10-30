#!/usr/bin/env python3
"""LFASTA dot-plot signatures: true indel, tandem repeat, periodic repeat. Section 12."""
from seqcore import dotplot_offsets
def demo():
    indel=dotplot_offsets("ABCDEFGHIJ","ABCXDEFGHIJ",1)   # X inserted after ABC
    assert 0 in indel and -1 in indel, indel                # diagonal shift 0 -> -1
    rep=dotplot_offsets("ABCDABCD","ABCDABCD",1,drop_main=True)
    assert 4 in rep and -4 in rep, rep                      # tandem repeat offsets
    print("lfasta: indel diagonals %s (0->-1 shift) ; tandem-repeat offsets %s"
          %(sorted(set(indel)&{0,-1}), [x for x in rep if abs(x)==4]))
    return indel,rep
if __name__=="__main__": demo()
