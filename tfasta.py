#!/usr/bin/env python3
"""TFASTA: six-frame translation, coordinate back-mapping, protein-space scoring. Section 10."""
from seqcore import revcomp, six_frames, codon_span, sw_local
def demo():
    dna="ATGGCTGATAAATGG"; L=len(dna)
    assert revcomp(dna)=="CCATTTATCAGCCAT"
    fr=six_frames(dna)
    assert fr['+1']=="MADKW" and fr['+2']=="WLIN" and fr['+3']=="GXXM"
    assert fr['-1']=="PFISH" and fr['-2']=="HLSA" and fr['-3']=="IYQP"
    assert fr['+3'].count('X')==2 and fr['+1'].count('X')==0
    assert codon_span(1,1,L)==(1,3,'+') and codon_span(1,5,L)==(13,15,'+')
    assert codon_span(-1,1,L)==(13,15,'-')
    scores={k:sw_local("MADKW",fr[k].replace("X","Z"))[0] for k in fr}
    assert scores['+1']==10 and all(v<=2 for k,v in scores.items() if k!='+1')
    print("tfasta: frames %s ; +1 query score %d others<=2 ; back-map +1 M->[1,3] -1 P->[13,15]rev"
          %({k:fr[k] for k in ['+1','+2','+3','-1','-2','-3']}, scores['+1']))
    return fr,scores
if __name__=="__main__": demo()
