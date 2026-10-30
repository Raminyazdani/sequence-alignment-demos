#!/usr/bin/env python3
"""FASTA Stage 3 compatible-region chaining (init1/initn) and Stage 4 band. Section 9."""
def chain_initn(regions, J):
    """regions: list of (qstart,qend,lstart,lend,score) sorted; returns (init1, initn, chain).
    Compatible = strictly ordered and non-overlapping in BOTH sequences."""
    R=sorted(regions, key=lambda r:(r[0],r[2]))
    init1=max(r[4] for r in R)
    C=[0]*len(R); back=[-1]*len(R)
    for i,r in enumerate(R):
        C[i]=r[4]
        for k in range(i):
            q=R[k]
            if q[1]<r[0] and q[3]<r[2] and C[k]-J+r[4]>C[i]:
                C[i]=C[k]-J+r[4]; back[i]=k
    bi=max(range(len(R)), key=lambda i:C[i]); chain=[]; i=bi
    while i!=-1: chain.append(R[i]); i=back[i]
    return init1, C[bi], list(reversed(chain))
def band_cells(m,n,d0,w):
    """Indices (i,j) inside the band |（i-j)-d0|<=w (Stage 4 region)."""
    return [(i,j) for i in range(m) for j in range(n) if abs((i-j)-d0)<=w]
def demo():
    regs=[(1,5,1,5,14),(8,14,10,16,25)]  # two compatible regions, scores 14 and 25
    init1,initn,chain=chain_initn(regs, J=8)
    assert init1==25 and initn==14+25-8==31, (init1,initn)
    nb=len(band_cells(14,13,1,3)); assert nb>0
    print("fasta: init1=%d initn=14+25-8=%d ; band cells(d0=1,w=3)=%d"%(init1,initn,nb))
    return init1,initn
if __name__=="__main__": demo()
