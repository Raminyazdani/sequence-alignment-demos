#!/usr/bin/env python3
"""RDF2 significance: global/local shuffle, z-score, empirical p-value. Section 11."""
import random
from seqcore import zstats, SEED
def demo():
    random.seed(SEED)
    S,mean,sd,z,mx,p = zstats("ACTGACTG","ACTGACTG","global",200)
    assert S==16 and mx<S and round(z,2)==5.10 and abs(p-1/201)<1e-9, (S,round(mean,2),round(sd,2),round(z,2),mx,p)
    random.seed(SEED); g=zstats("LLVLLVKKRKRK","LLVLLVKKRKRK","global",300)
    random.seed(SEED); l=zstats("LLVLLVKKRKRK","LLVLLVKKRKRK","local",300,6)
    assert g[3]>l[3], (g[3],l[3])
    print("rdf2: toy z=%.2f p=(1+0)/201=%.4f max=%d ; global z=%.2f vs local z=%.2f (local stricter)"
          %(z,p,mx,g[3],l[3]))
    return z,p,g[3],l[3]
if __name__=="__main__": demo()
