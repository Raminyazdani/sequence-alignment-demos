#!/usr/bin/env python3
"""
Shared educational sequence-alignment primitives. Standard library only.

The source-only snapshot includes synthetic offline tests under tests/.
PAM250 below is a selected-cell teaching table, not a complete scoring matrix.
No clinical, experimental or external scientific validation is claimed.
"""
from __future__ import annotations
import random
from collections import Counter
from typing import Dict, List, Tuple

SEED = 20260628

# --- genetic code + reverse complement ------------------------------------------------
CODON = {
 'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
 'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
 'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
 'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
 'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
 'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
 'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
 'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}
_COMP = {'A':'T','T':'A','G':'C','C':'G'}

def revcomp(s: str) -> str:
    """Reverse complement of a DNA string."""
    return ''.join(_COMP[c] for c in reversed(s))

def translate(s: str, stop: str = 'X') -> str:
    """Translate a DNA string in frame 0; stop codons become `stop` (default 'X')."""
    out = []
    for i in range(0, len(s) - 2, 3):
        aa = CODON[s[i:i+3]]
        out.append(stop if aa == '*' else aa)
    return ''.join(out)

def six_frames(dna: str) -> Dict[str, str]:
    """All six reading-frame translations (+1/+2/+3 forward, -1/-2/-3 on revcomp)."""
    rc = revcomp(dna)
    return {'+1': translate(dna), '+2': translate(dna[1:]), '+3': translate(dna[2:]),
            '-1': translate(rc),  '-2': translate(rc[1:]),  '-3': translate(rc[2:])}

def codon_span(frame: int, p: int, L: int) -> Tuple[int, int, str]:
    """1-based inclusive DNA codon span for protein position p in `frame`.
    Reverse frames are returned as forward-strand coordinates with strand '-'."""
    if frame > 0:
        start = frame + 3*(p-1)
        return (start, start+2, '+')
    f = -frame
    rcs = f + 3*(p-1)
    return (L - (rcs+2) + 1, L - rcs + 1, '-')

# --- dynamic programming: Needleman-Wunsch (global), Smith-Waterman (local) ------------
def nw_global(A: str, B: str, match: int = 1, mismatch: int = -1, gap: int = 2) -> int:
    """Optimal global-alignment score (linear gap). Tie order diag>up>left."""
    m, n = len(A), len(B)
    H = [[0]*(n+1) for _ in range(m+1)]
    for i in range(m+1): H[i][0] = -i*gap
    for j in range(n+1): H[0][j] = -j*gap
    for i in range(1, m+1):
        for j in range(1, n+1):
            s = match if A[i-1] == B[j-1] else mismatch
            H[i][j] = max(H[i-1][j-1]+s, H[i-1][j]-gap, H[i][j-1]-gap)
    return H[m][n]

def sw_local(A: str, B: str, match: int = 2, mismatch: int = -1, gap: int = 2
             ) -> Tuple[int, str, str]:
    """Best local score and the gapped local alignment (zero floor; tie diag>up>left)."""
    m, n = len(A), len(B)
    H = [[0]*(n+1) for _ in range(m+1)]
    ptr = [[None]*(n+1) for _ in range(m+1)]
    best = bi = bj = 0
    for i in range(1, m+1):
        for j in range(1, n+1):
            s = match if A[i-1] == B[j-1] else mismatch
            diag, up, left = H[i-1][j-1]+s, H[i-1][j]-gap, H[i][j-1]-gap
            val = max(0, diag, up, left)
            H[i][j] = val
            if val > 0:
                ptr[i][j] = 'D' if val == diag else ('U' if val == up else 'L')
            if val > best:
                best, bi, bj = val, i, j
    a, b = [], []
    i, j = bi, bj
    while i > 0 and j > 0 and H[i][j] > 0:
        d = ptr[i][j]
        if d == 'D':   a.append(A[i-1]); b.append(B[j-1]); i -= 1; j -= 1
        elif d == 'U': a.append(A[i-1]); b.append('-');    i -= 1
        elif d == 'L': a.append('-');    b.append(B[j-1]); j -= 1
        else: break
    return best, ''.join(reversed(a)), ''.join(reversed(b))

# --- PAM250 (selected verified Dayhoff MDM78 log-odds cells) ---------------------------
PAM250 = {('W','W'):17, ('C','C'):12, ('L','I'):2, ('K','R'):3, ('A','G'):1,
          ('Y','Y'):10, ('F','F'):9, ('C','W'):-8}

def pam(a: str, b: str) -> int:
    """Symmetric PAM250 log-odds lookup for the verified cell set."""
    return PAM250.get((a, b), PAM250.get((b, a)))

# --- k-tuple diagonals (Wilbur-Lipman) ------------------------------------------------
def ktup_diagonals(Q: str, L: str, ktup: int = 2) -> Dict[int, int]:
    """Map diagonal d=i-j to the number of length-`ktup` word hits on it."""
    pos: Dict[str, List[int]] = {}
    for j in range(len(L)-ktup+1):
        pos.setdefault(L[j:j+ktup], []).append(j)
    c: Counter = Counter()
    for i in range(len(Q)-ktup+1):
        for j in pos.get(Q[i:i+ktup], []):
            c[i-j] += 1
    return dict(c)

# --- RDF/RDF2 shuffling -----------------------------------------------------------------
def shuffle_whole(s: str) -> str:
    l = list(s); random.shuffle(l); return ''.join(l)

def shuffle_local(s: str, w: int = 6) -> str:
    l = list(s)
    for k in range(0, len(l), w):
        blk = l[k:k+w]; random.shuffle(blk); l[k:k+w] = blk
    return ''.join(l)

def zstats(Q: str, T: str, mode: str, N: int, w: int = 6):
    """Return (S_obs, mean, sd, z, max_shuffle, p_hat) for N shuffles of Q vs T.
    p_hat = (1 + exceedances)/(N+1); exceedances = #{shuffle score >= S_obs}."""
    S = sw_local(Q, T)[0]
    sc = []
    for _ in range(N):
        Qr = shuffle_whole(Q) if mode == 'global' else shuffle_local(Q, w)
        sc.append(sw_local(Qr, T)[0])
    mean = sum(sc)/len(sc)
    sd = (sum((x-mean)**2 for x in sc)/len(sc))**0.5
    z = (S-mean)/sd if sd > 0 else float('inf')
    exc = sum(1 for x in sc if x >= S)
    return S, mean, sd, z, max(sc), (1+exc)/(N+1)

# --- dot-plot diagonals -----------------------------------------------------------------
def dotplot_offsets(A: str, B: str, ktup: int = 1, drop_main: bool = False):
    """Sorted set of diagonals d=i-j carrying a length-`ktup` match."""
    offs = set()
    for i in range(len(A)-ktup+1):
        for j in range(len(B)-ktup+1):
            if A[i:i+ktup] == B[j:j+ktup]:
                if not (drop_main and i == j):
                    offs.add(i-j)
    return sorted(offs)
