"""Independent small synthetic checks, not experimental-data validation."""
import ast
from collections import Counter
from functools import lru_cache
import io
import itertools
from pathlib import Path
import random
import unittest
from contextlib import redirect_stdout

import seqcore
import fasta
import lfasta
import rdf2
import tfasta


def exhaustive_global(a, b, match=1, mismatch=-1, gap=2):
    @lru_cache(None)
    def score(i, j):
        if i == len(a): return -(len(b) - j) * gap
        if j == len(b): return -(len(a) - i) * gap
        return max(score(i + 1, j + 1) + (match if a[i] == b[j] else mismatch),
                   score(i + 1, j) - gap, score(i, j + 1) - gap)
    return score(0, 0)


class AlgorithmTests(unittest.TestCase):
    def test_source_compiles(self):
        for path in Path(__file__).resolve().parents[1].rglob("*.py"):
            compile(ast.parse(path.read_text(encoding="utf-8")), str(path), "exec")

    def test_global_scores_match_independent_exhaustive_paths(self):
        words = [""] + ["".join(t) for size in (1, 2) for t in itertools.product("AC", repeat=size)]
        for a in words:
            for b in words:
                self.assertEqual(seqcore.nw_global(a, b), exhaustive_global(a, b), (a, b))

    def test_local_alignment_score_and_substrings(self):
        score, a, b = seqcore.sw_local("TTACGT", "ACGA")
        self.assertEqual(score, 6)
        self.assertEqual(a, "ACG")
        self.assertEqual(b, "ACG")
        self.assertEqual(seqcore.sw_local("AAA", "CCC")[0], 0)

    def test_translation_and_coordinates(self):
        dna = "ATGGCTTAA"
        self.assertEqual(seqcore.translate(dna), "MAX")
        self.assertEqual(seqcore.revcomp(seqcore.revcomp(dna)), dna)
        self.assertEqual(len(seqcore.six_frames(dna)), 6)
        self.assertEqual(seqcore.codon_span(-1, 1, len(dna)), (7, 9, "-"))

    def test_chaining_band_and_shuffles(self):
        regions = [(1, 3, 1, 3, 7), (5, 7, 5, 7, 8)]
        self.assertEqual(fasta.chain_initn(regions, 2)[:2], (8, 13))
        self.assertTrue(all(abs(i-j) <= 1 for i, j in fasta.band_cells(4, 4, 0, 1)))
        random.seed(7)
        self.assertEqual(Counter(seqcore.shuffle_whole("AACCGGTT")), Counter("AACCGGTT"))
        self.assertEqual(seqcore.pam("L", "I"), seqcore.pam("I", "L"))
        self.assertIsNone(seqcore.pam("A", "A"))  # selected teaching cells, not a full matrix

    def test_existing_synthetic_demo_entrypoints(self):
        with redirect_stdout(io.StringIO()):
            fasta.demo()
            lfasta.demo()
            tfasta.demo()
            rdf2.demo()


if __name__ == "__main__":
    unittest.main()

