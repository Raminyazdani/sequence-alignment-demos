# Sequence Alignment Teaching Demos

Five small Python files preserved from a personal teaching workspace. They demonstrate global/local sequence-alignment scoring, selected PAM250 teaching cells, k-tuple diagonals, compatible-region chaining, six-frame translation, coordinate mapping and shuffle-based toy statistics.

These are educational demonstrations, not the complete FASTA/LFASTA/TFASTA/RDF2 software distributions. The PAM250 table contains only selected cells; unlisted pairs return None. DNA helpers expect uppercase A/C/G/T. No clinical, experimental-data or publication benchmark validation is claimed.

## Run

Python 3.10 or newer; standard library only:

```sh
python fasta.py
python lfasta.py
python tfasta.py
python rdf2.py
python -m unittest discover -s tests -v
```

### A one-letter insertion

Run `python lfasta.py` to compare the invented strings `ABCDEFGHIJ` and
`ABCXDEFGHIJ`. Each matching letter before the inserted `X` has the same
position in both strings. After the insertion, the position in the first
string minus the position in the second is −1. The command also compares
`ABCDABCD` with itself. Matches between its two copies of `ABCD` lie four
positions apart in either direction.

```text
lfasta: indel diagonals [-1, 0] (0->-1 shift) ; tandem-repeat offsets [-4, 4]
```

All shipped examples use small fictional/toy strings. Tests compare global scores with an independent exhaustive-path calculation for small inputs and cover local scores, translation, coordinates, chaining, shuffling and the existing demo entrypoints. No original biological dataset, paper PDF, slide deck, report, figure or student record is included.

The source was retained from one canonical five-file code export; repeated report exports were not turned into separate projects. Parent-directory import manipulation was removed so demos use the adjacent seqcore.py. Original code remains proprietary, and third-party intellectual-property rights are retained; see LICENSE and IMPORT_NOTES.md.

