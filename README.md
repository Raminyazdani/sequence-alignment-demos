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

All shipped examples use small fictional/toy strings. Tests compare global scores with an independent exhaustive-path calculation for small inputs and cover local scores, translation, coordinates, chaining, shuffling and the existing demo entrypoints. No original biological dataset, paper PDF, slide deck, report, figure or student record is included.

The source was retained from one canonical five-file code export; repeated report exports were not turned into separate projects. Parent-directory import manipulation was removed so demos use the adjacent seqcore.py. Original code remains proprietary, and third-party intellectual-property rights are retained; see LICENSE and IMPORT_NOTES.md.

