# Changelog

## 1.2.1

Input handling only; the statements and scores of 1.2 are reproduced exactly (byte-identical tables on 400 mutated reference proteins and on the 30 bundled examples, in `cutoff` and `match` mode, and `score` on a hit table written by 1.2).

* The search sees synthetic identifiers and the result carries yours: non-ASCII or emoji identifiers (they failed or crashed on a console with a narrow code page), NCBI-style `lcl|` and `gi|` identifiers (failed after the search) and two identifiers with the same accession such as `Q3YZV5` and `sp|Q3YZV5|NAME` (silently merged) now give the result of the same protein under a plain identifier. A hit table written with `--keep-hits` carries the identifiers of your FASTA.
* Refused before the search: an output path or `--keep-hits` folder that cannot be written, `--threads` below 1, a `--year` outside 1995–2100.
* `score`: a hit table whose lines are repeated exactly is read once, with a warning; hits to targets that are not reference enzymes are reported; two query identifiers with one accession are reported.
* On Windows a failed MMseqs2 search points to the BusyBox helpers that the official build needs.
* Python 3.9 or newer, as stated in the metadata (the program had still said 3.8).
* The README describes only what the repository contains; the tests of input handling and an end-to-end test with the real search were added, and the automatic checks now also run Python 3.14, the search itself with MMseqs2 18-8cc5c, and build, `twine check` and a wheel install.

## 1.2

* Requires Python 3.9 or newer (tested with 3.9, 3.11, 3.13 and 3.14).
* Experimental additions, kept apart from the validated rule and placed after the validated columns: `p_new_function` (ranking score for a function new to the record), orphan-activity suggestions (`orphan_candidates`, `orphan_names`; `--no-orphan` switches them off) and `--external NAME=FILE`, which reads the output files of CLEAN, DeepECtransformer, ProteInfer and DeepEC that you ran yourself and shows them next to the nearest enzyme.
* `selftest` also checks the new-function score, the orphan table and the external-file readers.
* Tested on Linux (Python 3.9 and 3.11) and on Windows 11 (Python 3.14, the official MMseqs2 Windows build).

## 1.1

* Fixed: FASTA identifiers of the form `sp|ACCESSION|NAME` (MMseqs2 prints only the accession in the hit table).
* Input validation with clear messages (malformed files, repeated identifiers, percentages instead of fractions, identifiers not in the FASTA); `--tmpdir`; the `reference` command; a `.params.json` file that records how a table was made.
* The risk of a three-digit error is capped at the risk of a four-digit error (no statement of the replay changed).
* A stronger `selftest`.
