# Changelog

## 1.2

* Experimental additions, kept apart from the validated rule and placed after the validated columns: `p_new_function` (ranking score for a function new to the record), orphan-activity suggestions (`orphan_candidates`, `orphan_names`; `--no-orphan` switches them off) and `--external NAME=FILE`, which reads the output files of CLEAN, DeepECtransformer, ProteInfer and DeepEC that you ran yourself and shows them next to the nearest enzyme.
* `selftest` also checks the new-function score, the orphan table and the external-file readers.
* Tested on Linux (Python 3.9 and 3.11) and on Windows 11 (Python 3.14, the official MMseqs2 Windows build).

## 1.1

* Fixed: FASTA identifiers of the form `sp|ACCESSION|NAME` (MMseqs2 prints only the accession in the hit table).
* Input validation with clear messages (malformed files, repeated identifiers, percentages instead of fractions, identifiers not in the FASTA); `--tmpdir`; the `reference` command; a `.params.json` file that records how a table was made.
* The risk of a three-digit error is capped at the risk of a four-digit error (no statement of the replay changed).
* A stronger `selftest`.
