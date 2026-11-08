# ecbackoff

**Inherit an enzyme's EC number only where the neighbourhood of its nearest relative supports it; otherwise say three digits, or nothing.**

Most enzyme annotation still inherits the EC number of the nearest characterized relative when sequence identity passes a threshold, by convention 60% for the complete EC number and 40% for its first three digits. `ecbackoff` is an alternative to that global rule. For every query protein it looks at the *neighbourhood* of the nearest characterized enzyme (how close it is, how many relatives there are, whether they agree, how many different functions they carry, how common the nearest function is), estimates the probability that inheriting would be wrong, and then says one of three things:

| statement | meaning |
|---|---|
| `EC4` | inherit the complete four-digit EC number(s) of the nearest enzyme |
| `EC3` | the sequence is probably in the same sub-subclass (first three digits) but the full number should not be inherited |
| `none` | the neighbourhood does not support a statement; do not inherit |

It is a rule for deciding whether to inherit, not a function predictor. It needs no deep-learning framework and no GPU: Python standard library plus one MMseqs2 search.

Version 1.2 adds three **experimental** features, kept apart from the validated rule (see "Experimental features"): a ranking score for a function that is new to the record, suggestions of orphan enzyme activities, and the import of the output of CLEAN, DeepECtransformer, ProteInfer or DeepEC that you have run yourself.

Version 1.2.1 changes only how queries are handled (identifiers of any form, early checks of options and output paths, repeated hit lines) and gives the same statements and scores as 1.2 (see "Checks").

The method and its evaluation are described in a manuscript that is in preparation.

## Install and run

```bash
pip install ecbackoff
ecbackoff selftest                       # checks the files and the logic, no MMseqs2 needed
ecbackoff run -q my_proteins.fasta -o result.tsv --threads 8
```

The development version: `pip install git+https://github.com/Eagan-lau/ecbackoff`. From a clone without installing: `PYTHONPATH=src python -m ecbackoff selftest`.

Requirements: Python 3.9 or newer (tested with 3.9, 3.11, 3.13 and 3.14); [MMseqs2](https://github.com/soedinglab/MMseqs2) on the `PATH` (or given with `--mmseqs`) for the `run` command (the model was built and checked with 18.8cc5c, commit 8cc5c; another version gives a warning). MMseqs2 is not installed by pip (for example `conda install -c bioconda mmseqs2`, or the release page). The package has no Python dependencies and contains the model and the reference data (about 6 MB).

Tested on Linux (`selftest` with Python 3.9 and 3.11; the end-to-end checks below with Python 3.11 and MMseqs2 18.8cc5c) and on Windows 11 (Python 3.14 and the official Windows build of MMseqs2 18-8cc5c, `mmseqs-win64.zip`; the program, its self-test and the command-line tests passed there, and the Windows search gave the same hit table and the same result table for the 30 bundled examples as Linux; paths with spaces and Chinese characters worked). Not tested on macOS.

**MMseqs2 on Windows.** Unpack `mmseqs-win64.zip`, then run `mmseqs\mmseqs.bat` once: it installs the BusyBox helper programs that MMseqs2's search workflows need (it may ask for administrator rights, or run `mmseqs\bin\busybox.exe --install -s mmseqs\bin` yourself). Without this step `run` stops with "Failed to execute ... easysearch.sh". Then put `mmseqs\bin` on the `PATH`, or pass `--mmseqs C:\path\to\mmseqs\bin\mmseqs.exe`.

Size and speed: a run needs room for the decompressed reference (about 10 MB) and the MMseqs2 work files (`--tmpdir` chooses where). A run on 400 proteins with 8 threads took about 7 seconds on a Linux cluster node and about 27 seconds on a Windows PC, search included.

If you prefer to run the search yourself:

```bash
ecbackoff reference -o reference.fasta   # the bundled reference sequences
mmseqs easy-search my_proteins.fasta reference.fasta hits.m8 tmp -s 7.5 -e 1e-3 --max-seqs 20000 --alignment-mode 3 \
      --format-output query,target,fident,qcov,tcov,evalue,bits
ecbackoff score -q my_proteins.fasta --hits hits.m8 -o result.tsv
```

The search settings matter: the model was trained on hit tables from exactly this search. Hit tables from BLAST, DIAMOND or other MMseqs2 settings have not been validated.

**Identifiers.** The identifier is the first word of each FASTA header, and any characters are accepted. With `run` the search itself sees synthetic identifiers (`q0000001`, ...) and the result carries yours, so non-ASCII characters, UniProt-style (`sp|P12345|NAME`) and NCBI-style (`gi|123|gb|AAA12345.1|`) identifiers, and two identifiers with the same accession all behave alike. With `score` the hit table comes from your own search: MMseqs2 prints only the accession of a UniProt-style header (`sp|P12345|NAME` appears as `P12345`), which the program accepts, but two identifiers sharing an accession cannot then be told apart and the program warns. Hits to the query itself (same identifier, or the accession inside a UniProt-style identifier) are dropped unless you give `--keep-self`; this matters only if your queries include enzymes that are in the reference.

**Input checks.** The program stops with a message, before searching, on a missing or empty file, a file that is not FASTA, a repeated identifier, an output path that could not be written, `--threads` below 1 or an implausible `--year`; and on a hit table with the wrong number of columns, percentages instead of fractions, or identifiers that are not in the FASTA. A hit table whose lines are repeated exactly is read once, with a warning. It tolerates a byte-order mark, Windows line ends, lower case and wrapped lines. Empty sequences are not searched and get `none` with a note.

**Record of a run.** With `-o result.tsv` the program also writes `result.tsv.params.json`: version, command line, mode and cutoffs, prediction year, MMseqs2 version and settings, checksums of the inputs, counts, and the absolute path of the query file. `--no-params` switches this off.

## Output

One row per query, in the order of the FASTA, tab-separated. Example (first rows of the output for 400 proteins sampled from proteins that have a neighbour, cutoff mode):

| query | nearest_enzyme | identity | p_wrong_4digit | p_wrong_3digit | statement | statement_ec | identity_rule |
|---|---|---|---|---|---|---|---|
| Q3YZV5 | P52006 | 0.985 | 0.0036 | 0.0002 | EC4 | 3.6.1.12;3.6.1.23;3.6.1.9 | EC4 |
| P68068 | Q8WWM9 | 0.256 | 0.1733 | 0.0525 | EC3 | 1.15.1 | none |
| P0DD38 | P0C0J0 | 1.000 | 0.0026 | 0.0017 | EC4 | 3.4.22.10 | EC4 |
| P41498 | P40347 | 0.444 | 0.0769 | 0.0501 | EC4 | 3.1.3.2;3.1.3.48 | EC3 |

(Columns `nearest_ec`, `n_neighbours` and `note` left out here. The 400 proteins of that file were sampled from proteins that have a neighbour, so every one of them has a statement or a probability; in a real batch, proteins without a neighbour get `none`. Rows 2 and 4 show where the rule differs from the identity rule: a protein at 26% identity receives the sub-subclass where the identity rule says nothing, and one at 44% identity receives the full number where the identity rule says three digits.)

| column | content |
|---|---|
| `query` | identifier |
| `nearest_enzyme` | accession of the nearest reference enzyme (neighbours need coverage of at least 70% of both sequences) |
| `identity` | identity to it |
| `nearest_ec` | its EC number(s) |
| `n_neighbours` | number of qualifying neighbours |
| `p_wrong_4digit` | estimated probability that none of the nearest enzyme's EC numbers is among the protein's |
| `p_wrong_3digit` | estimated probability that the first three digits of the nearest enzyme's EC numbers share nothing with the protein's (never above `p_wrong_4digit`, because a three-digit error is also a four-digit error) |
| `statement` | `EC4`, `EC3` or `none` |
| `statement_ec` | the EC number(s) for `EC4`; the sub-subclass(es) for `EC3` |
| `identity_rule` | what the plain 60%/40% identity rule would say, for comparison |
| `note` | why there is no statement, if there is none |
| `p_new_function` | *experimental*: ranking score for a function new to the record (see below); empty without a neighbour |
| `orphan_candidates` | *experimental*: suggested orphan activities, each as `EC [shared digits, similarity]`, separated by a vertical bar; empty unless `p_new_function` is at least 0.5 |
| `orphan_names` | *experimental*: their names and the year each EC number was created, in the same order |
| `<name>_ec`, `<name>_conf`, `<name>_vs_nearest` | *experimental*, only with `--external`: what another predictor said and how it relates to the nearest enzyme |

A protein with no neighbour covering 70% of both sequences gets `none` and empty probabilities.

## Modes

* `--mode cutoff` (default). `EC4` if `p_wrong_4digit` is at most `--cutoff4` (default 0.10); otherwise `EC3` if `p_wrong_3digit` is at most `--cutoff3` (default 0.10); otherwise `none`. Works for any batch size, including a single protein. The defaults were chosen after looking at the replay, so the numbers below for this mode are descriptive.
* `--mode match`. Gives exactly as many `EC4` and `EC3` statements as the identity rule gives on this batch (`--identity-hi 0.6`, `--identity-lo 0.4`) and chooses which proteins receive them by risk. This is the comparison made in the paper, and the headline numbers below are from it; it is meaningful only for batches of several hundred proteins, because for a very small batch it collapses to the identity rule.

Use the probabilities for ranking and for choosing your own cutoffs; see the limits before relying on their absolute values. `--year` changes the prediction year, which enters two features (years since the nearest enzyme was established, and the number of earlier proteins carrying its EC number); values other than the frozen 2027 are outside what was validated.

## Experimental features (version 1.2)

These additions are **not part of the validated rule**. They come after the validated columns, `--no-orphan` switches the suggestions off, and the evidence for each is weaker than for the statements above.

**1. `p_new_function`, a ranking score for a function that is new to the record.** The probability, from a logistic regression on the same 18 features, that no earlier characterized protein carries any of the protein's EC numbers. Use it to rank, not as a probability. In the replay it separated proteins with a new function from the rest with an AUC of 0.867 against 0.771 for identity, and 65.9% of the top-ranked tenth of each year carried a new function against 19.6% overall. In the external check on present-day proteins the gain over identity was small and not shown (0.679 against 0.659), and the probabilities were three to seven times too high.

**2. `orphan_candidates` and `orphan_names`, suggestions of orphan activities.** For proteins with `p_new_function` of at least 0.5 (`--orphan-min-p`), the program lists the five (`--orphan-top`) orphan enzyme activities that rank highest. Orphan activities are complete EC numbers in ExplorEnz that no characterized enzyme of the reference carries (1,687 at the time of the reference). The ranking is the rule of the paper: first by the EC digits shared with the nearest enzyme, then by the similarity of the reactions (Rhea reaction fingerprints; `NA` when the EC number has no reaction with a structure). Each suggestion shows `[shared digits, similarity]`; `orphan_names` gives the names and the year each EC number was created. In the forward test of the paper (640 new-function enzymes with an earlier neighbour whose function was an orphan activity) the true activity was among the first ten suggestions for 50.0% (45.2–55.1; first 16.9%, among the first five 38.1%), against 0.57% by chance. Over all proteins flagged by the score, only 7.9% had the true orphan activity among ten (first 2.1%), because most flagged proteins are not orphan-type. It is a way to generate hypotheses to test, not a prediction.

**3. `--external NAME=FILE`, the output of other predictors that you ran yourself.** The program does not run, install or include CLEAN, DeepECtransformer, ProteInfer or DeepEC: they need their own environments (PyTorch or TensorFlow of specific versions, downloaded weights, often a GPU) and have their own licences. Run a tool as its authors describe, then give its output file:

| tool | `--external` | file to give |
|---|---|---|
| CLEAN | `clean=` | the maximum-separation CSV (`results/<name>_maxsep.csv`: identifier, then `EC:x.x.x.x/distance` fields) |
| DeepECtransformer | `deepectransformer=` | `DeepECv2_result.txt` (identifier and prediction, tab-separated) |
| ProteInfer | `proteinfer=` | its tab-separated output (sequence name, predicted label, confidence) |
| DeepEC | `deepec=` | `DeepEC_Result.txt` |

```bash
ecbackoff run -q my_proteins.fasta -o result.tsv --external clean=my_maxsep.csv --external deepec=DeepEC_Result.txt
```

For every tool three columns are added: `<name>_ec` (the EC numbers it predicted), `<name>_conf` (its confidence for each, as the tool writes it; CLEAN writes distances, where smaller is closer) and `<name>_vs_nearest` (`same_4digits`, `same_3digits`, `differs`, `no_tool_answer` when the tool gave no EC number, or `no_neighbour`). JSON lines with `query_id`, `labels` and `confidence` are read too. Identifiers must be those of your FASTA (the accession of a UniProt-style identifier is accepted, and DeepEC's `_SEPARATED_SEQUENCE_` suffix is removed). The comparison is descriptive: agreement with the nearest enzyme has not been validated as a sign of correctness, and these tools were not better than the best sequence hit in the paper (none returned the correct complete EC number more often than the best hit on enzymes described after their training; CLEAN was not shown to differ from its own training set, the other three were clearly lower). **The readers were tested with small files written in these formats, not with output of the real tools, which were not run for this version**; if a tool writes a different format, the program says so or reads no EC numbers, and `--external` with a table whose first column is the identifier and whose other columns hold EC numbers usually works.

## Evidence

All numbers are retrospective. They come from replays of UniProtKB/Swiss-Prot in forward time in which each protein was annotated using only what was known before it, with models trained on earlier data only, and the statements were scored against the experimentally characterized function. Intervals are 95% intervals from resampling sequence clusters.

* **Main replay, 2000–2024, 10,559 proteins with a characterized neighbour (one 2026 snapshot dated by the papers cited).** At the same number of statements as the 60%/40% identity rule (4,013 four-digit and 3,255 three-digit statements; `match` mode): 3.1% (2.6–3.6) of the model's statements were wrong against 7.5% (5.5–9.9) for the identity rule (difference −4.4 points, −6.9 to −2.5; four digits 3.5% against 9.7%, three digits 2.6% against 4.8%); 3.0 more correct statements per 100 proteins. At the error rate of the 60% rule the model could give a four-digit statement to 61% of the proteins instead of 38%.
* **Release-clock replay with the labels recorded at the time (3,572 proteins, 2016–2026).** 7.8% (6.3–9.3) wrong against 15.4% (11.8–19.1) for the identity rule, same number of statements.
* **The 60% rule itself was about as accurate as intended:** 9.7% of its four-digit statements were wrong. The gain is fewer errors, or more statements at the same error rate, not a repair of a miscalibrated rule.
* **How much is family information.** An identity rule calibrated within Pfam families on earlier proteins (our own implementation, not the software of any earlier method) was wrong in 4.3% (3.7–4.9) in the main replay. Calibration within families gave 49–80% of the gain over the global thresholds, depending on the replay. The neighbourhood model's further gain over that rule was 1.2 points (0.5–1.8) in the main replay and 3.9 on the release clock, and was **not shown** (0.45 points, −0.96 to 0.03) after removing the 32.6% of proteins whose own or nearest enzyme's EC number was created after the test year. Treat the neighbourhood's advantage over a family-calibrated rule as the least secure result.
* **Default cutoffs in the main replay** (post hoc): `EC4` for 41.3% of proteins with 3.9% wrong; `EC3` for a further 41.9% with 4.5% wrong; `none` for 16.9%.
* **Default cutoffs on enzymes that entered Swiss-Prot after the training years** (174 enzymes that entered in 2025 and 2026, searched against a reference cut to the year before; program, model and cutoffs unchanged; plan fixed before the run). Of the 132 with a neighbour, `EC4` was given for 49, of which 2 were wrong (4.1%, 95% interval 0.5–14.0), and `EC3` for 53, of which 2 were wrong (3.8%, 0.5–13.0); the 60%/40% rule made 6 wrong statements out of 57 at four digits and 4 out of 41 at three. Three of the four wrong statements were on proteins whose function was new to the record. The sample is small: it shows no sign of loss of control, not that the error rate is 4%.
* **Proteins from sequence clusters with no earlier member** (11.9% of the replay; the identity rules make almost no statement there): the model ranked transfer failures with an AUC of 0.79 (identity 0.54), post hoc; the default cutoffs gave `EC4` to 57 of 1,256 proteins with 8.8% wrong and `EC3` to 347 with 9.8% wrong, so errors are more frequent than in the replay as a whole.

## Limits (please read)

* **It cannot name a function that is absent from the reference.** About one in five proteins that have a characterized homolog carry a function that no earlier protein in Swiss-Prot had; for those, every `EC4` statement is wrong, and the rule only tries to keep them out. About half of the errors it removes are of this kind.
* **The probabilities are approximate.** They followed observed frequencies in the replay but ran somewhat high in the upper deciles, and in an external check on present-day proteins with independent experimental annotations, transfer failed more often than predicted (16.7% against 13.0% for proteins with a propagated EC number; 62.8% against 39.6% for proteins with no catalytic annotation). Choose cutoffs on data of your own if you can.
* **Only proteins with a neighbour covering 70% of both sequences are scored**; multidomain proteins and fragments often get `none`.
* **The experimental features** (`p_new_function`, orphan suggestions, `--external`) are weaker than the rule: the score ranks but is not calibrated, the orphan suggestions are right for a minority of flagged proteins, the orphan table is a snapshot of ExplorEnz at the time of the reference, and the imported predictions are only shown next to the nearest enzyme.
* **The reference is Swiss-Prot 2026_04 only**: 18,517 experimentally characterized enzymes with a complete EC number (5,289 EC numbers). It is a static snapshot; a protein whose relatives are not in it gets `none`.
* **All evidence is retrospective** (including the test on 174 enzymes above, which used a reference cut to the year before, not the reference the program ships with), on labels, evidence tags and Pfam assignments as they stand in 2026; the main replay includes EC numbers created after the year replayed (see above). The only prospective test is the protocol for scoring future releases with the frozen model, which has not yet produced a result.
* **The model was fitted on enzymes that someone chose to study**, and the replay says nothing about proteins that a rule rejects.
* **The three-digit model and the cutoffs were added when this tool was built**: the three-digit model has the same features, scaler and training rows as the frozen four-digit model, and was not part of a pre-specified analysis. The risk of a three-digit error is capped at the risk of a four-digit error; in the replay this changed no statement.
* **Other search tools or settings are not validated**, and queries that are themselves in the reference are not what the model was evaluated on.
* Tested on Linux and Windows 11, not on macOS (see above).

## How it works

For each query the program takes all reference enzymes aligned over at least 70% of both sequences, ranked by bit score, and computes 18 features: identity, log bit score and the two coverages of the nearest one; the relative bit-score gap to the second; log number of neighbours; the fraction of the next nine that share an EC number, or its first three digits, with the nearest; whether there is only one neighbour; log number of distinct EC numbers and of distinct sub-subclasses among the ten nearest; the number of EC numbers of the nearest; log query length and absolute log length ratio; log number of earlier protein–EC pairs carrying the nearest enzyme's EC number; log number of distinct EC numbers among all neighbours; identity to the nearest enzyme that shares no EC number with the nearest; and years since the nearest enzyme was established. Two logistic regressions (standardized features, C = 1) turn them into `p_wrong_4digit` and `p_wrong_3digit`. All coefficients are in `data/model.json`.

Three of the features carry most of the signal: identity to the nearest enzyme, the number of distinct EC numbers among all neighbours, and how many earlier proteins carry the nearest enzyme's EC number.

## Files

| file | content |
|---|---|
| `src/ecbackoff/__init__.py` | the program (`ecbackoff` on the command line, `python -m ecbackoff`) |
| `src/ecbackoff/data/model.json` | scaler and coefficients of the two models, search settings, evidence |
| `src/ecbackoff/data/reference.tsv.gz` | reference enzymes: accession, length, year of entry, EC numbers with the year each was established |
| `src/ecbackoff/data/reference.fasta.gz` | their sequences (for the search) |
| `src/ecbackoff/data/orphan.tsv.gz` | experimental: for every distinct EC set of a reference enzyme, its ten highest-ranked orphan activities (shared digits, reaction similarity, year created, name) |
| `src/ecbackoff/data/CHECKSUMS.sha256` | checksums of the data and example files, checked by `selftest` |
| `src/ecbackoff/example/` | 30 example proteins, their MMseqs2 hits and the expected output of the frozen model; used by `selftest` |
| `tests/` | smoke tests of the installed package (version and `selftest`) |

The scripts that built the reference and the orphan table and ran the larger checks below read the project's own files and are not part of this repository; they are to be archived with the paper [DOI to be added].

## Checks

* **`selftest`** (passes with Python 3.9, 3.11, 3.13 and 3.14 on Linux and Windows 11; run for every commit by the repository's automatic checks): the checksums of the data and example files; the reference table and the reference FASTA list the same 18,517 enzymes; the 30 bundled example proteins reproduce the released scores of the frozen model from their stored hits (maximum difference 5e-11); the three-digit risk is never above the four-digit risk; `match` mode gives the numbers of statements of the identity rule; and edge cases built from the example: a query with no neighbour, an empty sequence, a `sp|ACCESSION|NAME` identifier whose hits carry the accession only, a hit to the query itself (dropped, and kept with `--keep-self`), a hit table with percentages (refused), identifiers of any form replaced by synthetic ones for the search and restored afterwards, a hit table pasted twice (counted once), two identifiers with one accession (reported); and, for the experimental features of 1.2: the new-function score reproduces the released values of the 30 examples (maximum difference 5e-11), the orphan table covers every EC set of the reference (5,817 sets, ten ranked suggestions each) and is ordered by shared digits, suggestions appear exactly for proteins at or above the score threshold, and small files in the formats of CLEAN, DeepECtransformer, ProteInfer, DeepEC and JSON lines are read and compared with the nearest enzyme.
* **Version 1.2.1 against 1.2** (Windows 11, Python 3.14, MMseqs2 18-8cc5c): on 400 mutated reference proteins (plain and UniProt-style identifiers) and on the 30 bundled examples, `run` gave tables byte-identical to those of 1.2 in `cutoff` mode; `match` mode gave identical tables; `score` with a hit table written by 1.2 gave the same table; non-ASCII, emoji, `lcl|` and `gi|` identifiers and two identifiers with the same accession each gave the result of the same protein under a plain identifier; an unwritable output path, `--threads 0`, a `--year` outside 1995–2100 and a missing `--keep-hits` folder were refused in under a second, before the search.
* **End-to-end on 400 proteins** from the released score table of the frozen model (run for 1.2 with Python 3.11 and MMseqs2 18.8cc5c): all 400 have the same nearest enzyme as in the released table, and the two probabilities agree with the released values to within the four decimals of the output table (the largest differences, 4.96e-5 and 5.00e-5, are the rounding of the fourth decimal); `p_new_function` agreed to 4.99e-5. In an earlier build the 18 features agreed with the stored features to 1e-9 and did not depend on the number of search threads.
* **`run` and `score`** (with the hit table written by `run --keep-hits`) give identical tables.
* **Command line with the real search** (18 checks on 1.2, all passed; the scripts are not part of this repository): UniProt-style identifiers; a file with a byte-order mark, Windows line ends, lower case and wrapped lines; empty and unmatched queries; the parameter file; refusal of repeated identifiers, a missing or non-FASTA file, a wrong `--mmseqs`, a missing `--tmpdir`, percentages, unknown identifiers and too few columns in a hit table; `--tmpdir` being used and cleaned up; a pipe closed early; the `reference` command.
* **Speed:** a run on 400 proteins took 7.4 seconds with 8 threads on a Linux cluster node and 27 seconds on the Windows PC.
* **Windows 11, Python 3.14, official MMseqs2 Windows build:** `selftest` passed; `score` on the 30 examples reproduced the released scores to the four printed decimals; `run` on the 30 examples gave a hit table identical to the Linux one (after sorting) and an identical result table; the 18 command-line checks passed with the real search; 28 further checks of `score`, the experimental features, input validation and error messages passed with paths containing spaces and Chinese characters; a run on 400 proteins took about 27 seconds with 8 threads.
* **Not checked:** macOS; other MMseqs2 versions; hit tables from other search tools; batches of many thousands of proteins; the readers of `--external` on the output of the real tools.

## Licence and data

The ecbackoff code is licensed under MIT. The reference and example sequences are subsets of [UniProtKB/Swiss-Prot](https://www.uniprot.org/uniprot) from the UniProt Consortium. The reference table combines selected EC annotations with data from the [ENZYME nomenclature database](https://enzyme.expasy.org/) maintained by the Swiss-Prot group at the SIB Swiss Institute of Bioinformatics. The sequences and annotations were selected and reformatted for ecbackoff. UniProtKB/Swiss-Prot and the ENZYME database use [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

The orphan table (`data/orphan.tsv.gz`) selects and ranks EC entries from [ExplorEnz](https://enzyme-database.org/), the IUBMB Enzyme Nomenclature database. ExplorEnz entries are copyright 2001-2026 IUBMB. The [database authors permit use in other applications with acknowledgement of the source](https://doi.org/10.1093/nar/gkn582). The reaction similarity values are calculated from [Rhea](https://www.rhea-db.org/), which is available under [CC BY 4.0](https://www.rhea-db.org/help/license-disclaimer). See `LICENSE`.

## Citation

A citation for the manuscript will be added here when it is public. Please also cite the resources used:

* The UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2025. *Nucleic Acids Res.* 53, D609–D617 (2025).
* Bairoch, A. The ENZYME database in 2000. *Nucleic Acids Res.* 28, 304-305 (2000). [DOI](https://doi.org/10.1093/nar/28.1.304).
* Rhea, the reaction knowledgebase in 2022. *Nucleic Acids Res.* 50, D693-D700 (2022). [DOI](https://doi.org/10.1093/nar/gkab1016).
* ExplorEnz: a MySQL database of the IUBMB enzyme nomenclature. *BMC Biochem.* 8, 14 (2007). [DOI](https://doi.org/10.1186/1471-2091-8-14).
* Steinegger, M. & Söding, J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. *Nat. Biotechnol.* 35, 1026–1028 (2017).
