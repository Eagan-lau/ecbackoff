#!/usr/bin/env python3
"""ecbackoff: inherit an enzyme's EC number only where the neighbourhood supports it; otherwise back off to three digits, or say nothing.

For every query protein the program reads the sequence-search hits against a bundled reference of experimentally
characterized enzymes (UniProtKB/Swiss-Prot 2026_04), computes 18 features of the neighbourhood of the nearest enzyme,
and estimates the probability that inheriting that enzyme's four-digit EC number, or its first three digits, would be
wrong. It then says

    EC4   the four-digit EC number(s) of the nearest enzyme,
    EC3   only the first three digits (the sub-subclass),
    none  nothing.

Experimental additions (version 1.2; see README, "Experimental features"): a ranking score for a function that is new to the record, suggestions of orphan
enzyme activities, and the import of the output of other predictors (CLEAN, DeepECtransformer, ProteInfer, DeepEC) that you have run yourself.

Commands
    run        FASTA in, table out (needs MMseqs2 on the PATH)
    score      MMseqs2 hit table + FASTA in, table out (no MMseqs2 needed)
    reference  write the reference sequences (FASTA) for running the search yourself
    selftest   checks the installation: file checksums, bundled example, edge cases
    info       model, reference and evidence summary

Standard library only (Python 3.9 or newer). See README.md.
"""
import argparse
import bisect
import contextlib
import datetime
import gzip
import hashlib
import io
import json
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile

VERSION = '1.2.1'
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data')
EXAMPLE = os.path.join(HERE, 'example')
SEARCH_ARGS = ['-s', '7.5', '-e', '1e-3', '--max-seqs', '20000', '--alignment-mode', '3']
FORMAT = 'query,target,fident,qcov,tcov,evalue,bits'
MMSEQS_TESTED = '18.8cc5c'
MIN_COV, TOPK = 0.7, 10
COLUMNS = ['query', 'nearest_enzyme', 'identity', 'nearest_ec', 'n_neighbours', 'p_wrong_4digit', 'p_wrong_3digit', 'statement', 'statement_ec', 'identity_rule', 'note']
EXTRA_COLUMNS = ['p_new_function', 'orphan_candidates', 'orphan_names']      # experimental, after the validated columns
ORPHAN_MIN_P, ORPHAN_TOP = 0.5, 5
EC_RE = re.compile(r'(\d+\.(?:\d+|-)\.(?:\d+|-)\.(?:\d+|-))')
EXTERNAL_TOOLS = ('clean', 'deepectransformer', 'proteinfer', 'deepec')


def die(msg):
    print(f'error: {msg}', file=sys.stderr)
    sys.exit(2)


def warn(msg):
    print(f'[ecbackoff] warning: {msg}', file=sys.stderr)


def lines(path):
    """Lines of a text or .gz file without line ends. Tolerates a UTF-8 byte-order mark, Windows line ends and stray bytes."""
    try:
        fh = gzip.open(path, 'rt', encoding='utf-8-sig', errors='replace') if path.endswith('.gz') else open(path, 'r', encoding='utf-8-sig', errors='replace')
    except OSError as e:
        die(f'cannot read {path}: {e.strerror or e}')
    try:
        for line in fh:
            yield line.rstrip('\r\n')
    except (OSError, EOFError) as e:
        die(f'cannot read {path}: {e}')
    finally:
        fh.close()


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def read_fasta(path):
    """Ordered {identifier: length}; the identifier is the first word of the header, as MMseqs2 uses it."""
    lengths, cur, odd = {}, None, 0
    for n, line in enumerate(lines(path), 1):
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            words = line[1:].split()
            if not words:
                die(f'{path}, line {n}: a header has no identifier')
            cur = words[0]
            if cur in lengths:
                die(f'identifier {cur!r} occurs twice in {path}')
            lengths[cur] = 0
        else:
            if cur is None:
                die(f'{path} does not look like a FASTA file: line {n} is sequence before any ">" header')
            seq = line.replace('*', '').replace('-', '')
            lengths[cur] += len(seq)
            if not seq.isalpha():
                odd += 1
    if not lengths:
        die(f'{path} contains no sequences')
    if odd:
        warn(f'{odd} sequence line(s) in {path} contain characters other than letters, "-" and "*"; they are counted in the length')
    return lengths


def write_clean_fasta(src, dst):
    """Copy of the query FASTA for the search: sequences under synthetic ASCII identifiers (q0000001, q0000002, ...), letters only, empty sequences left out.
    The search therefore never sees the identifiers of the user (non-ASCII characters, 'sp|ACC|NAME', 'gi|123|gb|AAA.1|' and two identifiers with the same
    accession all behave alike). Returns the identifiers left out and {synthetic identifier: original identifier}."""
    empty, mapping, cur, buf = [], {}, None, []

    def flush(fh):
        if cur is None:
            return
        seq = ''.join(buf)
        if seq:
            sid = 'q%07d' % (len(mapping) + 1)
            mapping[sid] = cur
            fh.write(f'>{sid}\n{seq}\n')
        else:
            empty.append(cur)

    with open(dst, 'w', encoding='utf-8', newline='\n') as fh:
        for line in lines(src):
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                flush(fh)
                cur, buf = line[1:].split()[0], []
            else:
                buf.append(line.replace('*', '').replace('-', '').replace(' ', '').upper())
        flush(fh)
    return empty, mapping


def translate_hits(src, dst, mapping):
    """Hit table of the search with the synthetic query identifiers replaced by the identifiers of the query FASTA."""
    with open(src, 'r', encoding='utf-8', errors='replace', newline='') as a, open(dst, 'w', encoding='utf-8', newline='\n') as b:
        for line in a:
            line = line.rstrip('\r\n')
            if not line:
                continue
            q, sep, rest = line.partition('\t')
            b.write(mapping.get(q, q) + sep + rest + '\n')


def check_writable(path, what):
    """Stop before a long search if the result could not be written afterwards."""
    if not path or path == '-':
        return
    full = os.path.abspath(path)
    if os.path.isdir(full):
        die(f'{what} {path!r} is a directory')
    folder = os.path.dirname(full)
    if not os.path.isdir(folder):
        die(f'{what}: the directory {folder!r} does not exist')
    if not os.access(folder, os.W_OK):
        die(f'{what}: the directory {folder!r} is not writable')


def warn_accession_clashes(qids):
    """Two query identifiers that share a UniProt accession cannot be told apart in a hit table that prints only the accession."""
    by = {}
    for q in qids:
        by.setdefault(accession_of(q), []).append(q)
    clash = {a: qs for a, qs in by.items() if len(qs) > 1}
    if clash:
        a, qs = next(iter(sorted(clash.items())))
        warn(f'{len(clash)} accession(s) are shared by several query identifiers (e.g. {a}: {qs[:3]}); if the hit table prints only the accession, their hits cannot be told apart. '
             'Give every protein a distinct identifier, or use the "run" command, which does not have this problem')


def pref3(ecs):
    return {'.'.join(e.split('.')[:3]) for e in ecs}


class Reference:
    """Experimentally characterized enzymes of the frozen reference: length, year of entry, EC numbers with the year each was established."""

    def __init__(self, data_dir=DATA):
        self.length, self.year, self.labels, self.by_ec = {}, {}, {}, {}
        self._p3 = {}
        path = os.path.join(data_dir, 'reference.tsv.gz')
        if not os.path.exists(path):
            die(f'reference data not found: {path} (is the installation complete? "ecbackoff selftest" checks it)')
        for line in lines(path):
            acc, length, year, ecs = line.split('\t')
            self.length[acc], self.year[acc] = int(length), int(year)
            labs = set()
            for item in ecs.split(';'):
                e, y = item.rsplit(':', 1)
                labs.add(e)
                self.by_ec.setdefault(e, []).append(int(y))
            self.labels[acc] = labs
        for e in self.by_ec:
            self.by_ec[e].sort()

    def p3(self, acc):
        if acc not in self._p3:
            self._p3[acc] = pref3(self.labels[acc])
        return self._p3[acc]


def load_model(data_dir=DATA):
    path = os.path.join(data_dir, 'model.json')
    if not os.path.exists(path):
        die(f'model not found: {path}')
    with open(path) as fh:
        return json.load(fh)


def load_orphans(data_dir=DATA):
    """Orphan enzyme activities (EC numbers with no characterized sequence among the reference enzymes), ranked for every distinct EC set of a reference enzyme:
    {ec_set: [(rank, orphan_ec, shared_ec_digits, reaction_similarity or None, year_created, name)]}."""
    path = os.path.join(data_dir, 'orphan.tsv.gz')
    out = {}
    if not os.path.exists(path):
        return out
    for i, line in enumerate(lines(path)):
        if i == 0 or not line.strip():
            continue
        key, rank, ec, shared, sim, created, name = line.split('\t')
        out.setdefault(key, []).append((int(rank), ec, int(shared), None if sim in ('', 'NA') else float(sim), int(created), name))
    for key in out:
        out[key].sort()
    return out


def read_external(path):
    """Predictions of another EC predictor, read from the file that tool wrote: {identifier: [(EC number, confidence or None), ...]}.
    Recognised: JSON lines (query_id with labels/confidence), CLEAN's maximum-separation CSV (id,EC:x.x.x.x/distance,...), and tab or comma separated tables
    with the identifier first (DeepECtransformer, DeepEC, ProteInfer). A protein the tool listed without an EC number is kept with an empty list."""
    recs = {}

    def key(q):
        return str(q).strip().split('_SEPARATED_SEQUENCE_')[0]

    def add(q, ec, conf):
        lst = recs.setdefault(key(q), [])
        if ec and ec not in [e for e, _ in lst]:
            lst.append((ec, conf))

    rows = [ln for ln in lines(path) if ln.strip()]
    if not rows:
        die(f'{path} is empty')
    if rows[0].lstrip().startswith('{'):
        for n, ln in enumerate(rows, 1):
            try:
                r = json.loads(ln)
            except ValueError:
                die(f'{path}, line {n}: not valid JSON')
            q = r.get('query_id') or r.get('sequence_name') or r.get('id') if isinstance(r, dict) else None
            if q is None:
                die(f'{path}, line {n}: no query_id')
            recs.setdefault(key(q), [])
            labs = r.get('labels') or r.get('native10') or r.get('prediction') or []
            labs = labs.split(';') if isinstance(labs, str) else labs
            conf = r.get('confidence')
            for i, lab in enumerate(labs):
                m = EC_RE.search(str(lab))
                if m:
                    add(q, m.group(1), float(conf[i]) if isinstance(conf, list) and i < len(conf) else None)
        return recs
    sep = '\t' if any('\t' in ln for ln in rows[:20]) else ','
    header = None
    if not EC_RE.search(rows[0]):
        header = [h.strip().lower() for h in rows[0].split(sep)]
        rows = rows[1:]
    conf_col = next((i for i, h in enumerate(header or []) if 'conf' in h or 'score' in h), None)
    pat = re.compile(r'(\d+\.(?:\d+|-)\.(?:\d+|-)\.(?:\d+|-))(?:/(\d+(?:\.\d+)?))?')
    for ln in rows:
        p = ln.split(sep)
        if not p[0].strip():
            continue
        recs.setdefault(key(p[0]), [])
        row_conf = None
        if conf_col is not None and conf_col < len(p):
            try:
                row_conf = float(p[conf_col])
            except ValueError:
                row_conf = None
        for cell in p[1:]:
            for m in pat.finditer(cell):
                add(p[0], m.group(1), float(m.group(2)) if m.group(2) else row_conf)
    return recs


def compare_external(ecs, nearest):
    """How an external prediction relates to the nearest reference enzyme's EC numbers."""
    if not nearest:
        return 'no_neighbour'
    if not ecs:
        return 'no_tool_answer'
    if {e for e in ecs if '-' not in e} & nearest:
        return 'same_4digits'
    t3 = {'.'.join(e.split('.')[:3]) for e in ecs if '-' not in e.split('.')[:3]}
    return 'same_3digits' if t3 & pref3(nearest) else 'differs'


def parse_external_args(items):
    out, seen = [], set()
    for item in items or []:
        name, eq, path = item.partition('=')
        name = name.strip().lower()
        if not eq or not re.fullmatch(r'[a-z0-9_]+', name) or not path:
            die(f'--external expects NAME=FILE (NAME of letters, digits and underscores, e.g. clean=maxsep.csv), got {item!r}')
        if name in seen:
            die(f'--external {name} is given twice')
        if not os.path.exists(path):
            die(f'cannot read {path}: no such file')
        if name not in EXTERNAL_TOOLS:
            warn(f'--external {name}: not one of {", ".join(EXTERNAL_TOOLS)}; the file is read in the generic formats')
        seen.add(name)
        out.append((name, path))
    return out


def attach_external(rows, qlens, ref, externals):
    """Adds to every row what each external predictor said, and whether it agrees with the nearest reference enzyme."""
    alias = make_alias(list(qlens))
    for name, path in externals:
        recs = read_external(path)
        by_query, unmatched = {}, 0
        for t, lst in recs.items():
            q = alias.get(t)
            if q is None:
                unmatched += 1
            else:
                by_query[q] = lst
        if not by_query:
            die(f'none of the {len(recs)} identifiers in {path} is in the query FASTA')
        if unmatched:
            warn(f'{unmatched} identifier(s) in {path} are not in the query FASTA and were ignored')
        for r in rows:
            lst = by_query.get(r['query'])
            ecs = [e for e, _ in lst] if lst else []
            nearest = ref.labels[r['nearest_enzyme']] if r['nearest_enzyme'] else set()
            r.setdefault('ext', {})[name] = {'ec': ';'.join(ecs), 'conf': ';'.join('' if c is None else f'{c:g}' for _, c in (lst or [])),
                                              'vs': 'no_tool_answer' if (lst is None and nearest) else compare_external(ecs, nearest)}


def attach_orphans(rows, ref, orphans, min_p, top):
    """Orphan-activity suggestions for proteins whose new-function score is at least min_p (the rule of the paper: ranked by shared EC digits with the nearest
    enzyme, then by reaction similarity)."""
    for r in rows:
        r['orphan'] = r['orphan_names'] = ''
        if top <= 0 or r['pn'] is None or r['pn'] < min_p or not r['nearest_enzyme']:
            continue
        cands = orphans.get(';'.join(sorted(ref.labels[r['nearest_enzyme']])), [])[:top]
        r['orphan'] = ' | '.join('%s [shared %d, similarity %s]' % (ec, sh, 'NA' if sim is None else '%.2f' % sim) for _, ec, sh, sim, _, _ in cands)
        r['orphan_names'] = ' | '.join(('%s (EC created %d)' % (nm, cr)) if nm else ('(EC created %d)' % cr) for _, _, _, _, cr, nm in cands)


def accession_of(identifier):
    """UniProt-style 'sp|P12345|NAME' -> 'P12345'; anything else unchanged."""
    parts = identifier.split('|')
    return parts[1] if len(parts) >= 3 else identifier


def make_alias(qids):
    """Identifier as it may appear in a hit table -> identifier in the FASTA. MMseqs2 prints the accession of a UniProt-style header
    ('sp|P12345|NAME' -> 'P12345'); an accession shared by two queries stays unresolved."""
    alias = {q: q for q in qids}
    by_acc = {}
    for q in qids:
        a = accession_of(q)
        if a != q:
            by_acc.setdefault(a, []).append(q)
    for a, qs in by_acc.items():
        if a not in alias and len(qs) == 1:
            alias[a] = qs[0]
    return alias


def read_hits(path, ref, qids, keep_self):
    """Hits that qualify as neighbours (both sequences covered to at least 70%), best bit score first. Stops with a clear message on a malformed table."""
    alias = make_alias(qids)
    per, unknown, n_hits, n_ref, n_badbits, n_foreign, n_dup, seen = {}, {}, 0, 0, 0, 0, 0, set()
    for n, line in enumerate(lines(path), 1):
        if not line.strip():
            continue
        p = line.split('\t')
        if len(p) < 7:
            die(f'{path}, line {n}: expected 7 tab-separated columns ({FORMAT}), found {len(p)}')
        try:
            fident, qcov, tcov, bits = float(p[2]), float(p[3]), float(p[4]), float(p[6])
        except ValueError:
            if n == 1:
                continue   # a header line
            die(f'{path}, line {n}: identity, coverages and bit score must be numbers ({FORMAT})')
        if not all(math.isfinite(x) for x in (fident, qcov, tcov, bits)):
            die(f'{path}, line {n}: a value is not finite')
        if min(fident, qcov, tcov) < 0 or max(fident, qcov, tcov) > 1.0001:
            die(f'{path}, line {n}: identity and coverages must be fractions between 0 and 1; was the table made with --format-output {FORMAT}?')
        n_hits += 1
        q0, t = p[0], p[1]
        q = alias.get(q0)
        if q is None:
            unknown[q0] = unknown.get(q0, 0) + 1
            continue
        if t not in ref.labels:
            n_foreign += 1
            continue
        key = (q,) + tuple(p[1:7])
        if key in seen:      # a line repeated exactly (a hit table pasted twice); a genuine search prints each alignment once
            n_dup += 1
            continue
        seen.add(key)
        n_ref += 1
        if not keep_self and (t == q or t == accession_of(q)):
            continue
        if bits <= 0:
            n_badbits += 1
            continue
        if qcov >= MIN_COV and tcov >= MIN_COV:
            per.setdefault(q, []).append((t, fident, qcov, tcov, bits))
    if unknown:
        ex = sorted(unknown)[:3]
        die(f'the hit table has {len(unknown)} query identifier(s) that are not in the FASTA, e.g. {ex}; the identifier is the first word of the FASTA header '
            '(an accession shared by several UniProt-style identifiers cannot be resolved)')
    if n_hits and not n_ref:
        die('none of the hits points to a reference enzyme; the search must be run against the bundled reference (ecbackoff reference -o reference.fasta)')
    if n_badbits:
        warn(f'{n_badbits} hit(s) with a bit score of 0 or less were ignored')
    if n_dup:
        warn(f'{n_dup} repeated line(s) of the hit table were ignored (was the table concatenated twice?)')
    if n_foreign:
        warn(f'{n_foreign} hit(s) point to targets that are not reference enzymes and were ignored (was the search run against the reference from "ecbackoff reference"?)')
    for q in per:
        per[q].sort(key=lambda x: -x[4])   # stable: ties keep file order
    return per


def features(g, qlen, ref, py):
    t = [x[0] for x in g]
    f = [x[1] for x in g]
    b = [x[4] for x in g]
    nn = t[0]
    ln = ref.labels[nn]
    top = t[:TOPK]
    others = top[1:]
    nn3 = ref.p3(nn)
    return {
        'ident': f[0], 'log_bits': math.log(b[0]), 'qcov': g[0][2], 'tcov': g[0][3],
        'margin': (b[0] - b[1]) / b[0] if len(b) > 1 else 1.0,
        'log_nnb': math.log1p(len(t)),
        'k_agree': (sum(1 for x in others if ref.labels[x] & ln) / len(others)) if others else 1.0,
        'k_single': 1 if not others else 0,
        'log_k_distinct': math.log1p(len(set().union(*[ref.labels[x] for x in top]))),
        'nn_nlabels': len(ln),
        'log_qlen': math.log(qlen), 'abs_log_lenratio': abs(math.log(qlen / ref.length[nn])),
        'log_label_freq': math.log1p(max(bisect.bisect_right(ref.by_ec[e], py - 1) for e in ln)),
        'l3_agree_k': (sum(1 for x in others if ref.p3(x) & nn3) / len(others)) if others else 1.0,
        'log_k_l3_distinct': math.log1p(len(set().union(*[ref.p3(x) for x in top]))),
        'log_fam_n_ec': math.log1p(len(set().union(*[ref.labels[x] for x in t]))),
        'ident_other': next((f[i] for i, x in enumerate(t) if not (ref.labels[x] & ln)), 0.0),
        'nn_age': float(py - ref.year[nn]),
    }


def p_wrong(model, key, feats):
    s, m = model['scaler'], model[key]
    z = m['intercept'] + sum(((feats[n] - mu) / sc) * c for n, mu, sc, c in zip(model['features'], s['mean'], s['scale'], m['coef']))
    if not math.isfinite(z):
        return None
    z = max(-60.0, min(60.0, z))
    return 1.0 / (1.0 + math.exp(-z))


def fmt(x, nd=4):
    return '' if x is None else f'{x:.{nd}f}'


def analyse(qlens, per, ref, model, py):
    """One record per query, in the order of the FASTA. The three-digit risk is capped at the four-digit risk: a three-digit error is also a four-digit error."""
    rows = []
    for q, qlen in qlens.items():
        r = {'query': q, 'nearest_enzyme': '', 'identity': None, 'nearest_ec': '', 'n_neighbours': 0, 'p4': None, 'p3': None, 'statement': 'none',
             'statement_ec': '', 'identity_rule': 'none', 'note': '', 'pn': None, 'orphan': '', 'orphan_names': ''}
        g = per.get(q)
        if qlen <= 0:
            r['note'] = 'empty_sequence'
        elif not g:
            r['note'] = 'no_neighbour_covering_70pct_of_both_sequences'
        else:
            ft = features(g, qlen, ref, py)
            p4, p3 = p_wrong(model, 'wrong_4digit', ft), p_wrong(model, 'wrong_3digit', ft)
            r.update({'nearest_enzyme': g[0][0], 'identity': ft['ident'], 'nearest_ec': ';'.join(sorted(ref.labels[g[0][0]])), 'n_neighbours': len(g)})
            if p4 is None or p3 is None:
                r['note'] = 'score_not_finite'
            else:
                r.update({'p4': p4, 'p3': min(p3, p4)})
                if 'new_function' in model:
                    r['pn'] = p_wrong(model, 'new_function', ft)
        rows.append(r)
    return rows


def decide_match(rows, hi, lo):
    """As many four-digit and three-digit statements as the identity rule gives on this batch; the model chooses which proteins get them."""
    ok = [r for r in rows if r['p4'] is not None]
    n4 = sum(1 for r in ok if r['identity'] >= hi)
    n3 = sum(1 for r in ok if lo <= r['identity'] < hi)
    first = sorted(ok, key=lambda r: (r['p4'], r['query']))[:n4]
    taken = {r['query'] for r in first}
    second = sorted([r for r in ok if r['query'] not in taken], key=lambda r: (r['p3'], r['query']))[:n3]
    for r in first:
        r['statement'] = 'EC4'
    for r in second:
        r['statement'] = 'EC3'


def decide_cutoff(rows, c4, c3):
    for r in rows:
        if r['p4'] is None:
            continue
        if r['p4'] <= c4:
            r['statement'] = 'EC4'
        elif r['p3'] <= c3:
            r['statement'] = 'EC3'


def finish(rows, ref, hi, lo):
    for r in rows:
        if r['identity'] is not None:
            r['identity_rule'] = 'EC4' if r['identity'] >= hi else ('EC3' if r['identity'] >= lo else 'none')
            labs = ref.labels[r['nearest_enzyme']]
            if r['statement'] == 'EC4':
                r['statement_ec'] = ';'.join(sorted(labs))
            elif r['statement'] == 'EC3':
                r['statement_ec'] = ';'.join(sorted(pref3(labs)))
            if r['statement'] == 'none' and not r['note']:
                r['note'] = 'neighbourhood does not support a statement'


def classify(qlens, per, ref, model, py, mode, hi, lo, c4, c3):
    rows = analyse(qlens, per, ref, model, py)
    if mode == 'match':
        decide_match(rows, hi, lo)
    else:
        decide_cutoff(rows, c4, c3)
    finish(rows, ref, hi, lo)
    return rows


def write_table(rows, out, ext_names=()):
    if out and out != '-':
        try:
            fh = open(out, 'w', encoding='utf-8', newline='\n')
        except OSError as e:
            die(f'cannot write {out}: {e.strerror or e}')
    else:
        fh = sys.stdout
    ext_cols = [f'{n}_{s}' for n in ext_names for s in ('ec', 'conf', 'vs_nearest')]
    try:
        fh.write('\t'.join(COLUMNS + EXTRA_COLUMNS + ext_cols) + '\n')
        for r in rows:
            cells = [r['query'], r['nearest_enzyme'], fmt(r['identity'], 3), r['nearest_ec'], str(r['n_neighbours']), fmt(r['p4']), fmt(r['p3']),
                     r['statement'], r['statement_ec'], r['identity_rule'], r['note'], fmt(r['pn']), r['orphan'], r['orphan_names']]
            for n in ext_names:
                e = r.get('ext', {}).get(n, {'ec': '', 'conf': '', 'vs': ''})
                cells += [e['ec'], e['conf'], e['vs']]
            fh.write('\t'.join(cells) + '\n')
    finally:
        if fh is not sys.stdout:
            fh.close()


def counts(rows):
    c = {'queries': len(rows), 'with_neighbour': sum(1 for r in rows if r['p4'] is not None),
         'statements': {k: sum(1 for r in rows if r['statement'] == k) for k in ('EC4', 'EC3', 'none')},
         'identity_rule': {k: sum(1 for r in rows if r['identity_rule'] == k) for k in ('EC4', 'EC3', 'none')},
         'with_orphan_suggestions': sum(1 for r in rows if r['orphan'])}
    names = sorted({n for r in rows for n in r.get('ext', {})})
    if names:
        c['external'] = {n: {k: sum(1 for r in rows if r['ext'][n]['vs'] == k) for k in ('same_4digits', 'same_3digits', 'differs', 'no_tool_answer', 'no_neighbour')} for n in names}
    return c


def summary(c, mode):
    print(f'[ecbackoff] {c["queries"]} queries, {c["with_neighbour"]} with a qualifying neighbour ({mode} mode)', file=sys.stderr)
    s, i = c['statements'], c['identity_rule']
    print(f'[ecbackoff] statements: four digits {s["EC4"]}, three digits {s["EC3"]}, none {s["none"]}', file=sys.stderr)
    print(f'[ecbackoff] identity rule (for comparison): four digits {i["EC4"]}, three digits {i["EC3"]}, none {i["none"]}', file=sys.stderr)
    if c['with_orphan_suggestions']:
        print(f'[ecbackoff] experimental: orphan-activity suggestions for {c["with_orphan_suggestions"]} proteins with a new-function score at or above the threshold', file=sys.stderr)
    for n, v in c.get('external', {}).items():
        print(f'[ecbackoff] {n} against the nearest enzyme: same four digits {v["same_4digits"]}, same three digits {v["same_3digits"]}, different {v["differs"]}, '
              f'no answer from the tool {v["no_tool_answer"]}, no neighbour {v["no_neighbour"]}', file=sys.stderr)
    if mode == 'match' and c['with_neighbour'] < 200:
        warn('fewer than 200 queries with a neighbour; the matched numbers are noisy, consider --mode cutoff')


def check_options(args):
    if not (0 <= args.identity_lo <= args.identity_hi <= 1):
        die('--identity-lo and --identity-hi must satisfy 0 <= lo <= hi <= 1')
    for name in ('cutoff4', 'cutoff3'):
        if not (0 <= getattr(args, name) <= 1):
            die(f'--{name} must be a probability between 0 and 1')
    if not (0 <= args.orphan_min_p <= 1):
        die('--orphan-min-p must be a probability between 0 and 1')
    if args.orphan_top < 0:
        die('--orphan-top must be 0 or more')
    if args.year and not (1995 <= args.year <= 2100):
        die('--year must be a year between 1995 and 2100 (0 means the frozen default)')
    if not getattr(args, 'externals', None):
        args.externals = parse_external_args(args.external)


def write_params(args, c, extra):
    """Record next to the output table how it was made."""
    if args.out in ('', '-') or args.no_params:
        return
    p = {'ecbackoff_version': VERSION, 'date_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds'), 'command': shlex.join(sys.argv),
         'mode': args.mode, 'cutoff4': args.cutoff4, 'cutoff3': args.cutoff3, 'identity_hi': args.identity_hi, 'identity_lo': args.identity_lo,
         'prediction_year': args.py, 'prediction_year_is_frozen': args.py == args.model['prediction_year'],
         'reference': args.model['reference_release'], 'query_fasta': os.path.abspath(args.query), 'query_fasta_sha256': sha256_of(args.query), 'counts': c,
         'experimental': {'orphan_suggestions': not args.no_orphan, 'orphan_min_p': args.orphan_min_p, 'orphan_top': args.orphan_top,
                          'external': {n: {'file': os.path.abspath(f), 'sha256': sha256_of(f)} for n, f in args.externals}}}
    p.update(extra)
    try:
        with open(args.out + '.params.json', 'w') as fh:
            json.dump(p, fh, indent=1)
    except OSError as e:
        warn(f'could not write the parameter file: {e.strerror or e}')


def process(args, hits_path, qlens=None):
    check_options(args)
    ref, model = Reference(), load_model()
    scoring_only = qlens is None
    qlens = qlens if qlens is not None else read_fasta(args.query)
    if scoring_only:
        warn_accession_clashes(list(qlens))
    per = read_hits(hits_path, ref, list(qlens), args.keep_self)
    py = args.year if args.year else model['prediction_year']
    if py != model['prediction_year']:
        warn(f'--year {py} differs from the frozen prediction year {model["prediction_year"]}; it changes two features (years since the nearest enzyme was established, '
             'and the number of earlier proteins carrying its EC number) and is outside what was validated')
    args.py, args.model = py, model
    rows = classify(qlens, per, ref, model, py, args.mode, args.identity_hi, args.identity_lo, args.cutoff4, args.cutoff3)
    if not args.no_orphan:
        attach_orphans(rows, ref, load_orphans(), args.orphan_min_p, args.orphan_top)
    if args.externals:
        attach_external(rows, qlens, ref, args.externals)
    return rows


def finish_command(args, rows, extra):
    c = counts(rows)
    summary(c, args.mode)
    write_table(rows, args.out, [n for n, _ in args.externals])
    write_params(args, c, extra)


def cmd_score(args):
    check_writable(args.out, 'output table')
    rows = process(args, args.hits)
    finish_command(args, rows, {'hit_table': os.path.abspath(args.hits), 'hit_table_sha256': sha256_of(args.hits)})


def find_mmseqs(arg):
    mm = shutil.which(arg) if arg else shutil.which('mmseqs')
    if not mm:
        if arg:
            die(f'--mmseqs {arg!r} is not an executable file')
        die('MMseqs2 was not found on the PATH; install it (https://github.com/soedinglab/MMseqs2) or give --mmseqs, '
            'or run the search yourself ("ecbackoff reference" writes the reference) and use the "score" command')
    try:
        v = subprocess.run([mm, 'version'], capture_output=True, text=True, timeout=120)
    except (OSError, subprocess.SubprocessError) as e:
        die(f'{mm} could not be run: {e}')
    if v.returncode != 0:
        die(f'{mm} version failed; is it an MMseqs2 executable?')
    version = v.stdout.strip().splitlines()[0] if v.stdout.strip() else 'unknown'
    # the conda build prints "18.8cc5c", the official Windows build prints the full commit hash, which starts with 8cc5c
    if version != MMSEQS_TESTED and not version.startswith(MMSEQS_TESTED.split('.')[-1]):
        warn(f'MMseqs2 {version} found; the model was built and checked with {MMSEQS_TESTED}')
    return mm, version


def cmd_run(args):
    check_options(args)
    qlens = read_fasta(args.query)   # validate before searching
    if args.threads < 1:
        die('--threads must be 1 or more')
    check_writable(args.out, 'output table')
    check_writable(args.keep_hits, '--keep-hits')
    mm, mm_version = find_mmseqs(args.mmseqs)
    if args.tmpdir and not os.path.isdir(args.tmpdir):
        die(f'--tmpdir {args.tmpdir!r} is not a directory')
    with tempfile.TemporaryDirectory(prefix='ecbackoff_', dir=args.tmpdir or None) as tmp:
        refq, clean = os.path.join(tmp, 'reference.fasta'), os.path.join(tmp, 'queries.fasta')
        with gzip.open(os.path.join(DATA, 'reference.fasta.gz'), 'rb') as src, open(refq, 'wb') as dst:
            shutil.copyfileobj(src, dst)
        empty, mapping = write_clean_fasta(args.query, clean)
        if empty:
            warn(f'{len(empty)} query sequence(s) are empty and are not searched (e.g. {empty[:3]})')
        if len(empty) == len(qlens):
            hits = os.path.join(tmp, 'hits.m8')
            open(hits, 'w').close()
        else:
            raw = os.path.join(tmp, 'hits_raw.m8')
            cmd = [mm, 'easy-search', clean, refq, raw, os.path.join(tmp, 'tmp')] + SEARCH_ARGS + ['--threads', str(args.threads), '-v', '1', '--format-output', FORMAT]
            print('[ecbackoff] ' + shlex.join(cmd), file=sys.stderr)
            r = subprocess.run(cmd)
            if r.returncode != 0:
                hint = ('; on Windows the official MMseqs2 build needs its BusyBox helpers installed once (run mmseqs.bat once, or "busybox --install -s bin" '
                        'in its bin folder)') if os.name == 'nt' else ''
                die(f'MMseqs2 failed with exit code {r.returncode}; the command is printed above{hint}')
            if not os.path.exists(raw):
                die('MMseqs2 finished without writing a hit table')
            hits = args.keep_hits or os.path.join(tmp, 'hits.m8')
            translate_hits(raw, hits, mapping)   # the search saw synthetic identifiers; the table kept or scored carries the identifiers of the query FASTA
        rows = process(args, hits, qlens)
        finish_command(args, rows, {'mmseqs': mm_version, 'search': ' '.join(SEARCH_ARGS), 'threads': args.threads})


def cmd_reference(args):
    src_path = os.path.join(DATA, 'reference.fasta.gz')
    if not os.path.exists(src_path):
        die(f'reference data not found: {src_path}')
    with gzip.open(src_path, 'rb') as src:
        if args.out == '-':
            shutil.copyfileobj(src, sys.stdout.buffer)
        else:
            try:
                with open(args.out, 'wb') as dst:
                    shutil.copyfileobj(src, dst)
            except OSError as e:
                die(f'cannot write {args.out}: {e.strerror or e}')
    print('[ecbackoff] search with:', file=sys.stderr)
    print(f'  mmseqs easy-search QUERIES.fasta {"REFERENCE.fasta" if args.out == "-" else args.out} hits.m8 tmp {" ".join(SEARCH_ARGS)} --format-output {FORMAT}', file=sys.stderr)
    print('then: ecbackoff score -q QUERIES.fasta --hits hits.m8 -o result.tsv', file=sys.stderr)


# ---------------------------------------------------------------- selftest

def read_checksums():
    path = os.path.join(HERE, 'data', 'CHECKSUMS.sha256')
    if not os.path.exists(path):
        return None
    out = {}
    for line in lines(path):
        if line.strip():
            digest, name = line.split(None, 1)
            out[name.strip()] = digest
    return out


def cmd_selftest(args):
    results = []

    def check(name, ok, detail=''):
        results.append(bool(ok))
        print(f'[selftest] {"ok  " if ok else "FAIL"} {name}' + (f' ({detail})' if detail else ''))

    # 1. files
    sums = read_checksums()
    if sums is None:
        check('data/CHECKSUMS.sha256 present', False, 'the checksum file is missing')
    else:
        bad = [n for n, d in sums.items() if not os.path.exists(os.path.join(HERE, n)) or sha256_of(os.path.join(HERE, n)) != d]
        check('file checksums', not bad and len(sums) >= 6, f'{len(sums)} files' if not bad else f'changed or missing: {bad}')
    ref, model = Reference(), load_model()
    n_fa, ids = 0, set()
    for line in lines(os.path.join(DATA, 'reference.fasta.gz')):
        if line.startswith('>'):
            n_fa += 1
            ids.add(line[1:].split()[0])
    check('reference sequences match the reference table', ids == set(ref.labels) and n_fa == len(ids), f'{n_fa} sequences, {len(ref.labels)} table rows')

    # 2. the bundled example against the released scores
    q, h, e = (os.path.join(EXAMPLE, n) for n in ('queries.fasta', 'hits.m8', 'expected.tsv'))
    qlens = read_fasta(q)
    per = read_hits(h, ref, list(qlens), False)
    rows = classify(qlens, per, ref, model, model['prediction_year'], 'cutoff', 0.6, 0.4, 0.10, 0.10)
    got = {r['query']: r for r in rows}
    exp = {}
    for i, line in enumerate(lines(e)):
        if i == 0:
            head = line.split('\t')
        elif line.strip():
            d = dict(zip(head, line.split('\t')))
            exp[d['query']] = d
    worst4 = worst3 = 0.0
    bad = 0
    for k, d in exp.items():
        r = got.get(k)
        if r is None or r['p4'] is None or r['nearest_enzyme'] != d['nearest_enzyme']:
            bad += 1
            continue
        worst4 = max(worst4, abs(r['p4'] - float(d['p_wrong_4digit'])))
        worst3 = max(worst3, abs(r['p3'] - float(d['p_wrong_3digit'])))
    check('example reproduces the released scores', bad == 0 and worst4 < 1e-6 and worst3 < 1e-6 and exp,
          f'{len(exp)} proteins, wrong or missing {bad}, max difference {worst4:.1e} and {worst3:.1e}')
    scored = [r for r in rows if r['p4'] is not None]
    check('three-digit risk never above four-digit risk', all(r['p3'] <= r['p4'] for r in scored))

    # 3. matched mode gives as many statements as the identity rule
    rm = classify(qlens, per, ref, model, model['prediction_year'], 'match', 0.6, 0.4, 0.10, 0.10)
    n4 = sum(1 for r in rm if r['p4'] is not None and r['identity'] >= 0.6)
    n3 = sum(1 for r in rm if r['p4'] is not None and 0.4 <= r['identity'] < 0.6)
    check('match mode keeps the numbers of the identity rule', sum(r['statement'] == 'EC4' for r in rm) == n4 and sum(r['statement'] == 'EC3' for r in rm) == n3, f'{n4} and {n3}')

    # 4. edge cases built from the example
    hit_lines = [x for x in lines(h) if x.strip()]
    first = next(iter(exp))                       # a query with a neighbour
    orig = got[first]
    with tempfile.TemporaryDirectory(prefix='ecbackoff_selftest_') as tmp:
        # a query without any hit
        ql = dict(qlens)
        ql['no_hit_query'] = 120
        pr = read_hits(h, ref, list(ql), False)
        r0 = classify(ql, pr, ref, model, model['prediction_year'], 'cutoff', 0.6, 0.4, 0.10, 0.10)[-1]
        check('query without a neighbour gets none and a note', r0['statement'] == 'none' and r0['p4'] is None and r0['note'] != '')
        # an empty sequence
        ql2 = {'empty_query': 0}
        r1 = classify(ql2, {}, ref, model, model['prediction_year'], 'cutoff', 0.6, 0.4, 0.10, 0.10)[0]
        check('empty sequence gets none and a note', r1['statement'] == 'none' and r1['note'] == 'empty_sequence')
        # accession in the hit table, UniProt-style identifier in the FASTA
        qa = {f'sp|{first}|TEST_HUMAN': qlens[first]}
        hp = os.path.join(tmp, 'alias.m8')
        with open(hp, 'w') as fh:
            fh.write('\n'.join(x for x in hit_lines if x.split('\t')[0] == first) + '\n')
        ra = classify(qa, read_hits(hp, ref, list(qa), False), ref, model, model['prediction_year'], 'cutoff', 0.6, 0.4, 0.10, 0.10)[0]
        check('sp|ACCESSION|NAME identifier is matched to the accession in the hit table', ra['p4'] is not None and abs(ra['p4'] - orig['p4']) < 1e-12)
        # a hit to the query itself (named after a reference accession that is not among its neighbours) is dropped unless --keep-self
        targets = {x.split('\t')[1] for x in hit_lines if x.split('\t')[0] == first}
        acc = next(a for a in ref.labels if a not in targets)
        qs = {f'sp|{acc}|SELF_HUMAN': qlens[first]}
        hp2 = os.path.join(tmp, 'self.m8')
        with open(hp2, 'w') as fh:
            fh.write(f'{acc}\t{acc}\t1.000\t1.000\t1.000\t0\t5000\n')
            fh.write('\n'.join(x.replace(first, acc, 1) for x in hit_lines if x.split('\t')[0] == first) + '\n')
        r_drop = classify(qs, read_hits(hp2, ref, list(qs), False), ref, model, model['prediction_year'], 'cutoff', 0.6, 0.4, 0.10, 0.10)[0]
        r_keep = classify(qs, read_hits(hp2, ref, list(qs), True), ref, model, model['prediction_year'], 'cutoff', 0.6, 0.4, 0.10, 0.10)[0]
        check('self hit is dropped', r_drop['p4'] is not None and abs(r_drop['p4'] - orig['p4']) < 1e-12 and r_drop['nearest_enzyme'] == orig['nearest_enzyme'])
        check('self hit is kept with --keep-self', r_keep['nearest_enzyme'] == acc and r_keep['identity'] == 1.0)
        # malformed input stops with a message
        bad_hits = os.path.join(tmp, 'bad.m8')
        with open(bad_hits, 'w') as fh:
            fh.write(f'{first}\t{orig["nearest_enzyme"]}\t88.0\t1.0\t1.0\t0\t100\n')
        try:
            with contextlib.redirect_stderr(io.StringIO()):   # the refusal message is expected here
                read_hits(bad_hits, ref, list(qlens), False)
            check('percent identity in the hit table is refused', False)
        except SystemExit:
            check('percent identity in the hit table is refused', True)
        # identifiers of any form are hidden from the search, and a hit table pasted twice is not counted twice
        odd = os.path.join(tmp, 'odd.fasta')
        with open(odd, 'w', encoding='utf-8', newline='\n') as fh:
            fh.write('>蛋白_one comment\nMKV\n>sp|P12345|A_HUMAN\nMKVL\n>P12345\nMKVLA\n>emptyone\n')
        empty_ids, mp = write_clean_fasta(odd, os.path.join(tmp, 'odd_clean.fasta'))
        check('query identifiers of any form are replaced by synthetic ASCII ones for the search',
              list(mp.values()) == ['蛋白_one', 'sp|P12345|A_HUMAN', 'P12345'] and empty_ids == ['emptyone'] and all(k.isascii() for k in mp))
        raw_hits, tr_hits = os.path.join(tmp, 'raw.m8'), os.path.join(tmp, 'translated.m8')
        with open(raw_hits, 'w', newline='\n') as fh:
            fh.write('q0000002\tT1\t0.5\t1\t1\t0\t10\n')
        translate_hits(raw_hits, tr_hits, mp)
        with open(tr_hits, encoding='utf-8') as fh:
            check('the hit table gets the identifiers of the query FASTA back', fh.read().startswith('sp|P12345|A_HUMAN\tT1\t'))
        hp3 = os.path.join(tmp, 'twice.m8')
        with open(hp3, 'w') as fh:
            fh.write(('\n'.join(x for x in hit_lines if x.split('\t')[0] == first) + '\n') * 2)
        with contextlib.redirect_stderr(io.StringIO()):
            once = read_hits(h, ref, list(qlens), False)[first]
            twice = read_hits(hp3, ref, list(qlens), False)[first]
        check('a hit table pasted twice is not counted twice', len(once) == len(twice))
        buf_err = io.StringIO()
        with contextlib.redirect_stderr(buf_err):
            warn_accession_clashes(['sp|P12345|A_HUMAN', 'P12345'])
        check('two query identifiers with the same accession are reported', 'shared by several query identifiers' in buf_err.getvalue())

    # 5. experimental features: new-function score, orphan suggestions, import of other predictors' output
    expn = {k: d for k, d in exp.items() if 'p_new_function' in d and d['p_new_function'] != ''}
    worst_n = max((abs(got[k]['pn'] - float(d['p_new_function'])) for k, d in expn.items() if got[k]['pn'] is not None), default=1.0)
    check('new-function score reproduces the released values', bool(expn) and worst_n < 1e-6, f'{len(expn)} proteins, max difference {worst_n:.1e}')
    check('new-function score is a probability for every scored protein', all(r['pn'] is not None and 0 < r['pn'] < 1 for r in scored))
    orphans = load_orphans()
    keys = {';'.join(sorted(ref.labels[a])) for a in ref.labels}
    ranked_ok = all(v and [x[0] for x in v] == list(range(1, len(v) + 1)) and all(0 <= x[2] <= 3 for x in v) and all(v[i][2] >= v[i + 1][2] for i in range(len(v) - 1)) for v in orphans.values())
    check('orphan table covers every EC set of the reference and is ranked by shared digits', set(orphans) == keys and ranked_ok, f'{len(orphans)} EC sets')
    ro = classify(qlens, per, ref, model, model['prediction_year'], 'cutoff', 0.6, 0.4, 0.10, 0.10)
    attach_orphans(ro, ref, orphans, ORPHAN_MIN_P, ORPHAN_TOP)
    flagged = [r for r in ro if r['pn'] is not None and r['pn'] >= ORPHAN_MIN_P]
    check('orphan suggestions only for proteins at or above the new-function threshold', all(bool(r['orphan']) for r in flagged) and all(not r['orphan'] for r in ro if r not in flagged),
          f'{len(flagged)} of {len(ro)} proteins')
    qids = list(qlens)
    q1, q2, q3 = qids[0], qids[1], qids[2]
    e1 = sorted(ref.labels[got[q1]['nearest_enzyme']])[0]
    e1_3 = '.'.join(e1.split('.')[:3] + ['9999'])
    with tempfile.TemporaryDirectory(prefix='ecbackoff_selftest_ext_') as tmp:
        def write(name, text):
            path = os.path.join(tmp, name)
            with open(path, 'w', newline='\n') as fh:
                fh.write(text)
            return path
        f_dect = write('dect.tsv', f'sequence_ID\tprediction\n{q1}\t{e1}\n{q2}\t9.9.9.9\n{q3}\tNone\n')
        f_deepec = write('deepec.tsv', f'Query ID\tPredicted EC number\n{q1}_SEPARATED_SEQUENCE_1\tEC:{e1_3}\n{q2}\tEC number not predicted\n')
        f_clean = write('clean.csv', f'{q1},EC:{e1}/5.25,EC:9.9.9.9/7.1\n{q2},EC:9.9.9.9/3.0\n')
        f_pi = write('proteinfer.tsv', f'sequence_name\tpredicted_label\tconfidence\tdescription\n{q1}\tEC:{e1}\t0.91\tsome enzyme\n{q1}\tEC:9.9.9.9\t0.60\tother\n')
        f_json = write('pred.jsonl', json.dumps({'query_id': q1, 'labels': ['EC:' + e1], 'confidence': [0.8]}) + '\n')
        rd, rc_, rcl, rp, rj = (read_external(f) for f in (f_dect, f_deepec, f_clean, f_pi, f_json))
        check('external formats are read (DeepECtransformer, DeepEC, CLEAN, ProteInfer, JSON lines)',
              rd[q1][0][0] == e1 and rd[q3] == [] and rc_[q1][0][0] == e1_3 and rc_[q2] == [] and rcl[q1][0] == (e1, 5.25) and len(rcl[q1]) == 2
              and rp[q1][0] == (e1, 0.91) and len(rp[q1]) == 2 and rj[q1] == [(e1, 0.8)])
        r6 = classify(ql, pr, ref, model, model['prediction_year'], 'cutoff', 0.6, 0.4, 0.10, 0.10)   # the 30 examples plus a query without neighbour
        attach_external(r6, ql, ref, [('dectransformer', f_dect), ('deepec', f_deepec)])
        v = {r['query']: (r['ext']['dectransformer']['vs'], r['ext']['deepec']['vs']) for r in r6}
        check('external predictions are compared with the nearest enzyme',
              v[q1] == ('same_4digits', 'same_3digits') and v[q2] == ('differs', 'no_tool_answer') and v[q3][0] == 'no_tool_answer' and v['no_hit_query'][0] == 'no_neighbour',
              f'{v[q1]} {v[q2]} {v[q3][0]} {v["no_hit_query"][0]}')
        f_other = write('other.tsv', 'sequence_ID\tprediction\nsome_other_protein\t1.1.1.1\n')
        try:
            with contextlib.redirect_stderr(io.StringIO()):
                attach_external(r6, ql, ref, [('dectransformer', f_other)])
            check('an external file with none of the query identifiers is refused', False)
        except SystemExit:
            check('an external file with none of the query identifiers is refused', True)

    ok = all(results)
    print('[selftest] ' + ('PASS' if ok else 'FAIL'))
    sys.exit(0 if ok else 1)


def cmd_info(args):
    m = load_model()
    ref = Reference()
    print(f'ecbackoff {VERSION}')
    print(f'reference: {m["reference_release"]}, {len(ref.labels)} experimentally characterized enzymes, {len(ref.by_ec)} EC numbers')
    print(f'prediction year of the frozen model: {m["prediction_year"]}')
    print(f'search: {m["search"]["program"]} {" ".join(SEARCH_ARGS)}; neighbours need coverage >= {MIN_COV} of both sequences')
    print(f'model: logistic regression on {len(m["features"])} standardized features, trained on {m["wrong_4digit"]["training_rows"]} proteins (years {m["training_years"][0]}-{m["training_years"][1]})')
    print('evidence (retrospective; see README):')
    for line in m['evidence']:
        print('  ' + line)
    if 'new_function' in m:
        n = sum(len(v) for v in load_orphans().values())
        print(f'experimental (version 1.2): ranking score for a new function (logistic regression on the same 18 features, trained on {m["new_function"]["training_rows"]} proteins), '
              f'orphan-activity table with {n} ranked suggestions, import of CLEAN, DeepECtransformer, ProteInfer and DeepEC output')


def build_parser():
    ap = argparse.ArgumentParser(prog='ecbackoff', description=__doc__.split('\n\n')[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--version', action='version', version=f'ecbackoff {VERSION}')
    sub = ap.add_subparsers(dest='cmd', required=True)

    def common(p):
        p.add_argument('-q', '--query', required=True, help='FASTA of the query proteins (identifier = first word of the header)')
        p.add_argument('-o', '--out', default='-', help='output table (TSV); default standard output. A file <out>.params.json records how the table was made')
        p.add_argument('--no-params', action='store_true', help='do not write the .params.json file')
        p.add_argument('--mode', choices=['match', 'cutoff'], default='cutoff',
                       help='cutoff: four digits if P(wrong) <= --cutoff4, else three digits if P(wrong) <= --cutoff3, else nothing; works for any batch size (default). '
                            'match: as many four-digit and three-digit statements as the identity rule gives on this batch, chosen by model risk (the comparison made in the paper); '
                            'meaningful only for batches of several hundred proteins')
        p.add_argument('--identity-hi', type=float, default=0.6, help='identity rule: four digits at or above this (default 0.6)')
        p.add_argument('--identity-lo', type=float, default=0.4, help='identity rule: three digits at or above this (default 0.4)')
        p.add_argument('--cutoff4', type=float, default=0.10, help='cutoff mode: largest P(four-digit wrong) for a four-digit statement (default 0.10)')
        p.add_argument('--cutoff3', type=float, default=0.10, help='cutoff mode: largest P(three-digit wrong) for a three-digit statement (default 0.10)')
        p.add_argument('--year', type=int, default=0, help='prediction year (default: the frozen 2027). It enters two features, years since the nearest enzyme was '
                                                            'established and the number of earlier proteins carrying its EC number; other values are outside what was validated')
        p.add_argument('--keep-self', action='store_true', help='do not drop a hit that is the query itself (by identifier or UniProt accession)')
        p.add_argument('--no-orphan', action='store_true', help='experimental: do not suggest orphan enzyme activities')
        p.add_argument('--orphan-min-p', type=float, default=ORPHAN_MIN_P, help=f'experimental: suggest orphan activities for proteins with p_new_function at or above this (default {ORPHAN_MIN_P})')
        p.add_argument('--orphan-top', type=int, default=ORPHAN_TOP, help=f'experimental: number of orphan activities suggested (default {ORPHAN_TOP})')
        p.add_argument('--external', action='append', default=[], metavar='NAME=FILE',
                       help='experimental: add the output of another predictor that you ran yourself (NAME one of ' + ', '.join(EXTERNAL_TOOLS) +
                            '); can be repeated, e.g. --external clean=maxsep.csv --external deepec=DeepEC_Result.txt')

    r = sub.add_parser('run', help='FASTA in, table out (needs MMseqs2)')
    common(r)
    r.add_argument('--threads', type=int, default=4)
    r.add_argument('--mmseqs', default='', help='path to the mmseqs executable (default: from PATH)')
    r.add_argument('--keep-hits', default='', help='write the MMseqs2 hit table here instead of a temporary file')
    r.add_argument('--tmpdir', default='', help='directory for the temporary files (decompressed reference, MMseqs2 work files); default: the system temporary directory')
    r.set_defaults(func=cmd_run)

    s = sub.add_parser('score', help='MMseqs2 hit table in, table out')
    common(s)
    s.add_argument('--hits', required=True, help='hit table from: mmseqs easy-search QUERY reference.fasta HITS tmp -s 7.5 -e 1e-3 --max-seqs 20000 --alignment-mode 3 '
                                                  f'--format-output {FORMAT}  (reference.fasta: "ecbackoff reference -o reference.fasta")')
    s.set_defaults(func=cmd_score)

    f = sub.add_parser('reference', help='write the reference sequences for running the search yourself')
    f.add_argument('-o', '--out', default='reference.fasta', help='output FASTA (default reference.fasta; - for standard output)')
    f.set_defaults(func=cmd_reference)

    t = sub.add_parser('selftest', help='check the installation')
    t.set_defaults(func=cmd_selftest)
    i = sub.add_parser('info', help='model, reference and evidence summary')
    i.set_defaults(func=cmd_info)
    return ap


def main(argv=None):
    if sys.version_info < (3, 9):
        die('Python 3.9 or newer is needed')
    for stream in (sys.stdout, sys.stderr):   # a console with a narrow code page must not crash the program on an unusual identifier
        try:
            stream.reconfigure(errors='replace')
        except (AttributeError, ValueError, OSError):
            pass
    args = build_parser().parse_args(argv)
    try:
        args.func(args)
    except BrokenPipeError:
        try:
            os.dup2(os.open(os.devnull, os.O_WRONLY), sys.stdout.fileno())
        except OSError:
            pass
        sys.exit(0)
    except KeyboardInterrupt:
        sys.exit(130)


if __name__ == '__main__':
    main()
