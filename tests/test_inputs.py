"""Input handling of ecbackoff. Most tests need no search; the end-to-end tests run only if MMseqs2 is on the PATH."""
import csv
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import ecbackoff

EXAMPLE = Path(ecbackoff.__file__).parent / 'example'


def cli(*args):
    return subprocess.run([sys.executable, '-m', 'ecbackoff', *args], capture_output=True, text=True, encoding='utf-8', errors='replace')


def table(path):
    with open(path, encoding='utf-8', newline='') as fh:
        return {r['query']: r for r in csv.DictReader(fh, delimiter='\t')}


def first_record():
    text = (EXAMPLE / 'queries.fasta').read_text(encoding='utf-8').split('>')[1]
    head, _, seq = text.partition('\n')
    return head.split()[0], seq.replace('\n', '')


class InputHandling(unittest.TestCase):
    def test_identifiers_of_any_form_are_hidden_from_the_search(self):
        with tempfile.TemporaryDirectory() as tmp:
            src = os.path.join(tmp, 'odd.fasta')
            with open(src, 'w', encoding='utf-8', newline='\n') as fh:
                fh.write('>蛋白_one x\nMKV\n>sp|P12345|A_HUMAN\nMKVL\n>P12345\nMKVLA\n>gi|12345|gb|AAA12345.1|\nMK\n>nothing\n')
            empty, mapping = ecbackoff.write_clean_fasta(src, os.path.join(tmp, 'clean.fasta'))
            self.assertEqual(empty, ['nothing'])
            self.assertEqual(list(mapping.values()), ['蛋白_one', 'sp|P12345|A_HUMAN', 'P12345', 'gi|12345|gb|AAA12345.1|'])
            self.assertTrue(all(k.isascii() and k.startswith('q') for k in mapping))
            with open(os.path.join(tmp, 'clean.fasta'), encoding='ascii') as fh:
                self.assertEqual(fh.read().count('>'), 4)

    def test_unwritable_output_is_refused_before_the_search(self):
        missing = os.path.join(tempfile.gettempdir(), 'ecbackoff_no_such_dir', 'x.tsv')
        r = cli('run', '-q', str(EXAMPLE / 'queries.fasta'), '-o', missing)
        self.assertEqual(r.returncode, 2, r.stderr)
        self.assertIn('does not exist', r.stderr)

    def test_threads_and_year_are_checked(self):
        r = cli('run', '-q', str(EXAMPLE / 'queries.fasta'), '--threads', '0')
        self.assertEqual(r.returncode, 2)
        self.assertIn('--threads', r.stderr)
        r = cli('score', '-q', str(EXAMPLE / 'queries.fasta'), '--hits', str(EXAMPLE / 'hits.m8'), '--year', '1900')
        self.assertEqual(r.returncode, 2)
        self.assertIn('--year', r.stderr)

    def test_a_hit_table_pasted_twice_is_counted_once(self):
        with tempfile.TemporaryDirectory() as tmp:
            hits = (EXAMPLE / 'hits.m8').read_text(encoding='utf-8')
            doubled = os.path.join(tmp, 'twice.m8')
            Path(doubled).write_text(hits + hits, encoding='utf-8', newline='\n')
            one, two = os.path.join(tmp, 'one.tsv'), os.path.join(tmp, 'two.tsv')
            a = cli('score', '-q', str(EXAMPLE / 'queries.fasta'), '--hits', str(EXAMPLE / 'hits.m8'), '-o', one, '--no-params')
            b = cli('score', '-q', str(EXAMPLE / 'queries.fasta'), '--hits', doubled, '-o', two, '--no-params')
            self.assertEqual((a.returncode, b.returncode), (0, 0), a.stderr + b.stderr)
            self.assertIn('repeated line', b.stderr)
            self.assertEqual(Path(one).read_bytes(), Path(two).read_bytes())

    def test_shared_accession_is_reported_when_scoring(self):
        with tempfile.TemporaryDirectory() as tmp:
            qid, seq = first_record()
            fa = os.path.join(tmp, 'clash.fasta')
            Path(fa).write_text(f'>{qid}\n{seq}\n>sp|{qid}|ONE_HUMAN\n{seq}\n', encoding='utf-8', newline='\n')
            r = cli('score', '-q', fa, '--hits', str(EXAMPLE / 'hits.m8'), '-o', os.path.join(tmp, 'o.tsv'), '--no-params')
            self.assertIn('shared by several query identifiers', r.stderr)


@unittest.skipUnless(shutil.which('mmseqs'), 'MMseqs2 is not on the PATH')
class EndToEnd(unittest.TestCase):
    def test_example_run_reproduces_the_released_scores(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = os.path.join(tmp, 'out.tsv')
            r = cli('run', '-q', str(EXAMPLE / 'queries.fasta'), '-o', out, '--no-params', '--threads', '2', '--tmpdir', tmp)
            self.assertEqual(r.returncode, 0, r.stderr)
            got, exp = table(out), table(str(EXAMPLE / 'expected.tsv'))
            self.assertEqual(set(got), set(exp))
            for q, e in exp.items():
                self.assertEqual(got[q]['nearest_enzyme'], e['nearest_enzyme'], q)
                self.assertLess(abs(float(got[q]['p_wrong_4digit']) - float(e['p_wrong_4digit'])), 1e-4, q)
                self.assertLess(abs(float(got[q]['p_wrong_3digit']) - float(e['p_wrong_3digit'])), 1e-4, q)

    def test_identifier_form_does_not_change_the_result(self):
        qid, seq = first_record()
        ids = [qid, '蛋白_one', f'sp|{qid}|ONE_HUMAN', 'lcl|myid_5', 'gi|12345|gb|AAA12345.1|']
        with tempfile.TemporaryDirectory() as tmp:
            fa, out = os.path.join(tmp, 'ids.fasta'), os.path.join(tmp, 'out.tsv')
            Path(fa).write_text(''.join(f'>{i}\n{seq}\n' for i in ids), encoding='utf-8', newline='\n')
            r = cli('run', '-q', fa, '-o', out, '--no-params', '--threads', '2', '--tmpdir', tmp)
            self.assertEqual(r.returncode, 0, r.stderr)
            got = table(out)
            self.assertEqual(list(got), ids)
            self.assertEqual({got[i]['p_wrong_4digit'] for i in ids}, {got[qid]['p_wrong_4digit']})
            self.assertEqual({got[i]['nearest_enzyme'] for i in ids}, {got[qid]['nearest_enzyme']})


if __name__ == '__main__':
    unittest.main()
