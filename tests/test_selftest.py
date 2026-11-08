"""Smoke tests of the installed package: version and the built-in self-test (no MMseqs2 needed)."""
import subprocess
import sys
import unittest


def run(*args):
    r = subprocess.run([sys.executable, '-m', 'ecbackoff', *args], capture_output=True, text=True)
    return r.returncode, r.stdout + r.stderr


class SmokeTests(unittest.TestCase):
    def test_version(self):
        code, out = run('--version')
        self.assertEqual(code, 0, out)
        self.assertIn('ecbackoff', out)

    def test_selftest_passes(self):
        code, out = run('selftest')
        self.assertEqual(code, 0, out)
        self.assertIn('[selftest] PASS', out)


if __name__ == '__main__':
    unittest.main()
