"""Re-run the session's census script UNCHANGED on the pre-registered corpus (HEAD minus THM-4839, which is the corpus of
commit 593105aedd: the checkpoint 4e32454073 only added THM-4839 inside the window), under a chosen PYTHONHASHSEED.

The session script is executed from its committed source; the only intervention is that glob.glob drops THM-4839
(otherwise the run includes the session's own theorem and, on Windows with a long checkout path, crashes on the
260-character path limit at that file).

Usage (from the repository root):  PYTHONHASHSEED=<k> python rerun_session_census.py 4400 4839 2000
"""
import glob
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "attractor_census_canon_20261028.py")
_orig_glob = glob.glob


def _glob(pat, *a, **k):
    return [f for f in _orig_glob(pat, *a, **k) if not os.path.basename(f).startswith("THM-4839-")]


glob.glob = _glob
sys.argv = [SRC] + sys.argv[1:]
code = compile(open(SRC, encoding="utf-8").read(), SRC, "exec")
exec(code, {"__name__": "__main__", "__file__": SRC})
