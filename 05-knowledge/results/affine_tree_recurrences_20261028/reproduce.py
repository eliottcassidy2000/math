#!/usr/bin/env python3
"""Replay the finite research receipts, without assertions or third-party tools."""
from pathlib import Path
import contextlib
import importlib.util
import io
import json
import subprocess
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
expected=json.loads((HERE/'checks.json').read_text())
for options in ([],['-O']):
    raw=subprocess.check_output([sys.executable,'-B',*options,str(HERE/'check_recurrences.py')],text=True)
    if json.loads(raw)!=expected:
        raise RuntimeError('changed recurrence receipt: '+str(options))
    print('PASS', 'optimized' if options else 'normal', expected['checks'], 'checks')
spec=importlib.util.spec_from_file_location('historical',ROOT/'04-computation/experiments/collatz_catalan_ramsey_20260927.py')
module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
buffer=io.StringIO()
with contextlib.redirect_stdout(buffer):
    module.part1(); module.part2(); module.part5()
if buffer.getvalue()!=(HERE/'historical_checks.txt').read_text():
    raise RuntimeError('changed historical correction receipt')
print('PASS corrected historical sections 1, 2, 5')
