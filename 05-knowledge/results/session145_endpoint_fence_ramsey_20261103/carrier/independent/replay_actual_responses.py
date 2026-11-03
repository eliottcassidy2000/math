#!/usr/bin/env python3
"""Rebuild and replay the actual-graph counter in a temporary directory."""
from pathlib import Path
import subprocess,sys,tempfile
HERE=Path(__file__).resolve().parent
from audit_carrier import groups
_,signatures=groups(6)
expected={m:(signature,sign) for signature,rows in signatures.items() for m,sign in rows}
actual=[tuple(map(int,line.split())) for line in (HERE/'coefficient_map.txt').read_text().splitlines()]
if len(actual)!=len(expected) or len({m for m,key,sign in actual})!=len(expected):
 raise RuntimeError('coefficient map must contain every forest exactly once')
keys={}
for mask,key,sign in actual:
 if mask not in expected or expected[mask][1]!=sign:raise RuntimeError('coefficient map mask/sign mismatch')
 signature=expected[mask][0]
 if key in keys and keys[key]!=signature:raise RuntimeError('coefficient map merges different signatures')
 keys[key]=signature
if set(keys)!=set(range(571)) or len(set(keys.values()))!=571:
 raise RuntimeError('coefficient map must index all 571 signatures bijectively')
with tempfile.TemporaryDirectory(prefix='s145-actual-graphs-') as d:
 executable=Path(d)/'actual_responses'
 subprocess.run(['g++','-O2','-std=c++17',str(HERE/'actual_responses.cpp'),'-o',str(executable)],check=True)
 result=subprocess.run([str(executable),str(HERE/'coefficient_map.txt')],capture_output=True,check=True)
 sys.stdout.buffer.write(result.stdout)
