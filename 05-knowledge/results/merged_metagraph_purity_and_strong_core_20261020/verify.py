#!/usr/bin/env python3
"""Portable replay: standard library only; no repository or census dependencies."""
from pathlib import Path
import subprocess
import sys

base=Path(__file__).resolve().parent
for stem in ('probe_n5_rewiring','verify_transitive_mixed_landing','verify_strong_core_counts'):
    result=subprocess.run([sys.executable,str(base/(stem+'.py'))],capture_output=True,check=True)
    expected=(base/(stem+'.out')).read_bytes()
    if result.stdout!=expected:
        raise AssertionError(f'{stem}: stdout differs from saved exact output')
    print(f'PASS {stem}: exact stdout replay')
print('PASS: all three verification programs; no repository imports or external packages.')
