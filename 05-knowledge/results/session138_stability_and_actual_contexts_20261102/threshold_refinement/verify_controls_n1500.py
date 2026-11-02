#!/usr/bin/env python3
"""Check normal/-O agreement and explicit failure after a false-check injection."""
from pathlib import Path
import subprocess
import sys

root = Path(__file__).resolve().parent
verifier = root / 'verify_n1500.py'
source = verifier.read_text()
needle = 'failed = [name for name, passed in checks.items() if not passed]'
if source.count(needle) != 1:
    raise RuntimeError('Expected exactly one explicit failure guard')
bad = source.replace(needle, "checks['negative_control_injected_false'] = False\n" + needle)
outputs = []
for optimized in (False, True):
    flags = ['-O'] if optimized else []
    label = 'optimized' if optimized else 'normal'
    good = subprocess.run([sys.executable, *flags, str(verifier)],
                          capture_output=True, text=True, check=False)
    if good.returncode != 0 or 'exact_endpoint_checks: 11 / 11 PASS' not in good.stdout:
        raise RuntimeError(f'{label}: positive verifier failed')
    outputs.append(good.stdout)
    print(f'{label}_positive: PASS (11/11 exact checks)')
    control = subprocess.run([sys.executable, *flags, '-c', bad],
                             capture_output=True, text=True, check=False)
    expected = 'RuntimeError: Exact endpoint verification failed: negative_control_injected_false'
    if control.returncode == 0 or expected not in control.stderr:
        raise RuntimeError(f'{label}: false-check negative control did not fail explicitly')
    print(f'{label}_negative_control: PASS (explicit RuntimeError, exit {control.returncode})')
if outputs[0] != outputs[1]:
    raise RuntimeError('Normal and optimized positive outputs differ')
print('normal_optimized_output_identity: PASS')
print('controls: 5 / 5 PASS')
