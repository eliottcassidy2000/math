#!/usr/bin/env python3
"""Check fail-closed behavior in disposable copies, both normal and optimized."""
from pathlib import Path
from tempfile import TemporaryDirectory
import json
import shutil
import subprocess
import sys


ROOT = Path(__file__).resolve().parent
MUTATIONS = [
    ('bad_field_reduction', 'c[k-3] -= c[k]', 'c[k-3] -= 2*c[k]', 'z^3 reduction'),
    ('wrong_extension_formula', 'n125 == 126-t**3+15*t', 'n125 == 126-t**3+14*t', 'Frobenius-cube formula'),
    ('wrong_magic_entry', 'witness = [[108,148,122]', 'witness = [[109,148,122]', 'magic witness entries'),
    ('confused_hecke_cube', 'hecke3 == (0,-10,0,1)', 'hecke3 == (0,-15,0,1)', 'Hecke / symmetric-cube coefficient'),
    ('wrong_fiber_claim', 'tfibers == [3,2,1,1,1,1,1,1,1]', 'tfibers == [2,2,2,1,1,1,1,1,1]', 'tournament spectrum fibers'),
]


def main():
    source = (ROOT/'check_bridge.py').read_text()
    expected = json.loads((ROOT/'result.json').read_text())
    outcomes = []
    for optimized in (False, True):
        command = [sys.executable] + (['-O'] if optimized else [])
        baseline = subprocess.run(command + [str(ROOT/'check_bridge.py')], capture_output=True, text=True)
        if baseline.returncode or json.loads(baseline.stdout) != expected:
            raise RuntimeError('positive control did not reproduce result')
        outcomes.append({'case':'unmodified_positive','optimized':optimized,'accepted':True})
        for name, old, new, message in MUTATIONS:
            if source.count(old) != 1:
                raise RuntimeError('mutation target must be unique: '+name)
            with TemporaryDirectory(prefix='frobenius-hostile-') as temp:
                root = Path(temp)
                shutil.copytree(ROOT/'sources', root/'sources')
                (root/'check_bridge.py').write_text(source.replace(old,new))
                run = subprocess.run(command + [str(root/'check_bridge.py')], capture_output=True, text=True)
                if run.returncode == 0 or message not in run.stderr:
                    raise RuntimeError('mutation not rejected at intended gate: '+name)
                outcomes.append({'case':name,'optimized':optimized,'rejected':True,
                                 'expected_gate':message})
    print(json.dumps({'status':'PASS','runs':len(outcomes),'outcomes':outcomes},indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
