"""Replay the six exact proof-derived engines and the synthesis consequence."""
from pathlib import Path
import json
import subprocess
import sys

ROOT = Path(__file__).resolve().parent


def main():
    scripts = [
        'geometry/check_switching.py', 'geometry/digraph_candidate.py',
        'orbits/check_affine.py', 'constructive/examples.py',
        'sequences/apex_words.py', 'sequences/moment_recurrence.py',
    ]
    for script in scripts:
        p = ROOT / script
        print('\nREPLAY', script, flush=True)
        subprocess.run([sys.executable, str(p)], cwd=p.parent, check=True)
    cert = json.loads((ROOT/'sequences/moment_certificate.json').read_text())
    m = [(int(n == 0)+3*4**n+3*8**n+12**n)//8 for n in range(6)]
    c = [(-s+6*a-5*4**n)//4 for n,(s,a) in enumerate(zip(cert['initial_terms'],m))]
    assert c[:4] == [0, 1, 4, -442]
    sys.path.insert(0, str(ROOT/'sequences'))
    from apex_words import h
    assert all(h([0]*(n-1)+[2]) == 3 for n in range(1, 21))
    print('SYNTHESIS support-bound initial values:', c)
    print('Persistent H=3 construction checked through spine length20')
    print('ALL SIX ENGINES AND SYNTHESIS CONTROLS PASSED')


if __name__ == '__main__':
    main()
