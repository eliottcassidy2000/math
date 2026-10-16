#!/usr/bin/env python3
"""S37 audit 2, script b6b: run the SESSION's collatz_pole_certificate_other_zeros_20261015.py (unchanged, output to a
temporary file outside the repository) at a zero of D with Re > 1, where F(s) = sum t(n) n^-s is holomorphic and the
residue must vanish.  A correct pipeline must return an interval for |r| that contains 0 (no certification).
D has no zeros with 1 < Re s < 2 (2^-sigma + (3/2)^sigma/3 < 1 there), so the test zero has Re = 2.3029.
Usage: python b6b_session_script_at_re_gt_1.py moments37.npz moments37.json TMPDIR
"""
import sys, os, json, subprocess
HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.normpath(os.path.join(HERE, '..', 'collatz_pole_certificate_other_zeros_20261015.py'))
mom, cen, tmp = sys.argv[1:4]
for re_, im_ in (('2.30288', '23.6037'), ('0.99820922', '240.2478301')):
    out = os.path.join(tmp, f'b6b_{re_}_{im_}.json')
    subprocess.run([sys.executable, SCRIPT, '--moments', mom, '--census', cen, '--output', out, '--re', re_, '--im', im_],
                   check=True, stdout=subprocess.DEVNULL)
    R = json.load(open(out))
    for row in R['rows']:
        print(f"zero {R['zero'][0][:10]} + {R['zero'][1][:10]} i, budget ({row['a']}, {row['b']}): |r| lower "
              f"{row['abs_residue_lower'][:40]}, upper {row['abs_residue_upper'][:40]}, certified {row['nonzero_certified']}")
    print(f"  residue before tail {R['residue_before_future_tail']}")
