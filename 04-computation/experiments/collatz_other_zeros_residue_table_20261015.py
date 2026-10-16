#!/usr/bin/env python3
"""Residues of F(s) = sum t(n) n^(-s) at the 28 zeros of D with Re > 0.65, Im < 3000: two independent estimates.

(1) CONDITIONAL interval (Arb): collatz_pole_certificate_other_zeros_20261015.py at N = 2^37 under the energy budget
    (a, b) = (2, 0) for dyadic X >= 2^37 (finite boundary numerator of Theorem B3 plus a rigorous tail bound).
(2) Least-squares: the coefficient c_rho of the explicit-formula fit S(N) ~ Re sum_rho c_rho N^rho on [2^22, 2^37)
    (collatz_mercator_waveform_prediction_20261015.json); a pole of F with residue r gives S(N) ~ Re(2 r N^rho / rho),
    so |r_fit| = |c_rho| |rho| / 2.
The table lists both, and whether the conditional interval certifies a nonzero residue.
Session opus-2026-10-15-S37.  Usage: python collatz_other_zeros_residue_table_20261015.py moments.npz census.json pred.json out.json
"""
import sys, json, subprocess, os, tempfile


def arbf(text, side):
    # parse an Arb string '[m +/- r]' (or a plain number) to a float bound on the given side
    t = text.strip()
    if t.startswith('['):
        m, r = t.strip('[]').split('+/-')
        m, r = float(m), float(r)
        return m - r if side == 'lo' else m + r
    return float(t)

MOM, CEN, PRED, OUT = sys.argv[1:5]
HERE = os.path.dirname(os.path.abspath(__file__))
P = json.load(open(PRED))
rows = []
for (a, b), c in zip(P['zeros'], P['coefficients']):
    re_, im_ = float(a), float(b)
    cabs = abs(complex(*c))
    rfit = cabs * abs(complex(re_, im_)) / 2
    tmp = os.path.join(tempfile.gettempdir(), f'oz_{re_:.4f}_{im_:.2f}.json')
    subprocess.run([sys.executable, os.path.join(HERE, 'collatz_pole_certificate_other_zeros_20261015.py'), '--moments', MOM,
                    '--census', CEN, '--output', tmp, '--re', repr(re_), '--im', repr(im_)], check=True,
                   stdout=subprocess.DEVNULL)
    R = json.load(open(tmp))
    r20 = R['rows'][0]
    lo, hi = arbf(r20['abs_residue_lower'], 'lo'), arbf(r20['abs_residue_upper'], 'hi')
    rows.append(dict(zero=[re_, im_], abs_r_fit=rfit, abs_r_interval_2_0=[lo, hi], nonzero_certified=r20['nonzero_certified'],
                     fit_inside_interval=bool(lo <= rfit <= hi)))
    print(f'rho = {re_:.6f} + {im_:9.4f} i : |r_fit| = {rfit:.5g}; conditional |r| in [{lo:.5g}, {hi:.5g}] nonzero {r20["nonzero_certified"]}; '
          f'fit inside {lo <= rfit <= hi}', flush=True)
json.dump(dict(N=2 ** 37, budget='(2, 0) for dyadic X >= 2^37', rows=rows), open(OUT, 'w'), indent=1)
