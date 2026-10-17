#!/usr/bin/env python3
"""Conditional Arb certificates (collatz_pole_certificate_other_zeros_20261015.py, budget rows (2, 0) and (64, 16)) at the
zeros of D whose boundary residue estimate most exceeds a rough radius model c |rho| N^(-(sigma - 1/2)), with exact block
moments through N = 2^42 (S38).  CONDITIONAL: the energy budgets are unproved; the intervals are rigorous given them.
Session opus-2026-10-16-S38.  Usage: python collatz_pole_certificates_many_zeros_2p42_20261016.py moments.npz census.json
       candidates.json out.json [k]
"""
import sys, json, subprocess, os, tempfile

MOM, CEN, CAND, OUT = sys.argv[1:5]
K = int(sys.argv[5]) if len(sys.argv) > 5 else 16
here = os.path.dirname(os.path.abspath(__file__))
cands = json.load(open(CAND))[:K]
rows = []
for c in cands:
    tmp = os.path.join(tempfile.gettempdir(), f"cert_{c['re']:.6f}_{c['im']:.4f}.json")
    subprocess.run([sys.executable, '-B', os.path.join(here, 'collatz_pole_certificate_other_zeros_20261015.py'),
                    '--moments', MOM, '--census', CEN, '--output', tmp, '--re', repr(c['re']), '--im', repr(c['im'])],
                   check=True, stdout=subprocess.DEVNULL)
    R = json.load(open(tmp))
    r20 = next(r for r in R['rows'] if int(r['a']) == 2 and int(r['b']) == 0)
    r64 = next(r for r in R['rows'] if int(r['a']) == 64 and int(r['b']) == 16)
    rows.append(dict(zero=R['zero'][:2] if isinstance(R['zero'], list) else R['zero'], re=c['re'], im=c['im'],
                     residue_before_tail=R['residue_before_future_tail'],
                     abs_r_2_0=[r20['abs_residue_lower'], r20['abs_residue_upper']], certified_2_0=r20['nonzero_certified'],
                     abs_r_64_16=[r64['abs_residue_lower'], r64['abs_residue_upper']], certified_64_16=r64['nonzero_certified']))
    print(f"{c['re']:.6f} + {c['im']:.4f}i: (2,0) {'CERTIFIED' if r20['nonzero_certified'] else 'no'}, "
          f"(64,16) {'CERTIFIED' if r64['nonzero_certified'] else 'no'}; residue {R['residue_before_future_tail']}", flush=True)
json.dump(dict(N=json.load(open(CEN))['N'], rows=rows,
               certified_2_0=sum(1 for r in rows if r['certified_2_0']),
               certified_64_16=sum(1 for r in rows if r['certified_64_16'])), open(OUT, 'w'), indent=1)
print(f"certified under (2,0): {sum(1 for r in rows if r['certified_2_0'])} of {len(rows)}; "
      f"under (64,16): {sum(1 for r in rows if r['certified_64_16'])}")
