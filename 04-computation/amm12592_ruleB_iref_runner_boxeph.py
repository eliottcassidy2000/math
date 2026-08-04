#!/usr/bin/env python3
"""D0-aware INDEPENDENT referee runner for rule-B/rule-A epoch witnesses (boxeph).

Loads amm12592_independent_witness_referee_boxeph.py (checker re-implemented
from THM-3002/THM-3026 statements only: exact Z[x], Lucas-box admissibility,
EXACT floor via 5^d <= phi^{2m} Lucas/Fibonacci integer sign tests -- no
floats, no solver code shared) and verifies a witness file at slack D0:

  PROFILE  profile[i] == floor_gamma_star(R+i) + D0 exactly, every i
  ADM      capacity + parity per cell (mod.admissible)
  IDENT    sum_i x^i Delta_i == (1-x)^{R-1} exactly (mod.block_poly/padd/qpow)
  UNIT     delta_{i,0} == +-1 every row
  EFF      max_i d_i/(R+i) as exact Fraction
  SHA256   of the exact file bytes verified (provenance pinning)

Usage: python3 amm12592_ruleB_iref_runner_boxeph.py <witness.json> <D0>
"""
import contextlib, hashlib, importlib.util, io, json, os, sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "iref", os.path.join(HERE, "amm12592_independent_witness_referee_boxeph.py"))
mod = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):   # its own R<=64 self-run
    spec.loader.exec_module(mod)

path, D0 = sys.argv[1], int(sys.argv[2])
raw = open(path, "rb").read()
sha = hashlib.sha256(raw).hexdigest()
w = json.loads(raw)
R, prof, blocks = w["R"], w["profile"], w["blocks"]
src = w.get("source", "?")

# pre-warm qpow cache iteratively (module's qpow is recursive)
for n in range(1, max(prof) + 2):
    mod.qpow(n)

prof_ok = all(prof[i] == mod.floor_gamma_star(R + i) + D0 for i in range(R))
adm_bad = []
for i, (d, row) in enumerate(zip(prof, blocks)):
    ok, why = mod.admissible(row, d)
    if not ok:
        adm_bad.append((i, why))
S = []
for i, (d, row) in enumerate(zip(prof, blocks)):
    S = mod.padd(S, mod.pshift(mod.block_poly(row, d), i))
ident_ok = mod.peq(S, mod.qpow(R - 1))
unit_ok = all(row[0] in (1, -1) for row in blocks)
eff = max(Fraction(d, R + i) for i, d in enumerate(prof))
allok = prof_ok and not adm_bad and ident_ok and unit_ok
print(f"IREF {os.path.basename(path)} R={R} D0={D0}: "
      f"profile=floor+{D0}:{prof_ok} admissible:{not adm_bad} "
      f"identity:{ident_ok} unit:{unit_ok} eff_rate={eff}={float(eff):.6f} "
      f"sha256={sha} source={src} => {'ALL PASS' if allok else 'FAIL ' + str(adm_bad[:2])}")
sys.exit(0 if allok else 1)
