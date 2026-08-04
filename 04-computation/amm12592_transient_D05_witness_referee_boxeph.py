"""Hostile audit of the NEW R=512 D0=5 rule-A witness using ONLY the
independent referee's primitives (fresh implementation lineage: exact
Lucas/Fibonacci floor, THM-3026 admissibility, THM-3002 identity).
No solver code imported."""
import json, os, sys, io, contextlib, importlib.util
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "iref", os.path.join(HERE, "amm12592_independent_witness_referee_boxeph.py"))
iref = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(iref)

W = json.load(open(os.path.join(HERE, "amm12592_witness_R512_ruleA_D0_5_boxeph.json")))
R, D0 = W["R"], W["D0"]
prof, blocks = W["profile"], W["blocks"]
assert R == 512 and D0 == 5 and len(blocks) == R

ok_prof = all(prof[i] == iref.floor_gamma_star(R + i) + 5 for i in range(R))
bad = []
for i, (d, row) in enumerate(zip(prof, blocks)):
    ok, why = iref.admissible(row, d)
    if not ok:
        bad.append((i, why))
S = []
for i, (d, row) in enumerate(zip(prof, blocks)):
    S = iref.padd(S, iref.pshift(iref.block_poly(row, d), i))
ok_id = iref.peq(S, iref.qpow(R - 1))
ok_unit = all(row[0] in (1, -1) for row in blocks)
from fractions import Fraction
eff = max(Fraction(d, R + i) for i, d in enumerate(prof))
print(f"profile == exact_floor + 5 everywhere: {ok_prof}")
print(f"all {R} blocks admissible (Lucas box + parity): {len(bad) == 0} {bad[:3]}")
print(f"epoch identity sum x^i Delta_i == q^{R-1}: {ok_id}")
print(f"forced unit cells delta_0 = +-1 every row: {ok_unit}")
print(f"effective rate max d_i/(R+i) = {eff} = {float(eff):.6f} "
      f"(gamma* = 0.597987...) -> {'BELOW' if eff < Fraction(598, 1000) else 'CHECK'}")
verdict = ok_prof and not bad and ok_id and ok_unit
print("VERDICT:", "PASS -- witness verified by independent referee machinery"
      if verdict else "FAIL")
