#!/usr/bin/env python3
"""Audit B: decomposition of the HYP-9241 ratios from the two independent km_audit runs.
c_2 := q2/q1^3 = [q2/(q01 q12 q02)] * (q12/q01) * (q02/q01): independence factor times lag ratios.
Benchmark: for three independent Brownian walkers killed at the first meeting (Karlin-McGregor / vicious walkers)
P(no collision)/prod_pairs P(pair no collision) -> pi/4 = 0.785 as T -> infinity (computed below from Mehta's integral)."""
import json, math, glob
OUT = []
def say(s): print(s); OUT.append(s)
# Karlin-McGregor constant for 3 standard BMs: P_x(no collision by t) ~ h(x) t^{-3/2} * (1/prod_{j<3} j!) * (1/3!) E|Vandermonde(Z)|
EabsV = math.gamma(1.5)/math.gamma(1.5) * math.gamma(2)/math.gamma(1.5) * math.gamma(2.5)/math.gamma(1.5)   # Mehta, beta = 1, n = 3
c3 = (1/2) * (1/6) * EabsV
pair = 1/math.sqrt(math.pi)          # P(two BMs at distance x do not meet by t) ~ x/sqrt(pi t)
say(f"vicious walkers (3 BMs): P(no collision) ~ {c3:.6f} h(x) t^-1.5 ; product of pair probabilities ~ {pair**3:.6f} h(x) t^-1.5 ;"
    f" ratio = {c3/pair**3:.4f} (= pi/4 = {math.pi/4:.4f})")
for fn in sorted(glob.glob("km_audit_N*_s*.json")):
    tab = json.load(open(fn))
    say(f"-- {fn}")
    say("   T     q2/q1^3   indep2=q2/(q01q12q02)  q12/q01  q02/q01 | q3/q1^6  indep3=q3/prod6  prod6/q1^6  q03/q01")
    for r in tab:
        if r['q2'] == 0: continue
        prod6 = r['q01']*r['q12']*r['q23']*r['q02']*r['q13']*r['q03']
        s3 = f"{r['r3']:.3f}   {r['ind3']:.3f}            {prod6/r['q1']**6:.3f}" if r['q3'] > 0 else "-"
        say(f"  {r['T']:5d}   {r['r2']:.3f}     {r['ind2']:.3f}                  {r['q12']/r['q01']:.3f}    {r['q02']/r['q01']:.3f}  | {s3}       {r['q03']/r['q01']:.3f}")
with open("km_audit_summary.out", "w") as f: f.write("\n".join(OUT) + "\n")
