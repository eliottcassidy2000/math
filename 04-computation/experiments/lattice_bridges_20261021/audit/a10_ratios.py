#!/usr/bin/env python3
"""AUDIT of HYP-9362 (b) "1 - margin decays roughly like 2.5/q^2" and (c) the rotational/Paley ratios and the
normalised values H 2^(q-1)/q! quoted as "Paley normalised value (2.55, 2.60, 2.63, 2.66)".
Inputs: own computations (a3_classes_*.txt, a4_q31.out, a10_rot.out); H(P_27) = 418425945688463432031 is taken from
the repo (two earlier independent methods: paley_extremal_20261021/P27.out and paley_h31_audit_20261014), since P_27 is
not a circulant and a3_hcount cannot compute it.
Output: a10_ratios.out
"""
import os, json, re
from math import factorial
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
out = open(os.path.join(HERE, 'a10_ratios.out'), 'w')
def say(*a):
    s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n')
def classes(n):
    return {tuple(json.loads(l)['S']): json.loads(l)['H'] for l in open(os.path.join(HERE, f'a3_classes_{n}.txt'))}
HP, HF, HR = {}, {}, {}
for q in (7, 11, 19, 23):
    C = classes(q); m = (q - 1) // 2
    QR = sorted({x * x % q for x in range(1, q)})
    from math import gcd
    U = [a for a in range(1, q) if gcd(a, q) == 1]
    canon = lambda S: min(tuple(sorted((a * s) % q for s in S)) for a in U)
    HP[q] = C[canon(QR)]
    HF[q] = C[canon(sorted((set(QR) - {1}) | {q - 1}))]
    HR[q] = C[canon(range(1, m + 1))]
for line in open(os.path.join(HERE, 'a4_q31.out')):
    mm = re.match(r'31 ([\d ]+) : (\d+)', line)
    if mm:
        S = tuple(int(x) for x in mm.group(1).split()); h = int(mm.group(2))
        if S == tuple(sorted({x * x % 31 for x in range(1, 31)})): HP[31] = h
        else: HF[31] = h
rot = os.path.join(HERE, 'a10_rot.out')
if os.path.exists(rot):
    for line in open(rot):
        mm = re.match(r'(\d+) ([\d ]+) : (\d+)', line)
        if mm: HR[int(mm.group(1))] = int(mm.group(3))
HP[27] = 418425945688463432031
say("q   H(P_q)                          H(flip)/H(P)   (1 - margin) q^2")
for q in sorted(HF):
    mg = Fraction(HF[q], HP[q])
    say(f"{q:2d}  {HP[q]:<32d} {float(mg):.6f}      {float((1 - mg) * q * q):.3f}")
say("")
say("q   H(R_q)/H(P_q)   rotational H 2^(q-1)/q!   Paley H 2^(q-1)/q!")
for q in sorted(set(HP) & set(HR)):
    if q < 19: continue
    nr = Fraction(HR[q] * 2 ** (q - 1), factorial(q)); npl = Fraction(HP[q] * 2 ** (q - 1), factorial(q))
    say(f"{q:2d}  {HR[q]/HP[q]:.6f}        {float(nr):.5f}                  {float(npl):.5f}")
for n in sorted(HR):
    say(f"rotational n={n}: H = {HR[n]}, normalised {float(Fraction(HR[n] * 2 ** (n - 1), factorial(n))):.5f}")
out.close()
