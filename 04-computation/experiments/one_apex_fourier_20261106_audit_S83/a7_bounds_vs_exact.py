# Audit S83, step 7: THM-4886 (A) and the audit's two-block bound against EXACT |S_n(theta)| for every admissible q:
# all reduced r/q with |2 theta - r/q| <= 1/q^2 and q <= 2^min(n, 18), enumerated exhaustively in exact integer arithmetic.
import math, time
from fractions import Fraction
import numpy as np
from oa_engine import g_block, S_exact, IRR
from oa_bounds import thmA, twoblock

def admissible(beta, qmax):
    # beta: Fraction; returns list of q (reduced fractions r/q with |beta - r/q| <= 1/q^2)
    num, den = beta.numerator, beta.denominator
    out = []
    for q in range(1, qmax + 1):
        # r nearest to beta*q; |beta - r/q| <= 1/q^2  <=>  |num*q - r*den| * q <= den
        r = (2 * num * q + den) // (2 * den)
        for rr in (r - 1, r, r + 1):
            if abs(num * q - rr * den) * q <= den and math.gcd(rr, q) == 1:
                out.append(q); break
    return out

thetas = dict(IRR)
for fr in (Fraction(1, 7), Fraction(3, 37), Fraction(5, 101), Fraction(11, 93), Fraction(1, 127), Fraction(1, 4095),
           Fraction(1234, 65537), Fraction(7, 64), Fraction(1, 12), Fraction(1, 3), Fraction(3, 16)):
    thetas[str(fr)] = fr

tot = 0; vA = 0; vB = 0; worstA = 0.0; worstB = 0.0; nontrivA = 0; nontrivB = 0; betterB = 0
rows = []
for n in range(12, 27, 2):
    t0 = time.time()
    vals = g_block(n)
    for name, th in thetas.items():
        s = abs(S_exact(vals, th)) / 2 ** (n - 1)
        qs = admissible(2 * th, 1 << min(n, 18))
        bestA = (math.inf, None); bestB = (math.inf, None)
        for q in qs:
            bA, _ = thmA(n, q); bB, _ = twoblock(n, q)
            tot += 1
            vA += s > bA * (1 + 1e-12); vB += s > bB * (1 + 1e-12)
            worstA = max(worstA, s / bA); worstB = max(worstB, s / bB)
            nontrivA += bA < 1; nontrivB += bB < 1; betterB += bB < bA
            if bA < bestA[0]: bestA = (bA, q)
            if bB < bestB[0]: bestB = (bB, q)
        rows.append((n, name, s, len(qs), bestA, bestB))
    print(f'n = {n} done [{time.time() - t0:.1f}s]', flush=True)
print(f'\n{tot} (n, theta, q) instances (n = 12, 14, ..., 26; 16 theta; every admissible q <= 2^min(n,18))')
print(f'  THM-4886 (A):   violations {vA}; max exact/bound {worstA:.4f}; nontrivial (< 1) in {nontrivA}')
print(f'  two-block (U1): violations {vB}; max exact/bound {worstB:.4f}; nontrivial (< 1) in {nontrivB}; smaller than (A) in {betterB}')
print('\nbest bound over admissible q (n, theta, exact |S|/2^(n-1), #q, (A) best (bound, q), two-block best (bound, q)):')
for r in rows:
    if r[0] in (16, 20, 24, 26):
        print(f'  n={r[0]:2d} {r[1]:>11s} exact={r[2]:.5f} #q={r[3]:3d}  A={r[4][0]:.4f} (q={r[4][1]})  U1={r[5][0]:.4f} (q={r[5][1]})')
