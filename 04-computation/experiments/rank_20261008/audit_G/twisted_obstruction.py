#!/usr/bin/env python3
"""Audit G: a 'twisted' congruence obstruction missed by HYP-9244's obstruction lemma.
Z_3 map m = (1,5,7), r = (0,1,1) (HYP-9244 evidence row 'Z_3: 1, 5, 7 ... 1.000 at T = 4096'; audit F's obstruction scan
reported no obstruction for l < 2000).  Claim: e_t = e_0 + A_t + B_t (mod 2) along the pair chain, where
M_t = 5^A_t 7^B_t.  Reason: d = 3 and all m_i are odd, so m_i/d = 1 and M = 1 (mod 2); hence e' = e + r_i - r_j (mod 2),
and r_i = psi(m_i) (mod 2) for the character psi(5^a 7^b) = a + b (mod 2) (psi(1) = 0 = r_0, psi(5) = 1 = r_1,
psi(7) = 1 = r_2), while M' = M m_i/m_j shifts psi(M) by psi(m_i) - psi(m_j).  At M = 1: e = e_0 (mod 2), so odd offsets
never merge.  The HYP-9244 lemma needs (m_i - d)c + r_i = 0 (mod l) for a CONSTANT c: impossible here for l = 2.
General form (audit suggestion): if l does not divide d, all m_i = d (mod l), and r_i = d psi(m_i) (mod l) for a
homomorphism psi: <m_0..m_(d-1)> -> Z/l, then e_n - psi(M_n) = e_0 (mod l)."""
from fractions import Fraction as Fr
import random
d = 3; m = [1, 5, 7]; r = [0, 1, 1]; expo = {1: (0, 0), 5: (1, 0), 7: (0, 1)}
def res(x, q): return (x.numerator * pow(x.denominator, -1, q)) % q
rnd = random.Random(5); bad = 0; checked = 0; m1_odd_zero = 0
for _ in range(3000):
    e0 = rnd.choice([1, 3, 5, -1, 7, 2, 4, 6]); M = Fr(1); e = Fr(e0); A = B = 0
    for t in range(150):
        j = rnd.randrange(d); i = (res(M, d) * j + res(e, d)) % d
        e = (m[i] * e + r[i] - Fr(m[i], m[j]) * M * r[j]) / d; M = M * Fr(m[i], m[j])
        A += expo[m[i]][0] - expo[m[j]][0]; B += expo[m[i]][1] - expo[m[j]][1]
        checked += 1
        if (res(e, 2) - (e0 + A + B)) % 2 != 0: bad += 1
        if M == 1 and e == 0 and e0 % 2: m1_odd_zero += 1
print(f"invariant e_t = e_0 + A_t + B_t (mod 2): {checked} chain steps checked, violations {bad}; "
      f"absorptions from odd e_0: {m1_odd_zero}")
# HYP-9244 lemma (constant c) at l = 2: is there c with (m_i - d) c + r_i = 0 mod 2 for all i?
print("constant-c lemma at l = 2 applicable:", any(all(((mi - d) * c + ri) % 2 == 0 for mi, ri in zip(m, r)) for c in (0, 1)))
# which other HYP-9244 / audit maps have a mod-2 twisted obstruction?  (d odd, all m_i odd, psi consistent)
def twisted2(d, m, r):
    if d % 2 == 0 or any(x % 2 == 0 for x in m): return None
    # psi must satisfy psi(m_i) = r_i mod 2; consistency only matters for repeated or dependent multipliers;
    # check: same multiplier => same parity of r; multiplier 1 => r even
    seen = {}
    for mi, ri in zip(m, r):
        if mi == 1 and ri % 2: return False
        if mi in seen and seen[mi] != ri % 2: return False
        seen[mi] = ri % 2
    return True    # (multiplicative relations among distinct multipliers would need checking; none below)
for name, d_, m_, r_ in [('Z3 (1,5,7) r=(0,1,1)', 3, [1, 5, 7], [0, 1, 1]), ('Z3 (1,1,5) r=(0,2,2)', 3, [1, 1, 5], [0, 2, 2]),
                         ('Z3 (1,5,7) r=(0,4,1)', 3, [1, 5, 7], [0, 4, 1]), ('Z5 (3,4,6,7,9)', 5, [3, 4, 6, 7, 9], [0, 1, 3, 4, 4]),
                         ('3x+1 (d=2)', 2, [1, 3], [0, 1])]:
    print(f"   {name}: mod-2 twisted obstruction for odd offsets: {twisted2(d_, m_, r_)}")
