"""
mac-mini-2026-07-07-S44 (HYP-4917, part 2) -- EXACT test of kps-S62's palindromic-
extremizer conjecture, in its sharpest decidable form:

    Fix a step MULTISET; over all arrangements (words) of it, which word minimizes
    mu_{1/7}?  Conjecture predicts: a palindrome.

Uses opus-S136's exact order-cell engine (mu_exact) -- kps's own descents were noisy
(reported honestly in S62); exact rationals settle each multiset completely.
Also reports WHICH palindrome wins (valley / mountain / block shape) to seed the
symmetrization-inequality statement, and the exact spread of mu over arrangements
(how much the ARRANGEMENT matters at fixed multiset = the pure 'word' effect,
separating kps's step-alphabet stratification from word order).

NOTE (proved en route, part-1 script): U(E, 1-x) = U(E, x) identically for every E
(e(1-x) = -ex mod 1 + THM-639-A pointwise), so mu/E[U] computations may restrict to
[0, 1/2]; and 'reversal-symmetrized moments' are vacuously plain moments.
"""
import sys, time
sys.path.insert(0, '04-computation')
from fractions import Fraction as F
from itertools import permutations
from math import comb

from lrc_exact_mu_ordercells_opus_S136 import mu_exact  # credit opus-S136

THETA = F(1,7)

def word_to_E(word):
    E=[0]
    for s in word: E.append(E[-1]+s)
    return E

def distinct_words(multiset):
    """all distinct arrangements, modulo reversal (canonical = min(word, reversed))."""
    seen=set()
    for p in permutations(multiset):
        c = min(p, p[::-1])
        seen.add(c)
    return sorted(seen)

def is_palin(w): return tuple(w)==tuple(w[::-1])

CASES = [
    ("k=13 {2^8,1^4} (record's multiset)", (2,)*8+(1,)*4, None),
    ("k=13 {1^11,2}",                     (1,)*11+(2,),   None),
    ("k=13 {1^10,2^2}",                   (1,)*10+(2,)*2, None),
    ("k=13 {1^9,2^3}",                    (1,)*9+(2,)*3,  None),
    ("k=8  {1^6,2} (binding leg)",        (1,)*6+(2,),    None),
    ("k=8  {1^5,2^2}",                    (1,)*5+(2,)*2,  None),
    ("k=8  {1^4,2^3}",                    (1,)*4+(2,)*3,  None),
    ("k=8  {3,2,1^5}",                    (3,2)+(1,)*5,   None),
]

for label, ms, _ in CASES:
    t0=time.time()
    words = distinct_words(ms)
    vals = []
    for w in words:
        E = word_to_E(list(w))
        vals.append((mu_exact(E), w))
    vals.sort()
    lo_mu, lo_w = vals[0]
    hi_mu, hi_w = vals[-1]
    palin_vals = [(m,w) for m,w in vals if is_palin(w)]
    lo_p_mu, lo_p_w = palin_vals[0] if palin_vals else (None,None)
    winner_palin = is_palin(lo_w)
    # ties at the min?
    ties = [w for m,w in vals if m == lo_mu]
    print(f"== {label}: {len(words)} words (mod reversal), {time.time()-t0:.1f}s")
    print(f"   min mu = {lo_mu} = {float(lo_mu):.6f} at word {lo_w}  "
          f"PALINDROME? {winner_palin}  (#tied minimizers: {len(ties)})")
    if not winner_palin and lo_p_mu is not None:
        print(f"   best PALINDROME: {lo_p_mu} = {float(lo_p_mu):.6f} at {lo_p_w} "
              f"(excess over min: {float(lo_p_mu-lo_mu):+.6f})")
    print(f"   max mu = {float(hi_mu):.6f} at {hi_w}; arrangement spread = {float(hi_mu-lo_mu):.6f}")
    if len(ties)>1:
        for w in ties: print(f"      tie: {w} palin={is_palin(w)}")
    sys.stdout.flush()
