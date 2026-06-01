#!/usr/bin/env python3
"""
n6_inside_debt_three_channel_s533.py    oracle-2026-06-01-S533

THE n=6 INSIDE DEBT AS A FUNCTION OF THE 3-PAIR JOINT STATE, and the mod-6
three-channel analogue of the n=4 parity law 'a+b+c odd'.

Framework (S529): observer-loneliness measure
   |LONELY(v)| = sum_{ (m_i): sum m_i v_i = 0 } prod_i ghat(m_i),
   ghat(0) = 1 - 2/n,  ghat(m) = -sin(2 pi m / n)/(pi m).
ZEROS: ghat(m)=0 iff n* | m, where n* = n/2 (even n) or n (odd n).
   n=4 -> n*=2: characters on ODD m. m_i v_i = v_i (mod 2). Resonance => sum v_i = 0
                (mod 2). DEBT-FREE iff sum v_i ODD. UNITS mod 2 = {1}: NO sign freedom
                -> a single fixed sum (one channel).
   n=6 -> n*=3: characters on m != 0 (mod 3), i.e. m = +-1 (mod 3). m_i v_i = +- v_i
                (mod 3). UNITS mod 3 = {+1,-1}: a SIGN on each runner. Resonance =>
                sum eps_i v_i = 0 (mod 3) for SOME eps in {+-1}. The sign freedom is
                the new mod-2 channel; mod 6 = (mod 3 residue) (x) (mod 2 sign).

CLAIM (the three-channel parity law, n=6):
   inside debt (genuine multi-channel, order >= 3) VANISHES  <=>
   no sign pattern eps in {+-1}^5 gives  sum eps_i v_i = 0  (mod 3)
   <=>  AT MOST ONE runner has v_i != 0 (mod 3)   (k<=1; for primitive sets k=1).
The 3 'channels' are the residue classes mod 3 (= the 3 antipodal/diameter axes of
the hexagon = the 3 independent pairs). 'a+b+c odd' (n=4) is the degenerate
no-sign-freedom case; n=6 is the first GENUINELY multi-channel parity law.

This script: (A) prints ghat and its mod-6 sign table; (B) computes the order-graded
resonance sums (incl. the inside debt = order>=3) for curated speed sets and
validates the total against direct integration; (C) verifies the law on many
primitive 5-sets: inside-debt-zero <=> 0 not sign-reachable mod 3 <=> k<=1.
"""
from math import sin, pi, gcd
from itertools import product, combinations
from functools import reduce
import random

N = 6
def ghat(m):
    if m == 0: return 1.0 - 2.0/N
    return -sin(2*pi*m/N)/(pi*m)

# ----------------------------------------------------------------------
# (A) the character and its mod-6 sign structure
# ----------------------------------------------------------------------
def part_A():
    print("="*70); print("(A) n=6 character ghat(m): zeros at 3|m; sign by m mod 6"); print("="*70)
    print("   m : ghat(m)      |ghat|*pi   m%6   sign   (zero iff 3|m)")
    for m in range(1, 9):
        g = ghat(m)
        sgn = '0' if abs(g) < 1e-12 else ('+' if g > 0 else '-')
        print(f"   {m} : {g:+.5f}   {abs(g)*pi:.5f}    {m%6}     {sgn}")
    print("   => |ghat(m)| = sqrt(3)/(2 pi m) for m not= 0 mod 3; sign: m%6 in {1,2}-> -,")
    print("      m%6 in {4,5}-> +. The mod-6 pattern = (mod-3 channel)(x)(mod-2 sign).")

# ----------------------------------------------------------------------
# (B) order-graded resonance sums incl. inside debt (order>=3)
# ----------------------------------------------------------------------
def resonance_by_order(speeds, M=6):
    by = {}
    rng = range(-M, M+1)
    for ms in product(rng, repeat=len(speeds)):
        if sum(a*b for a, b in zip(ms, speeds)) != 0:
            continue
        r = sum(1 for a in ms if a != 0)
        t = 1.0
        for a in ms: t *= ghat(a)
        by[r] = by.get(r, 0.0) + t
    return by

def lonely_direct(speeds, G=300000):
    lo, hi = 1.0/N, 1.0 - 1.0/N
    c = 0
    for i in range(G):
        t = (i+0.5)/G
        if all(lo < (s*t) % 1.0 < hi for s in speeds): c += 1
    return c/G

def part_B():
    print("\n"+"="*70); print("(B) order-graded resonance sums; inside debt = sum of order>=3"); print("="*70)
    sets = [
        (1,2,3,4,5),     # AP = regular hexagon (tight)
        (1,2,4,5,3),     # same multiset, perm (sanity)
        (3,6,9,12,1),    # k=1: exactly ONE runner !=0 mod 3 -> predicted DEBT-FREE
        (3,6,9,12,15),   # all =0 mod3 (non-primitive); k=0
        (1,2,3,6,9),     # k=2 active {1,2}
        (1,4,7,2,5),     # residues mod3: 1,1,1,2,2 -> k=5
    ]
    for v in sets:
        by = resonance_by_order(v, M=6)
        inside = sum(val for r, val in by.items() if r >= 3)
        tot = sum(by.values())
        k = sum(1 for s in v if s % 3 != 0)
        direct = lonely_direct(v)
        ordstr = " ".join(f"r{r}:{by[r]:+.4f}" for r in sorted(by))
        res3 = mod3_zero_reachable(v)
        print(f"  v={v} (res mod3={[s%3 for s in v]}, k={k}, 0-reach mod3={res3})")
        print(f"     {ordstr}")
        print(f"     INSIDE DEBT(order>=3) = {inside:+.5f} ; total resonance-sum={tot:.5f} ; direct|LONELY|={direct:.5f}")

# ----------------------------------------------------------------------
# (C) verify the law on many primitive 5-sets
# ----------------------------------------------------------------------
def mod3_zero_reachable(speeds):
    """is there eps in {+-1}^len with sum eps_i v_i = 0 (mod 3)?"""
    for eps in product((1,-1), repeat=len(speeds)):
        if sum(e*s for e, s in zip(eps, speeds)) % 3 == 0:
            return True
    return False

def has_order3_resonance(speeds, M=7):
    """does an order>=3 resonance exist with all nonzero coeffs != 0 mod 3?"""
    rng = [x for x in range(-M, M+1) if x % 3 != 0]
    n = len(speeds)
    # search over 3- and 4- and 5-subsets with small coeffs
    for r in (3, 4, 5):
        for S in combinations(range(n), r):
            for coeffs in product(rng, repeat=r):
                if sum(coeffs[i]*speeds[S[i]] for i in range(r)) == 0:
                    return True
    return False

def part_C():
    print("\n"+"="*70); print("(C) LAW CHECK on primitive 5-sets: debt-zero <=> 0 not sign-reachable mod3 <=> k<=1"); print("="*70)
    rnd = random.Random(533)
    tested = 0; law_holds = 0; debtfree = []
    pool = list(combinations(range(1, 13), 5))
    rnd.shuffle(pool)
    for v in pool:
        if reduce(gcd, v) != 1:
            continue
        tested += 1
        if tested > 200: break
        k = sum(1 for s in v if s % 3 != 0)
        reach = mod3_zero_reachable(v)
        has = has_order3_resonance(v, M=6)
        # law: NO order-3 resonance  <=>  not reach  <=>  k<=1
        law = ((not has) == (not reach)) and ((not reach) == (k <= 1))
        if law: law_holds += 1
        if not has:
            debtfree.append((v, k))
    print(f"  tested {tested} primitive 5-sets; law (no-debt <=> 0-not-reachable-mod3 <=> k<=1) holds: {law_holds}/{tested}")
    print(f"  debt-FREE sets found: {debtfree[:12]}{' ...' if len(debtfree)>12 else ''}")
    print(f"  (debt-free require exactly one runner != 0 mod 3, the rest divisible by 3)")

def part_D_n4_parallel():
    print("\n"+"="*70); print("(D) n=4 parallel: units mod 2 = {1} -> single sum 'a+b+c odd' (no sign freedom)"); print("="*70)
    # n=4 inside debt-free iff a+b+c odd; verify via resonance existence with odd coeffs
    def n4_debtfree(a,b,c,M=7):
        rng=[x for x in range(-M,M+1) if x%2!=0]  # odd m (chars live on odd for n=4)
        for ma,mb,mc in product(rng,repeat=3):
            if ma*a+mb*b+mc*c==0: return False
        return True
    print("   (a,b,c) : sum  parity   debt-free(no odd-coeff resonance)?")
    for v in [(1,2,3),(1,2,4),(1,3,5),(2,3,4),(1,1,1),(1,2,5)]:
        df=n4_debtfree(*v); s=sum(v)
        print(f"   {v} : {s:2d}  {'odd ' if s%2 else 'even'}   debt-free={df}  (predict {bool(s%2)})")
    print("   => n=4: ONE fixed sum mod 2 (units {1}). n=6: SIGNED sum mod 3 (units {+-1})")
    print("      = the first genuinely multi-channel parity law; mod6 = mod3 (x) mod2-sign.")

def main():
    part_A(); part_B(); part_C(); part_D_n4_parallel()

if __name__ == "__main__":
    main()
