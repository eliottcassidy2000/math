#!/usr/bin/env python3
"""monad-explorer-2026-06-06-S708: STRUCTURAL REDUCTION signed-LRC -> unsigned LRC.

Core thesis (rigorous, gauge-invariant reframing of S699 T3 / HYP-2262):
  "carries a signed zero pair-clock (shell-partner v_i+v_j = 2n-1)"
   <=>  "S mod C=2n-1 is NOT a shell-transversal" (some shell {a,-a} occupied twice),
  a purely UNSIGNED arithmetic property of S, invariant under the full sign gauge group
  AND under the THM-407 group G=<2,-1>.

Findings to verify:
 (A) Exact M(S) via arrangement-vertex method; reproduce floor M=1/n.
 (B) n=8 ALREADY carries a shell-partner tight config (1,2,3,4,5,7,12)=AP_8 with 6->12,
     shell-partner (3,12), 3+12=15=2*8-1.  => the worry-set split is at n=8, NOT n=14.
 (C) Family II "double the (n-2) speed": AP_n with (n-2)->2(n-2). Characterize EXACTLY
     when it is tight (M=1/n) and a shell-partner config, across BOTH parities of n.
 (D) Gauge invariance of M over the sign-orbit; transversality is the gauge invariant.
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations

def ndist(x):
    x %= 1
    return min(x, 1 - x)

def M_exact(V):
    """Exact M(S)=max_t min_i ||v_i t||, V a list/tuple of positive ints.
    Candidate optimal t are arrangement vertices: peaks (2k+1)/(2 v_i) and
    edge-crossings k/(v_i+v_j), k/|v_i-v_j|.  Evaluate min-gap exactly, take max."""
    V = list(V)
    dens = set()
    for v in V:
        dens.add(2 * v)
    for i in range(len(V)):
        for j in range(i + 1, len(V)):
            dens.add(V[i] + V[j])
            d = abs(V[i] - V[j])
            if d: dens.add(d)
    best = F(0)
    seen = set()
    for q in dens:
        for p in range(1, q):  # t=p/q in (0,1)
            t = F(p, q)
            if t > F(1, 2):     # symmetry t<->1-t
                continue
            key = (t.numerator, t.denominator)
            if key in seen:
                continue
            seen.add(key)
            g = min(ndist(v * t) for v in V)
            if g > best:
                best = g
    return best

def shells(V, C):
    """multiset of shells min(r,C-r) for residues r=v mod C; report partners."""
    occ = {}
    partners = []
    for v in V:
        r = v % C
        s = min(r, C - r)
        occ.setdefault(s, []).append(v)
    for s, vs in occ.items():
        if len(vs) >= 2:
            partners.append((s, sorted(vs)))
    return occ, partners

def is_transversal(V, n):
    C = 2 * n - 1
    sset = set()
    for v in V:
        r = v % C
        if r == 0:
            return False
        sset.add(min(r, C - r))
    return len(sset) == len(V)

def report(label, V, n):
    C = 2 * n - 1
    m = M_exact(V)
    occ, partners = shells(V, C)
    tv = is_transversal(V, n)
    tight = (m == F(1, n))
    missing = [a for a in range(1, n) if a not in occ]
    print(f"  {label}: V={tuple(V)}")
    print(f"      n={n} C={C}  M={m}  (1/n={F(1,n)})  tight={tight}  transversal={tv}")
    print(f"      shell-partners (sum=C): {partners}   missing shells: {missing}")
    return tight, partners, tv

def main():
    print("="*70)
    print("(A)+(B): the n=8 shell-partner tight config vs the AP")
    print("="*70)
    report("AP_8        ", tuple(range(1, 8)), 8)
    report("AP_8 6->12  ", (1, 2, 3, 4, 5, 7, 12), 8)       # Family II at n=8
    report("AP_8 sporad2", (1, 4, 5, 6, 7, 11, 13), 8)

    print()
    print("="*70)
    print("(C): Family II  AP_n with (n-2)->2(n-2), both parities, n=5..29")
    print("     (2(n-2) mod (2n-1) = -3, so shell-partner with 3 whenever 3<=n-2)")
    print("="*70)
    fam = []
    for n in range(5, 30):
        a = n - 2
        a2 = 2 * a
        V = sorted([x for x in range(1, n) if x != a] + [a2])
        C = 2 * n - 1
        m = M_exact(V)
        tight = (m == F(1, n))
        _, partners = shells(V, C)
        sp = bool(partners)
        cond3 = (a2 % C == (C - 3) % C)  # is 2(n-2) == -3 mod C ?
        div3 = (C % 3 == 0)
        flag = "TIGHT+shellpartner" if (tight and sp) else ("tight" if tight else "")
        print(f"  n={n:2d} C={C:2d} (n-2={a}->{a2}) M={str(m):8s} tight={int(tight)} "
              f"shellpartner={int(sp)} 3|C={int(div3)} parity={'even' if n%2==0 else 'odd '} {flag}")
        if tight and sp:
            fam.append(n)
    print(f"\n  Family II tight+shell-partner at n = {fam}")
    print(f"  (n mod 6 for these: {[x%6 for x in fam]})")

    print()
    print("="*70)
    print("(D): gauge invariance of M over the SIGN ORBIT (T1) -- transversality invariant")
    print("="*70)
    base = (1, 2, 3, 4, 5, 7, 12)
    n = 8
    Mbase = M_exact(base)
    print(f"  base {base}  M={Mbase}")
    import itertools
    alleq = True
    multisets = set()
    for signs in itertools.product([1, -1], repeat=len(base)):
        signed = [s * v for s, v in zip(signs, base)]
        # M of signed (negative speeds): ||v t|| sign-blind, so use abs for M
        m = M_exact([abs(x) for x in signed])
        if m != Mbase:
            alleq = False
        # shell multiset (residue shells) -- should be invariant
        C = 2 * n - 1
        ms = tuple(sorted(min(x % C, (C - x) % C) for x in signed))
        multisets.add(ms)
    print(f"  M constant over all 2^{len(base)} sign patterns: {alleq}")
    print(f"  distinct shell-multisets over sign orbit: {len(multisets)} (should be 1)")
    print(f"  -> shell-multiset (hence 'has shell-partner') is a GAUGE INVARIANT")

if __name__ == '__main__':
    main()
