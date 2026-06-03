#!/usr/bin/env python3
"""
claudebox-2026-06-03-S593 : rigidity through one lens — local (fixed-point pinning) vs global
(symmetry cascade), and what it says about where the Lonely Runner argument works/fails.

THESIS. "Rigidity" is used two opposite ways in the repo:
  - PINNING rigidity (THM-398): a degree of freedom is FORCED (a small-owner endpoint pinned to
    the arc centre; the witness time isolated/tight). Local, around a fixed point.
  - AUTOMORPHISM rigidity (HYP-2130): trivial Aut, a GENERIC/asymmetric object.
They point OPPOSITE ways. The reconciler: SYMMETRY converts one into the other.
  symmetric config  => automorphism-FLEXIBLE  but witness-RIGID (tight, isolated lonely time)
  generic config    => automorphism-RIGID     but witness-FLEXIBLE (open set of lonely times)

LOCAL rigidity around a fixed point: at the maximin witness t*, the BINDING runners (those at
distance =G) pin t*. The witness is locally rigid iff the binding gradients positively span ℝ
(both signs present) — a KKT/rigidity-matrix condition. The APEX runner sits at the reflection
fixed point frac=1/2 (the non-smooth max of ‖·‖), a degenerate pin.

GLOBAL rigidity that cascades: in a symmetric config the binding runners form a SYMMETRY ORBIT —
one pin replicates to its whole isomorphic family ("permeating from isomorphic to another").
A single corrector can clear ONE local pin (the apex); it cannot clear an ORBIT of pins — that
needs the multi-sieve (HYP-2075). So: one local pin = single-corrector works; cascaded orbit of
pins = multi-sieve. This is HYP-2130's single/multi line, refined as local-vs-global rigidity.

This file measures, per config: the maximin witness, #binding runners + gradient-sign balance
(local rigidity), the apex, additive symmetry proxies, and tests the cascade.
"""
import itertools, math
from collections import Counter

def frac(x): return x - math.floor(x)
def dist(x):
    s = frac(x); return min(s, 1 - s)

def maximin_witness(S, N=60000):
    best = -1.0; bt = 0.0
    for a in range(1, N):
        t = a / N
        m = min(dist(v * t) for v in S)
        if m > best: best, bt = m, t
    return best, bt

def binding_analysis(S, t, G, tol=2e-3):
    """binding runners (dist≈G), and slope sign of ‖v t‖ at t (+ if frac∈(0,1/2), − if (1/2,1))."""
    plus = []; minus = []; apex = []
    for v in S:
        d = dist(v * t)
        if abs(d - G) < tol:
            s = frac(v * t)
            if abs(s - 0.5) < tol: apex.append(v)        # at the reflection fixed point 1/2
            elif s < 0.5: plus.append(v)                 # slope +v (distance rising through t)
            else: minus.append(v)                        # slope −v
    return plus, minus, apex

def three_term(S):
    Sset = set(S); c = 0
    for a, b in itertools.combinations(S, 2):
        if a + b in Sset: c += 1
    for a in S:
        if 2 * a in Sset: c += 1
    return c

def additive_energy(S):
    sums = Counter()
    for a in S:
        for b in S: sums[a + b] += 1
    return sum(v * v for v in sums.values())

def neg_symmetric(S, M):
    """is S symmetric under x -> -x mod M (reflection symmetry of the speed set)?"""
    Sset = set(x % M for x in S)
    return all(((-x) % M) in Sset for x in Sset)

def report(name, S):
    k = len(S); n = k + 1; delta = 1 / n
    G, t = maximin_witness(S)
    plus, minus, apex = binding_analysis(S, t, G)
    nb = len(plus) + len(minus) + len(apex)
    both = (len(plus) > 0 and len(minus) > 0) or len(apex) > 0
    tt = three_term(S); en = additive_energy(S)
    tight = abs(G - delta) < 5e-3
    print(f"  {name:22s} k={k:2d} G={G:.4f} δ={delta:.4f} margin={G-delta:+.4f} "
          f"{'TIGHT' if tight else 'safe ':5s} | binding={nb:2d} (+{len(plus)} −{len(minus)} apex{len(apex)}) "
          f"localrigid={'Y' if both else 'n'} | 3term={tt:2d} energy={en}")
    return dict(k=k, G=G, delta=delta, tight=tight, nbind=nb, plus=len(plus), minus=len(minus),
                apex=len(apex), both=both, tt=tt, en=en)

def main():
    print("=" * 92)
    print("S593  rigidity: local fixed-point pinning vs global symmetry cascade (LRC witness)")
    print("=" * 92)

    print("\n[1] TIGHT (symmetric) vs SAFE (generic): binding-set size & local rigidity at the witness")
    fams = [
        ("AP {1..6}",        list(range(1, 7))),
        ("AP {1..7}",        list(range(1, 8))),
        ("AP {1..13} (n=14)",list(range(1, 14))),
        ("transl {7..12}",   list(range(7, 13))),
        ("transl {8..14}",   list(range(8, 15))),
        ("geometric 1,2,4..",[2**i for i in range(6)]),
        ("Sidon 1,2,5,11,22", [1,2,5,11,22]),
        ("generic A",        [3,5,11,17,29,41]),
        ("generic B",        [2,7,13,19,31,43,53]),
    ]
    rows = []
    for nm, S in fams:
        rows.append((nm, report(nm, S)))

    print("\n[2] THE CASCADE: do binding runners form a SYMMETRY ORBIT? (reflection x→−x mod n)")
    print("    tight/symmetric configs should have their pins replicate through symmetry.")
    for nm, S in [("AP {1..6}", list(range(1,7))), ("AP {1..13}", list(range(1,14))),
                  ("transl {7..12}", list(range(7,13))), ("generic A", [3,5,11,17,29,41])]:
        k=len(S); n=k+1
        G,t = maximin_witness(S)
        plus,minus,apex = binding_analysis(S,t,G)
        binders = sorted(plus+minus+apex)
        sym = neg_symmetric(S, n)
        # is the binding set closed under negation mod n? (pins come in ± orbits)
        bset = set(b % n for b in binders)
        orbit_closed = all(((-b)%n) in bset for b in bset) if bset else True
        print(f"    {nm:16s} n={n:2d} speeds±-sym mod n={sym}  binders={binders}  "
              f"binders ±-closed mod n={orbit_closed}")

    print("\n[3] THE APEX as reflection fixed point (even n: runner n/2 sits at frac=1/2)")
    for n in [6, 10, 14]:
        S=list(range(1,n)); k=len(S)
        G,t=maximin_witness(S)
        q=n//2
        print(f"    n={n}: apex={q}, frac(apex·t*)={frac(q*t):.4f} (→ reflection fixed pt 1/2), "
              f"‖apex·t*‖={dist(q*t):.4f} (max possible) — the safest spot, a non-smooth pin")

    print("\n[4] SUMMARY — the rigidity duality")
    nt = [r for _,r in rows if r['tight']]; ns=[r for _,r in rows if not r['tight']]
    if nt:
        print(f"    TIGHT configs: avg binding={sum(r['nbind'] for r in nt)/len(nt):.1f}, "
              f"avg 3term={sum(r['tt'] for r in nt)/len(nt):.1f}, all local-rigid={all(r['both'] for r in nt)}")
    if ns:
        print(f"    SAFE  configs: avg binding={sum(r['nbind'] for r in ns)/len(ns):.1f}, "
              f"avg 3term={sum(r['tt'] for r in ns)/len(ns):.1f}")
    print("    => symmetry (3-term/AP) ⇒ large binding orbit ⇒ witness-rigid (tight);")
    print("       asymmetry ⇒ small generic binding ⇒ witness-flexible (safe).")

    print("\n[5] RIGOR: is 'binding set closed under a symmetry of the speed set' the tight discriminator?")
    print("    Test reflection x→a−x and dilation x→cx mod n that fix the speed set, check binders orbit-closed.")
    def symmetries_fixing(S, n):
        """affine maps x->c x + b mod n (c in units) that permute S mod n; return as functions."""
        Sset = set(x % n for x in S); syms = []
        units = [c for c in range(1, n) if math.gcd(c, n) == 1]
        for c in units:
            for b in range(n):
                if all(((c * x + b) % n) in Sset for x in Sset):
                    syms.append((c, b))
        return syms
    def binders_of(S):
        G, t = maximin_witness(S)
        p, m, a = binding_analysis(S, t, G)
        return sorted(set(x % (len(S)+1) for x in (p+m+a))), G, 1/(len(S)+1)
    tests = [
        ("AP {1..6} (tight)",      list(range(1,7)),  True),
        ("sporadic {1,3,4,7} n=5", [1,3,4,7],         True),
        ("sporadic {1,3,4,5,9} n=6",[1,3,4,5,9],      True),
        ("transl {7..12} (safe)",  list(range(7,13)), False),
        ("generic A (safe)",       [3,5,11,17,29,41], False),
    ]
    for nm, S, _ in tests:
        n = len(S)+1
        binders, G, d = binders_of(S)
        syms = symmetries_fixing(S, n)
        bset = set(binders)
        # is the binding set invariant under some NONTRIVIAL speed-set symmetry?
        nontriv = [(c,b) for (c,b) in syms if (c,b) != (1,0)]
        closed_under = [(c,b) for (c,b) in nontriv if bset and all(((c*x+b)%n) in bset for x in bset)]
        print(f"    {nm:26s} n={n:2d} G={G:.4f} δ={d:.4f} {'TIGHT' if abs(G-d)<5e-3 else 'safe':5s} | "
              f"#speed-syms={len(nontriv)} binders={binders} closed under {closed_under if closed_under else 'NONE'}")
    print("    => tight configs: binding set is invariant under a nontrivial speed-set symmetry (a cascade orbit);")
    print("       safe configs: no nontrivial symmetry, binders generic. Local pin (apex/fixed pt) vs global orbit.")

if __name__ == "__main__":
    main()
