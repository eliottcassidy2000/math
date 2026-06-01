#!/usr/bin/env python3
"""
lrc_nowhere_zero_flow_s537.py    oracle-2026-06-01-S537o

Attacking LRC from the NOWHERE-ZERO FLOW mindset.

THE DICTIONARY. The covering inside debt (S529/S533) is a sum over resonances
sum_i m_i v_i = 0. A FULL-SUPPORT resonance (all m_i != 0) is exactly a
NOWHERE-ZERO FLOW on the speed-weighted dipole graph G_v (2 vertices, edge i
carrying weight v_i): conservation = sum_i v_i m_i = 0, nowhere-zero = m_i != 0.
Mod n* (n*=n/2 even, n odd; the character modulus, S533/S534):
   full-support inside debt is FEASIBLE  <=>  G_v has a nowhere-zero Z_{n*}-flow
   <=>  exists m in (Z_{n*}\0)^{n-1} with sum m_i v_i = 0 (mod n*).
So the THREE-CHANNEL PARITY LAW (S533) and its n=18 VACUITY (S534) are
nowhere-zero-flow existence statements, where Tutte/Seymour theory applies.

KEY STRUCTURES we verify:
 (1) FACTORIZATION: a runner with v_i = 0 (mod n*) is an INERT/FREE edge (any of the
     n*-1 nonzero values works, contributes 0 to the sum) -> factor (n*-1). An ACTIVE
     runner (v_i coprime-nonzero mod n*) is CONSTRAINED. So
        NZ-flow count = (n*-1)^{#inert} * C_active.
 (2) BRIDGE CHARACTERIZATION: debt-free (no NZ Z_{n*}-flow) <=> the active edges form
     a 'Z_{n*}-bridge' = EXACTLY ONE active runner (the rest inert). [matches S533 k=1]
     With >=2 active edges the (bridgeless) dipole always carries a NZ flow -- the
     flow-theoretic reason debt is present, echoing Seymour's 6-flow theorem.
 (3) FLOW-POLYNOMIAL form: the unweighted m-dipole flow polynomial is
        F(D_m;k) = ((k-1)^m + (-1)^m (k-1))/k ;
     the inside-debt VALUE is its ghat-weighted analogue
        debt ~ sum over NZ integer flows of prod ghat(m_i)  (a flow enumerator).
 (4) TENSION DUAL: runner positions x_i=frac(v_i t) are a TENSION (circular coloring)
     on the observer-star; loneliness = a circular-n-coloring (all star-edges at
     circular distance >= 1 on circumference n). Flow<->tension(coloring) duality.

We verify (1)-(3) against the S533/S534 inside-debt facts, and sanity-check (4).
"""
from itertools import product
from functools import reduce
from math import gcd, sin, pi
import random

def nstar(n): return n//2 if n % 2 == 0 else n

# ---------- nowhere-zero Z_k flow count on the speed-weighted dipole ----------
def nz_flow_count(weights, k):
    """#{m in (Z_k\0)^E : sum m_i w_i = 0 mod k} -- NZ Z_k-flows of G_v."""
    # DP over residues mod k
    from collections import defaultdict
    dist = defaultdict(int); dist[0] = 1
    for w in weights:
        nd = defaultdict(int)
        for r, cnt in dist.items():
            for m in range(1, k):              # nonzero flow value
                nd[(r + m*w) % k] += cnt
        dist = nd
    return dist[0]

def n_active(weights, k):
    return sum(1 for w in weights if w % k != 0)

# ---------- inside-debt enumerator (ghat-weighted NZ integer flows) ----------
def ghat(m, n):
    if m == 0: return 1.0 - 2.0/n
    return -sin(2*pi*m/n)/(pi*m)

def inside_debt_flow_enum(weights, n, M=6):
    """sum over NZ integer flows (all m_i != 0, sum m_i w_i = 0) of prod ghat(m_i),
    truncated |m_i|<=M. (the full-support inside debt as a flow enumerator)"""
    tot = 0.0
    rng = [m for m in range(-M, M+1) if m != 0]
    for ms in product(rng, repeat=len(weights)):
        if sum(a*b for a, b in zip(ms, weights)) == 0:
            t = 1.0
            for a in ms: t *= ghat(a, n)
            tot += t
    return tot

def main():
    print("="*74)
    print("LRC via NOWHERE-ZERO FLOWS: full-support resonance = NZ flow on speed dipole")
    print("="*74)

    print("\n(1)+(2) FACTORIZATION & BRIDGE: NZ Z_{n*}-flow count, inert factor, active count")
    print("   n=6 (n*=3) and n=18 (n*=9). debt-free <=> NZ-flow-count=0 <=> exactly 1 active.")
    tests6 = [(1,2,3,4,5),(3,6,9,12,1),(3,6,9,15,2),(1,2,3,6,9),(6,1,12,18,24)]
    for v in tests6:
        k = 3; cnt = nz_flow_count(v, k); a = n_active(v, k)
        free = (k-1)**(len(v)-a)
        cact = cnt // free if free else 0
        print(f"   n=6 v={v}: active={a}, inert factor (k-1)^inert={free}, "
              f"NZ-Z3-flow count={cnt} (active part {cact}) -> {'DEBT-FREE' if cnt==0 else 'debt present'}")
    print()
    tests18 = [tuple(range(1,18)), (3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,1)]
    for v in tests18:
        k = 9; cnt = nz_flow_count(v, k); a = n_active(v, k)
        tag = "DEBT-FREE (1 active = Z9-bridge)" if cnt == 0 else "debt present"
        print(f"   n=18 v={'AP 1..17' if v[0]==1 and v[-1]==17 else v}: active(coprime+v3=1 nonzero)={a}, "
              f"NZ-Z9-flow count={cnt} -> {tag}")
    print()

    print("(2') BRIDGE LAW verified on random primitive sets (debt-free <=> 1 active):")
    for n in (6, 8, 18):
        k = nstar(n); rnd = random.Random(900+n); tot = 0; law = 0; df = 0
        for _ in range(6000):
            v = tuple(sorted(rnd.sample(range(1, 7*n), n-1)))
            if reduce(gcd, v) != 1: continue
            tot += 1
            if tot > 400: break
            cnt = nz_flow_count(v, k); a = n_active(v, k)
            free = (cnt == 0)
            if free == (a == 1): law += 1
            if free: df += 1
        print(f"   n={n} (n*={k}): bridge law (NZ-flow=0 <=> 1 active) holds {law}/{tot}; debt-free count={df}")
    print()

    print("(3) FLOW POLYNOMIAL: unweighted m-dipole F(D_m;k)=((k-1)^m+(-1)^m(k-1))/k")
    for m in (3,4,5):
        for k in (3,4,5):
            # all weights coprime to k -> 'unweighted-like'? use weights all =1
            cnt = nz_flow_count(tuple([1]*m), k)
            formula = ((k-1)**m + (-1)**m * (k-1))//k
            print(f"   D_{m}, k={k}: NZ-flow count(weights all 1)={cnt}, formula={formula}, "
                  f"match={cnt==formula}")
    print()

    print("(3') INSIDE-DEBT as ghat-weighted NZ-integer-flow enumerator (n=6):")
    for v in [(1,2,3,4,5),(3,6,9,12,1),(1,2,3,6,9)]:
        d = inside_debt_flow_enum(v, 6, M=6)
        cnt = nz_flow_count(v, 3)
        print(f"   v={v}: full-support inside debt (flow enumerator) = {d:+.5f}; "
              f"NZ-Z3-flow count={cnt} ({'0 => debt exactly 0' if cnt==0 else 'nonzero'})")
    print()

    print("(4) TENSION DUAL sanity: runner positions = circular coloring (circumference n)")
    print("    loneliness <=> observer at circular distance >= 1 from every runner (scaled).")
    n = 7; v = (1,2,3,4,6)   # non-AP primitive set: genuine open lonely window
    found = None
    for s in range(40000):
        t = (s+0.5)/40000
        y = [ (n*((vi*t) % 1.0)) for vi in v ]   # positions on circumference-n circle
        if all(min(yi % n, n - (yi % n)) >= 1.0 for yi in y):
            found = (t, [round(yi,2) for yi in y]); break
    if found:
        print(f"    n={n} v={v}: circular-{n}-coloring with observer far found at t={found[0]:.4f}: "
              f"positions {found[1]} (all >= 1 from 0). LRC = orbit-constrained circular coloring.")
    else:
        print(f"    n={n} v={v}: no strict circular coloring on the grid (tight set -- "
              f"colorable only at measure-zero t=k/n, the boundary extremal).")

    print("\n" + "="*74)
    print("SYNTHESIS: inside debt / parity = NZ Z_{n*}-flow on the speed dipole; debt-free")
    print("= a Z_{n*}-bridge (1 active runner). Tutte/Seymour: bridgeless => NZ flow exists")
    print("=> debt present once >=2 active -- the flow reason large n* is parity-vacuous.")
    print("Runner positions are the dual TENSION (circular coloring). " )
    print("="*74)

if __name__ == "__main__":
    main()
