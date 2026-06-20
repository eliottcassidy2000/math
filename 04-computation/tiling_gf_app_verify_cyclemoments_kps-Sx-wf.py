#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ADVERSARIAL VERIFICATION of claimed tiling-model cycle-moments (THM-554 application).
kps-Sx-wf. EXACT rationals via fractions.Fraction. No git.

Claims under test:
  E[c3] = (C(n,3)+(n-2))/4 = (n^3-3n^2+8n-12)/24       [PROVED upstream]
  Var[c3] = (n^3-7n^2+20n-16)/32
  E[c5] = (n^5-10n^4+45n^3-140n^2+294n-280)/160
  E[c7] = (n^7-21n^6+189n^5-1015n^4+3836n^3-10514n^2+18458n-15204)/896
  E[c_k] leading term (k-1)!/2^k * C(n,k)
  E[H] OCF brute: n=3..7 -> 2, 4, 79/8, 29, 3175/32
  n=5 identity 1+2E[c3]+2E[c5] = E[H]

Strategy:
 (1) Z-engine for the score census (exact, no 2^F enumeration) -> exact E,Var of c3.
 (2) Independent BRUTE enumeration over all 2^{C(n-1,2)} tilings (n<=7) for c3,c5,c7,H.
 (3) Per-subset linearity recomputation of E[c_k] independent of brute & polynomial.
 (4) Check claimed closed-form polynomials at each n vs the exact measured values.
"""
import sys
from fractions import Fraction as Fr
from collections import defaultdict, Counter
from itertools import product, combinations, permutations
from math import comb, factorial
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

# -------------------------------------------------------------- Z-engine
def beta_step(dist, n):
    nd = defaultdict(int)
    for vec, cnt in dist.items():
        l = list(vec) + [0]; l[n-1] += 1; nd[tuple(l)] += cnt
    dist = nd
    for b in range(1, n-1):
        nd = defaultdict(int)
        for vec, cnt in dist.items():
            l0 = list(vec); l0[n-1] += 1; nd[tuple(l0)] += cnt
            l1 = list(vec); l1[b-1] += 1; nd[tuple(l1)] += cnt
        dist = nd
    return dist

def build_Z(N):
    dist = {(0,): 1}
    for n in range(2, N+1):
        dist = beta_step(dist, n)
    return dist

def c3_moments_from_Z(distZ, n):
    """Exact E[c3], Var[c3] from full score census."""
    tot = 0; s1 = 0; s2 = 0
    for vec, cnt in distZ.items():
        c3 = comb(n,3) - sum(comb(s,2) for s in vec)
        tot += cnt; s1 += cnt*c3; s2 += cnt*c3*c3
    E = Fr(s1, tot); E2 = Fr(s2, tot)
    return E, E2 - E*E, tot

# -------------------------------------------------------------- brute tiling
def tiles(n):
    return [(a, b) for a in range(3, n+1) for b in range(1, a-1)]

def build_adj(n, T, bv):
    adj = [[0]*(n+1) for _ in range(n+1)]
    for k in range(n, 1, -1):
        adj[k][k-1] = 1
    for (a, b), bit in zip(T, bv):
        if bit == 0: adj[a][b] = 1
        else: adj[b][a] = 1
    return adj

def count_dicycles_k(adj, n, k):
    """Count directed k-cycles (as cyclic vertex sequences / distinct cycles)."""
    # count distinct directed cycles of length k. Use combinations of vertex sets
    # then count Hamiltonian directed cycles within each k-subset.
    total = 0
    for S in combinations(range(1, n+1), k):
        # count directed Hamiltonian cycles in tournament induced on S
        # fix first vertex S[0], permute the rest, check directed cycle, divide by nothing
        # (each directed cycle counted once by fixing smallest vertex as start and
        #  taking the orientation; to avoid double count of reverse, we count both
        #  orientations as the cycle is directed -> a directed cycle has ONE orientation)
        s0 = S[0]
        rest = S[1:]
        cnt = 0
        for perm in permutations(rest):
            seq = (s0,) + perm
            ok = all(adj[seq[i]][seq[(i+1) % k]] == 1 for i in range(k))
            if ok:
                cnt += 1
        total += cnt
    return total

def OCF_H(adj, n):
    """H = sum over odd-cycle vertex-subsets ... use OCF = I(Omega(T),2).
    Simplest exact: H = number of Hamiltonian paths = ham(T). We compute via
    independence polynomial of the odd-cycle conflict graph? Use the canonical
    definition H = #Hamiltonian paths in the tournament (Redei). Count directly."""
    # Redei H(T) = number of Hamiltonian directed paths.
    cnt = 0
    verts = list(range(1, n+1))
    for perm in permutations(verts):
        if all(adj[perm[i]][perm[i+1]] == 1 for i in range(n-1)):
            cnt += 1
    return cnt

def brute_moments(n, ks=(3,5,7), do_H=False):
    T = tiles(n); F = len(T)
    accum = {k: [0,0] for k in ks}  # [sum, sumsq]
    Hsum = 0
    tot = 0
    for bv in product((0,1), repeat=F):
        adj = build_adj(n, T, bv)
        tot += 1
        for k in ks:
            if k <= n:
                c = count_dicycles_k(adj, n, k)
                accum[k][0] += c; accum[k][1] += c*c
        if do_H:
            Hsum += OCF_H(adj, n)
    res = {}
    for k in ks:
        if k <= n:
            E = Fr(accum[k][0], tot); E2 = Fr(accum[k][1], tot)
            res[k] = (E, E2 - E*E)
    H = Fr(Hsum, tot) if do_H else None
    return res, H, tot

# -------------------------------------------------------------- per-subset linearity
def per_subset_Eck(n, k):
    """E[c_k] = sum over k-subsets S of P(S is a directed k-cycle).
    For each k-subset, the induced tournament's arcs: base-path arcs are FORCED
    (i>j among consecutive => arc i->j fixed by base path n->..->1; actually base path
    gives arc v->v-1 for all v, and tiles between non-consecutive are free fair coins).
    Within subset S, an arc (i,j) i>j is forced to i->j IFF j=i-1 (consecutive in the
    FULL path, i.e. they are adjacent integers), else it's a fair independent coin.
    Compute P(S induces some directed k-cycle) exactly by enumerating the free arcs."""
    total = Fr(0)
    for S in combinations(range(1, n+1), k):
        Ssort = sorted(S)
        # arcs among S: for each pair i>j, forced if i-j==1 else free coin
        pairs = [(Ssort[b], Ssort[a]) for a in range(k) for b in range(a)]  # (low,high)
        free = []
        forced = {}
        for (lo, hi) in pairs:
            if hi - lo == 1:
                forced[(hi, lo)] = 1  # arc hi->lo forced (base path)
            else:
                free.append((hi, lo))
        nf = len(free)
        good = 0
        for bits in product((0,1), repeat=nf):
            # build orientation
            arc = set(forced.keys())
            for (hi, lo), bit in zip(free, bits):
                if bit == 0: arc.add((hi, lo))
                else: arc.add((lo, hi))
            # adjacency among S
            adjset = arc
            # count directed Hamiltonian cycles on S (length k)
            s0 = Ssort[0]; rest = Ssort[1:]
            has = 0
            for perm in permutations(rest):
                seq = (s0,)+perm
                if all((seq[i], seq[(i+1) % k]) in adjset for i in range(k)):
                    has += 1
            good += has  # E[c_k] counts number of directed k-cycles, so add count
        total += Fr(good, 1<<nf)
    return total

# -------------------------------------------------------------- claimed polys
def poly_Ec3(n): return Fr(n**3 - 3*n**2 + 8*n - 12, 24)
def poly_Varc3(n): return Fr(n**3 - 7*n**2 + 20*n - 16, 32)
def poly_Ec5(n): return Fr(n**5 - 10*n**4 + 45*n**3 - 140*n**2 + 294*n - 280, 160)
def poly_Ec7(n):
    return Fr(n**7 - 21*n**6 + 189*n**5 - 1015*n**4 + 3836*n**3
              - 10514*n**2 + 18458*n - 15204, 896)

def main():
    print("="*78)
    print("PART 1: E[c3], Var[c3] from Z-engine vs claimed closed forms (n=3..10)")
    print("="*78)
    dist = {(0,): 1}
    ok_c3 = True
    for n in range(2, 11):
        dist = beta_step(dist, n)
        if n < 3: continue
        E, V, tot = c3_moments_from_Z(dist, n)
        pE, pV = poly_Ec3(n), poly_Varc3(n)
        mE = (E == pE); mV = (V == pV)
        ok_c3 &= mE and mV
        print(f" n={n:2d} tot=2^{comb(n-1,2):<2d}  E[c3]={str(E):>8s} poly={str(pE):>8s} {mE}"
              f"   Var={str(V):>8s} poly={str(pV):>8s} {mV}")
    print(f"  --> c3 closed forms match Z-engine for all n<=10: {ok_c3}")

    print("\n" + "="*78)
    print("PART 2: BRUTE tiling enumeration of c3,c5,c7,H (n<=7) vs Z-engine & polys")
    print("="*78)
    ok_brute = True
    for n in range(3, 8):
        ks = tuple(k for k in (3,5,7) if k <= n)
        res, H, tot = brute_moments(n, ks=ks, do_H=True)
        # c3 from Z
        Ez, Vz, _ = c3_moments_from_Z(build_Z(n), n)
        line = f" n={n}: tot={tot}"
        print(line)
        # c3
        E3b, V3b = res[3]
        c3match = (E3b == Ez == poly_Ec3(n)) and (V3b == Vz == poly_Varc3(n))
        ok_brute &= c3match
        print(f"   c3: brute E={E3b} V={V3b} | Z E={Ez} V={Vz} | poly E={poly_Ec3(n)} V={poly_Varc3(n)}  match={c3match}")
        # c5
        if 5 in res:
            E5b, V5b = res[5]
            m5 = (E5b == poly_Ec5(n))
            ok_brute &= m5
            print(f"   c5: brute E={E5b} | poly E={poly_Ec5(n)}  match={m5}  (brute Var={V5b})")
        # c7
        if 7 in res:
            E7b, V7b = res[7]
            m7 = (E7b == poly_Ec7(n))
            ok_brute &= m7
            print(f"   c7: brute E={E7b} | poly E={poly_Ec7(n)}  match={m7}")
        # H
        claimed_H = {3:Fr(2),4:Fr(4),5:Fr(79,8),6:Fr(29),7:Fr(3175,32)}
        mh = (H == claimed_H[n])
        print(f"   E[H](Redei ham paths): brute={H} | claimed={claimed_H[n]}  match={mh}")
        # n=5 identity
        if n == 5:
            ident = Fr(1) + 2*poly_Ec3(5) + 2*poly_Ec5(5)
            print(f"   n=5 identity 1+2E[c3]+2E[c5] = {ident}  vs E[H]={H}  match={ident==H}")
        ok_brute &= mh
    print(f"  --> brute matches all: {ok_brute}")

    print("\n" + "="*78)
    print("PART 3: independent PER-SUBSET linearity E[c_k] vs polys (no brute, no Z)")
    print("="*78)
    ok_ps = True
    for n in range(3, 9):
        for k in (3,5,7):
            if k > n: continue
            ps = per_subset_Eck(n, k)
            poly = {3:poly_Ec3, 5:poly_Ec5, 7:poly_Ec7}[k](n)
            m = (ps == poly)
            ok_ps &= m
            print(f" n={n} k={k}: per-subset E[c{k}]={str(ps):>10s} poly={str(poly):>10s} match={m}")
    print(f"  --> per-subset matches polys: {ok_ps}")

    print("\n" + "="*78)
    print("PART 4: leading-term check E[c_k] ~ (k-1)!/2^k * C(n,k)")
    print("="*78)
    for k in (3,5,7):
        lead_coeff = Fr(factorial(k-1), 2**k)  # coeff of C(n,k)
        print(f" k={k}: per-subset baseline (k-1)!/2^k = {lead_coeff} "
              f"(leading n^k coeff = {lead_coeff}/{factorial(k)} = {lead_coeff/factorial(k)})")

    print("\nDONE.")

if __name__ == "__main__":
    main()
