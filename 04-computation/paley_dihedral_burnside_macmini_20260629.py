#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Paley tournament <-> prime: the DIHEDRAL BURNSIDE identity for Hamiltonian paths.
(mac-mini-2026-06-29-S9)

T_p (p=3 mod4) has full symmetry D_{2p} = <r: v->v+1 (automorphism, order p),
s: v->-v (anti-automorphism, involution)>.  It acts on the set of directed
Hamiltonian PATHS: rotations act directly r^k.P = (v_1+k,...,v_n+k); reflections
g (anti-automorphisms) act by reversal-conjugation rho_g(P) = g(reverse(P))
(THM-582).  Burnside:
   #orbits = (1/2p) ( Fix(id) + sum_{k!=0} Fix(r^k) + sum_reflections Fix(g) ).
Fix(id)=H, Fix(r^k)=0 (free rotation: no path is rotation-invariant), Fix(g)=f_g
(the g-palindromic paths).  We VERIFY: all p reflections have the SAME f_g = f,
so   #orbits = (H + p f)/(2p) in Z,   whence H = p*(odd) and H/p == f (mod 2).
"""
from __future__ import annotations
import functools, itertools
print = functools.partial(print, flush=True)


def paley(p):
    qr = set((x * x) % p for x in range(1, p))
    return [[(i != j and ((j - i) % p) in qr) for j in range(p)] for i in range(p)]


def ham_paths(arc, n):
    out = []
    for perm in itertools.permutations(range(n)):
        if all(arc[perm[k]][perm[k + 1]] for k in range(n - 1)):
            out.append(perm)
    return out


def main():
    print("=" * 78)
    print("Paley dihedral Burnside identity for Hamiltonian paths (mac-mini-S9)")
    print("=" * 78)
    for p in (3, 7, 11):
        arc = paley(p)
        HP = ham_paths(arc, p)
        HPset = set(HP)
        H = len(HP)

        # rotations r^k: v -> v+k ; Fix = paths P with (v_1+k,...)=P  -> only k=0
        rot_fix = []
        for k in range(p):
            cnt = sum(1 for P in HP if tuple((v + k) % p for v in P) in HPset and
                      tuple((v + k) % p for v in P) == P)
            rot_fix.append(cnt)

        # reflections g = (v -> -v - c) for c in Z_p (the p involutory anti-automorphisms).
        # check each is an involutory anti-automorphism, then Fix(rho_g) = g-palindromic count.
        refl_fix = []
        for c in range(p):
            g = [(-v - c) % p for v in range(p)]
            # involution?  g[g[v]] = -(-v-c)-c = v  -> always involution
            # anti-automorphism?  arc(u,w) <=> arc(g(w),g(u))
            anti = all(arc[u][w] == arc[g[w]][g[u]] for u in range(p) for w in range(p) if u != w)
            if not anti:
                refl_fix.append(None); continue
            # rho_g(P) = g(reverse(P))
            cnt = 0
            for P in HP:
                if tuple(g[v] for v in reversed(P)) == P:
                    cnt += 1
            refl_fix.append(cnt)

        f_values = [x for x in refl_fix if x is not None]
        all_same = len(set(f_values)) == 1 and len(f_values) == p
        f = f_values[0] if f_values else None
        total_fix = sum(rot_fix) + sum(x for x in refl_fix if x)
        orbits = total_fix / (2 * p)
        print(f"\n--- Paley T_{p}: H={H} ---")
        print(f"  rotation fixed counts (k=0..p-1): {rot_fix}  (only k=0 nonzero: {sum(rot_fix[1:])==0})")
        print(f"  reflection fixed counts: {refl_fix}")
        print(f"  all p reflections give the SAME f? {all_same}  (f={f})")
        print(f"  Burnside #orbits = (H + sum_refl_fix)/(2p) = ({H} + {sum(f_values)})/{2*p} "
              f"= {orbits}  (integer: {orbits == int(orbits)})")
        if f:
            print(f"  => (H + p f)/(2p) = ({H} + {p*f})/{2*p} = {(H+p*f)/(2*p)}")
            print(f"  => H/p = {H//p} (odd: {(H//p)%2==1}),  H/p == f mod2: {(H//p)%2 == f%2}")

    print("\n" + "=" * 78)
    print("RESULT: #orbits of Ham paths under D_{2p} = (H + p f)/(2p) in Z, with the")
    print("free rotation giving p|H and the p reflections each fixing f paths (THM-582).")
    print("Hence H(T_p) = p * (odd), and H(T_p)/p == f (mod 2) -- a Paley refinement of Redei.")
    print("=" * 78)


if __name__ == "__main__":
    main()
