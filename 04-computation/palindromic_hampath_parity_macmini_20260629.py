#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""The PALINDROMIC HAMILTONIAN PATH parity / odd index (mac-mini-2026-06-29-S6).

CLAIM (the "odd index" the project sought; repo lacks an explicit statement):
Let T be a SELF-CONVERSE tournament with an INVOLUTORY anti-automorphism phi
(phi^2 = id, and u->v in T  iff  phi(v)->phi(u) in T, i.e. phi: T ~= T^op).
Define the reversal-conjugation on directed Hamiltonian paths
   rho(P) = phi(reverse(P)),   P = (v_1,...,v_n) -> (phi(v_n),...,phi(v_1)).
Then:
  (1) rho maps Ham paths of T to Ham paths of T (well-defined),
  (2) rho is an INVOLUTION (rho^2 = id),
  (3) hence H(T) = #fixed(rho) + 2*#(free 2-orbits), so H(T) == #fixed(rho) (mod 2),
  (4) the fixed points of rho are the phi-PALINDROMIC Ham paths,
  (5) since H(T) is ODD (Redei), the number of phi-palindromic Ham paths is ODD (>=1).

This is the Hamiltonian-path-level twin of THM-281 (SC tiling fibers are odd because
the grid-symmetric/palindromic tilings are odd).  It is the sigma-ODD index, the
WITNESS/Borsuk-Ulam side -- distinct from the sigma-EVEN lonely-measure (existence).

We VERIFY (1)-(5) on Paley tournaments T_p (p=3 mod4, self-converse via phi(x)=-x)
and on small self-converse tournaments by exhaustive enumeration.
"""
from __future__ import annotations
import functools, itertools
print = functools.partial(print, flush=True)


def paley(p):
    """Paley tournament on Z_p: i->j iff (j-i) is a nonzero QR mod p. (p=3 mod4.)"""
    qr = set((x * x) % p for x in range(1, p))
    arc = [[False] * p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and ((j - i) % p) in qr:
                arc[i][j] = True
    return arc


def ham_paths(arc, n):
    """All directed Hamiltonian paths (as vertex tuples) of tournament `arc`."""
    out = []
    for perm in itertools.permutations(range(n)):
        if all(arc[perm[k]][perm[k + 1]] for k in range(n - 1)):
            out.append(perm)
    return out


def is_anti_aut_involution(arc, phi, n):
    """phi^2=id and (u->v) iff (phi(v)->phi(u))."""
    if any(phi[phi[x]] != x for x in range(n)):
        return False
    for u in range(n):
        for v in range(n):
            if u == v:
                continue
            if arc[u][v] != arc[phi[v]][phi[u]]:
                return False
    return True


def rho(P, phi):
    """rho(P) = phi(reverse(P))."""
    return tuple(phi[v] for v in reversed(P))


def analyze(name, arc, phi, n):
    if not is_anti_aut_involution(arc, phi, n):
        print(f"  {name}: phi is NOT an involutory anti-automorphism -- skip")
        return
    HP = ham_paths(arc, n)
    HPset = set(HP)
    H = len(HP)
    # (1) rho maps HP -> HP
    maps_in = all(rho(P, phi) in HPset for P in HP)
    # (2) involution
    invol = all(rho(rho(P, phi), phi) == P for P in HP)
    # (4) fixed points = palindromic
    fixed = [P for P in HP if rho(P, phi) == P]
    nf = len(fixed)
    print(f"  {name} (n={n}): H(T)={H} ({'ODD' if H%2 else 'EVEN'}), "
          f"rho:HP->HP={maps_in}, involution={invol}, "
          f"#palindromic={nf} ({'ODD' if nf%2 else 'EVEN'}), "
          f"H==#palin mod2: {H%2==nf%2}")
    if fixed:
        print(f"      example palindromic path: {fixed[0]}  (rho = {rho(fixed[0],phi)})")
    return H, nf


def main():
    print("=" * 80)
    print("Palindromic Hamiltonian path parity -- the odd index (mac-mini-S6)")
    print("=" * 80)

    print("\n[Paley tournaments T_p, p=3 mod4, self-converse via phi(x) = -x mod p]")
    for p in (3, 7, 11):
        arc = paley(p)
        phi = [(-x) % p for x in range(p)]
        analyze(f"Paley_{p}", arc, phi, p)

    print("\n[Small self-converse tournaments by search, n=4,5]")
    for n in (4, 5):
        # enumerate tournaments via upper-triangular arc bits
        pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        found = 0
        for bits in range(1 << len(pairs)):
            arc = [[False] * n for _ in range(n)]
            for b, (i, j) in enumerate(pairs):
                if (bits >> b) & 1:
                    arc[i][j] = True
                else:
                    arc[j][i] = True
            # find an involutory anti-automorphism phi
            phi_found = None
            for perm in itertools.permutations(range(n)):
                if all(perm[perm[x]] == x for x in range(n)) and is_anti_aut_involution(arc, perm, n):
                    phi_found = perm
                    break
            if phi_found is not None:
                found += 1
                if found <= 3:
                    analyze(f"SC_n{n}_#{found}", arc, list(phi_found), n)
        print(f"   n={n}: {found} self-converse tournaments (with involutory phi) found")

    print("\n" + "=" * 80)
    print("Conclusion: for self-converse T with involutory anti-aut phi, rho=phi.reverse")
    print("is an involution on Ham paths => H(T) == #palindromic (mod2); Redei (H odd)")
    print("=> ODD number of palindromic Ham paths >=1.  This is the sigma-ODD index")
    print("(witness side); the LRC lonely measure is the sigma-EVEN index (existence).")
    print("=" * 80)


if __name__ == "__main__":
    main()
