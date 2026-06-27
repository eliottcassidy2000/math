#!/usr/bin/env python3
"""Exact audit for the THM-079-style LRC14 proof template.

This is not a proof of LRC(14).  It checks the local facts used by the
template and prints the remaining theorem-level obligation:

    STAR: M(S)=1/14 forces the denom-14 apex/tight-locus row.

The important guardrail is that "covering blocks denom 14" is only useful
after STAR.  Covering rows can still be safely lonely off the apex.

Post-pull integration: THM-568 proves the apex-denominator half of STAR0 for
tight rows.  The remaining theorem target is the covering-strictness residual,
especially rows with at least seven multiples of 14.
"""

from fractions import Fraction as F
from itertools import combinations
from math import gcd


TARGET = F(1, 14)


def frac_norm(x):
    r = x % 1
    return r if r <= F(1, 2) else 1 - r


def g_exact(speeds, tau):
    return min(frac_norm(v * tau) for v in speeds)


def candidate_taus(speeds):
    speeds = sorted(set(speeds))
    cands = {F(1, 2)}
    for v in speeds:
        k = 0
        while True:
            t = F(2 * k + 1, 2 * v)
            if t > F(1, 2):
                break
            cands.add(t)
            k += 1
    for a, b in combinations(speeds, 2):
        for d in (a + b, b - a):
            if d <= 0:
                continue
            k = 1
            while True:
                t = F(k, d)
                if t > F(1, 2):
                    break
                cands.add(t)
                k += 1
    return cands


def m_exact(speeds):
    best = F(0)
    arg = None
    for tau in candidate_taus(speeds):
        value = g_exact(speeds, tau)
        if value > best:
            best = value
            arg = tau
    return best, arg


def primitive(speeds):
    g = 0
    for v in speeds:
        g = gcd(g, v)
    return g == 1


def covering(speeds):
    return all(any(v % q == 0 for v in speeds) for q in range(2, 15))


def blocked_at_denom14(speeds):
    units = [1, 3, 5, 9, 11, 13]
    return all(g_exact(speeds, F(a, 14)) < TARGET for a in units)


def relation(value):
    if value < TARGET:
        return "<"
    if value == TARGET:
        return "="
    return ">"


def independent_polynomial_at_two(n, edges):
    edges = {tuple(sorted(e)) for e in edges}
    total = 0
    for mask in range(1 << n):
        verts = [i for i in range(n) if (mask >> i) & 1]
        independent = True
        for a, b in combinations(verts, 2):
            if tuple(sorted((a, b))) in edges:
                independent = False
                break
        if independent:
            total += 2 ** len(verts)
    return total


def connected(n, edges):
    if n == 0:
        return False
    adj = [[] for _ in range(n)]
    for a, b in edges:
        adj[a].append(b)
        adj[b].append(a)
    seen = {0}
    stack = [0]
    while stack:
        v = stack.pop()
        for w in adj[v]:
            if w not in seen:
                seen.add(w)
                stack.append(w)
    return len(seen) == n


def connected_i7_atoms(max_n=5):
    atoms = []
    for n in range(1, max_n + 1):
        pairs = list(combinations(range(n), 2))
        for mask in range(1 << len(pairs)):
            edges = [pairs[i] for i in range(len(pairs)) if (mask >> i) & 1]
            if connected(n, edges) and independent_polynomial_at_two(n, edges) == 7:
                atoms.append((n, tuple(edges)))
    return atoms


def graph_name(n, edges):
    if n == 3 and len(edges) == 3:
        return "K3"
    return f"n={n}, edges={edges}"


def row_report(name, speeds):
    m_value, tau = m_exact(speeds)
    print(
        f"{name:24s} primitive={primitive(speeds)!s:5s} "
        f"covering={covering(speeds)!s:5s} "
        f"has14={any(v % 14 == 0 for v in speeds)!s:5s} "
        f"denom14_blocked={blocked_at_denom14(speeds)!s:5s} "
        f"M={m_value} {relation(m_value)} 1/14 tau={tau}"
    )


def main():
    print("LRC14 THM-079-template STAR audit (codex S119)")
    print("=" * 72)
    print("Exact M checks")
    rows = [
        ("AP {1..13}", list(range(1, 14))),
        ("GW {1..11,13,24}", [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 24]),
        ("cover C_1", [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 84]),
        ("cover C_2", [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 168]),
        ("cover C_6", [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 504]),
    ]
    for name, speeds in rows:
        row_report(name, speeds)

    print()
    print("Apex exclusion is necessary but not sufficient")
    print("- AP and GW are tight, non-covering, and not denom-14 blocked.")
    print("- Covering rows are denom-14 blocked because they contain a multiple of 14.")
    print("- The covering rows above still have M>1/14 at off-apex binding-pair times.")
    print("This was the pre-THM568 gap: apex blocking needed STAR to close.")

    print()
    print("Connected graph preimage of I(G,2)=7")
    atoms = connected_i7_atoms(5)
    for n, edges in atoms:
        print(f"  {graph_name(n, edges)}")
    print(f"total connected atoms through n<=5: {len(atoms)}")

    print()
    print("STAR theorem target")
    print("  Original STAR+: every primitive reduced atom has M(S)>=1/14,")
    print("  equality only on the AP/Goddyn-Wong apex tight locus.")
    print("  Incoming THM-568 proves the structural STAR0 half: primitive")
    print("  tight rows have apex denominator D=14.  The residual is")
    print("  covering strictness: a 14-covering primitive atom cannot be")
    print("  tight, hence M(S)>1/14.")

    print()
    print("Tournament Analysis / assumption challenge")
    print("  Vertices used here are proof obligations and graph atoms, not runners.")
    print("  Dependency orientation: peel_exhaustive -> THM568 apex denominator")
    print("  -> covering_strictness -> state_lift.  This preserves the LRC predicate 'counterexample")
    print("  impossible' but destroys raw speed identity; the exact-M rows above")
    print("  show why that destruction is acceptable only after THM568 plus covering strictness.")


if __name__ == "__main__":
    main()
