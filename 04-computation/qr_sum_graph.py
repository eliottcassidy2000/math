#!/usr/bin/env python3
"""
Investigate the "QR sum graph" Γ_p: vertices = QR_p, edges = {a,b} with a+b ∈ QR.
Hypothesis: β_m^(0) = m(m-3)/2 = number of non-edges of cycle C_m.
Does Γ_p have cycle structure related to C_m?

Also: what is the graph Γ'_p with a~b iff a+b ∈ NQR?
For p≡3 mod 4: edges of Γ + edges of Γ' = C(m,2) - 0 (no pair sums to 0).

opus-2026-03-13-S71c
"""

import sys
sys.stdout.reconfigure(line_buffering=True)

def qr(p):
    return sorted(set((a*a)%p for a in range(1,p)))

def analyze_sum_graph(p):
    Q = qr(p)
    m = len(Q)
    QR_set = set(Q)

    print(f"\nP_{p}: m={m}, QR={Q}")

    # Build sum graph Γ: a~b iff a+b mod p ∈ QR
    edges_qr = []
    edges_nqr = []
    for i in range(m):
        for j in range(i+1, m):
            s = (Q[i] + Q[j]) % p
            if s in QR_set:
                edges_qr.append((Q[i], Q[j]))
            elif s != 0:
                edges_nqr.append((Q[i], Q[j]))

    print(f"  Γ (sum in QR):  |E| = {len(edges_qr)}")
    print(f"  Γ' (sum in NQR): |E| = {len(edges_nqr)}")
    print(f"  Total pairs: C(m,2) = {m*(m-1)//2}")
    print(f"  m(m-3)/2 = {m*(m-3)//2} (predicted β_m^(0))")

    # Degree sequence of Γ
    deg = {q: 0 for q in Q}
    for a, b in edges_qr:
        deg[a] += 1
        deg[b] += 1
    deg_seq = sorted(deg.values())
    print(f"  Γ degree sequence: {deg_seq}")

    # Is Γ vertex-transitive under QR-scaling? (a ↦ ca for c ∈ QR)
    # Check: is the degree sequence constant?
    is_regular = len(set(deg_seq)) == 1
    print(f"  Γ regular: {is_regular}" + (f" (degree {deg_seq[0]})" if is_regular else ""))

    # Check if Γ is a cycle C_m
    if len(edges_qr) == m:
        print(f"  Γ has m edges — could be C_m!")
        # Check if it's connected and all degrees = 2
        if all(d == 2 for d in deg_seq):
            print(f"  Γ IS a cycle C_m! ✓")

    # Complement: non-edges of Γ that are also not in Γ'
    # Since sum ≠ 0 for all pairs (since -1 ∉ QR):
    # |E(Γ)| + |E(Γ')| = C(m,2)
    print(f"  |E(Γ)| + |E(Γ')| = {len(edges_qr) + len(edges_nqr)} = C(m,2) = {m*(m-1)//2}")

    # Non-edges of C_m = diagonals of m-gon = m(m-3)/2
    # Does |E(Γ')| = m(m-3)/2? Or |E(Γ)| - something?
    print(f"  |E(Γ')| = {len(edges_nqr)}, m(m-3)/2 = {m*(m-3)//2}")
    print(f"  C(m,2) - |E(Γ')| = {m*(m-1)//2 - len(edges_nqr)}")

    # Cycle structure: find all cycles of Γ
    # First find connected components
    adj = {q: set() for q in Q}
    for a, b in edges_qr:
        adj[a].add(b)
        adj[b].add(a)

    visited = set()
    components = []
    for q in Q:
        if q not in visited:
            comp = set()
            stack = [q]
            while stack:
                v = stack.pop()
                if v in comp: continue
                comp.add(v)
                for u in adj[v]:
                    if u not in comp:
                        stack.append(u)
            visited |= comp
            components.append(comp)

    print(f"  Γ connected components: {len(components)}, sizes: {sorted(len(c) for c in components)}")

    # Check if Γ' is also interesting
    deg_nqr = {q: 0 for q in Q}
    for a, b in edges_nqr:
        deg_nqr[a] += 1
        deg_nqr[b] += 1
    deg_seq_nqr = sorted(deg_nqr.values())
    print(f"  Γ' degree sequence: {deg_seq_nqr}")
    is_regular_nqr = len(set(deg_seq_nqr)) == 1
    print(f"  Γ' regular: {is_regular_nqr}" + (f" (degree {deg_seq_nqr[0]})" if is_regular_nqr else ""))

    # Key question: β_m^(0) = m(m-3)/2 = C(m,2) - m
    # This equals number of non-edges of ANY m-cycle
    # But does it relate to Γ or Γ' specifically?

    # First Betti number of Γ: β_1(Γ) = |E| - |V| + #components
    b1_gamma = len(edges_qr) - m + len(components)
    print(f"  β_1(Γ) = {len(edges_qr)} - {m} + {len(components)} = {b1_gamma}")

    b1_gamma_nqr = len(edges_nqr) - m + 1  # assuming connected
    print(f"  β_1(Γ') ≈ {b1_gamma_nqr} (if connected)")

    return edges_qr, edges_nqr

print("=" * 60)
print("QR SUM GRAPH ANALYSIS")
print("=" * 60)

for p in [7, 11, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83]:
    try:
        analyze_sum_graph(p)
    except Exception as e:
        print(f"  Error: {e}")

print("\n" + "=" * 60)
