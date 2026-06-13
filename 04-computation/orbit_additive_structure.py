#!/usr/bin/env python3
"""
orbit_additive_structure.py -- What distinguishes the 4 extended orbits at p=11?

Each orbit has identical H. The orbits are under (Z/pZ)* action on S.
What additive invariant distinguishes them?

Candidates:
1. Number of "sum-free" pairs: #{(i,j): s_i + s_j not in S}
2. Number of "sum-closed" triples: #{(i,j,k): s_i + s_j = s_k}
3. Schur number / Rado-type invariants
4. The "3-cycle overlap weight" (how many vertex pairs share a 3-cycle)

Author: kind-pasteur-2026-03-12-S60
"""

from itertools import combinations


def legendre(a, p):
    if a % p == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def additive_invariants(S, p):
    """Compute additive structure invariants of S subset Z_p."""
    S_set = set(S)
    m = len(S)

    # Sum-closed pairs: s_i + s_j in S (mod p), i <= j
    sum_in = 0
    sum_out = 0
    for i in range(m):
        for j in range(i, m):
            s = (S[i] + S[j]) % p
            if s in S_set:
                sum_in += 1
            elif s != 0:  # s=0 means s_i + s_j = p
                sum_out += 1

    # Difference pairs: s_i - s_j in S (mod p), i != j
    diff_in = 0
    for i in range(m):
        for j in range(m):
            if i == j:
                continue
            d = (S[i] - S[j]) % p
            if d in S_set:
                diff_in += 1

    # Schur triples: s_i + s_j = s_k (i,j,k not necessarily distinct)
    schur = 0
    for i in range(m):
        for j in range(i, m):
            s = (S[i] + S[j]) % p
            if s in S_set:
                schur += 1

    # 3-AP count: s_i + s_k = 2*s_j (arithmetic progressions)
    ap3 = 0
    for j in range(m):
        for d_val in range(1, p):
            a = (S[j] - d_val) % p
            c = (S[j] + d_val) % p
            if a in S_set and c in S_set:
                ap3 += 1

    # Energy: E(S) = #{(a,b,c,d) in S^4: a+b=c+d mod p}
    # = sum_t |hat(1_S)(t)|^4
    # For our half-set, this is a multiplicative invariant
    from collections import Counter
    sumset = Counter()
    for i in range(m):
        for j in range(m):
            s = (S[i] + S[j]) % p
            sumset[s] += 1
    energy = sum(v**2 for v in sumset.values())

    # Product-set: |S*S| where s*s' mod p
    prodset = set()
    for i in range(m):
        for j in range(m):
            prodset.add((S[i] * S[j]) % p)
    prod_size = len(prodset)

    return {
        'sum_in': sum_in,
        'schur': schur,
        'diff_in': diff_in,
        'ap3': ap3,
        'energy': energy,
        'prod_size': prod_size,
    }


def main():
    for p in [7, 11, 13]:
        m = (p - 1) // 2
        half = 1 << m
        pairs = [(s, p - s) for s in range(1, m + 1)]
        qr = set(a for a in range(1, p) if legendre(a, p) == 1)

        print("=" * 70)
        print(f"p={p}, m={m}, QR={sorted(qr)}")
        print("=" * 70)

        # Build orientations and extended orbits
        orientations = {}
        for bits in range(half):
            S = sorted(pairs[i][0] if bits & (1 << i) else pairs[i][1]
                       for i in range(m))
            orientations[bits] = S

        # Extended orbits under (Z/pZ)*
        ext_orbits = {}
        ext_map = {}
        for bits in range(half):
            S = set(orientations[bits])
            orbit_members = set()
            for q in range(1, p):
                S_q = frozenset((s * q) % p for s in S)
                for b2 in range(half):
                    if frozenset(orientations[b2]) == S_q:
                        orbit_members.add(b2)
                        break
            orbit_key = frozenset(orbit_members)
            if orbit_key not in ext_orbits:
                ext_orbits[orbit_key] = []
            if bits not in ext_map:
                ext_map[bits] = orbit_key
                ext_orbits[orbit_key].append(bits)

        print(f"  {len(ext_orbits)} extended orbits")

        # For each orbit, compute additive invariants
        for i, (key, members) in enumerate(sorted(ext_orbits.items(), key=lambda x: min(x[1]))):
            rep = min(members)
            S = orientations[rep]
            inv = additive_invariants(S, p)
            is_paley = frozenset(S) == frozenset(qr)
            is_anti = frozenset(S) == frozenset(set(range(1, p)) - qr)
            label = ""
            if any(frozenset(orientations[b]) == frozenset(qr) for b in members):
                label = " [PALEY]"

            print(f"\n  Orbit {i} ({len(members)} members){label}:")
            print(f"    Rep S = {S}")
            print(f"    sum_in={inv['sum_in']}, schur={inv['schur']}, "
                  f"diff_in={inv['diff_in']}")
            print(f"    ap3={inv['ap3']}, energy={inv['energy']}, "
                  f"prod_size={inv['prod_size']}")

            # Check if all members have same invariants
            all_same = True
            for b in members:
                inv2 = additive_invariants(orientations[b], p)
                if inv2 != inv:
                    all_same = False
                    print(f"    WARNING: bits={b:05b} has different invariants!")
                    break
            if all_same and len(members) > 1:
                print(f"    (all {len(members)} members have identical invariants)")

        # Which invariant DISTINGUISHES the orbits?
        print(f"\n  Distinguishing power:")
        all_invs = {}
        for key, members in ext_orbits.items():
            rep = min(members)
            all_invs[key] = additive_invariants(orientations[rep], p)

        for inv_name in ['sum_in', 'schur', 'diff_in', 'ap3', 'energy', 'prod_size']:
            vals = [all_invs[k][inv_name] for k in sorted(all_invs.keys(), key=lambda x: min(x))]
            n_distinct = len(set(vals))
            distinguishes = n_distinct == len(ext_orbits)
            print(f"    {inv_name}: values={vals}, "
                  f"{'DISTINGUISHES' if distinguishes else f'only {n_distinct}/{len(ext_orbits)}'}")


if __name__ == '__main__':
    main()
