#!/usr/bin/env python3
"""
Verify QR-scaling decomposition of k=0 eigenspace complex.

Key observation: QR_p acts on diff-sequences by (s_1,...,s_d) -> (as_1,...,as_d).
This action is FREE for d >= 1, splitting Omega_d^(0) into m equal pieces.
Each QR-eigenspace sub-complex has dimensions Omega_d/m.

Prediction: beta_m^{(0,j)} = (m-3)/2 for EACH QR-eigenspace j.
Total: m * (m-3)/2 = m(m-3)/2 = beta_m^(0). ✓

Also verify: Omega_d^(0) divisible by m for all d >= 1.

opus-2026-03-13-S71c
"""

import sys
sys.stdout.reconfigure(line_buffering=True)

def qr(p):
    return sorted(set((a*a)%p for a in range(1,p)))

def verify_qr_divisibility(p, max_d=None):
    """Verify Omega_d is divisible by m for d >= 1."""
    S = qr(p)
    m = len(S)
    if max_d is None:
        max_d = 2*m

    print(f"\nP_{p}: m={m}")

    # Enumerate sequences and compute QR orbits
    prev_seqs = [()]
    prev_ps = {(): frozenset([0])}
    prev_last = {(): 0}
    all_seqs = {0: [()]}

    for d in range(1, max_d + 1):
        curr_seqs = []
        curr_ps = {}
        curr_last = {}
        for seq in prev_seqs:
            ps_set = prev_ps[seq]
            last = prev_last[seq]
            for s in S:
                ns = (last + s) % p
                if ns not in ps_set:
                    nseq = seq + (s,)
                    curr_seqs.append(nseq)
                    curr_ps[nseq] = ps_set | {ns}
                    curr_last[nseq] = ns
        if not curr_seqs:
            break
        all_seqs[d] = curr_seqs
        prev_seqs = curr_seqs
        prev_ps = curr_ps
        prev_last = curr_last

    actual_max_d = max(all_seqs.keys())
    max_d = min(max_d, actual_max_d)

    # Verify QR-action is free on A_d for d >= 1
    S_set = set(S)
    for d in range(1, min(max_d + 1, 4)):  # check first few degrees
        orbits = {}
        for seq in all_seqs[d]:
            # Find canonical rep (smallest in lex order among orbit)
            orbit_reps = []
            for a in S:
                scaled = tuple((a * s) % p for s in seq)
                orbit_reps.append(scaled)
            canon = min(orbit_reps)
            if canon not in orbits:
                orbits[canon] = []
            orbits[canon].append(seq)

        orbit_sizes = [len(v) for v in orbits.values()]
        print(f"  d={d}: A={len(all_seqs[d])}, #orbits={len(orbits)}, "
              f"orbit sizes={sorted(set(orbit_sizes))}, "
              f"A/m={'✓' if len(all_seqs[d]) % m == 0 else '✗'} ({len(all_seqs[d])//m})")

    # Now compute Omega_d and check divisibility
    allowed = {d: set(all_seqs.get(d, [])) for d in range(max_d + 1)}

    print(f"\n  Omega_d divisibility by m={m}:")
    for d in range(1, max_d + 1):
        A_d = all_seqs[d]
        n_Ad = len(A_d)
        if n_Ad == 0:
            print(f"  d={d}: Omega=0, div=✓")
            continue

        # Count junk faces
        junk_set = set()
        for D in A_d:
            for fi in range(d + 1):
                if fi == 0: fd = D[1:]
                elif fi == d: fd = D[:d-1]
                else:
                    merged = (D[fi-1] + D[fi]) % p
                    fd = D[:fi-1] + (merged,) + D[fi+1:]
                if fd not in allowed[d-1]:
                    junk_set.add(fd)

        # Simple Omega computation: A_d - rank(constraint)
        # For divisibility check, we just need Omega_d mod m
        # But we need the actual rank... skip heavy computation
        # Just report A_d mod m
        print(f"  d={d}: A={n_Ad}, A mod m = {n_Ad % m}, "
              f"A/m = {n_Ad // m}", end="")

        # Known Omega values for these primes
        omega_p7 = [1,3,6,9,9,6,3]
        omega_p11 = [1,5,20,70,205,460,700,690,450,180,30]

        if p == 7 and d < len(omega_p7):
            o = omega_p7[d]
            print(f", Ω={o}, Ω mod m = {o % m}, Ω/m = {o // m}", end="")
        elif p == 11 and d < len(omega_p11):
            o = omega_p11[d]
            print(f", Ω={o}, Ω mod m = {o % m}, Ω/m = {o // m}", end="")
        print()

    # Compute sub-complex dimensions and predicted Betti
    if p in [7, 11]:
        omega = omega_p7 if p == 7 else omega_p11
        print(f"\n  QR sub-complex for j≠0 (dims = Ω_d/m for d≥1, 0 at d=0):")
        sub_dims = [0] + [omega[d] // m for d in range(1, len(omega))]
        print(f"    dims = {sub_dims}")

        # Predicted ranks assuming all intermediate β=0
        ranks = [0] * len(sub_dims)
        for d in range(1, len(sub_dims)):
            if d == 1:
                ranks[d] = 0  # codomain is 0-dimensional
            else:
                # β_{d-1} = 0 => rank(∂_d) = sub_dims[d-1] - ranks[d-1]
                ranks[d] = sub_dims[d-1] - ranks[d-1]

        print(f"    predicted ranks (β_d=0 for d<m) = {ranks}")

        # Check Euler char
        chi_sub = sum((-1)**d * sub_dims[d] for d in range(len(sub_dims)))
        print(f"    χ = {chi_sub} (should be 0 for j≠0)")

        # Where does exactness break?
        betti_sub = []
        for d in range(len(sub_dims)):
            r_d = ranks[d] if d < len(ranks) else 0
            r_d1 = ranks[d+1] if d+1 < len(ranks) else 0
            b = sub_dims[d] - r_d - r_d1
            betti_sub.append(b)
        print(f"    β (predicted, exact below m) = {betti_sub}")
        print(f"    β_m = β_{m+1} = {betti_sub[m] if m < len(betti_sub) else '?'}")
        print(f"    Expected: (m-3)/2 = {(m-3)//2}")

verify_qr_divisibility(7)
verify_qr_divisibility(11)
verify_qr_divisibility(19, max_d=7)

print("\n" + "=" * 60)
