#!/usr/bin/env python3
"""
orientation_symmetry_p11.py -- Analyze which orientation flips preserve H at p=11

The QR-orbit structure on chord pairs creates equivalence classes.
At p=7: only 2 H-values (175, 189). What's the structure at p=11?

Also: test whether the overlap weight matrix W determines the H-class.

Author: kind-pasteur-2026-03-12-S60
"""

import sys

def legendre(a, p):
    if a % p == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def qr_orbit_class(bits, m, p):
    """Which QR-orbit class does this orientation belong to?

    Two orientations are QR-equivalent if one can be obtained from
    the other by multiplying all gaps by a QR element.

    At p=11, QR = {1,3,4,5,9}. Multiplying gap i by q permutes the gaps.
    """
    pairs = [(s, p - s) for s in range(1, m + 1)]
    S = set()
    for i in range(m):
        if bits & (1 << i):
            S.add(pairs[i][0])
        else:
            S.add(pairs[i][1])

    # Compute all QR-equivalent S sets
    qr = [a for a in range(1, p) if legendre(a, p) == 1]
    orbit = set()
    for q in qr:
        S_q = frozenset((s * q) % p for s in S)
        # Convert to bits
        b = 0
        for i in range(m):
            if pairs[i][0] in S_q:
                b |= (1 << i)
            elif pairs[i][1] in S_q:
                pass
            else:
                # This shouldn't happen for valid QR orbit
                b = -1
                break
        if b >= 0:
            orbit.add(b)

    return frozenset(orbit)


def main():
    p = 11
    m = (p - 1) // 2  # 5
    half = 1 << m  # 32

    pairs = [(s, p - s) for s in range(1, m + 1)]
    qr = set(a for a in range(1, p) if legendre(a, p) == 1)

    print("=" * 70)
    print(f"ORIENTATION SYMMETRY ANALYSIS at p={p}")
    print(f"QR = {sorted(qr)}")
    print("=" * 70)

    # Map each orientation to its S set
    orientations = {}
    for bits in range(half):
        S = set()
        for i in range(m):
            if bits & (1 << i):
                S.add(pairs[i][0])
            else:
                S.add(pairs[i][1])
        orientations[bits] = frozenset(S)

    # Find QR orbits
    orbits = {}
    orbit_map = {}
    for bits in range(half):
        S = orientations[bits]
        # Apply all QR multipliers
        orbit_members = set()
        for q in qr:
            S_q = frozenset((s * q) % p for s in S)
            # Find which bits this corresponds to
            for b2 in range(half):
                if orientations[b2] == S_q:
                    orbit_members.add(b2)
                    break
        orbit_key = frozenset(orbit_members)
        if orbit_key not in orbits:
            orbits[orbit_key] = []
        if bits not in orbit_map:
            orbit_map[bits] = orbit_key
            orbits[orbit_key].append(bits)

    print(f"\n  {len(orbits)} QR-orbits:")
    for i, (key, members) in enumerate(sorted(orbits.items(), key=lambda x: min(x[1]))):
        S_list = [sorted(orientations[b]) for b in sorted(members)]
        print(f"\n    Orbit {i}: {len(members)} members")
        for b in sorted(members):
            is_paley = orientations[b] == qr
            is_anti = orientations[b] == set(range(1,p)) - qr
            label = " <-- PALEY" if is_paley else (" <-- anti-Paley" if is_anti else "")
            print(f"      bits={b:05b}, S={sorted(orientations[b])}{label}")

    # Which orbit is Paley in?
    for bits in range(half):
        if orientations[bits] == qr:
            print(f"\n  Paley: bits={bits:05b}, orbit has {len(orbits[orbit_map[bits]])} members")
            break

    # NQR orbit analysis: multiply by NQR
    nqr = set(a for a in range(1, p) if legendre(a, p) == -1)
    print(f"\n  NQR = {sorted(nqr)}")

    # For each orbit, what does NQR multiplication do?
    print(f"\n  NQR orbit pairing:")
    for key, members in sorted(orbits.items(), key=lambda x: min(x[1])):
        rep = min(members)
        S = orientations[rep]
        nq = min(nqr)
        S_nq = frozenset((s * nq) % p for s in S)
        for b2 in range(half):
            if orientations[b2] == S_nq:
                target_orbit = orbit_map[b2]
                same = (target_orbit == key)
                print(f"    Orbit of bits={rep:05b} "
                      f"-> NQR*S = bits={b2:05b} "
                      f"{'(SELF)' if same else '(different orbit)'}")
                break

    # The key prediction: QR-equivalent orientations have IDENTICAL H
    # (this is because multiplication by QR gives a graph automorphism)
    # NQR-equivalent orientations ALSO have identical H
    # (because -1 is NQR at p=11, so NQR mult = complement = T^op, and H(T)=H(T^op))

    # Extended orbits under full (Z/pZ)* action
    print(f"\n\n  Extended orbits under FULL (Z/pZ)* action:")
    ext_orbits = {}
    ext_map = {}
    for bits in range(half):
        S = orientations[bits]
        orbit_members = set()
        for q in range(1, p):
            S_q = frozenset((s * q) % p for s in S)
            for b2 in range(half):
                if orientations[b2] == S_q:
                    orbit_members.add(b2)
                    break
        orbit_key = frozenset(orbit_members)
        if orbit_key not in ext_orbits:
            ext_orbits[orbit_key] = []
        if bits not in ext_map:
            ext_map[bits] = orbit_key
            ext_orbits[orbit_key].append(bits)

    for i, (key, members) in enumerate(sorted(ext_orbits.items(), key=lambda x: min(x[1]))):
        print(f"\n    Extended Orbit {i}: {len(members)} members")
        for b in sorted(members):
            is_paley = orientations[b] == qr
            is_anti = orientations[b] == set(range(1,p)) - qr
            label = " <-- PALEY" if is_paley else (" <-- anti-Paley" if is_anti else "")
            print(f"      bits={b:05b}, S={sorted(orientations[b])}{label}")

    print(f"\n  {len(ext_orbits)} extended orbits (should partition into H-classes)")

    # Prediction: each extended orbit is an H-class
    # At p=7, there are 2 extended orbits -> 2 H-values
    # At p=11, the number of extended orbits = number of distinct H-values

    # Also check: what's the complement structure?
    print(f"\n\n  Complement structure (S -> p-S):")
    for bits in range(half):
        comp_bits = ((1 << m) - 1) ^ bits  # flip all bits
        same_orbit = ext_map.get(bits) == ext_map.get(comp_bits)
        if bits <= comp_bits:
            print(f"    bits={bits:05b} <-> comp={comp_bits:05b} "
                  f"{'same ext orbit' if same_orbit else 'DIFFERENT ext orbit'}")


if __name__ == '__main__':
    main()
