#!/usr/bin/env python3
"""
paley19_packed_omega.py — opus-2026-03-13-S71b

Compute Omega dims for P_19 using packed integer diff-seq encoding.
Each diff-seq step (from 9 possible values) is encoded in 4 bits.
A d-step diff-seq fits in 4*d bits (at most 72 bits for d=18).

This reduces memory from ~120 bytes/tuple to 8 bytes/int.
Target: push past degree 8 to compute remaining Omega_d for chi(P_19).
"""

import sys, time, struct
sys.path.insert(0, '04-computation')
from circulant_homology import sparse_rank_mod

# P_19 connection set
P = 19
S_LIST = sorted([1, 4, 5, 6, 7, 9, 11, 16, 17])  # QR_19
M = len(S_LIST)  # = 9

# Map S values to 4-bit codes and back
S_TO_CODE = {s: i for i, s in enumerate(S_LIST)}
CODE_TO_S = {i: s for i, s in enumerate(S_LIST)}

PRIME = 191  # Working prime for P_19

def pack_seq(tup):
    """Pack a diff-seq tuple into a single integer."""
    val = 0
    for i, s in enumerate(tup):
        val |= (S_TO_CODE[s] << (4 * i))
    return val

def unpack_seq(val, length):
    """Unpack an integer back to a diff-seq tuple."""
    tup = []
    for i in range(length):
        code = (val >> (4 * i)) & 0xF
        tup.append(CODE_TO_S[code])
    return tuple(tup)

def compute_face_packed(packed, length, face_idx, n):
    """Compute face of packed diff-seq, return as tuple (for hashing into lower set)."""
    # Unpack, compute face, return tuple
    D = unpack_seq(packed, length)
    m = length
    if face_idx == 0:
        return D[1:]
    elif face_idx == m:
        return D[:m - 1]
    else:
        merged = (D[face_idx - 1] + D[face_idx]) % n
        return D[:face_idx - 1] + (merged,) + D[face_idx + 1:]

def enumerate_degree_packed(prev_level, n, S_sorted):
    """
    Given prev_level = dict {packed_seq: (visited_as_int, last_vertex)},
    extend by one step.

    visited_as_int: bitmask of visited vertices (bit i = vertex i visited).
    Returns new dict of same format.
    """
    next_level = {}
    for packed, (visited, last) in prev_level.items():
        deg = 0
        tmp = packed
        while tmp > 0 or deg == 0:
            if packed == 0 and deg == 0:
                break
            deg += 1
            tmp >>= 4
        # Actually, we need to track degree separately
        for s in S_sorted:
            nv = (last + s) % n
            if not (visited & (1 << nv)):
                code = S_TO_CODE[s]
                # Append code at position deg
                # But we don't know deg from packed alone...
                # Need to pass deg
                pass
    return next_level


def enumerate_packed_with_deg(prev_dict, deg, n, S_sorted):
    """
    prev_dict: {packed_seq: (visited_bitmask, last_vertex)} for degree (deg-1) seqs.
    Returns: {packed_seq: (visited_bitmask, last_vertex)} for degree deg seqs.
    """
    next_dict = {}
    count = 0
    t0 = time.time()
    total = len(prev_dict)
    for packed, (visited, last) in prev_dict.items():
        for s in S_sorted:
            nv = (last + s) % n
            if not (visited & (1 << nv)):
                code = S_TO_CODE[s]
                new_packed = packed | (code << (4 * (deg - 1)))
                new_visited = visited | (1 << nv)
                next_dict[new_packed] = (new_visited, nv)
        count += 1
        if count % 5000000 == 0:
            elapsed = time.time() - t0
            rate = count / elapsed
            eta = (total - count) / rate
            print(f"    Extending: {count}/{total} ({100*count/total:.1f}%), "
                  f"{len(next_dict)} new, {elapsed:.0f}s, ETA {eta:.0f}s", flush=True)
    return next_dict


def compute_omega_degree(target_deg, verbose=True):
    """
    Compute Omega_{target_deg}/p for P_19.

    Strategy:
    1. Enumerate degree (target_deg-1) diff-seqs as packed ints → lower_set
    2. Enumerate degree target_deg diff-seqs, build constraint columns
    3. Compute sparse rank
    """
    n = P

    if verbose:
        print(f"\n{'='*60}")
        print(f"Computing Omega_{target_deg}/p for P_{n}")
        print(f"{'='*60}")

    # Build up from degree 0
    # degree 0: just the empty seq, starting at vertex 0
    current = {0: (1 << 0, 0)}  # packed=0 (empty), visited={0}, last=0

    for deg in range(1, target_deg + 1):
        t0 = time.time()
        if verbose:
            print(f"  Enumerating degree {deg}...", end=" ", flush=True)

        next_level = enumerate_packed_with_deg(current, deg, n, S_LIST)
        elapsed = time.time() - t0

        if verbose:
            print(f"|A_{deg}|/p = {len(next_level)}, time = {elapsed:.1f}s", flush=True)

        if deg == target_deg - 1:
            # Save as lower set (unpack to tuples for face checking)
            if verbose:
                print(f"  Building lower set from degree {deg}...", end=" ", flush=True)
            t1 = time.time()
            lower_set = set()
            for packed, _ in next_level.items():
                lower_set.add(unpack_seq(packed, deg))
            elapsed = time.time() - t1
            if verbose:
                print(f"done, {len(lower_set)} entries, {elapsed:.1f}s", flush=True)

        if deg < target_deg:
            current = next_level
        else:
            # This is the target degree — build constraints
            if verbose:
                print(f"  Building constraint columns for degree {deg}...", flush=True)

            t1 = time.time()
            junk_idx = {}
            columns = []
            n_Am = 0

            for packed, (visited, last) in next_level.items():
                D = unpack_seq(packed, deg)
                n_Am += 1

                col = {}
                for fi in range(deg + 1):
                    fd = compute_face_packed(packed, deg, fi, n)
                    if fd not in lower_set:
                        if fd not in junk_idx:
                            junk_idx[fd] = len(junk_idx)
                        row = junk_idx[fd]
                        sign = 1 if fi % 2 == 0 else -1
                        col[row] = (col.get(row, 0) + sign) % PRIME
                col = {r: v for r, v in col.items() if v != 0}
                columns.append(col)

                if verbose and n_Am % 2000000 == 0:
                    elapsed = time.time() - t1
                    print(f"    {n_Am/1e6:.1f}M processed, {len(junk_idx)} junk rows, "
                          f"{elapsed:.0f}s", flush=True)

            elapsed = time.time() - t1
            n_junk = len(junk_idx)
            if verbose:
                print(f"    Done: |A_{deg}|/p = {n_Am}, junk rows = {n_junk}, "
                      f"constraint build time = {elapsed:.1f}s", flush=True)

            # Free memory before rank computation
            del next_level
            del lower_set
            del current

            # Compute rank
            if n_junk == 0:
                rank = 0
            else:
                if verbose:
                    print(f"  Computing rank of {n_junk} x {n_Am} matrix...", flush=True)
                t2 = time.time()
                rank, mpd = sparse_rank_mod(columns, PRIME,
                                           verbose=verbose and n_Am > 500000,
                                           verbose_interval=200000)
                elapsed = time.time() - t2
                if verbose:
                    print(f"    Rank = {rank}, mpd = {mpd}, time = {elapsed:.1f}s", flush=True)

            omega = n_Am - rank
            return omega, n_Am


def main():
    print("=" * 70)
    print(f"P_{P} PACKED OMEGA COMPUTATION")
    print(f"S = {S_LIST}, working mod {PRIME}")
    print("=" * 70)

    known = [1, 9, 72, 540, 3753, 23832, 136260, 688266, 2987622]
    print(f"Known Omega/p (d=0..{len(known)-1}): {known}")
    chi_partial = sum((-1)**d * o for d, o in enumerate(known))
    print(f"Partial chi/p (d=0..{len(known)-1}) = {chi_partial}")

    # Try degree 9
    for d in range(9, 19):
        t0 = time.time()
        try:
            omega, n_Am = compute_omega_degree(d, verbose=True)
            elapsed = time.time() - t0

            known.append(omega)
            chi_partial = sum((-1)**dd * o for dd, o in enumerate(known))

            print(f"\n  RESULT: Omega_{d}/p = {omega}")
            print(f"  |A_{d}|/p = {n_Am}")
            print(f"  Omega/A = {omega/n_Am:.6f}")
            print(f"  chi/p (d=0..{d}) = {chi_partial}")
            print(f"  Total time: {elapsed:.1f}s")

        except MemoryError:
            print(f"\n  OOM at degree {d}")
            break
        except Exception as e:
            import traceback
            print(f"\n  Error at degree {d}: {e}")
            traceback.print_exc()
            break

    print(f"\nFinal Omega/p: {known}")
    chi = sum((-1)**d * o for d, o in enumerate(known))
    print(f"chi/p = {chi}")
    print(f"chi = {chi * P}")
    print("\nDONE.")


if __name__ == "__main__":
    main()
