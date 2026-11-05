"""aud_family.py -- audit of the S77 explicit digit-1 family (ternary_parent_digit_families_20261104, section 3 of the note).

Family (as in the S77 program header): for t >= 1 and bits b_1..b_t, C(b) has vertices r, w, z_1..z_t with r -> w,
w -> every z_i, z_i -> z_j (i < j), r -> z_i iff b_i = 1.  B = C^op plus a sole exit a -> r.  Then sigma_a(B) = 1 and
H(B) = h(C) + F_r(C) - 1 (one-exit law on the converse).

Independent closed form derived for this audit (C is the transitive order w > z_1 > ... > z_t with r inserted):
    F_r(C) = 1 + sum_{i: b_i = 0} 2^i
    h(C)   = 1 + [b_t = 0] + sum_{i: b_i = 0} 2^(i-1)
               + sum_{j: b_j = 1} ( [j >= 2 and b_(j-1) = 0] + sum_{i > j, b_i = 0} 2^(i-j-1) )
It is checked against (1) the C dynamic programme of aud_c2.c (`aud_c2 family 11 0`, every b, t <= 11) and
(2) brute-force permutation counting of H(B) in this script (every b, t <= 5).

Usage: python3 aud_family.py <dp_output_file> [TMAX]
"""
import sys
from itertools import permutations


def closed_H(t, b):
    bit = [None] + [(b >> (i - 1)) & 1 for i in range(1, t + 1)]
    F = 1 + sum(2 ** i for i in range(1, t + 1) if bit[i] == 0)
    h = 1 + (1 if bit[t] == 0 else 0) + sum(2 ** (i - 1) for i in range(1, t + 1) if bit[i] == 0)
    for j in range(1, t + 1):
        if bit[j] == 1:
            if j >= 2 and bit[j - 1] == 0:
                h += 1
            h += sum(2 ** (i - j - 1) for i in range(j + 1, t + 1) if bit[i] == 0)
    return h + F - 1


def family_tournament(t, b):
    """vertices: r=0, w=1, z_i=1+i (i=1..t), a=t+2; returns adjacency set of arcs (u, v) meaning u -> v in B."""
    n = t + 2
    C = set()
    C.add((0, 1))
    for i in range(1, t + 1):
        C.add((1, 1 + i))
        for j in range(i + 1, t + 1):
            C.add((1 + i, 1 + j))
        if (b >> (i - 1)) & 1:
            C.add((0, 1 + i))
        else:
            C.add((1 + i, 0))
    arcs = set((v, u) for (u, v) in C)  # converse
    a = n
    arcs.add((a, 0))
    for v in range(1, n):
        arcs.add((v, a))
    return n + 1, arcs


def brute_H(nv, arcs):
    cnt = 0
    for p in permutations(range(nv)):
        if all((p[i], p[i + 1]) in arcs for i in range(nv - 1)):
            cnt += 1
    return cnt


def main():
    dpfile = sys.argv[1]
    tmax = int(sys.argv[2]) if len(sys.argv) > 2 else 14
    dp = {}
    with open(dpfile) as fh:
        for line in fh:
            if line.startswith("F "):
                _, t, b, H = line.split()
                dp[(int(t), int(b))] = int(H)
    mism = sum(1 for (t, b), H in dp.items() if closed_H(t, b) != H)
    print(f"closed form vs C DP: {len(dp)} family members (t <= {max(t for t, _ in dp)}), mismatches {mism}")
    bad = 0
    nb = 0
    for t in range(1, 6):
        for b in range(2 ** t):
            nv, arcs = family_tournament(t, b)
            if brute_H(nv, arcs) != closed_H(t, b):
                bad += 1
            nb += 1
    print(f"closed form vs brute-force permutations: {nb} members (t <= 5), mismatches {bad}")

    vals = set()
    per_t_max = {}
    for t in range(1, tmax + 1):
        for b in range(2 ** t):
            H = closed_H(t, b)
            vals.add(H)
        per_t_max[t] = max(closed_H(t, b) for b in range(2 ** t))
    print(f"family values for t <= {tmax}: {len(vals)} distinct; largest {max(vals)}; per-t maxima {per_t_max}")
    # all values are odd?
    print("all values odd:", all(v % 2 == 1 for v in vals))
    print("dyadic-window coverage of odd numbers by the family (t <= %d):" % tmax)
    k = 1
    while 2 ** k <= max(vals):
        lo, hi = 2 ** k, 2 ** (k + 1)
        odd = [v for v in range(lo + 1, hi, 2)]
        hit = sum(1 for v in odd if v in vals)
        print(f"  [{lo}, {hi}): {hit}/{len(odd)} = {hit / len(odd):.3f}")
        k += 1
    gaps = [v for v in range(1, 401, 2) if v not in vals]
    print(f"odd values <= 400 never attained (t <= {tmax}): {len(gaps)}:", " ".join(map(str, gaps)))
    # what is the minimal t needed for each value <= 400? (stability check: are the small-window gaps final?)
    first_t = {}
    for t in range(1, tmax + 1):
        for b in range(2 ** t):
            H = closed_H(t, b)
            if H not in first_t:
                first_t[H] = t
    print("largest t at which a new value <= 400 first appears:", max((first_t[v] for v in first_t if v <= 400), default=None))
    # strong members only (b_t = 0 is necessary for strongness of C; report coverage restricted to b_t = 0)
    vals_s = set(closed_H(t, b) for t in range(1, tmax + 1) for b in range(2 ** t) if (b >> (t - 1)) & 1 == 0)
    print(f"members with b_t = 0 (strong C): {len(vals_s)} distinct values")
    for k in range(8, 16):
        lo, hi = 2 ** k, 2 ** (k + 1)
        odd = range(lo + 1, hi, 2)
        hit = sum(1 for v in odd if v in vals_s)
        print(f"  b_t = 0 only, [{lo}, {hi}): {hit}/{len(odd)} = {hit / len(odd):.3f}")


if __name__ == "__main__":
    main()
