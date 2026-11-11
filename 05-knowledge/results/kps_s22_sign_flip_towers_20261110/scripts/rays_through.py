"""rays_through.py - all rays r + 2^t M of pairs (B, O) of order <= 7 through given numbers (kind-pasteur-2026-11-10-S22).

Input: the file produced by
    for m in 1 2 3 4 5 6 7; do ./portrait $m 100000000 < cls_$m.txt; done > pt_all_le7.txt
(portrait.c of this folder; cls_m.txt = class lists of classgen.c, 05-knowledge/results/first_creators_genealogy_20261107/;
line format of the input: m code Omask r h M E).
Control: 7 and 21 are not Hamiltonian path counts (THM-1370), so no pair may have a ray through them.
usage: python rays_through.py pt_all_le7.txt"""
import sys, collections
rows = [tuple(map(int, l.split())) for l in open(sys.argv[1])]
targets = [7, 21, 13, 19, 27, 1039]
rays = collections.defaultdict(list)
for m, code, O, r, h, M, E in rows:
    if M <= 0:
        continue
    for N in targets:
        d = N - r
        if d > 0 and d % M == 0:
            q = d // M
            if q >= 2 and q & (q - 1) == 0:
                rays[N].append((m, code, O, r, h, M, q.bit_length() - 1))
print("pairs read:", len(rows))
for N in targets:
    L = rays[N]
    print("N = %d: %d rays r + 2^t M (t >= 1) through it" % (N, len(L)))
    seen = set()
    for m, code, O, r, h, M, t in sorted(L, key=lambda z: (z[0], abs(z[3]))):
        if (r, M) in seen:
            continue
        seen.add((r, M))
        if len(seen) <= 6:
            print("     %d = %d + 2^%d * %d   (base of order %d, code %d, out-set mask %d, H(B) = %d)" % (N, r, t, M, m, code, O, h))
print("control (no ray through 7 or 21):", len(rays[7]) == 0 and len(rays[21]) == 0)
