#!/usr/bin/env python3
"""glue_features.py file...: rebuild the printed gluing examples and report structural features of the bound-attained
cases: the roles of u, v (source/sink in A, in B), the cross pattern (all A'->B', all B'->A', mixed), the order statistics."""
import sys, re
def parse(s, n):
    adj = [[0]*n for _ in range(n)]; k = 0
    for i in range(n):
        for j in range(i+1, n):
            if s[k] == '1': adj[i][j] = 1
            else: adj[j][i] = 1
            k += 1
    return adj
for fn in sys.argv[1:]:
    print("==", fn)
    for l in open(fn):
        m = re.match(r'diff=(-?\d+) H\(U\)=(\d+) H\(A\)=(\d+) H\(B\)=(\d+) bound=(\d+) attained=(\d) A=(\d+) arc=\((\d+),(\d+)\) B=(\d+) arc=\((\d+),(\d+)\) cross=(\d+) hA\(uv\)=(\d+) hB\(uv\)=(\d+) strongU=(\d)', l)
        if not m: continue
        diff, HU, HA, HB, bound, att, As, u, v, Bs, u2, v2, cross, hAuv, hBuv, stU = m.groups()
        nA = int((1 + (1 + 8*len(As))**0.5)/2); nB = int((1 + (1 + 8*len(Bs))**0.5)/2)
        A = parse(As, nA); B = parse(Bs, nB); u, v, u2, v2, cross = int(u), int(v), int(u2), int(v2), int(cross)
        def role(adj, n, x):
            out = sum(adj[x]); inn = sum(adj[y][x] for y in range(n))
            return 'source' if inn == 0 else 'sink' if out == 0 else f'out{out}/in{inn}'
        ca = nA - 2; cb = nB - 2; npairs = ca * cb
        ones = bin(cross).count('1')
        pattern = "all A'->B'" if ones == npairs else "all B'->A'" if ones == 0 else f"mixed({ones}/{npairs})"
        print(f"  diff={diff} H(U)={HU} ({HA},{HB}) bound {bound} att={att} strongU={stU} | u: {role(A,nA,u)} in A, {role(B,nB,u2)} in B; v: {role(A,nA,v)} in A, {role(B,nB,v2)} in B | cross {pattern} | hA(uv)={hAuv}/{HA} hB(uv)={hBuv}/{HB}")
