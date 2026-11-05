"""H-values of seam words built only from two-port transitive gadgets TT_k (k <= K, all ports), THM-4879 (B) closed
forms. Seeds (1,1,1) and (1,0,2^z). Usage: python3 transitive_words_bfs.py X K"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from transitive_two_port import closed
X = int(sys.argv[1]); K = int(sys.argv[2])
maps = sorted(set(tuple(map(tuple, closed(k, u, z))) for k in range(1, K + 1) for u in range(k) for z in range(k)))
V = {(1, 1, 1)} | {(1, 0, 2 ** z) for z in range(1, K) if 2 ** z <= X}; fr = list(V)
while fr:
    nx = []
    for (H, S, R) in fr:
        for M in maps:
            h = M[0][0]*H + M[0][1]*S + M[0][2]*R
            r = M[2][0]*H + M[2][1]*S + M[2][2]*R
            if h > X or r > X: continue
            p = (h, M[1][0]*H + M[1][1]*S + M[1][2]*R, r)
            if p not in V: V.add(p); nx.append(p)
    fr = nx
Hs = {h for h, s, r in V}
miss = [n for n in range(1, X + 1, 2) if n not in Hs]
print(f'X={X} K={K} distinct maps={len(maps)} states={len(V)} missing odd H: {len(miss)}; first: {miss[:30]}')
