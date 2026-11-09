"""AUDIT D1 (THM-5055 Prop 3.3): for odd n in [5, NMAX], list EVERY graph G on n vertices whose complement has
maximum degree <= 2 (a disjoint union of paths P_k, k >= 1, and cycles C_k, k >= 3) - by the counting step of
Prop 3.3 these are the only candidates for Def(G) = n + 1 - and write them in aw_verify format with k = X(n).
aw_verify then reports d(G) and the parity; the check is that every G with d(G) = X(n) is ODD, and that the
realisers are exactly the shapes the proof lists (cG = P_3 u matching, or K_1 u H u matching, H in
{P_4, 2P_3, C_4}).
usage: python3 aw_prop33.py NMAX > file"""
import sys
NMAX = int(sys.argv[1])
def parts(n, minpart=(0, 1)):
    """multisets of components: ('P', k) k >= 1 or ('C', k) k >= 3, in nondecreasing order, total size n"""
    if n == 0: yield []; return
    comps = [('P', k) for k in range(1, n + 1)] + [('C', k) for k in range(3, n + 1)]
    comps.sort(key=lambda c: (c[1], c[0]))
    for c in comps:
        key = (c[1], c[0])
        if key < (minpart[1], minpart[0]) if isinstance(minpart[0], str) else False: continue
        if c[1] > n: continue
        for rest in parts(n - c[1], c):
            yield [c] + rest
def parts2(n):
    comps = sorted([('P', k) for k in range(1, n + 1)] + [('C', k) for k in range(3, n + 1)], key=lambda c: (c[1], c[0]))
    def rec(rem, start):
        if rem == 0: yield []; return
        for i in range(start, len(comps)):
            c = comps[i]
            if c[1] > rem: break
            for rest in rec(rem - c[1], i): yield [c] + rest
    yield from rec(n, 0)
for n in range(5, NMAX + 1, 2):
    X = 2**n - n - 2
    for shape in parts2(n):
        cedges = set(); v = 0
        for t, k in shape:
            for i in range(k - 1): cedges.add((v + i, v + i + 1))
            if t == 'C': cedges.add((v, v + k - 1))
            v += k
        g = [(a, b) for a in range(n) for b in range(a + 1, n) if (a, b) not in cedges]
        tag = "co[" + ",".join(f"{t}{k}" for t, k in shape) + "]"
        print(f"{tag} E {n} {X} | " + " ".join(f"{a}-{b}" for a, b in g))
