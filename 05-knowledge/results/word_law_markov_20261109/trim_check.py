"""M1 check (S88 audit): the trim block automaton of THM-5021.
Reachable and co-reachable block-boundary states, from the 18-state construction of markov_tight.py; the trim
transfer matrix, its characteristic polynomial and Jordan form, and the three generating functions.
usage: python trim_check.py
"""
import importlib.util, os
import sympy as sp
here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("mt", os.path.join(here, "markov_tight.py"))
mt = importlib.util.module_from_spec(spec); spec.loader.exec_module(mt)
states = sorted({(a, b, e) for a in (1, 2, 3) for b in (1, 2, 3) for e in (0, 1)})
def nxt(s, ph):
    out = []
    for col in (0, 1):
        t = mt.step(s, col, ph)
        if t[2] in (0, 1):
            out.append(t)
    return out
# states occurring at any position (phase tracked) starting from (3,3,0) at phase 0
occ = {(0, (3, 3, 0))}; frontier = [(0, (3, 3, 0))]
while frontier:
    new = []
    for ph, s in frontier:
        for t in nxt(s, ph):
            key = ((ph + 1) % 3, t)
            if key not in occ:
                occ.add(key); new.append(key)
    frontier = new
used = sorted({s for ph, s in occ})
boundary = sorted({s for ph, s in occ if ph == 0})
# co-reachable: can reach E = 0 at some later block boundary (or at the end of a word of any residue)
def block(s):
    res = {}
    for c1 in (0, 1):
        t1 = mt.step(s, c1, 0)
        if t1[2] not in (0, 1): continue
        for c2 in (0, 1):
            t2 = mt.step(t1, c2, 1)
            if t2[2] not in (0, 1): continue
            for c3 in (0, 1):
                t3 = mt.step(t2, c3, 2)
                if t3[2] not in (0, 1): continue
                res[t3] = res.get(t3, 0) + 1
    return res
def can_accept(s, depth=40):
    """co-reachable: some continuation (any number of letters, phases tracked from 0) ends with E = 0"""
    fr = {(0, s)}
    seen = set(fr)
    for _ in range(depth):
        if any(x[1][2] == 0 for x in fr):
            return True
        nf = set()
        for ph, x in fr:
            for y in nxt(x, ph):
                key = ((ph + 1) % 3, y)
                if key not in seen:
                    seen.add(key); nf.add(key)
        fr = nf
        if not fr:
            break
    return False
live = [s for s in boundary if can_accept(s)]
print("states occurring at some position:", len(used), "; at block boundaries:", boundary)
print("dead block-boundary states:", [s for s in boundary if s not in live])
idx = {s: i for i, s in enumerate(live)}
T = sp.zeros(len(live), len(live))
for s in live:
    for t, c in block(s).items():
        if t in idx: T[idx[s], idx[t]] += c
x = sp.symbols('x')
print("trim states:", live)
print("trim T =", T.tolist())
print("characteristic polynomial:", sp.factor(T.charpoly(x).as_expr()))
J = T.jordan_form(calc_transform=False)
print("Jordan form diagonal/superdiagonal:", [J[i, i] for i in range(len(live))], [J[i, i + 1] for i in range(len(live) - 1)])
