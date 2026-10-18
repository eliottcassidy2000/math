#!/usr/bin/env python3
"""HYP-9293 side claim: the symbol (top-length part) of the singular vector of V^k(sl2),
k = -2 + a/b, is proportional to e^(a-1) Cas^((a-1)(b-1)/2), Cas = e f + h^2/4.  Mod-P check.
Usage: python sl2_symbol_check.py a b [a b ...]"""
import sys
import flint
from alin import nullspace, rank
from sl2mod import AffSl2, colored_partitions, charge, pmul, E, H, F

P = 2147483647
one = flint.nmod(1, P)
args = list(map(int, sys.argv[1:])) or [3, 2, 2, 3, 3, 4, 5, 2, 4, 3, 2, 5]
for i in range(0, len(args), 2):
    a, b = args[i], args[i + 1]
    k = one * (a - 2 * b) / b
    V = AffSl2(k, one)
    D0, c0 = (a - 1) * b, 2 * (a - 1)
    Bs = [m for m in colored_partitions(D0) if charge(m) == c0]
    eqs = {}
    for (t, m) in [(E, 0), (E, 1), (H, 1), (F, 1)]:
        for bb in Bs:
            for mono, cc in V.act(t, m, bb).items():
                eqs.setdefault((t, m, mono), {})[bb] = cc
    ker = nullspace(list(eqs.values()), Bs, P)
    v = ker[0]
    sym = {m: c for m, c in v.items() if len(m) == D0}
    Cas = {((1, E), (1, F)): one, ((1, H), (1, H)): one / 4}
    target = {(): one}
    for _ in range(a - 1):
        target = pmul(target, {((1, E),): one})
    for _ in range((a - 1) * (b - 1) // 2):
        target = pmul(target, Cas)
    basis = sorted(set(sym) | set(target))
    ok = sym and rank([sym, target], basis, P) == 1
    print(f"k=-2+{a}/{b}: singular space dim {len(ker)}, #terms of v {len(v)}, "
          f"symbol proportional to e^{a-1} Cas^{(a-1)*(b-1)//2}: {bool(ok)}")
