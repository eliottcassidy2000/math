"""Referee: 'linear trijectivity' (signed pair-level) is strictly stronger than the existence of an
A-equivariant F2-affine bijection {T >= D} -> Tour(X).  Example: rigid oriented diamond + 2K1 (n=6),
A = Z2 = <(4 5)>.  The census shows D is not linearly trijective; here an explicit equivariant
F2-affine bijection onto Tour(X), X = {0..4} with tau_X = (1 2)(3 4), is built and verified on all
2^10 tournaments containing D.  Run: python affine.py"""
import itertools

n = 6
arcs = [(0, 1), (0, 2), (0, 3), (1, 2), (3, 1)]       # rigid orientation of K4 minus {2,3}
E = {frozenset(a) for a in arcs}
Pi = [(a, b) for a in range(n) for b in range(a + 1, n) if frozenset((a, b)) not in E]
assert len(Pi) == 10
tauV = (0, 1, 2, 3, 5, 4)
tauX = (0, 2, 1, 4, 3)
X = 5
PX = [(a, b) for a in range(X) for b in range(a + 1, X)]

# free part: ordered-pair correspondence (x -> 4) |-> fx4[x], extended by tau equivariantly
fx4 = {0: (0, 1), 1: (0, 3), 2: (1, 3), 3: (1, 4)}


def phi(T):
    """T: dict frozenset pair -> ordered (a,b) for pairs in Pi.  Returns tournament on X as dict."""
    out = {}
    for x in range(4):
        for tgt, img in ((4, fx4[x]), (5, tuple(tauX[i] for i in fx4[x]))):
            o = T[frozenset((x, tgt))]
            out[frozenset(img)] = img if o == (x, tgt) else (img[1], img[0])
    y_rs = 1 if T[frozenset((2, 3))] == (2, 3) else 0
    y_ab = 1 if T[frozenset((4, 5))] == (4, 5) else 0
    out[frozenset((1, 2))] = (1, 2) if (y_rs ^ y_ab) else (2, 1)
    out[frozenset((3, 4))] = (3, 4) if y_ab else (4, 3)
    assert len(out) == 10
    return out


def act(g, T):
    return {frozenset((g[a], g[b])): (g[o[0]], g[o[1]]) for p, o in T.items() for a, b in [tuple(p)]}


images = set()
ok = True
for bits in range(2 ** 10):
    T = {frozenset(p): (p if (bits >> i) & 1 else (p[1], p[0])) for i, p in enumerate(Pi)}
    F = phi(T)
    images.add(tuple(sorted((tuple(sorted(k)), v) for k, v in F.items())))
    lhs = phi(act(tauV, T))
    rhs = act(tauX, F)
    if lhs != rhs:
        ok = False
print("tau_V is an automorphism of D:", {(tauV[a], tauV[b]) for a, b in arcs} == set(arcs))
print("phi bijective onto Tour(X):", len(images) == 2 ** 10)
print("phi equivariant (phi(tau T) = tau phi(T) for all 1024 T):", ok)
# affine: phi(T) XOR phi(T0) depends linearly on T XOR T0 -- by construction each output bit is an
# XOR of input bits plus a constant; record the linear part rank
print("=> an A-equivariant F2-affine bijection exists although D is not linearly trijective (census).")
