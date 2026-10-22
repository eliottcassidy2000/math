#!/usr/bin/env python3
"""(d) Leech lattice facts: minimal vectors (explicit construction), 11-design test, inner products,
Voronoi facet count via Lambda/2Lambda, weighted 15-design ratio, Co0 order."""
import itertools, math
import numpy as np
import sympy as sp
import flint
# extended Golay code from the QR code of length 23
p = 23
f = flint.nmod_poly([1] + [0]*(p-1) + [1], 2)
g = [fa for fa, e in f.factor()[1] if fa.degree() == 11][0]
gc = [int(t) for t in g.coeffs()]
gens = []
for s in range(12):
    row = [0]*p
    for i, c in enumerate(gc): row[(i+s) % p] = c
    row.append(sum(row) % 2); gens.append(row)
G = np.array(gens, dtype=np.int64)
words = (np.array(list(itertools.product((0, 1), repeat=12))) @ G) % 2
octads = words[words.sum(1) == 8]
assert len(octads) == 759
vecs = []
# shape (+-4, +-4, 0^22): 1104
for i, j in itertools.combinations(range(24), 2):
    for si, sj in itertools.product((4, -4), repeat=2):
        v = np.zeros(24, np.int64); v[i] = si; v[j] = sj; vecs.append(v)
# shape (+-2^8 on an octad, even number of minus signs): 759 * 128 = 97152
for o in octads:
    idx = np.nonzero(o)[0]
    for signs in itertools.product((1, -1), repeat=8):
        if signs.count(-1) % 2 == 0:
            v = np.zeros(24, np.int64); v[idx] = 2*np.array(signs); vecs.append(v)
# shape (-3, 1^23) modified by Golay codewords: (1 - 2c) * ... standard: x = (-3,1,...,1) with signs flipped on a codeword, all 24 positions: 4096*24 = 98304
for c in words:
    for i in range(24):
        v = np.ones(24, np.int64); v[i] = -3
        v = v * (1 - 2*c)
        vecs.append(v)
V = np.array(vecs, dtype=np.int64)          # norms 32 = 8*4 (scaled by sqrt8)
print("number of vectors:", len(V), " all norm 32 (=4 after /8):", bool((np.einsum('ij,ij->i', V, V) == 32).all()))
# check lattice membership-type sanity: inner products / 8 integral and in {-4,-2,-1,0,1,2,4}
rng = np.random.default_rng(1)
sample = V[rng.choice(len(V), 300, replace=False)]
ips = set(((sample @ V.T) // 8).flatten().tolist()) if ((sample @ V.T) % 8 == 0).all() else None
print("inner products (norm-4 scaling) between minimal vectors (sample):", sorted(ips))
# spherical design test: for unit u, (1/N) sum (x.u)^{2k}/|x|^{2k} vs sphere moment (2k-1)!!/(24*26*...*(24+2k-2))
X = V / np.sqrt(32.0)
def sphere_moment(k, n=24):
    num = math.prod(range(1, 2*k, 2)); den = math.prod(n + 2*j for j in range(k))
    return num/den
for trial in range(3):
    u = rng.normal(size=24); u /= np.linalg.norm(u)
    t = X @ u
    rel = [ (np.mean(t**(2*k)) / sphere_moment(k)) - 1 for k in range(1, 8)]
    print(f"  random direction {trial}: relative deviation of moments k=1..7 (degree 2..14):", ["%.1e" % r for r in rel])
print("  -> degrees <= 10 match the sphere (11-design); degree 12 deviates; degree 14 (k=7) mixes deg-12 harmonic part")
# Voronoi facets via Lambda/2Lambda
N4, N6, N8 = 196560, 16773120, 398034000
print("Lambda/2Lambda classes: 1 + N4/2 + N6/2 + N8/48 =", 1 + N4//2 + N6//2 + N8//48, "= 2^24:", 1 + N4//2 + N6//2 + N8//48 == 2**24)
print("=> relevant vectors (Voronoi facets) = N4 + N6 =", N4 + N6)
# weighted 15-design ratio from theta_{P} = c Delta^2 (deg-12 harmonics): S3 sum = -48 * S2 sum
print("w3/w2 = (6/4)^6/48 =", sp.Rational(6, 4)**6/48, "= (3/4)^5:", sp.Rational(6, 4)**6/48 == sp.Rational(3, 4)**5)
co0 = 8315553613086720000
print("Co0 order factorization:", sp.factorint(co0), " = 2 * |Co1| with |Co1| = 4157776806543360000:", co0 == 2*4157776806543360000)
print("|Co0|/196560 = |Co2| =", co0 // 196560, " (Co2 order 42305421312000:", co0 // 196560 == 42305421312000, ")")
