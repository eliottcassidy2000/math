"""Nuclear (trace) norm of the skew-adjacency matrix of tournaments: exact checks and the parse obstruction numbers."""
import numpy as np, math, sys, itertools, random
def skew(out, n):
    S = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i != j: S[i, j] = 1.0 if (out[i] >> j) & 1 else -1.0
    return S
def nuc(S): return float(np.abs(np.linalg.eigvalsh(1j * S)).sum())
def tt(n): return [sum(1 << j for j in range(i + 1, n)) for i in range(n)]
def paley(q):
    qr = {(x * x) % q for x in range(1, q)}
    return [sum(1 << j for j in range(q) if j != i and ((j - i) % q) in qr) for i in range(q)]
def nu_tt_formula(n): return sum(abs(1 / math.tan((2 * k - 1) * math.pi / (2 * n))) for k in range(1, n + 1) if (2 * k - 1) != n)
print("TT_n: nuclear norm by SVD vs cot formula; ratio to (2/pi) n ln n")
for n in [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 16, 20, 32, 50, 100, 200, 400]:
    a = nuc(skew(tt(n), n)); b = nu_tt_formula(n)
    print(f"  n={n:4d}  svd={a:12.6f}  formula={b:12.6f}  diff={a-b:.2e}  nu/(n ln n)={a/(n*math.log(n)):.5f}  (2/pi={2/math.pi:.5f})")
print("formula only, large n:")
for n in [1000, 10000, 100000]:
    b = nu_tt_formula(n); print(f"  n={n:7d} nu={b:.3f} nu/(n ln n)={b/(n*math.log(n)):.5f}  nu - (2/pi) n ln n = {b-(2/math.pi)*n*math.log(n):.3f}  per n: {(b-(2/math.pi)*n*math.log(n))/n:.5f}")
print("Paley P_q: nuclear norm vs (q-1) sqrt q")
for q in [3, 7, 11, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83]:
    S = skew(paley(q), q); a = nuc(S); print(f"  q={q:3d} nu={a:10.5f} (q-1)sqrt(q)={(q-1)*math.sqrt(q):10.5f}  S S^T == qI-J: {np.allclose(S @ S.T, q*np.eye(q)-np.ones((q,q)))}")
print("Bordered Paley (source + P_q), n=q+1: nu vs n sqrt(n-1)")
for q in [3, 7, 11, 19, 23]:
    n = q + 1; out = [sum(1 << j for j in range(1, n))] + [(m << 1) for m in paley(q)]
    S = skew(out, n); a = nuc(S); print(f"  n={n:3d} nu={a:10.5f} n sqrt(n-1)={n*math.sqrt(n-1):10.5f}  S S^T==(n-1)I: {np.allclose(S@S.T,(n-1)*np.eye(n))}")
# caterpillar C3[v, v, C3[v, v, ...]] : vertices 0,1 then rest; arcs 0->1, 1->X, X->0
def caterpillar(n):
    out = [0] * n
    def build(verts):
        if len(verts) == 1: return
        if len(verts) == 2:
            out[verts[0]] |= 1 << verts[1]; return
        a, b, X = verts[0], verts[1], verts[2:]
        out[a] |= 1 << b
        for x in X: out[b] |= 1 << x; out[x] |= 1 << a
        build(X)
    build(list(range(n))); return out
print("Caterpillar C3[v,v,C3[v,v,...]]: nu vs n ln n and n^1.5")
for n in [10, 20, 40, 80, 160, 320, 640]:
    a = nuc(skew(caterpillar(n), n)); print(f"  n={n:4d} nu={a:10.3f}  nu/(n ln n)={a/(n*math.log(n)):.4f}  nu/n^1.5={a/n**1.5:.4f}")
# C3 tower
def tower(k):
    n = 3 ** k; out = [0] * n
    def build(verts):
        if len(verts) == 1: return
        m = len(verts) // 3; A, B, C = verts[:m], verts[m:2*m], verts[2*m:]
        for X, Y in ((A, B), (B, C), (C, A)):
            for x in X:
                for y in Y: out[x] |= 1 << y
        build(A); build(B); build(C)
    build(list(range(n))); return out
print("C3 tower n=3^k: nu vs (2/sqrt3) n log_3 n")
for k in [1, 2, 3, 4, 5]:
    n = 3 ** k; a = nuc(skew(tower(k), n)); print(f"  k={k} n={n:4d} nu={a:10.4f} formula={(2/math.sqrt(3))*n*k:10.4f}")
# random tournaments: nu ~ ?
random.seed(1)
print("random tournaments: nu / ((n-1) sqrt n)")
for n in [11, 23, 47, 95, 191]:
    vals = []
    for _ in range(5):
        out = [0] * n
        for i in range(n):
            for j in range(i + 1, n):
                if random.random() < 0.5: out[i] |= 1 << j
                else: out[j] |= 1 << i
        vals.append(nuc(skew(out, n)) / ((n - 1) * math.sqrt(n)))
    print(f"  n={n:4d} mean ratio={np.mean(vals):.4f}  (semicircle prediction 8/(3 pi)={8/(3*math.pi):.4f})")
