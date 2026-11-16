"""Middle-insertion Horner readers.  Digits = single sites (letters u) inserted between an arbitrary top
(past) and a fixed bottom reader rho.  Horner: f(w u) = lam f(w) + d_u  for every top and every w  <=>
rho (M_u - lam) (M_v - 1) r = 0 for letters u, v and r in the reachable space R0.  Search physical readers."""
import itertools, sys, numpy as np
from w45_kmid import KApex
P = 2147483629  # prime < 2^31

def modmat(M):
    return np.array([[int(x) % P for x in row] for row in M], dtype=np.int64)

def mm(A, B):  # mod-P matrix product, chunked to avoid overflow
    out = np.zeros((A.shape[0], B.shape[1]), dtype=np.int64)
    for j in range(A.shape[1]):
        out = (out + (A[:, j:j+1] * B[j:j+1, :]) % P) % P
    return out

def rowreduce(rows):
    """mod-P echelon basis of the row span"""
    basis = []; pivots = []
    for r in rows:
        r = r.copy() % P
        for b, pc in zip(basis, pivots):
            if r[pc]:
                r = (r - r[pc] * b) % P
        nz = np.nonzero(r)[0]
        if len(nz):
            pc = nz[0]; inv = pow(int(r[pc]), P - 2, P)
            r = (r * inv) % P
            for i in range(len(basis)):
                if basis[i][pc]:
                    basis[i] = (basis[i] - basis[i][pc] * r) % P
            basis.append(r); pivots.append(pc)
    return basis

def run(k, lam, maxbottom):
    K = KApex(k); N = K.N
    Ms = {b: modmat(K.M[b]) for b in K.M}
    I = np.eye(N, dtype=np.int64)
    # reachable space from init under all site types
    R = rowreduce([np.array([int(x) for x in K.init()], dtype=np.int64)])
    frontier = list(R)
    while frontier:
        new = []
        for r in frontier:
            for b in Ms:
                v = mm(Ms[b], r.reshape(-1, 1)).ravel()
                if len(rowreduce(R + [v])) > len(R):
                    R = rowreduce(R + [v]); new.append(v)
        frontier = new
    Rm = np.array(R, dtype=np.int64).T            # N x dimR
    print(f"k={k}: {N} states, reachable dim {Rm.shape[1]}")
    letters = list(Ms)
    cons = {}
    for u0, u1 in itertools.combinations(letters, 2):
        rows = []
        for u in (u0, u1):
            A = (Ms[u] - lam * I) % P
            for v in (u0, u1):
                B = mm((Ms[v] - I) % P, Rm)       # N x dimR
                C = mm(A, B)                       # N x dimR : constraint columns
                rows.extend(list(C.T))
        cons[(u0, u1)] = np.array(rowreduce(rows), dtype=np.int64)
    # nontriviality vectors: (M_u0 - M_u1) R
    diff = {pr: mm((Ms[pr[0]] - Ms[pr[1]]) % P, Rm) for pr in cons}
    hits = []
    arcsets = []
    pairs = [(X, Y) for X in range(k) for Y in range(X + 1, k)]
    for bits in range(1 << len(pairs)):
        arcs = {}
        for t, (X, Y) in enumerate(pairs):
            arcs[(X, Y)] = bits >> t & 1; arcs[(Y, X)] = 1 - (bits >> t & 1)
        arcsets.append(arcs)
    count = 0
    for arcs in arcsets:
        acc = np.array([int(x) for x in K.accept(lambda X, Y: arcs[(X, Y)])], dtype=np.int64)
        stack = [((), acc)]
        while stack:
            word, rho = stack.pop()
            count += 1
            for pr, Cm in cons.items():
                if len(Cm) == 0 or not ((Cm.dot(rho) % P) if False else mm(Cm, rho.reshape(-1, 1)).ravel()).any():
                    d = mm(rho.reshape(1, -1), diff[pr]).ravel()
                    if d.any():
                        hits.append((pr, word, tuple(sorted(arcs.items()))))
            if len(word) < maxbottom:
                for b in Ms:
                    stack.append(((b,) + word, mm(rho.reshape(1, -1), Ms[b]).ravel()))
    print(f"  lam={lam}: readers checked {count}; Horner readers with a nonzero digit: {len(hits)}")
    for h in hits[:8]:
        print("   letters", h[0], "bottom (top->bottom)", h[1], "apex arcs", h[2])
    return hits

if __name__ == '__main__':
    k = int(sys.argv[1]); lam = int(sys.argv[2]); mb = int(sys.argv[3])
    run(k, lam, mb)
