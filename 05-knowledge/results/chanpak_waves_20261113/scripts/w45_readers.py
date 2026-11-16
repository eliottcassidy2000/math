"""Top-specific Horner readers.  For letters (u0,u1), reader rho and top sigma0, the series
F(w) = rho M_w sigma0 is Horner (F(wu) = lam F(w) + d_u for all w) iff sigma0 is annihilated by
O = span{ rho (M_u - lam)(M_v - 1) M_w }.  Search physical readers with O deficient on the reachable space,
then physical tops in its annihilator with a nonzero digit."""
import itertools, sys, numpy as np
from w45_kmid import KApex
P = 67108859

def rref(rows, N):
    B = np.zeros((0, N), dtype=np.int64); piv = []
    for r in rows:
        r = np.array(r, dtype=np.int64) % P
        for i, pc in enumerate(piv):
            if r[pc]: r = (r - r[pc] * B[i]) % P
        nz = np.nonzero(r)[0]
        if len(nz):
            pc = int(nz[0]); r = (r * pow(int(r[pc]), P - 2, P)) % P
            for i in range(len(piv)):
                if B[i][pc]: B[i] = (B[i] - B[i][pc] * r) % P
            B = np.vstack([B, r]); piv.append(pc)
    return B, piv

def words(alpha, maxlen):
    for L in range(maxlen + 1):
        for w in itertools.product(alpha, repeat=L):
            yield w

def run(k, lam, maxbottom, maxtop):
    K = KApex(k); N = K.N
    Ms = {b: (np.array(K.M[b], dtype=np.int64) % P) for b in K.M}
    init = np.array([int(x) for x in K.init()], dtype=np.int64)
    # physical tops
    tops = []
    for w in words(list(Ms), maxtop):
        v = init.copy()
        for b in w: v = Ms[b].dot(v) % P
        tops.append((w, v))
    Tmat = np.array([v for _, v in tops], dtype=np.int64).T        # N x ntops
    pairs = [(X, Y) for X in range(k) for Y in range(X + 1, k)]
    arcsets = []
    for bits in range(1 << len(pairs)):
        arcs = {}
        for t, (X, Y) in enumerate(pairs):
            arcs[(X, Y)] = bits >> t & 1; arcs[(Y, X)] = 1 - (bits >> t & 1)
        arcsets.append(arcs)
    readers = []
    for arcs in arcsets:
        acc = np.array([int(x) for x in K.accept(lambda X, Y: arcs[(X, Y)])], dtype=np.int64)
        for w in words(list(Ms), maxbottom):
            r = acc.copy()
            for b in reversed(w): r = r.dot(Ms[b]) % P
            readers.append(((w, tuple(sorted(arcs.items()))), r))
    I = np.eye(N, dtype=np.int64)
    found = 0; deficient = 0
    for u0, u1 in itertools.combinations(list(Ms), 2):
        A0 = (Ms[u0] - lam * I) % P; A1 = (Ms[u1] - lam * I) % P
        B0 = (Ms[u0] - I) % P; B1 = (Ms[u1] - I) % P
        for rkey, rho in readers:
            psi = [rho.dot(A0) % P, rho.dot(A1) % P]
            gens = [p.dot(B) % P for p in psi for B in (B0, B1)]
            O, _ = rref(gens, N)
            frontier = list(O)
            while frontier and len(O) < N:
                new = []
                for g in frontier:
                    for M in (Ms[u0], Ms[u1]):
                        h = g.dot(M) % P
                        O2, _ = rref(list(O) + [h], N)
                        if len(O2) > len(O): O = O2; new.append(h)
                frontier = new
            if len(O) == 0:
                good = np.ones(Tmat.shape[1], dtype=bool)
            else:
                good = ~((O.dot(Tmat) % P).any(axis=0))
            if not good.any(): continue
            deficient += 1
            dig = ((psi[0] - psi[1]) % P).dot(Tmat) % P
            ok = np.nonzero(good & (dig != 0))[0]
            if len(ok):
                found += 1
                if found <= 6:
                    t = tops[ok[0]][0]
                    print(f"   HORNER: letters {u0},{u1} reader {rkey} top {t}")
    print(f"k={k} lam={lam}: readers {len(readers)}, tops {len(tops)}; (pair,reader) with some good top: {deficient}; with a nonzero digit: {found}")

if __name__ == '__main__':
    run(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))
