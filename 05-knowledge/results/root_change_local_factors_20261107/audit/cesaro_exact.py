"""cesaro_exact.py - auditor's EXACT (rational) computation of the long-block limits of the letter-count model.

Generic finite-Markov-chain method, no structural assumption: the O/Z block chain P (letters O, Z with probability
1/2 each) on (Z/q)^4 is decomposed into strongly connected components; for each closed class C the stationary law
pi_C is solved exactly (Fractions); absorption probabilities h_C(s) of the transient states are solved exactly in
reverse topological order. The Cesaro projector is then Pi delta_s = sum_C h_C(s) pi_C (valid for periodic classes
too). It outputs the class law of F = D - E + U mod q for
    two J:  Pi J Pi J Pi delta_t      one J:  Pi J Pi delta_t,      t = (1,0,2,1),
and compares with the closed forms proposed for primes p (coordinator addenda):
    two J:  L(0) = 1 + (3p-2)(p-1)/(2p^3),  L(r) = 1 - (3p-2)/(2p^3) (r != 0)
    one J:  L(0) = 1 + (3p-4)/(4p^2),  L(-1) = 1 - 1/p^2 + 1/(4p),  L(r) = 1 - 1/p^2 otherwise.
It also checks the block-mixing lemma (M1)-(M3) on the exact class structure:
    (M1) for d != 0 the states {E in {0,d}} with that d form ONE closed class, aperiodic (self-loop), with uniform
         stationary law (so (D,U) uniform and independent of E in {0,d});
    (M2) for d = 0 the closed classes with U != 0 are {(D,0,U,0): D any, U in a doubling orbit}, uniform law;
    (M3) for d = U = 0 the closed states are the fixed points (D,0,0,0), and from a post-J state (a = 0) the first
         letter sends D to D-E+U (Z) or keeps D = E (O).
usage: python3 cesaro_exact.py q [q ...]
"""
import sys
from fractions import Fraction as Fr

sys.setrecursionlimit(100000)


def run(q):
    S = q ** 4
    enc = lambda D, E, U, d: ((D * q + E) * q + U) * q + d
    dec = lambda s: (s // q ** 3, (s // q ** 2) % q, (s // q) % q, s % q)
    tO = [0] * S; tZ = [0] * S; tJ = [0] * S; Fv = [0] * S
    for s in range(S):
        D, E, U, d = dec(s)
        tO[s] = enc(D, 0, (2 * U) % q, d)
        tZ[s] = enc((D - E + U) % q, d, (2 * U - d) % q, d)
        tJ[s] = enc((D - E + U) % q, (D - E) % q, (D + U) % q, D)
        Fv[s] = (D - E + U) % q
    succ = lambda s: (tO[s], tZ[s])
    # iterative Tarjan SCC
    index = [-1] * S; low = [0] * S; onst = [False] * S; st = []; comp = [-1] * S; comps = []; idx = 0
    for root in range(S):
        if index[root] != -1:
            continue
        work = [(root, 0)]
        while work:
            v, i = work.pop()
            if i == 0:
                index[v] = low[v] = idx; idx += 1; st.append(v); onst[v] = True
            nbrs = succ(v)
            if i < 2:
                w = nbrs[i]
                work.append((v, i + 1))
                if index[w] == -1:
                    work.append((w, 0))
                elif onst[w]:
                    low[v] = min(low[v], index[w])
                continue
            for w in nbrs:
                if comp[w] == -1 and onst[w]:
                    low[v] = min(low[v], low[w])
            if low[v] == index[v]:
                c = []
                while True:
                    w = st.pop(); onst[w] = False; comp[w] = len(comps); c.append(w)
                    if w == v:
                        break
                comps.append(c)
    closed = [ci for ci, c in enumerate(comps) if all(comp[w] == ci for v in c for w in succ(v))]
    cset = set(closed)
    # stationary laws of closed classes, exact Gaussian elimination
    pi = {}
    for ci in closed:
        c = sorted(comps[ci]); n = len(c); pos = {v: k for k, v in enumerate(c)}
        if n == 1:
            pi[ci] = {c[0]: Fr(1)}; continue
        # equations: pi(x) = sum_y pi(y) P(y,x); replace last by sum = 1
        A = [[Fr(0)] * (n + 1) for _ in range(n)]
        for k in range(n):
            A[k][k] -= 1
        for y in c:
            for x in succ(y):
                A[pos[x]][pos[y]] += Fr(1, 2)
        A[n - 1] = [Fr(1)] * n + [Fr(1)]
        for col in range(n):
            piv = next(r for r in range(col, n) if A[r][col] != 0)
            A[col], A[piv] = A[piv], A[col]
            inv = 1 / A[col][col]
            A[col] = [x * inv for x in A[col]]
            for r in range(n):
                if r != col and A[r][col] != 0:
                    f = A[r][col]
                    A[r] = [a - f * b for a, b in zip(A[r], A[col])]
        pi[ci] = {c[k]: A[k][n] for k in range(n)}
    # absorption probabilities: comps are produced by Tarjan in reverse topological order (sinks first)
    h = [None] * S
    for ci, c in enumerate(comps):
        if ci in cset:
            for v in c:
                h[v] = {ci: Fr(1)}
            continue
        if len(c) == 1 and c[0] not in succ(c[0]):
            v = c[0]; acc = {}
            for w in succ(v):
                for k, x in h[w].items():
                    acc[k] = acc.get(k, 0) + x / 2
            h[v] = acc; continue
        raise RuntimeError('non-trivial transient SCC of size %d (not expected; solve generically)' % len(c))

    def Pi(mu):
        out = {}
        for s, m in mu.items():
            for ci, x in h[s].items():
                for v, pv in pi[ci].items():
                    out[v] = out.get(v, 0) + m * x * pv
        return out

    def Jop(mu):
        out = {}
        for s, m in mu.items():
            out[tJ[s]] = out.get(tJ[s], 0) + m
        return out

    def Flaw(mu):
        cls = [Fr(0)] * q
        for s, m in mu.items():
            cls[Fv[s]] += m
        tot = sum(cls)
        return [q * c / tot for c in cls]

    t = enc(1 % q, 0, 2 % q, 1 % q)
    mu0 = Pi({t: Fr(1)})
    one = Flaw(Pi(Jop(mu0)))
    two = Flaw(Pi(Jop(Pi(Jop(mu0)))))
    sizes = sorted(len(comps[ci]) for ci in closed)
    print('q=%d: %d SCCs, %d closed classes (sizes %s), transient SCCs all singletons without self-loop: %s' % (
        q, len(comps), len(closed), {sz: sizes.count(sz) for sz in sorted(set(sizes))},
        all(len(c) == 1 for ci, c in enumerate(comps) if ci not in cset)))
    print('  two-J limit L(r), r = 0..%d: %s' % (q - 1, ' '.join(str(x) for x in two)))
    print('  one-J limit L(r), r = 0..%d: %s' % (q - 1, ' '.join(str(x) for x in one)))
    # closed forms for primes
    p = q
    if all(q % k for k in range(2, q)):
        c2 = [1 + Fr((3 * p - 2) * (p - 1), 2 * p ** 3)] + [1 - Fr(3 * p - 2, 2 * p ** 3)] * (p - 1)
        c1 = [1 + Fr(3 * p - 4, 4 * p ** 2)] + [1 - Fr(1, p ** 2)] * (p - 1)
        c1[p - 1] = 1 - Fr(1, p ** 2) + Fr(1, 4 * p)
        print('  closed form two-J: %s  -> EXACT MATCH: %s' % (c2[0], two == c2))
        print('  closed form one-J: %s, L(-1) = %s -> EXACT MATCH: %s' % (c1[0], c1[p - 1], one == c1))
        # (M1)-(M3) on the exact structure
        okM1 = True; okM2 = True; okM3 = True
        for ci in closed:
            c = comps[ci]; ds = {dec(v)[3] for v in c}
            if len(ds) != 1:
                okM1 = False; continue
            d = ds.pop()
            if d:
                want = {enc(D, E, U, d) for D in range(q) for U in range(q) for E in (0, d)}
                okM1 &= set(c) == want and all(x == Fr(1, len(c)) for x in pi[ci].values())
                okM1 &= any(v in succ(v) for v in c)          # self-loop: aperiodic
            else:
                Us = {dec(v)[2] for v in c}
                if Us == {0}:
                    okM3 &= len(c) == 1 and dec(c[0])[1] == 0 and tO[c[0]] == c[0] and tZ[c[0]] == c[0]
                else:
                    okM2 &= all(dec(v)[1] == 0 for v in c) and len(c) == q * len(Us) and all(x == Fr(1, len(c)) for x in pi[ci].values())
                    okM2 &= {(2 * u) % q for u in Us} == Us and 0 not in Us
        # every d != 0 slice is a single closed class (no other closed class with that d)
        dcount = {}
        for ci in closed:
            d = dec(comps[ci][0])[3]; dcount[d] = dcount.get(d, 0) + 1
        okM1 &= all(dcount.get(d, 0) == 1 for d in range(1, q))
        # (M3) first letter from post-J confined states: J(p) with D_p = U_p = 0 mod p
        for s in range(S):
            D, E, U, d = dec(s)
            if D == 0 and U == 0:
                j = tJ[s]; Dj, Ej, Uj, dj = dec(j)
                okM3 &= (Dj - Ej - Uj + dj) % q == 0 and dec(tZ[j])[0] == (Dj - Ej + Uj) % q == 0 and dec(tO[j])[0] == Ej
        print('  (M1) exact: %s   (M2) exact: %s   (M3) exact: %s' % (okM1, okM2, okM3))


for a in sys.argv[1:]:
    run(int(a))
