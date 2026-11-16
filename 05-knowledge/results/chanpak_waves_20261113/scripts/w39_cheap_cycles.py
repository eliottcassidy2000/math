# Wave 39 (mac-mini-2026-11-15-frame). Run from a working directory; w39_lit_setup.py writes lit_base.pkl used by the others.
# Cheap-cycle dictionary: every primitive periodic literal word p of period <= PMAX whose per-letter growth (spectral
# radius of M_p in the minimal representation) is below RHO_CUT is a bottleneck (its k_p-th power has rank <= 1).
# Add all factors of the bi-infinite word p^Z of length < k_p |p| (the memory needed to cut the cycle).
# usage: python3 cheap_dict.py PMAX RHO_CUT
import pickle, sys, itertools, numpy as np
B = pickle.load(open('lit_base.pkl', 'rb'))
A, O, Rb = B['A'], B['O'], B['Rb']
PMAX, RHO = int(sys.argv[1]), float(sys.argv[2])
OR = O @ Rb
pinv = np.linalg.pinv(OR.astype(float))
Mm = [(O @ A[l] @ Rb).astype(float) @ pinv for l in range(8)]


def rank_le1(w):
    M = O.copy()
    for l in w[::-1]:
        M = M @ A[l]
    Mx = M @ Rb
    if Mx.max() > 2 ** 60:
        raise OverflowError
    cols = Mx[:, Mx.any(axis=0)]
    if cols.shape[1] <= 1:
        return True
    g = np.gcd.reduce(cols, axis=0)
    P = cols // g
    return bool((P == P[:, :1]).all())


D = set()
cyc = []
notbf = 0
for n in range(1, PMAX + 1):
    for w in itertools.product(range(8), repeat=n):
        if min(w[i:] + w[:i] for i in range(n)) != w or any(w == w[d:] + w[:d] for d in range(1, n) if n % d == 0):
            continue
        P = np.identity(Mm[0].shape[0])
        for l in w:
            P = Mm[l] @ P
        g = max(abs(np.linalg.eigvals(P))) ** (1.0 / n)
        if g >= RHO:
            continue
        # smallest length of a rank <= 1 factor of p^Z (check windows starting at every rotation)
        best = None
        for L in range(1, 40):
            ww = list(w) * (L // n + 2)
            if any(rank_le1(tuple(ww[i:i + L])) for i in range(n)):
                best = L
                break
        if best is None:
            notbf += 1
            print('CHEAP CYCLE WITHOUT A SHORT BOTTLENECK', w, g)
            continue
        cyc.append((g, w, best))
        ww = list(w) * (best // n + 2)
        for i in range(n):
            for L in range(1, best):
                D.add(tuple(ww[i:i + L]))
cyc.sort()
print('period <= %d, growth < %.3f: %d cheap primitive cycles; %d without a bottleneck of length < 40' % (PMAX, RHO, len(cyc), notbf))
print('rank-one lengths: max %d; distribution %s' % (max(c[2] for c in cyc), sorted({c[2] for c in cyc})))
print('dictionary strings %d (max length %d)' % (len(D), max(map(len, D))))
pickle.dump({'D': D, 'cycles': cyc}, open('cheapdict_P%d_%.3f.pkl' % (PMAX, RHO), 'wb'))
