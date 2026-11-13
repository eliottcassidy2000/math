"""Search the nested-funnel language for TWO-CLOCK ready families
     R(j,k) + 1 = 3^j * F(k)   (exact H -> 3H+2 under j -> j+1, for every k)
Construction: core prefix (nested words) -> inner funnel L1 b1^j R1 -> middle words -> outer funnel L2 b2^k R2 -> ready z.
Test: psi.(c1(j+1) - 3 c1(j)) = 0 for j=1..7, for psi in {phi_A M_mid, phi_B M_mid}."""
import sys, itertools, collections
import numpy as np
sys.path.insert(0, 'scripts')
from funnel import step6, word_polys, SINGLETON
from atlas import cores, words_upto

def contexts(maxctx):
    out = []
    for l in range(0, maxctx+1):
        for r in range(0, maxctx+1-l):
            for Lb in itertools.product('01', repeat=l):
                for Rb in itertools.product('01', repeat=r):
                    for b in '01':
                        Lw, Rw = ''.join(Lb), ''.join(Rb)
                        if not (Rw.endswith('1') or (Rw == '' and b == '1')): continue
                        out.append((Lw, b, Rw))
    return out

def outer_functionals(ctx):
    Lw, b, Rw = ctx
    rows = []
    for k in (1, 2):
        w = Lw + b*k + Rw
        L1, D1, B1 = word_polys(w, 1); L2, D2, B2 = word_polys(w, 2)
        rows.append(np.array([B1, -D1, L2, L2, -1, 0], dtype=object))
    # phi(k) = 2^k phiA + phiB  (exact for k>=1; checked at k=3)
    phiA = (rows[1] - rows[0]); phiA = np.array([x//2 for x in phiA], dtype=object)
    phiB = rows[0] - 2*phiA
    w = Lw + b*3 + Rw
    L1, D1, B1 = word_polys(w, 1); L2, D2, B2 = word_polys(w, 2)
    r3 = np.array([B1, -D1, L2, L2, -1, 0], dtype=object)
    assert all(r3 == 8*phiA + phiB)
    return phiA, phiB

def apply_words(state, words):
    for w in words: state = step6(state, w)
    return state

def mid_matrix(words):
    cols = []
    for i in range(6):
        e = [0]*6; e[i] = 1
        cols.append(apply_words(tuple(e), words))
    return np.array(cols, dtype=object).T   # M with M e_i = col i

if __name__ == '__main__':
    maxcore = int(sys.argv[1]); inctx = int(sys.argv[2]); midmax = int(sys.argv[3]); outctx = int(sys.argv[4])
    C = cores(maxcore)
    ICTX = contexts(inctx); OCTX = contexts(outctx)
    # middle word sequences: up to 2 words, total vertices <= midmax
    mids = [()]
    ws = words_upto(3)
    for w in ws:
        if len(w)+1 <= midmax: mids.append((w,))
    for w1 in ws:
        for w2 in ws:
            if len(w1)+len(w2)+2 <= midmax: mids.append((w1, w2))
    print('cores', len(C), 'inner ctx', len(ICTX), 'mids', len(mids), 'outer ctx', len(OCTX), flush=True)
    OUT = [(ctx,) + outer_functionals(ctx) for ctx in OCTX]
    MID = [(m, mid_matrix(m)) for m in mids]
    # psi list
    PSI = []
    for (m, Mm) in MID:
        for (ctx, pA, pB) in OUT:
            PSI.append((m, ctx, pA.dot(Mm), pB.dot(Mm)))
    PA = np.array([[int(x) for x in p[2]] for p in PSI], dtype=np.int64)
    PB = np.array([[int(x) for x in p[3]] for p in PSI], dtype=np.int64)
    print('psi pairs', len(PSI), flush=True)
    found = []
    for seq, (c0, nv) in C.items():
        for (L1, b1, R1) in ICTX:
            if seq == () and False: pass
            cs = [np.array(step6(c0, L1 + b1*j + R1), dtype=object) for j in range(1, 10)]
            D = np.array([[int(x) for x in (cs[j+1] - 3*cs[j])] for j in range(0, 8)], dtype=np.int64).T  # 6 x 8
            ok = (np.abs(PA.dot(D)).max(axis=1) == 0) & (np.abs(PB.dot(D)).max(axis=1) == 0)
            idx = np.nonzero(ok)[0]
            for i in idx:
                m, ctx, a, b = PSI[i]
                # exact recheck with integers
                vals = []
                for j in range(1, 6):
                    st = apply_words(c0, [L1 + b1*j + R1] + list(m))
                    row = []
                    for k in (1, 2, 3):
                        st2 = step6(st, ctx[0] + ctx[1]*k + ctx[2])
                        row.append(st2[0] + st2[4] - 1)
                    vals.append(row)
                good = all(vals[j+1][kk] == 3*vals[j][kk] + 2 for j in range(4) for kk in range(3))
                if good:
                    found.append((seq, (L1, b1, R1), m, ctx, vals[0]))
                    print('FOUND', seq, (L1, b1, R1), m, ctx, 'R(1,k)=', vals[0], flush=True)
    print('total found', len(found))
