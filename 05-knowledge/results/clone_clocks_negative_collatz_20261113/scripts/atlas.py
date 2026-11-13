"""Outer-tail atlas: for nested-funnel cores (prefix words) and outer words L b^k R,
the ready count is H(k) = A*2^k + B (k>=1).  Collect (A,B) and look for C2 chains
(A,B) -> (3A, 3B+2) -> ... and the C2-fixed line B = -1."""
import sys, itertools, collections
sys.path.insert(0, 'scripts')
from funnel import step6, word_polys, SINGLETON, ready_count

def words_upto(t, need_end1=True):
    out = []
    for L in range(0, t+1):
        for bits in itertools.product('01', repeat=L):
            w = ''.join(bits)
            if L and need_end1 and w[-1] != '1': continue
            out.append(w)
    return out

def cores(max_vertices, max_word=4):
    """nested word sequences from the singleton (each word adds len+1 vertices)."""
    res = {(): (SINGLETON, 1)}
    frontier = [((), SINGLETON, 1)]
    ws = words_upto(max_word)
    while frontier:
        nf = []
        for seq, st, nv in frontier:
            for w in ws:
                if w == '' and len(seq) == 0: continue   # empty word only on nontrivial core
                nv2 = nv + len(w) + 1
                if nv2 > max_vertices: continue
                st2 = step6(st, w)
                seq2 = seq + (w,)
                res[seq2] = (st2, nv2)
                nf.append((seq2, st2, nv2))
        frontier = nf
    return res

def AB(state, Lw, b, Rw):
    h, e, s2, p2, e2, s3 = state
    vals = []
    for k in (1, 2, 3):
        w = Lw + b*k + Rw
        L1, D1, B1 = word_polys(w, 1)
        L2, D2, B2 = word_polys(w, 2)
        vals.append(L2*(s2+p2) + B1*h - D1*e - e2 - 1)
    A = vals[1] - vals[0]; A //= 2
    B = vals[0] - 2*A
    assert A*8 + B == vals[2], (vals, A, B)
    return A, B

if __name__ == '__main__':
    maxv = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    maxctx = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    C = cores(maxv)
    print('cores', len(C))
    ctx = []
    for l in range(0, maxctx+1):
        for r in range(0, maxctx+1-l):
            for Lb in itertools.product('01', repeat=l):
                for Rb in itertools.product('01', repeat=r):
                    for b in '01':
                        Lw, Rw = ''.join(Lb), ''.join(Rb)
                        if not (Rw.endswith('1') or (Rw == '' and b == '1')): continue
                        ctx.append((Lw, b, Rw))
    print('contexts', len(ctx))
    table = collections.defaultdict(list)
    for seq, (st, nv) in C.items():
        if seq == ():
            pass
        for (Lw, b, Rw) in ctx:
            A, B = AB(st, Lw, b, Rw)
            key = (A, B)
            if key not in table or len(table[key]) < 3:
                table[key].append((nv + len(Lw) + len(Rw) + 2, seq, Lw, b, Rw))
    print('distinct (A,B)', len(table))
    import pickle
    pickle.dump(dict(table), open(f'atlas_{maxv}_{maxctx}.pkl', 'wb'))
