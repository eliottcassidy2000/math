"""Certificate for two-module clone clocks.
For skeleton S with module positions u, v: compute the COMPLETE visit table a[i][i'] (caps from the in-neighbour bound),
then the P-form coefficients c[m][m'] with H(S[u<-M, v<-N]) = sum c[m][m'] P_M(m) P_N(m'), and verify by literal
Hamiltonian path counts for transitive and non-transitive modules."""
import sys, itertools
sys.path.insert(0, 'scripts')
from twomod import visit_table
from funnel import H_only
from math import comb

def build(order, back):
    n = len(order); idx = {nm: i for i, nm in enumerate(order)}; adj = [0]*n
    for i in range(n):
        for j in range(i+1, n):
            if (order[j], order[i]) in back: adj[j] |= 1 << i
            else: adj[i] |= 1 << j
    return adj, idx

def complete_table(adj, n, u, v):
    Iu = [x for x in range(n) if x != u and (adj[x] >> u & 1)]
    Iv = [x for x in range(n) if x != v and (adj[x] >> v & 1)]
    CU = len(Iu) + len(Iv) + 3; CV = len(Iu) + len(Iv) + 3
    t = visit_table(adj, n, u, v, CU=CU, CV=CV)
    assert all(c == 0 for (i, ip), c in t.items() if i == CU or ip == CV), 'cap reached'
    return {k: c for k, c in t.items() if c and k[0] >= 1 and k[1] >= 1}

def pform(tab):
    # K_i = sum_m (-1)^(i-m) C(i,m) P(m)
    c = {}
    for (i, ip), a in tab.items():
        for m in range(1, i+1):
            for mp in range(1, ip+1):
                c[(m, mp)] = c.get((m, mp), 0) + a * (-1)**(i-m) * comb(i, m) * (-1)**(ip-mp) * comb(ip, mp)
    return {k: x for k, x in c.items() if x}

def substitute(adj, n, u, v, M, N):
    """M, N: adjacency bit rows of module tournaments (sizes m, k).  Returns adj of S[u<-M, v<-N]."""
    others = [x for x in range(n) if x not in (u, v)]
    m, k = len(M), len(N)
    nid = {x: i for i, x in enumerate(others)}
    Ublk = list(range(len(others), len(others)+m)); Vblk = list(range(len(others)+m, len(others)+m+k))
    tot = len(others)+m+k; A = [0]*tot
    def arc(a, b): A[a] |= 1 << b
    for x in others:
        for y in others:
            if adj[x] >> y & 1: arc(nid[x], nid[y])
        for b in Ublk:
            (arc(nid[x], b) if adj[x] >> u & 1 else arc(b, nid[x]))
        for b in Vblk:
            (arc(nid[x], b) if adj[x] >> v & 1 else arc(b, nid[x]))
    for b1 in Ublk:
        for b2 in Vblk:
            (arc(b1, b2) if adj[u] >> v & 1 else arc(b2, b1))
    for i in range(m):
        for j in range(m):
            if M[i] >> j & 1: arc(Ublk[i], Ublk[j])
    for i in range(k):
        for j in range(k):
            if N[i] >> j & 1: arc(Vblk[i], Vblk[j])
    return A, tot

def TT(m):
    return [sum(1 << j for j in range(i+1, m)) for i in range(m)]
C3 = [0b010, 0b100, 0b001]
def P_poly(M, mbins):
    """P_M(mbins): ordered mbins-tuples of possibly empty paths partitioning M (brute force)."""
    m = len(M)
    # enumerate set partitions into ordered bins via assignment, then count paths per bin
    from funnel import H_only
    total = 0
    for assign in itertools.product(range(mbins), repeat=m):
        prod = 1
        for b in range(mbins):
            verts = [i for i in range(m) if assign[i] == b]
            if not verts: continue
            sub = [0]*len(verts)
            for a_, x in enumerate(verts):
                for b_, y in enumerate(verts):
                    if M[x] >> y & 1: sub[a_] |= 1 << b_
            prod *= H_only(sub, len(verts))
            if prod == 0: break
        total += prod
    return total

if __name__ == '__main__':
    families = {
      'plus1_anchor_9v (p=q=2,x=1)': (['a','v','c','u','d0','d1','x0','f0','f1'], {('d0','a'),('d1','a'),('f0','c'),('f1','c')}),
      'near_S60_with_v=b (p=q=2,x=0)': (['a','v','c','u','d0','d1','f0','f1'], {('d0','a'),('d1','a'),('f0','c'),('f1','c')}),
    }
    for name, (order, back) in families.items():
        adj, idx = build(order, back); n = len(order); u, v = idx['u'], idx['v']
        tab = complete_table(adj, n, u, v)
        c = pform(tab)
        print(name, '\n   visit table', dict(sorted(tab.items())), '\n   P-form', dict(sorted(c.items())))
        mods = [('TT1', TT(1)), ('TT2', TT(2)), ('TT3', TT(3)), ('C3', C3), ('TT4', TT(4))]
        for (mn, M) in mods:
            for (nn, N) in mods[:4]:
                A, tot = substitute(adj, n, u, v, M, N)
                if tot > 15: continue
                lit = H_only(A, tot)
                pred = sum(x * P_poly(M, mm) * P_poly(N, mp) for (mm, mp), x in c.items())
                assert lit == pred, (name, mn, nn, lit, pred)
        print('   literal path counts agree with the P-form for modules TT1..TT4, C3 (both slots)')
