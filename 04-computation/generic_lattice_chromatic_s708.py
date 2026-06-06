#!/usr/bin/env python3
"""
generic_lattice_chromatic_s708.py  (monad-explorer-2026-06-06-S708)

QUESTION (open handoff from THM-418): does a GENERIC 2D lattice (w=2, only +-1
symmetry) attain chi=4 at some norm?  THM-418 settled square (chi=2) and
triangular/Eisenstein (chi=3); generic lattices were left open ("small patches
found only <=3").

ANGLE (dispatched): the richest w=2 lattices "between the triangular lattice and
the CM field" are ORDERS in imaginary quadratic fields Q(sqrt(-d)): CM-like
multiplicative norm forms (popular norms, big connection sets) but w=2 (THM-416:
density quantum 1).  We test their unit-distance graphs U(L,D) at POPULAR norms.

A unit-distance graph U(L,D) = Cayley graph on Z^2 with connection set
  S = { v in Z^2 : Q(v) = D },  Q the lattice's integral quadratic form.

THEORY USED (Cayley graph Cay(Z^2,S), S symmetric):
  * chi=2 (bipartite) IFF exists p in {(1,0),(0,1),(1,1)} with p.s ODD mod 2 for all s in S.
  * chi<=3 SUFFICIENT if exists q in {(0,1),(1,0),(1,1),(1,2)} with q.s != 0 mod 3 for all s.
    (linear coloring c(v)=q.v mod 3; generalises the triangular (a-b) mod 3 trick.)
  * RIGOROUS chi>=4: a finite patch (subgraph) that is NOT 3-colorable.
  * RIGOROUS chi<=3: a proper 3-coloring of a torus Z_N^2 (periodic coloring of Z^2).

We combine these to PIN chi for each (lattice, norm).
"""
import sys, itertools

# ---------- quadratic forms / connection sets ----------
def Qval(A,B,C,x,y):
    return A*x*x + B*x*y + C*y*y

def connection_set(A,B,C,D):
    """All integer (x,y)!=0 with A x^2 + B xy + C y^2 = D.  Form must be pos.def."""
    # bound: Q >= lam_min*(x^2+y^2); lam_min = smallest eigenvalue of [[A,B/2],[B/2,C]]
    tr = A + C
    det = A*C - (B*B)/4.0
    disc = (tr*tr - 4*det)
    lam_min = (tr - disc**0.5)/2.0
    if lam_min <= 0:
        raise ValueError("form not positive definite")
    bound = int((D/lam_min)**0.5) + 2
    S = []
    for x in range(-bound, bound+1):
        for y in range(-bound, bound+1):
            if x==0 and y==0: continue
            if Qval(A,B,C,x,y)==D:
                S.append((x,y))
    return S

# ---------- Cayley-graph chromatic tests ----------
def is_bipartite(S):
    for p in [(1,0),(0,1),(1,1)]:
        if all(((p[0]*s[0]+p[1]*s[1]) & 1)==1 for s in S):
            return True
    return False

def linear_3col(S):
    for q in [(0,1),(1,0),(1,1),(1,2)]:
        if all(((q[0]*s[0]+q[1]*s[1]) % 3)!=0 for s in S):
            return True
    return False

# ---------- generic graph 3-colorability (backtracking, DSATUR-ish) ----------
def three_colorable(adj):
    """adj: list of sets (neighbor indices). Return True if 3-colorable."""
    n = len(adj)
    if n==0: return True
    color = [0]*n      # 0 = uncolored, 1..3
    # order vertices by degree desc for better pruning
    order = sorted(range(n), key=lambda v: -len(adj[v]))
    pos_of = [0]*n
    # iterative-ish recursion via explicit stack to avoid recursion limits
    sys.setrecursionlimit(100000)
    def bt(idx):
        if idx==n: return True
        v = order[idx]
        used = set()
        for u in adj[v]:
            if color[u]: used.add(color[u])
        if len(used)==3: return False
        for c in (1,2,3):
            if c in used: continue
            color[v]=c
            if bt(idx+1): return True
            color[v]=0
        return False
    return bt(0)

def has_K4(S):
    """quick: 4 mutually adjacent lattice points 0,a,b,c with all pairwise diffs in S."""
    Sset=set(S)
    Slist=[s for s in S if s>(0,0)]  # half, by lex
    L=len(Slist)
    for i in range(L):
        a=Slist[i]
        for j in range(i+1,L):
            b=Slist[j]
            if (a[0]-b[0],a[1]-b[1]) not in Sset: continue
            for k in range(j+1,L):
                c=Slist[k]
                if (a[0]-c[0],a[1]-c[1]) in Sset and (b[0]-c[0],b[1]-c[1]) in Sset:
                    return True
    return False

# ---------- patch + torus builders ----------
def patch_adj(S, R):
    """vertices = lattice points in [-R,R]^2; edges from S."""
    verts = [(x,y) for x in range(-R,R+1) for y in range(-R,R+1)]
    idx = {v:i for i,v in enumerate(verts)}
    adj = [set() for _ in verts]
    Sset=set(S)
    for v in verts:
        i=idx[v]
        for s in S:
            w=(v[0]+s[0], v[1]+s[1])
            if w in idx:
                adj[i].add(idx[w])
    return adj

def torus_3colorable(S, N):
    """Cay(Z_N^2, S mod N) — if 3-colorable, gives periodic 3-coloring of Z^2."""
    verts=[(x,y) for x in range(N) for y in range(N)]
    idx={v:i for i,v in enumerate(verts)}
    adj=[set() for _ in verts]
    for v in verts:
        i=idx[v]
        for s in S:
            w=((v[0]+s[0])%N,(v[1]+s[1])%N)
            j=idx[w]
            if j!=i:
                adj[i].add(j); adj[j].add(i)
    return three_colorable(adj)

# ---------- main classifier ----------
def classify(A,B,C,D, patchR=6, max_torus=14):
    S = connection_set(A,B,C,D)
    if not S:
        return None
    res = {'D':D, 'rD':len(S), 'S':S}
    if is_bipartite(S):
        res['chi']=2; res['why']='bipartite (mod-2 functional)'; return res
    # not bipartite => chi>=3
    if has_K4(S):
        res['chi_lb']=4; res['why_lb']='contains K4'
    else:
        res['chi_lb']=3
    # try to certify chi<=3
    if linear_3col(S):
        res['chi']=3; res['why']='non-bip + linear mod-3 coloring'; return res
    # search periodic 3-coloring on tori
    for N in range(2, max_torus+1):
        if torus_3colorable(S,N):
            res['chi']=3; res['why']=f'non-bip + torus Z_{N}^2 3-coloring'; return res
    # no small periodic 3-coloring found: test a finite patch for 3-colorability
    adj = patch_adj(S, patchR)
    if three_colorable(adj):
        res['chi_note']=f'patch R={patchR} IS 3-colorable but no periodic <= {max_torus} found'
        res['chi_lb']=res.get('chi_lb',3)
        res['why']='UNDETERMINED (>=3, patch 3-col but periodic search incomplete)'
        return res
    else:
        res['chi_lb']=4
        res['chi']='>=4'
        res['why_lb']=f'patch R={patchR} NOT 3-colorable -> chi>=4 (RIGOROUS)'
        return res
    return res

# ---------- lattice catalog: imaginary quadratic orders ----------
def imag_quad_form(d):
    """Maximal order of Q(sqrt(-d)), d>0 squarefree.
       d=3 mod4: form x^2+xy+((1+d)/4)y^2 ; else x^2 + d y^2 (disc -4d)."""
    if d % 4 == 3:
        return (1,1,(1+d)//4), -d
    else:
        return (1,0,d), -4*d

def popular_norms(A,B,C, Dmax, top=6):
    counts={}
    tr=A+C; det=A*C-(B*B)/4.0
    lam_min=(tr-(tr*tr-4*det)**0.5)/2.0
    bound=int((Dmax/lam_min)**0.5)+2
    for x in range(-bound,bound+1):
        for y in range(-bound,bound+1):
            if x==0 and y==0: continue
            D=Qval(A,B,C,x,y)
            if 0<D<=Dmax:
                counts[D]=counts.get(D,0)+1
    # sort by representation count desc, then D asc
    ranked=sorted(counts.items(), key=lambda kv:(-kv[1], kv[0]))
    return ranked[:top]

if __name__=='__main__':
    print("="*78)
    print("GENERIC 2D LATTICE CHROMATIC NUMBERS — imaginary quadratic orders (w=2)")
    print("="*78)
    # squarefree d, w=2 means disc<-4 i.e. d>=2 (d=1 Gaussian w=4, d=3 Eisenstein w=6 excluded)
    test_d = [2,5,6,7,10,11,13,14,15,19,21,22,23,30,31]
    Dmax=400
    overall_max_chi=2
    found4=[]
    for d in test_d:
        (A,B,C),disc = imag_quad_form(d)
        w = 2  # all these have only +-1 units
        if d==1: w=4
        if d==3: w=6
        pops = popular_norms(A,B,C,Dmax,top=5)
        print(f"\n--- Q(sqrt(-{d})), disc={disc}, form ({A},{B},{C}), w={w} ---")
        print(f"    popular norms (D:r(D)): {pops}")
        for D,rD in pops:
            r=classify(A,B,C,D)
            chi = r.get('chi','?')
            note = r.get('why', r.get('why_lb',''))
            extra = r.get('chi_note','')
            print(f"    D={D:4d} r={rD:2d}: chi={chi}  [{note}] {extra}")
            if chi==3: overall_max_chi=max(overall_max_chi,3)
            if chi=='>=4' or r.get('chi_lb',0)>=4:
                found4.append((d,D,rD,r))
    print("\n"+"="*78)
    print(f"OVERALL max certified chi among tested orders/popular-norms: {overall_max_chi}")
    if found4:
        print(f"!!! chi>=4 candidates found: {[(d,D) for d,D,_,_ in found4]}")
    else:
        print("NO chi>=4 found: every tested generic-order popular norm is chi<=3.")
    print("grep RESULT:", "FOUND4" if found4 else "ALL_LE_3")
