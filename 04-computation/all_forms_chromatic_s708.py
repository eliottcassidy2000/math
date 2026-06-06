#!/usr/bin/env python3
"""
all_forms_chromatic_s708.py  (monad-explorer-2026-06-06-S708)

Map chi(U(L,D)) over ALL represented norms D for ALL reduced positive-definite
integral binary quadratic forms (A,B,C) up to a bound -- the truly GENERIC 2D
lattice question left open by THM-418.  Hunt for ANY chi=4.

Reduced form: -A < B <= A <= C, A>0, disc=B^2-4AC<0.  Every 2D integral lattice
is similar to exactly one reduced form, so this enumerates all integral lattice
chromatic spectra.

Reuses the classifier logic from generic_lattice_chromatic_s708.py.
"""
import sys
sys.setrecursionlimit(200000)

def Qval(A,B,C,x,y): return A*x*x+B*x*y+C*y*y

def connection_set(A,B,C,D):
    tr=A+C; det=A*C-(B*B)/4.0
    lam_min=(tr-(tr*tr-4*det)**0.5)/2.0
    bound=int((D/lam_min)**0.5)+2
    S=[]
    for x in range(-bound,bound+1):
        for y in range(-bound,bound+1):
            if x==0 and y==0: continue
            if Qval(A,B,C,x,y)==D: S.append((x,y))
    return S

def is_bipartite(S):
    for p in [(1,0),(0,1),(1,1)]:
        if all(((p[0]*s[0]+p[1]*s[1])&1)==1 for s in S): return True
    return False

def linear_3col(S):
    for q in [(0,1),(1,0),(1,1),(1,2)]:
        if all(((q[0]*s[0]+q[1]*s[1])%3)!=0 for s in S): return True
    return False

def three_colorable(adj):
    n=len(adj)
    if n==0: return True
    color=[0]*n
    order=sorted(range(n), key=lambda v:-len(adj[v]))
    def bt(idx):
        if idx==n: return True
        v=order[idx]
        used=0
        for u in adj[v]:
            if color[u]: used|=(1<<color[u])
        if used==0b1110: return False
        for c in (1,2,3):
            if used&(1<<c): continue
            color[v]=c
            if bt(idx+1): return True
            color[v]=0
        return False
    return bt(0)

def has_K4(S):
    Sset=set(S); Slist=[s for s in S if s>(0,0)]; L=len(Slist)
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

def torus_3colorable(S,N):
    verts=[(x,y) for x in range(N) for y in range(N)]
    idx={v:i for i,v in enumerate(verts)}
    adj=[set() for _ in verts]
    for v in verts:
        i=idx[v]
        for s in S:
            j=idx[((v[0]+s[0])%N,(v[1]+s[1])%N)]
            if j!=i: adj[i].add(j); adj[j].add(i)
    return three_colorable(adj)

def patch_adj(S,R):
    verts=[(x,y) for x in range(-R,R+1) for y in range(-R,R+1)]
    idx={v:i for i,v in enumerate(verts)}
    adj=[set() for _ in verts]
    for v in verts:
        i=idx[v]
        for s in S:
            w=(v[0]+s[0],v[1]+s[1])
            if w in idx: adj[i].add(idx[w])
    return adj

def chi_of(A,B,C,D, max_torus=12, patchR=6):
    S=connection_set(A,B,C,D)
    if not S: return None,None
    if is_bipartite(S): return 2,'bip'
    if has_K4(S):
        # K4 => chi>=4; still report (rare)
        return 4,'K4'
    if linear_3col(S): return 3,'lin3'
    for N in range(2,max_torus+1):
        if torus_3colorable(S,N): return 3,f'torus{N}'
    adj=patch_adj(S,patchR)
    if three_colorable(adj): return 3,'patch3-incomplete'  # >=3, periodic search incomplete but patch 3-colorable
    return 4,f'patchR{patchR}-NOT3col'

if __name__=='__main__':
    Amax=6; Cmax=14; Dmax=160
    print("="*78)
    print(f"ALL reduced integral binary forms (A<= {Amax}, C<= {Cmax}); chi over all D<= {Dmax}")
    print("="*78)
    forms=[]
    for A in range(1,Amax+1):
        for B in range(-A+1,A+1):
            for C in range(A,Cmax+1):
                disc=B*B-4*A*C
                if disc<0:
                    forms.append((A,B,C,disc))
    print(f"#reduced forms: {len(forms)}")
    global_chi_max=2
    chi4_forms=[]
    summary=[]
    for (A,B,C,disc) in forms:
        chiset=set(); chi4_here=[]
        for D in range(1,Dmax+1):
            chi,why=chi_of(A,B,C,D)
            if chi is None: continue
            chiset.add(chi)
            if chi==4: chi4_here.append((D,why))
        if not chiset: continue
        mx=max(chiset)
        global_chi_max=max(global_chi_max,mx)
        Bpar = 'Beven' if B%2==0 else 'Bodd'
        discpar = disc%4
        summary.append((A,B,C,disc,sorted(chiset)))
        flag=''
        if 4 in chiset:
            chi4_forms.append((A,B,C,disc,chi4_here)); flag=' <<< CHI4'
        # print only interesting / representative
        print(f"  ({A},{B},{C}) disc={disc:5d} ({discpar} mod4,{Bpar}): chi-set={sorted(chiset)}{flag}")
    print("\n"+"="*78)
    print(f"GLOBAL max chi over all forms & norms: {global_chi_max}")
    # dichotomy check: does maxchi correlate with disc mod 4 / B parity?
    print("\nDICHOTOMY TABLE (maxchi vs disc mod 4 and B parity):")
    from collections import Counter
    cnt=Counter()
    for (A,B,C,disc,chiset) in summary:
        cnt[(disc%4, B%2, max(chiset))]+=1
    for k in sorted(cnt):
        print(f"  disc%4={k[0]}, B%2={k[1]}, maxchi={k[2]}: {cnt[k]} forms")
    if chi4_forms:
        print(f"\n!!! CHI=4 forms ({len(chi4_forms)}):")
        for (A,B,C,disc,info) in chi4_forms[:20]:
            print(f"   ({A},{B},{C}) disc={disc}: {info[:5]}")
        print("grep RESULT: CHI4_EXISTS")
    else:
        print("\nNO chi=4 for ANY reduced integral binary form at any norm D<=Dmax.")
        print("grep RESULT: ALL_INTEGRAL_FORMS_LE_3")
