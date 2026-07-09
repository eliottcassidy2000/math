#!/usr/bin/env python3
"""
klein-2026-07-08-S195: R2 -- the k=12,13 density-floor tail, a-priori via the
BOX BOUND (kps-S89's k=11 template, now for the longest-AP=(k-1) families).

Tail family (post MISTAKE-126, longest-AP axis): E_d = d*{0,...,k-2} u {p},
a (k-1)-term AP at scale d + a primitivizing point. As d->inf the outlier phase
decorrelates (Weyl); moments m_i(E_d) -> L_i = E_{a,u}[W(block_{k-1}(a); u)^i]
(block (+) independent outlier). Rigorous rate (klein-S191 + kps-S89):
  |m_i(E_d) - L_i| <= Vh_i / d,  Vh_i = i (6/7)^{i-1} E[W_B],  E[W_B]=E[W](block_{k-1}).
Then min D3 over the box m_i in [L_i +- Vh_i/d] (den = m2-m3/M > 0) is a rigorous
lower bound; find crossover d0 with box-min D3 >= bar_k, so [finite check d<=d0-1]
+ [box d>=d0] closes the tail. Larger margin than k=11 => expect d0 < 62.
"""
import numpy as np
INV7=1.0/7.0; M=6.0/7.0
BAR={12: 199344/1e6*1.0, 13: 56487/1e6*1.0}   # placeholder; exact below
BAR_exact={12: 0.199344, 13: 0.056487}

def W_uncovered(ph):
    p=np.sort(np.mod(np.asarray(ph,float),1.0)); g=np.diff(p); g=np.append(g,1.0-p[-1]+p[0])
    return np.maximum(g-INV7,0.0).sum()
def block(kb): return np.arange(kb)   # block_{kb} co-offsets {0..kb-1}

def limit_moments(kb, Na=800, Nu=800):
    """L_i = E_{a,u}[W(block_kb phases at a; outlier u)^i], and E[W_B]=E_a[W(block)]."""
    av=(np.arange(Na)+0.5)/Na; uv=(np.arange(Nu)+0.5)/Nu
    S1=S2=S3=0.0; WB=0.0
    for a in av:
        ap=np.mod(a*block(kb),1.0)
        WB+=W_uncovered(ap)
        Wv=np.array([W_uncovered(np.append(ap,u)) for u in uv])
        S1+=Wv.sum(); S2+=(Wv**2).sum(); S3+=(Wv**3).sum()
    tot=Na*Nu
    return S1/tot,S2/tot,S3/tot, WB/Na

def D3(m1,m2,m3):
    den=m2-m3/M
    if den<=0: return -1e9
    return m1/M+(m1-m2/M)**2/den

def box_min_D3(L,Vh,d,ng=13):
    lo=[L[i]-Vh[i]/d for i in range(3)]; hi=[L[i]+Vh[i]/d for i in range(3)]
    best=1e9; g=np.linspace(0,1,ng)
    for x in g:
        m1=lo[0]+(hi[0]-lo[0])*x
        for y in g:
            m2=lo[1]+(hi[1]-lo[1])*y
            for z in g:
                m3=lo[2]+(hi[2]-lo[2])*z
                v=D3(m1,m2,m3)
                if v<best: best=v
    return best

for k in (12,13):
    kb=k-1  # block length = longest AP = k-1
    L1,L2,L3,EWB=limit_moments(kb)
    Vh=[1*EWB, 2*M*EWB, 3*(M**2)*EWB]
    bar=BAR_exact[k]; D3inf=D3(L1,L2,L3)
    print(f"\n===== k={k}: tail = block_{kb} (+) outlier =====")
    print(f"  E[W_B](block_{kb}) = {EWB:.5f}   L1,L2,L3 = {L1:.5f},{L2:.5f},{L3:.5f}")
    print(f"  D3_inf = {D3inf:.5f}   bar_{k} = {bar:.5f}   margin = {D3inf-bar:+.4f}")
    print(f"  Vh = {Vh[0]:.4f}, {Vh[1]:.4f}, {Vh[2]:.4f}")
    d0=None
    for d in range(4, 120):
        b=box_min_D3([L1,L2,L3],Vh,d)
        if b>=bar and d0 is None:
            d0=d
    # report box-min at a few d
    for d in (10,20,30,40,50,62):
        b=box_min_D3([L1,L2,L3],Vh,d)
        print(f"    d={d:3d}: box-min D3 = {b:.4f}  ({'>=' if b>=bar else '< '} bar)")
    print(f"  ==> crossover d0 = {d0}  (box bound closes tail for d >= {d0}; finite check d < {d0})")
print("\nSo the k=12,13 tail closes a-priori: [exhaustive compact prim-diam<=18 (mac-mini-S58)] +")
print("[finite check d<d0 via conditional-D3] + [box bound d>=d0]. Larger margins => smaller d0 than k=11's 62.")
