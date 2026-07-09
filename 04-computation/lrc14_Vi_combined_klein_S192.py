#!/usr/bin/env python3
"""
klein-2026-07-08-S192 (part 2): the CANCELLATION-aware a-priori constant.

The triangle-inequality C = sum_i |dD3/dm_i| V_i is loose by ~150x because the
three moment errors eps_i = m_i(E_d)-L_i share ONE grid-sampling defect and
cancel in D3. Capture it: to first order
   D3(E_d) - D3_inf = E_a[ (1/d) sum_k g_a(u_k) - integral g_a du ] + O(1/d^2),
   g_a(u) = c1 W(a;u) + c2 W(a;u)^2 + c3 W(a;u)^3,   c_i = dD3/dm_i |_L .
So the first-order constant is  C1 = E_a[ TV_u g_a ]  (Riemann error <= TV/d),
and g_a = phi(W) with phi(w)=c1 w + c2 w^2 + c3 w^3 a FIXED cubic. Since phi is
fixed, TV_u[phi(W)] <= Lip(phi on [0,6/7]) * TV_u[W], giving the a-priori
   C1 <= Lip(phi) * V1,   V1 = E_a[TV_u W] <= 2 E_a[W_full].
This script measures C1, the true deviation*d, E_a[W_full], Lip(phi), and checks
the a-priori chain against the required C <= (D3_inf - bar)*d0.
"""
import numpy as np
from math import gcd

INV7=1.0/7.0; M=6.0/7.0
bar=83549.0/252252.0

def W_from_phases(ph):
    p=np.sort(np.mod(ph,1.0)); g=np.diff(p); g=np.append(g,1.0-p[-1]+p[0])
    return np.maximum(g-INV7,0.0).sum()
def AP(a,K=10): return np.mod(a*np.arange(K),1.0)

# limit moments + partials (reuse S192 part1 values, recompute for self-containment)
Na,Nu=600,600
av=(np.arange(Na)+0.5)/Na; uv=(np.arange(Nu)+0.5)/Nu
S1=S2=S3=0.0; Wfull_sum=0.0
for a in av:
    ap=AP(a); Wfull_sum+=W_from_phases(ap)
    Wv=np.array([W_from_phases(np.append(ap,u)) for u in uv])
    S1+=Wv.sum();S2+=(Wv**2).sum();S3+=(Wv**3).sum()
tot=Na*Nu
L1,L2,L3=S1/tot,S2/tot,S3/tot
EWfull=Wfull_sum/Na
def D3_of(m1,m2,m3):
    Nn=m1-m2/M; D=m2-m3/M; return m1/M+Nn*Nn/D
D3inf=D3_of(L1,L2,L3)
Nn=L1-L2/M; D=L2-L3/M
c1=1.0/M+2*Nn/D; c2=-2*Nn/(M*D)-Nn*Nn/(D*D); c3=Nn*Nn/(M*D*D)
print(f"L1,L2,L3 = {L1:.5f},{L2:.5f},{L3:.5f}   D3inf={D3inf:.5f}  bar={bar:.5f}")
print(f"partials c1,c2,c3 = {c1:.3f},{c2:.3f},{c3:.3f}")
print(f"E_a[W_full] = {EWfull:.5f}   =>  a-priori V1 <= 2*E_a[W_full] = {2*EWfull:.5f}")

# phi(w) = c1 w + c2 w^2 + c3 w^3 ; Lip on [0,6/7]
ws=np.linspace(0,M,2000)
phip=c1+2*c2*ws+3*c3*ws*ws
Lip_phi=np.max(np.abs(phip))
print(f"Lip(phi) on [0,6/7] = {Lip_phi:.4f}  (attained at w={ws[np.argmax(np.abs(phip))]:.3f})")

# measured first-order constant C1 = E_a[ TV_u g_a ], g_a=phi(W(a;.))
def TV_gapWfull(a):
    ap=np.sort(AP(a)); g=np.diff(ap); g=np.append(g,1.0-ap[-1]+ap[0])
    return np.minimum(2*np.maximum(g-INV7,0.0),2*INV7).sum()
Na2=500; av2=(np.arange(Na2)+0.5)/Na2; Nu2=8000
TVg=[]; TVW=[]
for a in av2:
    ap=AP(a); us=(np.arange(Nu2)+0.5)/Nu2
    Wv=np.array([W_from_phases(np.append(ap,u)) for u in us])
    ga=c1*Wv+c2*Wv**2+c3*Wv**3
    tv=np.abs(np.diff(ga)).sum()+abs(ga[0]-ga[-1]); TVg.append(tv)
    TVW.append(TV_gapWfull(a))
C1=np.mean(TVg); V1=np.mean(TVW)
print(f"\nMEASURED first-order C1 = E_a[TV_u g_a]     = {C1:.4f}")
print(f"MEASURED V1 = E_a[TV_u W]                   = {V1:.4f}   (<= 2E[W_full]={2*EWfull:.4f} ok:{V1<=2*EWfull+1e-6})")
print(f"a-priori chain  C1 <= Lip(phi)*V1           = {Lip_phi*V1:.4f}")
print(f"a-priori chain  C1 <= Lip(phi)*2E[W_full]   = {Lip_phi*2*EWfull:.4f}")

# true deviation*d for real coprime families (second-order check)
print("\n--- true |D3(E_d)-D3inf|*d for coprime (d,p) [captures 1st+2nd order] ---")
def D3_family(elems,N=200000):
    E=np.array(sorted(set(elems)),dtype=float)
    x=(np.arange(N)+0.5)/N
    P=np.sort(np.mod(E[None,:]*x[:,None],1.0),axis=1)
    g=np.diff(P,axis=1); g=np.concatenate([g,(1.0-P[:,-1]+P[:,0])[:,None]],axis=1)
    W=np.maximum(g-INV7,0.0).sum(axis=1)
    m1,m2,m3=W.mean(),(W**2).mean(),(W**3).mean()
    Dl=m2-m3/M; return m1/M+(m1-m2/M)**2/Dl
import numpy as _np
for d in (30,50,80,120,200):
    # pick p coprime to d, interior-ish ratio ~ 8/3
    p=int(round(d*8/3));
    while gcd(p,d)!=1: p+=1
    val=D3_family(list(d*_np.arange(10))+[p])
    print(f"  d={d:3d} p={p:3d}: D3={val:.5f}  |dev|*d={abs(val-D3inf)*d:.4f}")

req26=(D3inf-bar)*26
print(f"\nrequired C for d>=26: (D3inf-bar)*26 = {req26:.3f}")
print(f"a-priori first-order C1 bound Lip(phi)*2E[W_full] = {Lip_phi*2*EWfull:.3f}  "
      f"clears(1st order):{Lip_phi*2*EWfull<=req26}")
