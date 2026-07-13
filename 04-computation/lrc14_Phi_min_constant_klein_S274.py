#!/usr/bin/env python3
"""
lrc14_Phi_min_constant_klein_S274.py
====================================
klein-2026-07-13-S274. The min-argument constant for the DEGREE-3 Phi functional.

Phi(E)=E_x[psi(N(x))], psi(N)=1-(2/3)N+(47/252)N(N-1)-(5/252)N(N-1)(N-2), N=#empty sectors.
Far element w: N -> N - 1{sec(wx) empty}. Error_Phi = -int Dpsi(N(x)) g_{Empty(x)}(wx) dx,
  Dpsi(N)=psi(N)-psi(N-1), g_{Empty}(y)=1{y in union of empty sectors}-N/7 (mean-zero in wx).
On each occupancy-constant interval I (N and empty-set fixed): |int_I| <= |Dpsi(N_I)|*min(osc(G_S)/w, ||g_S||*|I|),
osc(G_S)<=12/49, ||g_S||<=6/7. Min-argument => C_Phi <= (12/7)*max_N|Dpsi(N)| = (12/7)(2/3)=8/7~1.143 ABSOLUTE.

This script: (1) tabulate Dpsi(N); (2) for several clusters+w, compute the EXACT min-bound
Sum_I |Dpsi(N_I)|*min(osc(G_S_I)/w, ||g_S_I||*|I|) and confirm it bounds |Error_Phi|*w and is <= ~1.15.
"""
import math
NG=400009  # prime grid, independent of w
def occ_of(E,x):
    o=0
    for e in E:o|=1<<(int((e*x%1.0)*7.0)%7)
    return o
def psi(N):return 1-(2/3)*N+(47/252)*N*(N-1)-(5/252)*N*(N-1)*(N-2)
def Dpsi(N):return psi(N)-psi(N-1)
print("Dpsi(N) for N=1..7:", [round(Dpsi(N),4) for N in range(1,8)])
print(f"max_N |Dpsi(N)| = {max(abs(Dpsi(N)) for N in range(1,8)):.4f}  => universal C_Phi <= (12/7)*that = {(12/7)*max(abs(Dpsi(N)) for N in range(1,8)):.4f}")
print()

def Phi(E):
    s=0
    for k in range(1,NG):
        s+=psi(7-bin(occ_of(E,k/NG)).count("1"))
    return s/(NG-1)
def Phi_inf(C):
    # THM-710 transfer of moments then psi-combo (=Phi_from transferred m)
    s1=s2=s3=0
    for k in range(1,NG):
        N=7-bin(occ_of(C,k/NG)).count("1"); s1+=N;s2+=N*(N-1);s3+=N*(N-1)*(N-2)
    n=NG-1;m1,m2,m3=s1/n,s2/n,s3/n
    return 1-(2/3)*(6/7)*m1+(47/252)*(5/7)*m2-(5/252)*(4/7)*m3

def phi_minbound(C,w):
    """Sum over occupancy-constant intervals of |Dpsi(N_I)|*min(osc(G_S)/w, ||g_S||*|I|).
       osc(G_S) computed exactly per empty-set arrangement; ||g_S||=max(1-N/7,N/7)."""
    # precompute per empty-set-mask: osc(G_S), ||g_S||
    def osc_and_norm(mask):
        N=bin(mask).count("1")
        if N==0 or N==7: return 0.0,0.0
        a=N/7.0; gnorm=max(1-a,a)
        # antiderivative oscillation of 1_{union of empty sectors}-a over the circle:
        # walk 7 sectors, height += (1-a) if empty else -a? no: g=+ (1-a) on empty, -a on filled... wait
        # g_S = 1 on empty-union (measure a=N/7) - a. So on empty sector: 1-a; on filled: -a.
        vals=[]; h=0.0
        for s in range(7):
            slope=(1-a) if (mask>>s)&1 else (-a)
            h+=slope*(1/7); vals.append(h)
        return (max(vals)-min(vals)), gnorm
    cache={}
    minb=0.0; run=0; runmask=-1
    prevmask=None
    # scan intervals of constant (N, empty-set): empty-set mask = (~occ)&0x7F
    for k in range(1,NG):
        occ=occ_of(C,k/NG); emask=(~occ)&0x7F
        if emask==prevmask:
            run+=1
        else:
            if prevmask is not None and prevmask not in (0,) :
                N=bin(prevmask).count("1")
                if 1<=N<=6:
                    if prevmask not in cache: cache[prevmask]=osc_and_norm(prevmask)
                    osc,gn=cache[prevmask]
                    minb+=abs(Dpsi(N))*min(osc/w, gn*run/NG)
            prevmask=emask; run=1
    if prevmask is not None:
        N=bin(prevmask).count("1")
        if 1<=N<=6:
            if prevmask not in cache: cache[prevmask]=osc_and_norm(prevmask)
            osc,gn=cache[prevmask]; minb+=abs(Dpsi(N))*min(osc/w, gn*run/NG)
    return minb

print("cluster/w:  |Error_Phi|*w   Phi-min-bound   (bound>=err? bound<=1.15?)")
for C in [[0,1,2,3,4,5,6],[0,1,2,3,4,5,30],[0,2,4,6,8,10,12],[0,3,7,12,20,33,54],[0,1,2,118,119,120,60]]:
    pinf=Phi_inf(C)
    for w in [1009, 2003]:
        errw=abs(Phi(C+[w])-pinf)*w
        mb=phi_minbound(C,w)
        print(f"  {str(C):30s} w={w}: |Err|*w={errw:.4f}  min-bound={mb:.4f}  {'OK' if errw<=mb+1e-3 else 'CHK'}")
print("\ndone.")
