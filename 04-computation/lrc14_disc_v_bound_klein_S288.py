#!/usr/bin/env python3
"""
lrc14_disc_v_bound_klein_S288.py
================================
klein-2026-07-13-S288 (owner: prove the analytic disc_v bound).

disc_v = Sum_{m!=0}|chat_{mv}|^2 = v-grid discrepancy of the good-set autocorrelation A_{~v}, where
G'_{~v}=leave-one-out good set (speeds w!=v), chat_l=Fourier coeff of 1_{G'_{~v}}. THM-731 certificate
needs  disc_v < 6|G'_{~v}|^2.

Exact reformulations (this session):
  (I)  disc_v = Var(Phi_v),  Phi_v(t)=(1/v)#{j in [0,v): t+j/v in G'_{~v}}   [orbit-average variance]
  (II) disc_v = (1/2v^2) Sum_{p,q endpoints} eps_p eps_q B_2({v(p-q)})       [B_2 endpoint-pair disc
        = LITERALLY the density Q_s object of THM-729]

Candidate RIGOROUS upper bounds (want one that certifies, i.e. < 6|G'|^2):
  (A) all-energy:      disc_v <= |G'|(1-|G'|)                    [Parseval; loosest]
  (B) crude BV:        disc_v <= r^2/(3 v^2)                     [|U|<=2r, |B_2|<=1/6]
  (C) straddle L2:     disc_v <= 2 r |G'| / v                    [||Phi-|G'|||_inf<=r/v, L1<=2|G'|]
  (D) diagonal:        r/(6 v^2)  [the p=q part of (II); a LOWER-order reference, not an upper bound
        unless OffDiag<=0]
This script computes r (#arcs of G'_{~v}), the true disc_v, all bounds, and 6|G'|^2, for the best peel v.
"""
import numpy as np
NG=1<<21
THR=1.0/14.0
t=np.arange(NG,dtype=np.float64)/NG

def good_ind(S):
    g=np.ones(NG,dtype=np.float64)
    for w in S:
        fr=(w*t)%1.0; d=np.minimum(fr,1.0-fr); g*=(d>=THR)
    return g

def n_arcs(g):
    # number of maximal runs of 1 on the circle
    gi=g.astype(np.int8)
    diff=np.diff(np.concatenate([gi,gi[:1]]))
    ups=np.sum(diff==1)
    return int(ups)  # each up-edge starts one arc (circular)

def disc_true(g,v):
    G=np.fft.rfft(g); A=np.fft.irfft(G*np.conj(G),n=NG)/NG
    Abar=A.mean()
    idx=np.round((np.arange(v)/v)*NG).astype(np.int64)%NG
    return A[idx].mean()-Abar

def Dv_ind(v):
    fr=(v*t)%1.0; d=np.minimum(fr,1.0-fr); return (d<THR).astype(np.float64)

FAM=[([1,2,3,4,5,6,7,8,9,10,11,12,182],"deep well {1..12,182}",182),
     ([1,2,3,4,5,6,7,8,9,10,11,13,84], "residue {1..11,13,84}",84),
     ([2,3,4,5,6,7,8,9,10,11,12,13,14],"{2..14}",14),
     ([1,3,4,5,6,7,8,9,10,11,12,13,182],"variant {1,3..13,182}",182)]

print("="*112)
print("disc_v BOUNDS vs the certificate threshold 6|G'_{~v}|^2.  (need some RIGOROUS bound < threshold)")
print("="*112)
print("%-26s %4s %8s %6s %10s | %10s %10s %10s %10s %10s"%(
    "family(peel v)","v","|G'_~v|","r","6|G'|^2","disc_TRUE","(A)energy","(B)r^2/3v^2","(C)2r|G'|/v","(D)diag"))
for S,name,v in FAM:
    Snv=[w for w in S if w!=v]
    g=good_ind(Snv); L=g.mean(); r=n_arcs(g)
    dv=disc_true(g,v)
    thr=6*L*L
    A=L*(1-L); B=r*r/(3.0*v*v); C=2.0*r*L/v; D=r/(6.0*v*v)
    def mark(x): return ("<" if x<thr else " ")  # certifies?
    print("%-26s %4d %8.5f %6d %10.3e | %9.3e%s %9.3e%s %9.3e%s %9.3e%s %9.3e"%(
        name,v,L,r,thr, dv,mark(dv), A,mark(A), B,mark(B), C,mark(C), D))
print("-"*112)
print(" '<' after a value means it is BELOW the certificate threshold 6|G'|^2 (i.e. it certifies).")
print(" (A) all-energy, (B) crude-BV, (C) straddle-L2 are RIGOROUS upper bounds; (D) diagonal is the")
print(" p=q reference (rigorous upper bound only if OffDiag<=0). disc_TRUE is the actual value.")
print()
# also: does OffDiag<=0? i.e. is disc_TRUE <= diagonal r/(6v^2)?
print("OffDiag check (disc_TRUE vs diagonal r/(6v^2)):")
for S,name,v in FAM:
    Snv=[w for w in S if w!=v]; g=good_ind(Snv); r=n_arcs(g); dv=disc_true(g,v); D=r/(6.0*v*v)
    print("  %-26s disc_TRUE=%.3e  diag=%.3e  ratio disc/diag=%.3f  %s"%(
        name,dv,D,dv/D,"OffDiag<0 (disc<diag)" if dv<D else "OffDiag>0 (disc>diag)"))
print("\ndone.")
