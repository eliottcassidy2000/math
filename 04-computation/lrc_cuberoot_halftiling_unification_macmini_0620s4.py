#!/usr/bin/env python3
"""
lrc_cuberoot_halftiling_unification_macmini_0620s4.py  (mac-mini-2026-06-20-S4)

Ground the half-tiling <-> coverage unification via the shared 3+3+1 cube-root skeleton.

For 3 far runners {u,v,w} over a bounded core B, the Newton packets:
  one-far  A=Dlt_u, B=Dlt_v, C=Dlt_w   (Dlt_x = p0(B+x)-p0(B))
  two-far  D=I_uv,  E=I_uw,  F=I_vw     (I_xy = mixed 2nd diff)
  three-far G=Dlt_uvw                   (mixed 3rd diff)
Full correction H(1)=A+B+C+D+E+F+G = p0(B+u+v+w)-p0(B).
User's recursion R = A+B+C-D-E-F+G  (= H(1)-2(D+E+F), codex pair-tax shadow).
Eisenstein modes (cube root w3): S_w=A+w3 B+w3^2 C, P_w=D+w3 E+w3^2 F.

CHECK: (a) the 3+3+1 packet decomposition reproduces p0 exactly; (b) S_3 symmetry — permuting
(u,v,w) permutes (A,B,C) and (D,E,F), so |S_w| and |P_w| are S_3-INVARIANT magnitudes (the
cube-root modulus is the symmetric invariant); (c) the half-tiling parallel: 3 corner-sizes.
"""
import itertools, cmath, sys
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
def sector_of(p): return int((p%1)*7)
def measS7(E):
    E=sorted(set(E)); bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); tot=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        if len(set(sector_of(e*((x0+x1)/2)) for e in E))==7: tot+=x1-x0
    return tot
def packets(B,u,v,w):
    p0=measS7(B); 
    Du=measS7(B+[u])-p0; Dv=measS7(B+[v])-p0; Dw=measS7(B+[w])-p0
    Iuv=measS7(B+[u,v])-measS7(B+[u])-measS7(B+[v])+p0
    Iuw=measS7(B+[u,w])-measS7(B+[u])-measS7(B+[w])+p0
    Ivw=measS7(B+[v,w])-measS7(B+[v])-measS7(B+[w])+p0
    G=(measS7(B+[u,v,w])-measS7(B+[u,v])-measS7(B+[u,w])-measS7(B+[v,w])
       +measS7(B+[u])+measS7(B+[v])+measS7(B+[w])-p0)
    return (Du,Dv,Dw,Iuv,Iuw,Ivw,G,p0)

w3=cmath.exp(2j*cmath.pi/3)
print("3+3+1 cube-root packet structure (B=(0,1,2,3,4), far triples):")
for (u,v,w) in [(15,16,17),(15,23,41),(20,21,22),(15,30,45)]:
    A,Bp,C,D,E,Fp,G,p0=packets([0,1,2,3,4],u,v,w)
    full=A+Bp+C+D+E+Fp+G
    direct=measS7([0,1,2,3,4,u,v,w])-p0
    R=A+Bp+C-D-E-Fp+G
    Sw=complex(A)+w3*complex(Bp)+w3**2*complex(C)
    Pw=complex(D)+w3*complex(E)+w3**2*complex(Fp)
    print(f"  ({u},{v},{w}): full=H(1)={float(full):+.5f} (direct {float(direct):+.5f} {'OK' if full==direct else 'MISMATCH'})")
    print(f"     one-far(A,B,C)=({float(A):.4f},{float(Bp):.4f},{float(C):.4f})  two-far(D,E,F)=({float(D):.4f},{float(E):.4f},{float(Fp):.4f})  G={float(G):.5f}")
    print(f"     pair-tax R=A+B+C-D-E-F+G={float(R):+.5f}   |S_w|={abs(Sw):.5f}  |P_w|={abs(Pw):.5f}  (Eisenstein moduli, S3-invariant)")

# S_3 invariance check: permute (u,v,w), |S_w| of (A,B,C) should be permutation-invariant
print("\nS_3-invariance of the Eisenstein modulus |S_w| (permute far runners):")
B=[0,1,2,3,4]; base=(15,23,41)
mods=[]
for perm in itertools.permutations(base):
    A,Bp,C,_,_,_,_,_=packets(B,*perm)
    Sw=complex(A)+w3*complex(Bp)+w3**2*complex(C); mods.append(abs(Sw))
print(f"  |S_w| over 6 permutations of {base}: min={min(mods):.5f} max={max(mods):.5f} spread={max(mods)-min(mods):.2e}")
print("  (cyclic perms preserve |S_w| exactly; transpositions conjugate it -> same modulus)")
