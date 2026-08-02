#!/usr/bin/env python3
"""Infinite exact referee for the four primitive label-12 cell-90 lanes.

The proof combines a Dirichlet d<=8 reindexing, an actual-phase composite
trapezoid estimate, a uniform T>=5 estimate, and finitely many rigorously
bounded affine resonance channels when T<5.  It uses Fraction arithmetic
only; no floating-point comparison participates in the certificate.
"""
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as F
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

CORE=Path(__file__).with_name("lrc_label12_multiblock_core_thm3171.py")
EXPECTED_CORE_SHA256="7cff3a20e4d874f1143d096d3c9f4fa37448eb76c0ae0408513dd81780d8382b"
if sha256(CORE.read_bytes()).hexdigest()!=EXPECTED_CORE_SHA256:
 raise RuntimeError(("core hash",sha256(CORE.read_bytes()).hexdigest(),EXPECTED_CORE_SHA256))
SPEC=spec_from_file_location("multiblock", CORE)
MB=module_from_spec(SPEC); SPEC.loader.exec_module(MB)

LABELS=((1,12),(12,1),(2,12),(12,2))
# The finite run gives Q>10000.  Since Q<2P, the primitive smaller level
# satisfies P>=5001; this is the correct analytic-tail threshold.
P_MIN=5001
MAX_D=8
MAX_C=41
WORKERS=8
EXPECTED_CHANNELS=11608

def require(condition,detail):
 if not condition:raise RuntimeError(detail)
def ceildiv(a,b):return -((-a)//b)
def ceilf(x):return -((-x.numerator)//x.denominator)

def channel_tasks():
 for d in range(1,MAX_D+1):
  for a in range(d+1):
   D=d+a
   for c in range(-MAX_C,MAX_C+1):
    if c==0:continue
    for p0 in range(1,d+1):
     if (a*p0+c)%d:continue
     q0=p0+(a*p0+c)//d
     n0=ceildiv(P_MIN-p0,d)
     if a:
      n0=max(n0,ceildiv(1-(q0-p0),a))
     p=p0+d*n0;q=q0+D*n0
     if not p<q<2*p:continue
     for e,f in LABELS:
      A=168*c+e*D-d*f
      require(A!=0,("zero resonance step",d,a,c,p0,e,f))
      yield d,a,c,p0,q0,e,f,A,n0

def structural_certificate():
 label_perturb=max(abs(e*(d+a)-d*f)
                   for d in range(1,MAX_D+1) for a in range(d+1)
                   for e,f in LABELS)
 require(label_perturb==184,("label perturb",label_perturb))
 # Dirichlet gives |c|<=P/9, hence |A|<=(56/3)P+184 and rho<9/80.
 rho_cap=(F(56,3)*P_MIN+label_perturb)/(168*P_MIN-12)
 require(rho_cap<F(9,80),("rho cap",rho_cap))
 # If T<5, |A|<840*d*P/(P-7); integrality and the perturb cap force |c|<=41.
 A_strict=F(840*MAX_D*P_MIN,P_MIN-7)
 A_integer_ceiling=A_strict.numerator//A_strict.denominator
 require(F(A_integer_ceiling+label_perturb,168)<42,
         ("c cap",A_strict,A_integer_ceiling,label_perturb))
 return label_perturb,rho_cap,A_strict,A_integer_ceiling

def limit_and_error(task):
 d,a,c,p0,q0,e,f,A,n0=task
 P0=p0+d*n0
 y0=F(d+a,d); W0=168*y0; Z0=F(168)
 s=F(168*c,d)-f
 Zmin=Z0-F(e,P0); Wmin=W0-F(abs(s),P0)
 require(Zmin>0 and Wmin>0,("positive scaling",task,Zmin,Wmin))

 # Uniform scaled-tent perturbation constants: |F_P-F_inf|<=Kf/P.
 inv_const=abs(s)/(Wmin*W0)
 kH=24*inv_const; kL=168*inv_const
 kRatio=abs(168*s+W0*e)/(168*Zmin)
 ka=kRatio/14
 Wmax=W0+F(abs(s),P0)
 amax=max(F(1+y0,14),F(1+Wmax/Zmin,14))
 L0=168/W0
 kOuter=kL*amax+L0*ka
 kF=max(kH,kOuter)

 Kalpha=F(abs(A),d*Zmin)
 KG=d*kF+L0*Kalpha*d*(d-1)/2
 R=(90*e)%168; S=(90*f)%168
 phase0=R*y0-S
 Kx=abs(R*s+S*e+e*phase0)/(168*Zmin)
 Kpath=KG+d*L0*Kx

 r0=p0%d
 T0=F(abs(A),168*d)
 Kt=F(abs(A)*abs(168*r0-e),168*d)/Zmin
 Tmax=T0+Kt/P0
 Hmax=24/Wmin
 KJ=Tmax*Kpath+Kt*d*Hmax

 # Exact continuum limit along the canonical cell-90 phase.
 G0=MB.PeriodicPL(Z0,W0,F(a,d)%1,d)
 X=((R*y0-S)/168)%1
 m=T0.numerator//T0.denominator; theta=T0-m
 forward=A>0
 start=X if forward else X-theta
 J0=m*G0.mean+G0.interval(start,theta)
 continuum=F(168,abs(A))*J0

 Kmain=F(168,abs(A))*KJ+F(e,abs(A))*Tmax*d*Hmax
 Kend=d*Hmax/2
 Ktrap=F(abs(A)*84*d*ceilf(Tmax),1)/(Zmin*Wmin*P0)
 K=Kmain+Kend+Ktrap
 certified=continuum-K/P0
 return task,continuum,K,certified

def large_T_certificate():
 P=P_MIN; R=F(9,80)
 main=F(24*(P-7),1176*P)*F(4,5)
 endpoint=F(96,168*(P+1)-12)
 trap=F(84,168*(P+1)-12)*(P*R*R+8*R)
 # The positive main term increases; both negative terms decrease with P.
 require(R*R*156<8*R*168,("trap monotonicity",R))
 value=main-endpoint-trap
 require(value>F(1,105),("large-T bound",value))
 return value,value-F(1,105)

def main():
 structure=structural_certificate()
 large,large_margin=large_T_certificate()
 tasks=tuple(channel_tasks())
 require(len(tasks)==EXPECTED_CHANNELS,("channel count",len(tasks),EXPECTED_CHANNELS))
 with ProcessPoolExecutor(max_workers=WORKERS) as pool:
  rows=tuple(pool.map(limit_and_error,tasks,chunksize=1))
 require(all(row[3]>F(1,105) for row in rows),
         ("resonance failure",min(rows,key=lambda row:row[3])))
 worst=min(rows,key=lambda row:row[3])
 semantic=repr((structure,large,large_margin,rows)).encode()
 print("LRC LABEL-12 PRIMITIVE HORN -- DIRICHLET d-BLOCK REFEREE")
 print(f"large_T_floor={large};large_T_margin={large_margin}")
 print(f"label_perturb_cap={structure[0]};rho_cap={structure[1]};"
       f"Tsmall_A_strict_cap={structure[2]};Tsmall_A_integer_cap={structure[3]}")
 print(f"affine_resonance_channels={len(rows)};c_range=-{MAX_C}..-1,1..{MAX_C}")
 print(f"worst_channel={worst[0]};continuum={worst[1]};error_constant={worst[2]};"
       f"certified_floor={worst[3]};margin={worst[3]-F(1,105)}")
 print(f"semantic_sha256={sha256(semantic).hexdigest()}")

if __name__=="__main__":main()
