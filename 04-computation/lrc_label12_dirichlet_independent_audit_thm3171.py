#!/usr/bin/env python3
"""Independent exact audit of the label-12 Dirichlet tail referee.

This audit deliberately does not use the referee's PeriodicPL machinery for
integration or direct overlap masses.  It reconstructs the channel universe,
integrates individual periodic tents by polygonal quadrature, and evaluates
the original digital-strip samples by an integer recurrence.
"""
from concurrent.futures import ProcessPoolExecutor
import ast
from fractions import Fraction as F
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from math import gcd
from pathlib import Path

REFEREE=Path(__file__).with_name("lrc_label12_dirichlet_referee_thm3171.py")
CORE=Path(__file__).with_name("lrc_label12_multiblock_core_thm3171.py")
EXPECTED_REFEREE_SHA256="ae4669304e0d5a81ea6eaf75f385e2294fba6899ad72938c3ded8736e3ff5007"
EXPECTED_CORE_SHA256="7cff3a20e4d874f1143d096d3c9f4fa37448eb76c0ae0408513dd81780d8382b"
LABELS=((1,12),(12,1),(2,12),(12,2))
P_MIN=5001
MAX_D=8
MAX_C=41
WORKERS=8

def require(condition,detail):
 if not condition:raise RuntimeError(detail)

def file_hash(path):return sha256(path.read_bytes()).hexdigest()
require(file_hash(REFEREE)==EXPECTED_REFEREE_SHA256,
        ("referee hash",file_hash(REFEREE),EXPECTED_REFEREE_SHA256))
require(file_hash(CORE)==EXPECTED_CORE_SHA256,
        ("core hash",file_hash(CORE),EXPECTED_CORE_SHA256))
SPEC=spec_from_file_location("label12_referee",REFEREE)
REF=module_from_spec(SPEC);SPEC.loader.exec_module(REF)

def ceildiv(a,b):return -((-a)//b)
def ceilf(x):return -((-x.numerator)//x.denominator)

# ---- Independent one-tent geometry and integration. ----

def tent_breaks(z,w):
 z=F(z);w=F(w)
 outer=(w+z)/(14*z);inner=(w-z)/(14*z)
 require(0<=inner<outer<F(1,2),("tent geometry",z,w,inner,outer))
 return tuple(sorted({F(0),inner,outer,1-outer,1-inner,F(1)}))

def tent(x,z,w):
 x=F(x)%1;u=min(x,1-x);z=F(z);w=F(w)
 outer=(w+z)/(14*z);inner=(w-z)/(14*z)
 if u<=inner:return 24/w
 if u>=outer:return F(0)
 return F(168,w)*(outer-u)

def tent_mean(z,w):
 # Plateau plus two triangular shoulders, derived without PeriodicPL.
 z=F(z);w=F(w)
 outer=(w+z)/(14*z);inner=(w-z)/(14*z)
 return 2*(F(24,w)*inner+F(24,w)*(outer-inner)/2)

def tent_prefix(x,z,w):
 """Integral on [0,x], using exact trapezoids on the independent cuts."""
 x=F(x);require(0<=x<=1,("prefix",x))
 cuts=list(tent_breaks(z,w));total=F(0)
 for left,right in zip(cuts,cuts[1:]):
  if x<=left:break
  stop=min(x,right)
  total+=(tent(left,z,w)+tent(stop,z,w))*(stop-left)/2
  if stop==x:break
 return total

def tent_interval(start,length,z,w):
 start=F(start)%1;length=F(length)
 require(0<=length<=1,("interval",start,length))
 if start+length<=1:return tent_prefix(start+length,z,w)-tent_prefix(start,z,w)
 return tent_mean(z,w)-tent_prefix(start,z,w)+tent_prefix(start+length-1,z,w)

def g_value(x,z,w,alpha,d):
 return sum((tent(x+r*alpha,z,w) for r in range(d)),F(0))

def g_path_integral(x,z,w,alpha,d,forward,T):
 """Integral of G(x +/- t), independently as a sum of one-tent integrals."""
 m=T.numerator//T.denominator;theta=T-m
 total=m*d*tent_mean(z,w)
 start=x if forward else x-theta
 for r in range(d):total+=tent_interval(start+r*alpha,theta,z,w)
 return total

def tent_difference_sup(z,w,z0,w0):
 cuts={F(0)}
 for b in tent_breaks(z,w)[:-1]:cuts.add(b)
 for b in tent_breaks(z0,w0)[:-1]:cuts.add(b)
 return max(abs(tent(x,z,w)-tent(x,z0,w0)) for x in cuts)

def path_difference_sup(x,z,w,alpha,x0,z0,w0,alpha0,d,forward):
 """Exact sup over one period: a continuous PL difference peaks at a cut."""
 sign=1 if forward else -1;cuts={F(0)}
 for r in range(d):
  base=x+r*alpha
  for b in tent_breaks(z,w)[:-1]:cuts.add(((b-base)*sign)%1)
  base0=x0+r*alpha0
  for b in tent_breaks(z0,w0)[:-1]:cuts.add(((b-base0)*sign)%1)
 def difference(t):
  return g_value(x+sign*t,z,w,alpha,d)-g_value(x0+sign*t,z0,w0,alpha0,d)
 return max(abs(difference(t)) for t in cuts)

# ---- Independent channel reconstruction. ----

def independent_tasks():
 rows=[]
 for d in range(1,MAX_D+1):
  for a in range(d+1):
   for c in tuple(range(-MAX_C,0))+tuple(range(1,MAX_C+1)):
    for residue in range(d):
     if (a*residue+c)%d:continue
     p=P_MIN+((-P_MIN+residue)%d)
     gap=(a*p+c)//d
     if not (gap>=1 and p+gap<2*p):continue
     p0=d if residue==0 else residue
     n0=(p-p0)//d;q=p+gap;D=d+a;q0=q-D*n0
     for e,f in LABELS:
      A=168*c+e*D-d*f
      require(A!=0,("zero A",d,a,c,residue,e,f))
      rows.append((d,a,c,p0,q0,e,f,A,n0))
 return tuple(rows)

# ---- Independent digital-strip recurrence. ----

def literal_masses(e,p,f,q,complete_count):
 """Return complete-prefix and full masses from integer tent numerators."""
 z=168*p-e;w=168*q-f;R=(90*e)%168;S=(90*f)%168
 det=R*w-S*z;require(det%168==0,("phase divisibility",e,p,f,q,det))
 phase=(det//168)%z;prefix_num=0;full_num=0
 for k in range(p):
  distance=min(phase,z-phase);scaled=14*distance
  if scaled<=w-z:value_num=24*z
  elif scaled<w+z:value_num=12*(w+z-scaled)
  else:value_num=0
  full_num+=value_num
  if k<complete_count:prefix_num+=value_num
  phase=(phase+w)%z
 return F(prefix_num,w*z),F(full_num,w*z)

def independent_constants(task,p):
 """Re-derive every continuum-minus-1/P component at a channel point."""
 d,a,c,p0,q0,e,f,A,n0=task;D=d+a
 require((p-p0)%d==0,("channel residue",task,p))
 n=(p-p0)//d;q=q0+D*n
 require(p<q<2*p,("channel levels",task,p,q))
 Z0=F(168);W0=F(168*D,d);Z=F(168*p-e,p);W=F(168*q-f,p)
 s=F(168*c,d)-f;P0=p0+d*n0
 Zmin=Z0-F(e,P0);Wmin=W0-F(abs(s),P0)

 inv=abs(s)/(Wmin*W0);kH=24*inv;kL=168*inv
 kRatio=abs(168*s+W0*e)/(168*Zmin);ka=kRatio/14
 Wmax=W0+F(abs(s),P0);amax=max(F(1+F(D,d),14),F(1+Wmax/Zmin,14))
 L0=168/W0;kOuter=kL*amax+L0*ka;kF=max(kH,kOuter)
 require(tent_difference_sup(Z,W,Z0,W0)<=kF/p,("tent perturbation",task,p))

 Kalpha=F(abs(A),d*Zmin);KG=d*kF+L0*Kalpha*d*(d-1)/2
 R=(90*e)%168;S=(90*f)%168;phase0=R*F(D,d)-S
 Kx=abs(R*s+S*e+e*phase0)/(168*Zmin);Kpath=KG+d*L0*Kx

 z=168*p-e;w=168*q-f;alpha=F(w-z,z);alpha0=F(a,d)%1
 require(d*alpha-a==F(A,z),("alpha identity",task,p,d*alpha-a,F(A,z)))
 require(abs(alpha-F(a,d))<=Kalpha/p,("alpha perturbation",task,p))
 det=R*w-S*z;require(det%168==0,("phase",task,p))
 x=F(det//168,z)%1;X=(F(R*D,d)-S)/168%1
 require(path_difference_sup(x,Z,W,alpha,X,Z0,W0,alpha0,d,A>0)<=Kpath/p,
         ("path perturbation",task,p))

 residue=p0%d;blocks=p//d;rho=F(abs(A),z)
 require(rho<F(1,2),("rho",task,p,rho))
 beta=(d*alpha)%1
 require(beta==(rho if A>0 else 1-rho),("oriented beta",task,p,beta,rho))
 T=blocks*rho;T0=F(abs(A),168*d)
 Kt=F(abs(A)*abs(168*residue-e),168*d)/Zmin
 require(abs(T-T0)<=Kt/p,("T perturbation",task,p,T,T0,Kt/p))
 Tmax=T0+Kt/P0;Hmax=24/Wmin

 J=g_path_integral(x,Z,W,alpha,d,A>0,T)
 J0=g_path_integral(X,Z0,W0,alpha0,d,A>0,T0)
 KJ=Tmax*Kpath+Kt*d*Hmax
 require(abs(J-J0)<=KJ/p,("path integral perturbation",task,p,abs(J-J0),KJ/p))
 continuum=F(168,abs(A))*J0

 Kmain=F(168,abs(A))*KJ+F(e,abs(A))*Tmax*d*Hmax
 main=Z*J/abs(A)
 require(abs(main-continuum)<=Kmain/p,
         ("main perturbation",task,p,abs(main-continuum),Kmain/p))
 Kend=d*Hmax/2

 theta=T-(T.numerator//T.denominator);xn=x+theta if A>0 else x-theta
 endpoint=(g_value(x,Z,W,alpha,d)-g_value(xn,Z,W,alpha,d))/(2*p)
 require(abs(endpoint)<=Kend/p,("endpoint",task,p,abs(endpoint),Kend/p))

 Ktrap=F(abs(A)*84*d*ceilf(Tmax),1)/(Zmin*Wmin*P0)
 raw_trap=rho*ceilf(T)*F(672*d,W)/(8*p)
 require(raw_trap<=Ktrap/p,("trap envelope",task,p,raw_trap,Ktrap/p))
 K=Kmain+Kend+Ktrap
 return q,blocks*d,continuum,K,continuum-K/P0,main,endpoint,raw_trap

def audit_channel(task):
 d,a,c,p0,q0,e,f,A,n0=task;P0=p0+d*n0
 data=independent_constants(task,P0)
 q,complete_count,continuum,K,certified,main,endpoint,raw_trap=data
 target=REF.limit_and_error(task)
 require((continuum,K,certified)==target[1:],
         ("referee constant mismatch",task,(continuum,K,certified),target[1:]))
 prefix,full=literal_masses(e,P0,f,q,complete_count)
 require(prefix>=certified,("literal prefix below certificate",task,P0,q,prefix,certified))
 require(full>=prefix,("nonnegative residual",task,P0,q,prefix,full))
 require(prefix>F(1,105),("literal prefix threshold",task,P0,q,prefix))
 require(abs(prefix-main-endpoint)<=raw_trap,
         ("literal Peano envelope",task,P0,q,abs(prefix-main-endpoint),raw_trap))
 return (prefix-certified,full-F(1,105),task,P0,q,prefix,full,certified)

def dirichlet_choice(p,h):
 candidates=[]
 for d in range(1,MAX_D+1):
  for a in range(d+1):
   c=d*h-a*p
   if 9*abs(c)<=p:candidates.append((abs(c),d,a,c))
 require(candidates,("Dirichlet existence",p,h))
 return min(candidates)

def exhaustive_coverage():
 taskkeys={(t[0],t[1],t[2],t[3],t[5],t[6]) for t in independent_tasks()}
 pairs=0;small=0;large=0
 for p in range(P_MIN,P_MIN+48):
  for h in range(1,p):
   q=p+h
   if gcd(p,q)!=1:continue
   _,d,a,c=dirichlet_choice(p,h)
   require(c!=0,("primitive zero resonance",p,q,d,a))
   p0=d if p%d==0 else p%d
   for e,f in LABELS:
    A=168*c+e*(d+a)-d*f;z=168*p-e;n=p//d
    require(2*abs(A)<z,("rho orientation coverage",p,q,d,a,c,e,f,A))
    if n*abs(A)<5*z:
     small+=1
     require(abs(c)<=MAX_C,("small c coverage",p,q,d,a,c,e,f,A))
     require((d,a,c,p0,e,f) in taskkeys,
             ("missing channel",p,q,d,a,c,p0,e,f))
    else:large+=1
   pairs+=1
 return pairs,small,large

def hostile_large_integer_coverage():
 state=0x6A09E667F3BCC909;controls=0;small=0
 keys={(t[0],t[1],t[2],t[3],t[5],t[6]) for t in independent_tasks()}
 mask=(1<<64)-1
 for _ in range(20000):
  state=(6364136223846793005*state+1442695040888963407)&mask
  p=P_MIN+state%10**9
  state=(6364136223846793005*state+1442695040888963407)&mask
  h=1+state%(p-1)
  while gcd(p,h)!=1:
   h=1+(h%(p-1))
  _,d,a,c=dirichlet_choice(p,h)
  require(c!=0,("random primitive zero",p,h,d,a))
  p0=d if p%d==0 else p%d
  for e,f in LABELS:
   A=168*c+e*(d+a)-d*f;z=168*p-e;n=p//d
   if n*abs(A)<5*z:
    small+=1
    require(abs(c)<=MAX_C and (d,a,c,p0,e,f) in keys,
            ("random channel miss",p,h,d,a,c,p0,e,f,A))
  controls+=1
 return controls,small

def independent_large_floor():
 P=P_MIN;R=F(9,80)
 main=F(24*(P-7),1176*P)*F(4,5)
 endpoint=F(96,168*(P+1)-12)
 trap=F(84,168*(P+1)-12)*(P*R*R+8*R)
 value=main-endpoint-trap
 require(value>F(1,105),("independent large floor",value))
 require(REF.large_T_certificate()[0]==value,
         ("large floor mismatch",value,REF.large_T_certificate()[0]))
 return value

def large_boundary_controls():
 """Literal controls nearest T=5, generated without the channel list."""
 floor=independent_large_floor();candidates=[]
 for d in range(1,MAX_D+1):
  for a in range(d+1):
   for c in tuple(range(-64,0))+tuple(range(1,65)):
    for residue in range(d):
     if (a*residue+c)%d:continue
     p=P_MIN+((-P_MIN+residue)%d);h=(a*p+c)//d;q=p+h
     if not (p<q<2*p and gcd(p,q)==1):continue
     for e,f in LABELS:
      A=168*c+e*(d+a)-d*f;z=168*p-e;n=p//d
      rho=F(abs(A),z)
      if rho<=F(9,80) and n*rho>=5:
       candidates.append((n*rho-5,d,a,c,p,q,e,f,n*d))
 candidates.sort()
 chosen=tuple(candidates[:512]);require(len(chosen)==512,("large controls",len(chosen)))
 rows=[]
 for excess,d,a,c,p,q,e,f,complete_count in chosen:
  prefix,full=literal_masses(e,p,f,q,complete_count)
  require(prefix>=floor,("large literal floor",excess,d,a,c,p,q,e,f,prefix,floor))
  rows.append((prefix-floor,excess,d,a,c,p,q,e,f))
 return tuple(rows)

def tail_controls(tasks):
 rows=[]
 # Deterministic spread through the task list, with several much later points.
 for index in range(0,len(tasks),23):
  task=tasks[index];d=task[0];p0=task[3]+d*task[8]
  for jump in (1,17,257):
   p=p0+d*jump
   q,complete_count,continuum,K,certified,main,endpoint,raw=independent_constants(task,p)
   prefix,full=literal_masses(task[5],p,task[6],q,complete_count)
   require(prefix>=certified,("tail certificate",task,p,q,prefix,certified))
   require(abs(prefix-main-endpoint)<=raw,
           ("tail Peano envelope",task,p,q,abs(prefix-main-endpoint),raw))
   rows.append((prefix-certified,task,p,q))
 return tuple(rows)

def structural_independent():
 perturb=max(abs(e*(d+a)-d*f) for d in range(1,9) for a in range(d+1)
             for e,f in LABELS)
 require(perturb==184,("perturb cap",perturb))
 rho=F(F(168,9)*P_MIN+perturb,168*P_MIN-12)
 require(rho<F(9,80),("rho cap",rho))
 strict=F(840*8*P_MIN,P_MIN-7);integer_cap=(strict.numerator-1)//strict.denominator
 require(integer_cap==6729,("A cap",strict,integer_cap))
 require(integer_cap+perturb<42*168,("c cap",integer_cap,perturb))
 return perturb,rho,strict,integer_cap

def main():
 for path in (REFEREE,CORE,Path(__file__)):
  tree=ast.parse(path.read_text(),filename=str(path))
  require(not any(isinstance(node,ast.Assert) for node in ast.walk(tree)),
          ("optimization-sensitive assert",str(path)))
 structure=structural_independent()
 ours=independent_tasks();theirs=tuple(REF.channel_tasks())
 require(len(ours)==11608,("independent channel count",len(ours)))
 require(set(ours)==set(theirs),("channel universe mismatch",len(set(ours)^set(theirs))))
 coverage=exhaustive_coverage();random_coverage=hostile_large_integer_coverage()
 large_controls=large_boundary_controls()
 with ProcessPoolExecutor(max_workers=WORKERS) as pool:
  literal=tuple(pool.map(audit_channel,ours,chunksize=4))
 tails=tail_controls(ours)
 min_cert=min(literal,key=lambda x:x[0]);min_threshold=min(literal,key=lambda x:x[1])
 min_tail=min(tails,key=lambda x:x[0])
 min_large=min(large_controls,key=lambda x:x[0])
 semantic=repr((structure,coverage,random_coverage,min_cert,min_threshold,min_tail,
                min_large,len(literal),len(tails),len(large_controls))).encode()
 print("LRC LABEL-12 DIRICHLET REFEREE -- INDEPENDENT EXACT AUDIT")
 print(f"channels={len(ours)};literal_minimal_samples="
       f"{sum(task[3]+task[0]*task[8] for task in ours)}")
 print(f"coverage_pairs={coverage[0]};coverage_small_label_cases={coverage[1]};"
       f"coverage_large_label_cases={coverage[2]};random_pairs={random_coverage[0]};"
       f"random_small_label_cases={random_coverage[1]}")
 print(f"tail_controls={len(tails)}")
 print(f"near_T5_large_controls={len(large_controls)};"
       f"tightest_large_literal_minus_floor={min_large[0]};case={min_large[1:]}")
 print(f"tightest_literal_minus_certificate={min_cert[0]};case={min_cert[2:5]}")
 print(f"tightest_literal_minus_threshold={min_threshold[1]};case={min_threshold[2:5]}")
 print(f"tightest_tail_minus_certificate={min_tail[0]};case={min_tail[1:]}")
 print(f"semantic_sha256={sha256(semantic).hexdigest()}")

if __name__=="__main__":main()
