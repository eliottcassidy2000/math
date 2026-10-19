#!/usr/bin/env python3
"""Exact source controls for fixed-translation isolated Collatz collisions.

Reuses the finite atlas already frozen in integer_stabilizer_20261019.json;
this is not a new residue census. All assertions are explicit under -O.
"""
from fractions import Fraction as Q
from pathlib import Path
import json
CHECKS=0
def check(test,label):
 global CHECKS
 CHECKS+=1
 if not test:raise RuntimeError(label)
def T(x):return (3*x+1)//2 if x%2 else x//2

def v3(x):
 if x==0:raise ValueError('valuation of zero')
 x=Q(x);n=abs(x.numerator);d=x.denominator;v=0
 while n%3==0:n//=3;v+=1
 while d%3==0:d//=3;v-=1
 return v

def minimum_K(N,d):
 K=1
 while N+1>2*K*(1+((d+1)*3**K).bit_length()-1):K+=1
 return K

def analyse(x,d,N):
 u=x+d;v=x;a=b=Cu=Cv=0;Pu=Pv=Q(1)
 history=[];oddu=[];oddv=[];runs=[];run=0;active=0
 for t in range(N+1):
  k=a-b
  e=Q(u)-Q(3)**k*v
  E=3**max(0,-k)*e
  check(E.denominator==1,'normalized carry is an integer')
  check(2**t*u==3**a*(x+d)+Cu and 2**t*v==3**b*x+Cv,'actual affine numerators')
  if k!=0:check(v3(e)==min(0,k),'exact 3-adic translation valuation')
  history.append((t,u,v,a,b,e,int(E)))
  if t==N:break
  alpha=u%2;beta=v%2
  if alpha:oddu.append((t,u));Pu*=1+Q(1,3*u);a+=1;Cu=3*Cu+2**t
  if beta:oddv.append((t,v));Pv*=1+Q(1,3*v);b+=1;Cv=3*Cv+2**t
  if alpha or beta:active+=1;runs.append(run);run=0
  else:run+=1
  u=T(u);v=T(v)
 runs.append(run)
 check(u==v and a!=b,'requested isolated literal collision')
 K=max(a,b);m=min(a,b);B=(d+1)*3**K;L=B.bit_length()-1
 for t,U,V,aa,bb,e,E in history:
  check(E!=0,'isolated histories never hit zero translation')
  check(abs(E)<B,'uniform normalized carry height')
 check(active<=a+b<=2*K-1,'active step budget')
 check(len(runs)<=2*K,'common-even block budget')
 check(all(r<=L for r in runs),'bounded common-even run')
 check(N+1<=2*K*(1+L),'fixed-translation depth bound')
 check(x<Q(d,2)+Q(2**N,2**K),'odd-count-refined isolated height')
 Kmin=minimum_K(N,d)
 check(K>=Kmin,'exact integer inversion of depth bound')
 check(x<Q(d,2)+Q(2**N,2**Kmin),'fixed-d subdyadic head')
 check(Pv/Pu==Q(3)**(a-b)*Q(x+d,x),'exact multiplicative correction ratio')
 ratio=Pv/Pu
 certificate=None
 if ratio!=1:
  if ratio>1:path='lower';odd=oddv;count=b;r=ratio
  else:path='upper';odd=oddu;count=a;r=1/ratio
  bound=Q(count)*r/(3*(r-1))
  odd_time,odd_z=min(odd,key=lambda pair:pair[1])
  check(odd_z<bound,'rational product-forced small odd value')
  if bound<=x:
   column=2 if path=='lower' else 1
   first=next(h for h in history if h[column]<x)
   time,z=first[0],first[column]
   check(0<z<x and 0<time<N,'genuine premerge smaller source')
   endpoint=z
   for _ in range(N-time):endpoint=T(endpoint)
   check(endpoint==u,'clocked smaller-source receipt')
   certificate={'path':path,'smaller_source':z,'source_time':time,
                'first_descent_source':z,'forced_odd_witness':odd_z,'forced_odd_time':odd_time,
                'initial_clocks':[N,N],'smaller_source_clock':N-time,
                'common_endpoint':u,'forced_odd_strict_bound':str(bound)}
 if x>Q(d,2):
  if a<b:
   bound=Q(a*x,2*x-d);path=oddu
  else:
   bound=Q(b*(x+d),2*x+3*d);path=oddv
  check(min(z for _,z in path)<bound,'sign-specific rational source bound')
  if m<=2*x-d:
   check(bound<=x,'count-aware source descent threshold')
   check(certificate is not None,'count-certified common smaller source')
 if N<=2*x-d:
  check(certificate is not None,'linear-core source descent certificate')
 return {'lower_source':x,'translation':d,'depth':N,'odd_counts':[a,b],
         'common_endpoint':u,'normalized_carry_maximum':max(abs(h[-1]) for h in history),
         'common_even_runs':runs,'carry_height_upper_bound':B,
         'minimum_K_from_depth':Kmin,'odd_correction_ratio':str(ratio),
         'count_certificate_guaranteed':m<=2*x-d,
         'linear_core_certificate_guaranteed':N<=2*x-d,
         'certificate':certificate}

parent=Path(__file__).resolve().parent
atlas_path=parent.parent/'results'/'integer_stabilizer_20261019.json'
if not atlas_path.exists(): atlas_path=parent/'integer_stabilizer_20261019.json'
atlas=json.loads(atlas_path.read_text())
all_rows=[];births={}
for row in atlas['branch_valid_isolated_collision_atlas']:
 N=row['N'];d=row['translation']
 for x in row['isolated_positive_integer_roots']:
  result=analyse(x,d,N);all_rows.append(result)
  births.setdefault((x,d),result)

# The changing-d sharpness family remains admissible: no fixed-d claim is
# inferred from it. Its very small odd counts survive because log(d) grows.
changing_d=[]
for row in atlas['integer_height_sharpness_family']:
 changing_d.append(analyse(row['lower_source'],row['translation'],row['N']))

# Hostile: an identity cylinder does hit E=0 and cannot use the nonzero-carry
# proof. The out-of-phase pair (1,2) never joins at equal clocks.
u,v=5,4
for _ in range(3):u,v=T(u),T(v)
check(u==v==2,'identity-cylinder hostile has zero normalized carry')
u,v=2,1
for _ in range(20):
 check(u!=v,'out-of-phase hostile');u,v=T(u),T(v)

# Exact rational substitute for logarithms, checked on modest independent
# rational parameters. The proof in the note covers every m>=1, 0<t<1.
for m in range(1,21):
 for numerator in range(1,20):
  t=Q(numerator,20)
  check((1+t/m)**m<1/(1-t),'finite geometric domination of correction product')

# Expose exact integer thresholds at large depths without enumerating sources.
thresholds=[]
for d in (1,3,30,61):
 for N in (12,32,64,128,256,1024,4096,16384):
  K=minimum_K(N,d)
  check(N+1<=2*K*(1+((d+1)*3**K).bit_length()-1),'threshold K qualifies')
  if K>1:check(N+1>2*(K-1)*(1+((d+1)*3**(K-1)).bit_length()-1),'threshold K is minimal')
  thresholds.append({'d':d,'N':N,'minimum_maximum_odd_count':K,
                     'isolated_head_exponent_N_minus_K':N-K})

report={'status':'FINITE-EXACT source controls; fixed-d theorems proved in companion note',
        'checks':CHECKS,'reused_atlas_rows':len(all_rows),'distinct_first_births':len(births),
        'first_births':[births[k] for k in sorted(births)],
        'first_births_with_count_guarantee':sum(r['count_certificate_guaranteed'] for r in births.values()),
        'first_births_with_linear_guarantee':sum(r['linear_core_certificate_guaranteed'] for r in births.values()),
        'first_births_with_rational_product_certificate':sum(r['certificate'] is not None for r in births.values()),
        'changing_d_sharpness_controls':changing_d,'exact_depth_thresholds':thresholds}
print(json.dumps(report,indent=2,sort_keys=True))
