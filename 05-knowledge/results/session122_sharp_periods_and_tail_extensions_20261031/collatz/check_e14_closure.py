#!/usr/bin/env python3
from fractions import Fraction as F
from pathlib import Path
import json

def need(x,m):
 if not x:raise ValueError(m)
def lam(e):return F(5,8)-F(1,8*5**(e-1))
def ones(E,k):return sum(k//(2*5**(e-1)) for e in range(1,E+1))
def digits(n):
 s=0
 while n:s+=n%5;n//=5
 return s
controls=0
for E in range(2,11):
 for ell in range(1,E+1):
  T=2*5**(ell-1);M=5**(E-ell+1)
  for b in sorted(set(list(range(101))+[M-1,M,M+1,2*M,5*M+7]+[5**k for k in range(E-ell+2)])):
   for s in sorted(set([0,1,T//2,T-1])):
    k=b*T+s
    loss=(5*k-8*ones(E,k))-(5*s-8*ones(E,s))
    need(loss==2*digits(b%M)+10*(b//M),'phase digit identity')
    need((loss<=2)==(b==0 or b in [5**j for j in range(E-ell+1)]),'low-slack phase bank')
    controls+=1
receipt=json.loads((Path(__file__).parent/'near_tight10_cpp.json').read_text())
need((receipt['E'],receipt['R'],receipt['budget'],receipt['terminal'])==(10,2,2,-30),'receipt domain')
need(receipt['acyclic'],'DAG')
need(receipt['nodes']==3908743 and receipt['edges']==3908754,'complete graph census')
need(receipt['minimum_potential']==-128 and receipt['maximum_potential']==-10,'potential range')
need(-10+80-325<=-128-3 and 80-325<=-128-3,'strict omitted letters')
expected={30:[(342,326)],31:[],32:[(A,A-16) for A in [72,352,392,552,592,1592,6592,31592,156592,781592,3906592]]}
pairs=0
heights={}
for row in receipt['records']:
 actual=[tuple(p) for p in row['pairs']]
 need(actual==expected[row['score']],'exact clock pair list')
 for A,L in actual:
  need(1<=L<A,'positive word')
  need(pow(3,L,625)!=pow(2,A,625),'actual endpoint exclusion modulo625')
  h=next(e-1 for e in range(1,11) if pow(3,L,5**e)!=pow(2,A,5**e))
  heights[h]=heights.get(h,0)+1
  pairs+=1
need(heights=={2:10,3:2},'exact multiplier heights')
eta=F(183105426,48828125);u=F(33,8)
need(lam(14)==F(24,25)*F(5,8)+lam(12)/25,'slope interpolation')
eta14=F(24,25)*u+eta/25;alpha=eta14-4
need(eta14==F(5017089801,1220703125),'new terminal fee')
need(alpha==F(134277301,1220703125)>0,'positive margin')
for e,fee in [(11,F(31,125)),(12,F(6,25)),(13,F(1,5))]:
 need(lam(14)*(2*5**(e-1))-ones(e,2*5**(e-1))==fee,'all-one constituent fee')
for C in range(1,101):
 for B in range(1,C+1):need((eta14-1)*C-B-1>=eta14-3,'native core invoice')
need(lam(14)==F(762939453,1220703125),'sharp coefficient')
T=2*5**13;S=ones(14,T)
need(T==2441406250 and S==1525878906 and F(S,T)==lam(14),'all-one sharpness')
boundary_controls=0
for x in range(5):
 for y in range(5):
  phi=F((-x)%5,5)
  for dx in range(5):
   for dy in range(5):
    X=(x+dx)%5;Y=(y+dy)%5
    bill=int(X==Y==0)+F((-X)%5,5)-phi
    need(bill<=1,'arbitrary tangent invoice')
    if (dx,dy)==(4,4):need(bill<=F(1,5),'all-one tangent invoice')
    boundary_controls+=1
q=5**14
rho=3*pow(2,-1,q)%q
need(pow(rho,2*5**12,q)==1+4*5**13,'all-one13 tangent')
gamma12=F(61035156,48828125)
need(min(alpha,gamma12)==alpha,'unfinished suffix margin')
print(json.dumps(dict(status='PASS',scope='arithmetic and complete-pair-list controls; full potential and graph replays are separate',
 phase_identity_controls=controls,boundary_tangent_controls=boundary_controls,endpoint_pairs_excluded=pairs,multiplier_height_census=heights,actual_identity10_H_lower_bound=str(u),
 eta14=str(eta14),nonone_margin=str(alpha),sharp_S14_coefficient=str(lam(14)),sharp_all_one_clock=T,
 inherited_E10_R1_certificate_rerun=False,all_layer_bound='OPEN',Collatz_convergence='OPEN'),indent=2,sort_keys=True))
