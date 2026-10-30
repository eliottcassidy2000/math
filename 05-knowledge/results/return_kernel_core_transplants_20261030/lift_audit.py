#!/usr/bin/env python3
"""Exact universal small-kernel inequalities, ninth-to-tenth lift and next reduction."""
from fractions import Fraction as F
from pathlib import Path
import random,json
P=Path(__file__).resolve().parent
def need(ok,msg):
 if not ok:raise RuntimeError(str(msg))
def add(x,y):return 5*((x//5+y//5)%5)+(x%5+y%5)%5
def neg(x):return 5*(-(x//5)%5)+(-(x%5))%5
def conv(f,g):return [min(f[x]+g[add(z,neg(x))] for x in range(25)) for z in range(25)]
def shortest(g):
 c=g[:];c[24]=min(c[24],F(1,5));dist=[None]*25;dist[0]=F(0);done=set()
 for _ in range(25):
  u=min((x for x in range(25) if x not in done and dist[x] is not None),key=lambda x:dist[x]);done.add(u)
  for j in range(25):
   v=add(u,j);d=dist[u]+c[j]
   if dist[v] is None or d<dist[v]:dist[v]=d
 return min(dist[neg(x)]+c[x] for x in range(25))
def nine(g):
 pools={1:g,2:conv(g,g)};pools[3]=conv(pools[2],g);answers=[]
 for b in (1,2,3):
  for a in range(5-b):
   target=5*(a%5)+(a%5) # -a*(4,4)
   answers.append(dict(ones=a,other=b,cost=str(pools[b][target]+F(a,5))))
 return answers

def main():
 # Uniform potential: the raw five-adic state still has25 possibilities.
 phi=lambda x:F((-(x//5))%5,5);potential=0
 for u in range(25):
  v=add(u,24);need(F(v==0)+phi(v)-phi(u)<=F(1,5),'all-ones edge');potential+=1
  for v in range(25):
   need(F(v==0)+phi(v)-phi(u)<=1,'arbitrary non-one edge');potential+=1
 l9=F(488281,781250);l10=F(1220703,1953125);l11=F(5,8)*(1-F(1,5**11))
 delta=l10-l9;gamma10=2*l10-1
 need(2*l9>1 and delta==F(1,2*5**9),'one-fee lift hypothesis')
 need(gamma10==F(488281,1953125),'retained tenth-layer slack')
 need(gamma10+3*(l11-l10)>F(1,4),'eleventh non-one quarter-fee')
 T=2*5**9;charge=sum(T//(2*5**(e-1)) for e in range(1,11))
 need(charge==2441406==l10*T,'sharp all-one equality')
 # Boundary digit(4,4) at all levels, checked by independent modular powers.
 for e in range(1,16):
  q=5**(e+1);t=2*5**(e-1);m=pow(3,t,q)*pow(pow(2,t,q),-1,q)%q
  need((m-1)//5**e==4,'all-ones tangent lifting')
 rng=random.Random(20261030);tables=[[F(1)]*25,[F(1,4)]*25]
 for _ in range(70):tables.append([F(25+rng.randrange(101),100) for _ in range(25)])
 # Positive/hostile controls with high identity cost and a single off-identity cheap type.
 for x in range(1,25):
  g=[F(2)]*25;g[x]=F(3,10);tables.append(g)
 accepted=0
 for g in tables:
  ns=nine(g);verdict=all(F(r['cost'])>=1 for r in ns);need((shortest(g)>=1)==verdict,'nine-family equivalence');accepted+=verdict
 out=dict(status='exact algebraic controls; written universal lifting proof supplies quantifiers',paid_potential_edges=potential,
          all_ones_digits=15,abstract_kernels=len(tables),accepted_kernels=accepted,eleventh_short_families=nine([F(1)]*25),
          lambda9=str(l9),lambda10=str(l10),gamma10=str(gamma10),sharp_clock=T,sharp_charge=charge)
 (P/'lift_audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='eleventh_short_families'},indent=2))
if __name__=='__main__':main()
