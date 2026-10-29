#!/usr/bin/env python3
"""Independent finite CRT, clock-tree, and pointed SCC checks (stdlib only)."""
from itertools import product,permutations,combinations
from fractions import Fraction
from collections import Counter
from math import comb,factorial,gcd
import json
checks=0
def check(test,why):
 global checks
 checks+=1
 if not test:raise RuntimeError(why)
def fib(n):
 a,b=0,1
 for _ in range(n):a,b=b,a+b
 return a
def crt(n):
 a,b=fib(n-1)-1,fib(n);g=gcd(a,b);M=abs(a*a+a*b-b*b)//g**2;r=(-a//g*pow(b//g,-1,M))%M
 # Independent recursive independent-set generation by appending 0 or 10.
 def linear(k):
  if k==0:return [()]
  if k==1:return [(0,),(1,)]
  return [(0,)+w for w in linear(k-1)]+[(1,0)+w for w in linear(k-2)]
 words=[w for w in linear(n) if not(w[0] and w[-1])]
 def read(w):
  t=0;u=0
  for i,x in enumerate(w):
   if x:t^=(1,2,3)[i%3];u=(u+pow(r,i,M))%M
  return (t if g==2 else 0,u)
 keys=[read(w) for w in words];lookup={x:i for i,x in enumerate(keys)}
 check(len(lookup)==M*g*g,'CRT bijection')
 reverse=[lookup[read(w[:1]+w[:0:-1])] for w in words]
 trans=[lookup[(t^1 if g==2 else 0,(u+1)%M)] for t,u in keys]
 mods=[M] if n==9 else [31,11,M]
 out=[]
 for modulus in mods:
  obs=[u%modulus for t,u in keys]
  def refine(items):
   d={};return [d.setdefault(x,len(d)) for x in items]
  colors=refine(obs);sizes=[len(set(colors))]
  while True:
   nxt=refine([(colors[i],colors[reverse[i]],colors[trans[i]]) for i in range(len(keys))])
   if nxt==colors:break
   colors=nxt;sizes.append(len(set(colors)))
  triple={(obs[i],obs[reverse[i]],obs[reverse[trans[i]]]) for i in range(len(keys))}
  check(sizes[-1]==len(keys),'full TR refinement finite')
  out.append({'modulus':modulus,'sizes':sizes,'triple':len(triple)})
 return {'n':n,'size':len(keys),'charts':out}
def sparse_crt():
 N=167761;r=135722;n=25
 pairs=[(11,1,21495840),(1111,0,9049157),(1661,0,8692997),(101,0,2392656),(151,8388610,22185088),(15251,1,9471010)]
 check({d for d,_,_ in pairs}=={d for d in range(2,N) if N%d==0},'all proper divisor blocks')
 def value(mask):
  bits=tuple((mask>>i)&1 for i in range(n))
  check(all(not(bits[i] and bits[(i+1)%n]) for i in range(n)),'legal sparse word')
  x=0;y=0
  for bit in bits[::-1]:x=(r*x+bit)%N
  rev=bits[:1]+bits[:0:-1]
  for bit in rev[::-1]:y=(r*y+bit)%N
  return x,y
 for d,a,b in pairs:
  x,rx=value(a);y,ry=value(b)
  check((x-y)%d==0 and (rx-ry)%d!=0,'sparse reflection breaker')
 check(pow(37,55,121)==1 and all(pow(37,e,121)!=1 for e in (1,5,11)),'prime power clock55')
 check(pow(4,5,11)==1 and 11*11%121==0,'reduction loses clock and nilpotent')
 return {'period':25,'state_count':N,'divisor_blocks_excluded':len(pairs)}

def sparse55():
 N=312119004989;n=55;a,b=fib(54)-1,fib(55);r=(-a*pow(b,-1,N))%N
 check(abs(a*a+a*b-b*b)==N,'period55 norm')
 factors=(11,11,199,331,39161)
 from math import prod,isqrt
 check(prod(factors)==N and all(all(p%d for d in range(2,isqrt(p)+1)) for p in factors),'period55 prime factorization')
 pairs=[(11, 1, 12, 1, 23080948090273824), (121, 0, 121, 0, 18049582881571496), (199, 1, 200, 1, 24004537857607680), (331, 0, 331, 0, 11267795161452562), (2189, 0, 2189, 0, 11541023800994389), (3641, 0, 3641, 0, 20922744204034376), (24079, 1, 24080, 1, 9045304205984770), (39161, 0, 39161, 0, 5950565522248196), (40051, 2, 40053, 9007199254740994, 20310187382412296), (65869, 0, 65869, 0, 11580419396899412), (430771, 1, 430772, 1, 11964931256631634), (724559, 0, 724559, 0, 2851068548105297), (4738481, 0, 4738481, 0, 1482434439317520), (7793039, 0, 7793039, 0, 11971519802511701), (7970149, 0, 7970149, 0, 11968538407866709), (12962291, 0, 12962291, 0, 5790600041407621), (85723429, 0, 85723429, 0, 10417411502123541), (142585201, 0, 142585201, 0, 23679117078069888), (942957719, 0, 942957719, 0, 9183484182827266), (1568437211, 1, 1568437212, 1, 9024800034916354), (2579495909, 0, 2579495909, 0, 11417604999255058), (28374454999, 1, 28374455000, 1, 10169418513417250)]
 divisors={1}
 for p in factors:divisors|={d*p for d in tuple(divisors)}
 check({z[0] for z in pairs}==divisors-{1,N},'all22 divisor blocks')
 def val(mask):
  bits=[(mask>>i)&1 for i in range(n)]
  check(mask<1<<n and all(not(bits[i] and bits[(i+1)%n]) for i in range(n)),'period55 legal word')
  x=y=0
  for v in bits[::-1]:x=(r*x+v)%N
  for v in (bits[:1]+bits[:0:-1])[::-1]:y=(r*y+v)%N
  return x,y
 for d,x,y,wx,wy in pairs:
  u,ru=val(wx);v,rv=val(wy)
  check((u,v)==(x,y) and (u-v)%d==0 and (ru-rv)%d!=0,'period55 sparse witness')
 return {'states':N,'divisor_blocks':22,'enumeration':False}

def clock_tree():
 # Count tree shapes directly by slot evolution; valuation lift by weak compositions.
 out=[]
 for m in range(6):
  N=3*m+1;valid=0
  for internals in combinations(range(N),m):
   S=set(internals);slots=1
   for i in range(N):
    if slots==0:break
    slots+=2 if i in S else -1
   else:
    if slots==0:valid+=1
  check(valid==comb(3*m+1,m)//(3*m+1),'ternary cycle lemma')
  # For every specified shape, valuation-parity source has mass(1/3)^m(2/3)^(2m+1).
  mass=Fraction(valid*2**(2*m+1),3**N)
  mean=Fraction(4*m+1)+Fraction(2*N,3)
  check(mean==2*N-Fraction(1,3),'conditional clock mean')
  out.append({'N':N,'shapes':valid,'mass':str(mass),'mean_clock':str(mean)})
 def odd_step(x):
  y=3*x+1;a=0
  while y%2==0:y//=2;a+=1
  return y,a
 check(odd_step(-5)==(-7,1) and odd_step(-7)==(-5,2),'phase hostile')
 return out

def literal_scc():
 rows=[]
 for n in range(1,6):
  es=list(combinations(range(n),2)); perms=list(permutations(range(n)));total=[0]*4;strong=[0]*4;occupancy=[0]*3
  for mask in range(1<<len(es)):
   arcs={(a,b) if mask>>i&1 else (b,a) for i,(a,b) in enumerate(es)}
   reach=[[i==j or (i,j) in arcs for j in range(n)] for i in range(n)]
   for k in range(n):
    for i in range(n):
     for j in range(n):reach[i][j]|=reach[i][k] and reach[k][j]
   comp=[min(j for j in range(n) if reach[i][j] and reach[j][i]) for i in range(n)]
   for p in perms:
    if any((p[a],p[b]) not in arcs for a,b in arcs):continue
    fix=[i for i in range(n) if p[i]==i]
    for j in range(4):
     total[j]+=len(fix)**j
     if len(set(comp))==1:strong[j]+=len(fix)**j
    for a,b in product(fix,repeat=2):occupancy[1 if comp[a]==comp[b] else 2]+=1
  check(all(x%factorial(n)==0 for x in total+strong+occupancy),'labelled Burnside integral')
  total=[x//factorial(n) for x in total];strong=[x//factorial(n) for x in strong];occupancy=[x//factorial(n) for x in occupancy]
  check(total[0]==[1,1,2,4,12][n-1],'A000568 literal')
  check(sum(occupancy)==total[2],'joint occupancy conservation')
  rows.append({'n':n,'ports':total,'strong_ports':strong,'two_port_occupancy':occupancy})
 check(rows[1]['two_port_occupancy']==[0,2,2],'same different SCC hostile')
 return rows
if __name__=='__main__':
 report={'crt':[crt(9),crt(15)],'clock_tree':clock_tree(),'sparse_crt':sparse_crt(),'sparse55':sparse55(),'scc':literal_scc()};report['checks']=checks
 print(json.dumps(report,sort_keys=True,indent=2))
