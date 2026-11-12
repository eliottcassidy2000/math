#!/usr/bin/env python3
"""Standalone exact controls; explicit checks also execute under python -O."""
from itertools import combinations
from fractions import Fraction as F
from math import comb
import json

def check(x,msg):
 if not x:raise RuntimeError(msg)

def cyc(n,S):return [sum(1<<((i+d)%n)for d in S)for i in range(n)]
def attach(rows,C):
 n=len(rows);return [rows[i]|(0 if i in C else 1<<n)for i in range(n)]+[sum(1<<i for i in C)]
def hcount(rows):
 n=len(rows);N=1<<n;dp=[[0]*n for _ in range(N)]
 for i in range(n):dp[1<<i][i]=1
 for S in range(1,N):
  for i in range(n):
   if not S>>i&1:continue
   rem=((N-1)^S)&rows[i]
   while rem:
    bit=rem&-rem;rem-=bit;j=bit.bit_length()-1
    dp[S|bit][j]+=dp[S][i]
 return sum(dp[-1])
def coeff(rows):
 n=len(rows);N=1<<n;E=[[0]*n for _ in range(N)];S=[[0]*n for _ in range(N)]
 for i in range(n):E[1<<i][i]=S[1<<i][i]=1
 for mask in range(1,N):
  for v in range(n):
   if not mask>>v&1:continue
   rem=mask^(1<<v)
   if rem:
    E[mask][v]=sum(E[rem][u]for u in range(n)if rows[u]>>v&1)
    S[mask][v]=sum(S[rem][u]for u in range(n)if rows[v]>>u&1)
 H=[sum(row)for row in E];H[0]=1
 full=N-1;B=[[0]*n for _ in range(n)]
 for mask in range(1,N-1):
  other=full^mask
  for i in range(n):
   if E[mask][i]:
    for j in range(n):B[i][j]+=E[mask][i]*S[other][j]
 W=[[B[i][j]+B[j][i]if i!=j else 0 for j in range(n)]for i in range(n)]
 h=H[-1];Q=sum(H[m]*H[full^m]for m in range(N));return h,Q,W

def energy(C,h,r,W):return h+r*len(C)-sum(W[i][j]for i,j in combinations(C,2))
def corr(C,n):return tuple(sum((i+d)%n in C for i in C)for d in range(n))
def bc(n,k):return comb(n,k)if 0<=k<=n else 0

out={'small_layers':[],'literal_attachment_checks':0}
for n,S in [(5,[1,2]),(7,[1,2,3]),(7,[1,2,4]),(9,[1,7,6,5])]:
 rows=cyc(n,S);h,Q,W=coeff(rows);r=F(Q-2*h,n)
 check(r.denominator==1,'integral r');r=int(r)
 check(all(sum(row)==2*r for row in W),'row sum')
 for k in range(2,n-1):
  masks=[frozenset(C)for C in combinations(range(n),k)]
  vals={C:energy(C,h,r,W)for C in masks}
  if n<=7:
   for C in masks:
    check(vals[C]==hcount(attach(rows,C)),'literal attachment')
    out['literal_attachment_checks']+=1
  mu=F(h)+F(r*k*(n-k),n-1)
  check(sum(vals.values())==len(vals)*mu,'mean')
  b=[[F(W[i][j])-F(2*r,n-1)if i!=j else 0 for j in range(n)]for i in range(n)]
  for C in masks:
   drift=sum(vals[frozenset((C-{a})|{bb})]-vals[C]for a in C for bb in range(n)if bb not in C)
   check(drift==-2*(n-1)*(vals[C]-mu),'Johnson eigenlaw')
  for i,j in combinations(range(n),2):
   recovered=sum(mu-vals[C]for C in masks if i in C and j in C)/bc(n-4,k-2)
   check(recovered==b[i][j],'tomography')
  var=sum((v-mu)**2 for v in vals.values())/len(vals)
  rhs=F(k*(k-1)*(n-k)*(n-k-1),n*(n-1)*(n-2)*(n-3))*sum(b[i][j]**2 for i,j in combinations(range(n),2))
  check(var==rhs,'variance')
  out['small_layers'].append({'n':n,'increments':S,'k':k,'distinct_counts':len(set(vals.values())),'variance':str(var)})

n=13;rows=cyc(n,range(1,7));h,Q,W=coeff(rows);r=(Q-2*h)//n
A={0,1,3,9};B={1,2,5,7}
check(corr(A,n)==corr(B,n)==(4,)+(1,)*12,'old LRC homometry')
check(h==3711175 and Q==84542744,'carousel13 core')
check(W[0]==[0,1027338,1021076,1011596,991884,960716,919728,919728,960716,991884,1011596,1021076,1027338],'exact W')
base=energy(A,h,r,W);check(base==energy(B,h,r,W)==21508189,'twins')
records=[]
for v in range(n):
 if (v in A)!=(v in B):continue
 ap=A^{v};bp=B^{v};fa=energy(ap,h,r,W);fb=energy(bp,h,r,W)
 delta=sum(W[v][j]*((j in A)-(j in B))for j in range(n))
 check(fa-fb==-(1-2*(v in A))*delta,'physical phase breaker')
 records.append([v,fa,fb,fa-fb])
for C,expected in [(A,21508189),(B,21508189),(A|{8},23572029),(B|{8},23562137)]:
 check(hcount(attach(rows,C))==expected,'independent literal14')
# Carrier responses from direct path decomposition, independent of symmetry formula.
def direct_packet(rows,root):
 n=len(rows);N=1<<n;E=[[0]*n for _ in range(N)];S=[[0]*n for _ in range(N)]
 for i in range(n):E[1<<i][i]=S[1<<i][i]=1
 for mask in range(1,N):
  for v in range(n):
   if not mask>>v&1:continue
   rest=mask^(1<<v)
   if rest:
    E[mask][v]=sum(E[rest][u]for u in range(n)if rows[u]>>v&1)
    S[mask][v]=sum(S[rest][u]for u in range(n)if rows[v]>>u&1)
 H=[sum(row)for row in E];H[0]=1;full=N-1;D=full^(1<<root)
 u=v=d=e=0
 for mask in range(N):
  if not mask>>root&1:continue
  c=H[full^mask];sg=1 if mask.bit_count()%2 else -1
  u+=S[mask][root]*c;v+=E[mask][root]*c
  d+=sg*S[mask][root]*c;e+=sg*E[mask][root]*c
 q=sum(H[m]*H[D^m]for m in range(N)if m&D==m)
 return [H[full],H[D],q,u,v,d,e]
p0=direct_packet(attach(rows,A),13);p1=direct_packet(attach(rows,B),13)
p2=direct_packet(attach(rows,A|{8}),13);p3=direct_packet(attach(rows,B|{8}),13)
check(p0==p1,'full packet twin');check(p2[1:]==p3[1:],'full central released')
out['homometric13']={'A':sorted(A),'B':sorted(B),'core_h':h,'core_Q':Q,'twin_packet':p0,'released_packets':[p2,p3],'common_toggle_table':records}
# Exact common-toggle polynomial cancellation, tested on all128 choices.
J=[v for v in range(13)if (v in A)==(v in B)]
gaps={v:-(1-2*(v in A))*sum(W[v][j]*((j in A)-(j in B))for j in range(13))for v in J}
check(sum(gaps.values())==0,'particle-hole local-gap sum')
for mask in range(1<<len(J)):
 R={v for i,v in enumerate(J)if mask>>i&1}
 check(energy(A^R,h,r,W)-energy(B^R,h,r,W)==sum(gaps[v]for v in R),'replica cancellation')
rows2=cyc(13,[12,2,10,9,8,7]);hh,qq,ww=coeff(rows2);rr=(qq-2*hh)//13
AA={i for i in range(13)if 83>>i&1};BB={i for i in range(13)if 1128>>i&1}
check(corr(AA,13)==corr(BB,13),'second homometry')
R={2,11};check(not (R&AA) and not (R&BB),'common additions')
gg=[-sum(ww[v][j]*((j in AA)-(j in BB))for j in range(13))for v in sorted(R)]
check(sum(gg)==44,'gap44 cancellation')
pp=[direct_packet(attach(rows2,C),13)for C in [AA^R,BB^R]]
check(pp[0][0]-pp[1][0]==44 and pp[0][1:]==pp[1][1:],'literal full-central44')
out['replica_gap44']={'core_increments':[12,2,10,9,8,7],'original_cuts':[83,1128],'toggles':sorted(R),'new_cuts':[2135,3180],'local_gaps':gg,'packets':pp}
# Pair-transitive hostile uses exact Paley7 W; no assumption from the model.
paley=cyc(7,[1,2,4]);ph,pq,pw=coeff(paley);check(len(set(pw[0][1:]))==1,'Paley kernel constant')
for k in range(8):
 values={energy(C,ph,(pq-2*ph)//7,pw)for C in combinations(range(7),k)}
 check(len(values)==1,'Paley cardinality fiber rigid')
out['hostile']='Paley7: every cardinality layer is scalar-rigid'
# All-order unique-sum construction, tested at explicit odd modulus.
P={0,1,3};Qset={0,10,40};nn=101
AA={(p+q)%nn for p in P for q in Qset};BB={(p-q)%nn for p in P for q in Qset}
check(len(AA)==len(BB)==9,'Boolean unique-sum sidecar');check(corr(AA,nn)==corr(BB,nn),'group-ring homometry')
out['unique_sum_control']={'n':nn,'P':sorted(P),'Q':sorted(Qset),'A':sorted(AA),'B':sorted(BB)}
print(json.dumps(out,sort_keys=True,indent=2))
