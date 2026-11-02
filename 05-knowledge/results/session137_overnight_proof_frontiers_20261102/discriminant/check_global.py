"""Exact integer proof-identity checks, including weighted tournaments."""
from itertools import product, permutations, combinations
from math import factorial
from pathlib import Path
import hashlib,json,random
from probe import det

def need(ok,msg):
 if not ok:raise RuntimeError(msg)

def matrix(n,d,weights):
 a=[[0]*n for _ in range(n)]
 for z,(i,j) in zip(weights,combinations(range(n),2)):a[i][j]=z;a[j][i]=d-z
 return a

def sub(a,ids):return [[a[i][j]for j in ids]for i in ids]

def hp(a):
 n=len(a)
 if not n:return 1
 dp=[[0]*n for _ in range(1<<n)]
 for i in range(n):dp[1<<i][i]=1
 for s in range(1,1<<n):
  for i,c in enumerate(dp[s]):
   if c:
    for j,w in enumerate(a[i]):
     if w and not(s>>j&1):dp[s|1<<j][j]+=c*w
 return sum(dp[-1])

def data(a,d):
 n=len(a)
 if not n:return (1,1,1)
 h=d*hp(a);u=det([[a[i][j]+d*(i==j)for j in range(n)]for i in range(n)]);v=(-1)**n*det(a)
 need(u>abs(v),'algebraic positivity')
 need(det([[d*(i==j)+a[i][j]-a[j][i]for j in range(n)]for i in range(n)])==2**(n-1)*(u+v),'disc identity')
 return h,u,v

def transitive(a,d):
 return all(a[i][j]in(0,d)for i in range(len(a))for j in range(len(a))) and hp(a)==d**(len(a)-1)

def check(a,d):
 n=len(a);h,u,v=data(a,d)
 need(h>=u+abs(v),'GLOBAL inequality')
 binary=all(a[i][j]in(0,d)for i in range(n)for j in range(n))
 eq=n<=2 or (binary and (transitive(a,d) or (n==3 and hp(a)==3*d*d)))
 need((h==u+abs(v))==eq,'weighted equality classification')
 need((h==u+v)==transitive(a,d),'disc equality')
 return h,u,v

def insertion(a,d,x):
 ids=[i for i in range(len(a))if i!=x];old=sub(a,ids);n=len(old)
 dh,du,dv=data(a,d);ph,pu,pv=data(old,d)
 U=d*pu;V=0;H=d*ph;paths=empty=one=0;gaps=[d*(ph-pu),d*(ph-pu)]
 for ell in range(1,n+1):
  for P in permutations(range(n),ell):
   w=a[x][ids[P[0]]]*a[ids[P[-1]]][x]
   for i,j in zip(P,P[1:]):w*=old[i][j]
   if not w:continue
   rem=[i for i in range(n)if i not in P];hr,ur,vr=data(sub(old,rem),d)
   paths+=1;empty+=not rem;one+=ell==1
   U+=(1 if ell%2==0 else -1)*w*ur;V-=w*vr
   if ell%2==0:H+=2*w*hr
   for z,eps in enumerate((-1,1)):
    term=(2*(hr-ur)if ell%2==0 else 0)+ur-eps*vr
    need(term>=0,'signed summand positivity');gaps[z]+=w*term
 need((H,U,V)==(dh,du,dv),'vertex cycle identities')
 need(gaps==[dh-du-dv,dh-du+dv],'signed gaps')
 return paths,empty,one

def cycles(a):
 n=len(a);cs=[]
 for ell in range(3,n+1):
  for S in combinations(range(n),ell):
   for rest in permutations(S[1:]):
    P=(S[0],)+rest
    if all(a[P[i]][P[(i+1)%ell]] for i in range(ell)):cs.append(frozenset(P))
 return cs

def main():
 total=rec=paths=empty=one=defects=0
 for n in range(1,7):
  for bits in product(range(2),repeat=n*(n-1)//2):
   a=matrix(n,1,bits);h,u,v=check(a,1);total+=1
   if n<=4:
    for x in range(n):
     p,e,o=insertion(a,1,x);rec+=1;paths+=p;empty+=e;one+=o
   if n<=5:
    cs=cycles(a);proper=sum(len(c)<n for c in cs);disjoint=any(not c&b for i,c in enumerate(cs)for b in cs[i+1:])
    need(h-u>=len(cs),'all cycle bound')
    need(h-u-abs(v)>=proper,'proper cycle strengthening')
    need((h-u-abs(v)==proper)==(not disjoint),'proper cycle equality')
    if disjoint:need(h-u-abs(v)>=proper+2,'packing2 gap')
    defects+=1
 weighted=0
 for n in range(1,5):
  for w in product(range(4),repeat=n*(n-1)//2):
   a=matrix(n,3,w);check(a,3);weighted+=1
   if n<=3 or weighted%31==0:
    for x in range(n):
     p,e,o=insertion(a,3,x);rec+=1;paths+=p;empty+=e;one+=o
 rng=random.Random(137)
 for n in (5,6,7,8):
  for _ in range(40):
   a=matrix(n,5,[rng.randrange(6)for _ in range(n*(n-1)//2)]);check(a,5);weighted+=1
   if n<=6:
    p,e,o=insertion(a,5,0);rec+=1;paths+=p;empty+=e;one+=o
 need(data([],1)==(1,1,1) and 1<1+abs(1),'empty exception')
 c3=matrix(3,1,(1,0,1));h,u,v=data(c3,1)
 need((h,u,v)==(3,2,-1),'C3')
 # Deterministic controlled rejection of the named false implications.
 need(v!=1,'wrong v sign rejected')
 need(h!=2,'missing OCF factor2 rejected')
 w2=matrix(2,2,(1,));wh,wu,wv=data(w2,2)
 need((wh,wu,wv)==(4,3,-1) and wu!=4,'omitted weighted one-vertex path rejected')
 return {'status':'PASS','ordinary_tournaments':total,'weighted_tournaments':weighted,'root_recurrences':rec,'positive_weight_old_paths':paths,'empty_complements':empty,'weighted_one_vertex_paths':one,'ordinary_proper_cycle_checks':defects,'negative_controls':['empty theorem domain','wrong distinguished-cycle v sign','missing OCF factor2','omitted weighted one-vertex path'],'scope':'Finite exact identity/equality controls; the universal theorem is the written induction, not this census.'}
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
