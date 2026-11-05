from core import inventories,require
from itertools import permutations
from pathlib import Path
import json
P=Path(__file__).parent

def induced(T,vs):return [sum(((T[a]>>b)&1)<<j for j,b in enumerate(vs))for a in vs]
def reverse(T):return [sum(((T[b]>>a)&1)<<b for b in range(len(T)))for a in range(len(T))]
def clone(T,r):
 n=len(T);U=T[:]+[T[r]];U[r]|=1<<n
 for a in range(n):
  if a!=r and T[a]>>r&1:U[a]|=1<<n
 return U

def components(T,vs):
 reach={a:{a}|{b for b in vs if T[a]>>b&1}for a in vs}
 for z in vs:
  for a in vs:
   if z in reach[a]:reach[a]|=reach[z]
 groups=[];left=set(vs)
 while left:
  a=min(left);C={b for b in left if b in reach[a]and a in reach[b]};groups.append(C);left-=C
 return groups

def main():
 bank=json.loads((P/'controls.json').read_text());roots=zeros=units=literal=0;non=[]
 for n in range(2,8):
  for z in bank[str(n)]:
   T0=z['rows'];p=inventories(T0);F=(1<<n)-1
   for r in range(n):
    O=[a for a in range(n)if T0[r]>>a&1];I=[a for a in range(n)if T0[a]>>r&1];E=[(a,b)for a in O for b in I if T0[a]>>b&1]
    if not E:continue
    centers=[c for c in range(n)if all(c in edge for edge in E)]
    if not centers:continue
    c=centers[0];T=T0 if c in O else reverse(T0)
    O=[a for a in range(n)if T[r]>>a&1];I=[a for a in range(n)if T[a]>>r&1]
    Oi=inventories(induced(T,O));Ii=inventories(induced(T,I));X=sum(Ii['e'][j]for j,b in enumerate(I)if T[c]>>b&1);Y=Oi['s'][O.index(c)]
    h=p['H'];a2=(inventories(clone(T0,r))['H']-h)//2;d=p['allH'][F^(1<<r)];beta=h-2*a2
    require(d-h+a2==-X*Y,'rank-one boundary product');roots+=1
    # Small literal outside-word/cut skeleton checks do not use the clone recurrence.
    if n<=5:
     outside=[a for a in range(n)if a!=r];val=0
     for w in permutations(outside):
      z0=1-int(bool(T0[r]>>w[0]&1));z1=1-int(bool(T0[w[-1]]>>r&1));term=z0*z1
      for a,b in zip(w,w[1:]):term*=int(bool(T0[a]>>b&1))-int(bool(T0[a]>>r&1 and T0[r]>>b&1))
      val+=term
     require(val==d-h+a2,'literal gap evaluation');literal+=1
    twin=[a for a in range(n)if a!=r and all(((T0[a]>>b)&1)==((T0[r]>>b)&1)for b in range(n)if b not in [a,r])]
    if d-h+a2:continue
    zeros+=1
    if Y==0:
     C=next(C for C in components(T,O)if all(T[a]>>b&1 for a in C for b in set(O)-C));require(c not in C,'source block avoids crossing center')
    else:
     C=next(C for C in components(T,I)if all(T[b]>>a&1 for a in C for b in set(I)-C));require(all(T[a]>>c&1 for a in C),'terminal block beats center')
    require(all(((T[a]>>b)&1)==((T[r]>>b)&1)for a in C for b in range(n)if b not in C and b!=r),'adjacent homogeneous block')
    G=inventories(induced(T,sorted(C)))['H'];require(beta%G==0,'module offset divisibility')
    if abs(beta)==1:
     require(len(C)==1 and twin,'unit boundary implies actual twin');units+=1
    elif not twin:non.append(dict(rows=T0,root=r,H=h,a2=a2,beta=beta,deletion=d,block=sorted(C),block_H=G))
 # The smallest nonunit counterexample forbids dropping the unit-offset hypothesis.
 H=[2,32,11,19,7,29];r=0;p=inventories(H);a2=(inventories(clone(H,r))['H']-p['H'])//2
 require((p['allH'][62],p['H'],a2)==(15,27,12),'nonunit hostile')
 require(not any(all(((H[a]>>b)&1)==((H[r]>>b)&1)for b in range(6)if b not in [a,r])for a in range(1,6)),'hostile has no twin')
 print(json.dumps(dict(status='PASS',physical_rankone_roots=roots,literal_small_gap_evaluations=literal,boundary_zero_roots=zeros,unit_boundary_zero_roots=units,nonunit_nontwin_cases=non),sort_keys=True,indent=2))
if __name__=='__main__':main()
