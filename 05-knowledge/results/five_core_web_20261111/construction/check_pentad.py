"""Independent small module orbits and literal pentagon-cut polynomial control."""
import itertools,json
from pathlib import Path
HERE=Path(__file__).resolve().parent

def need(q,label):
 if not q:raise RuntimeError(label)
def apply(A,x):
 y=0
 for i,a in enumerate(A):
  if x>>i&1:y^=a
 return y
def mul(A,B):return tuple(apply(A,b) for b in B)
def rank(A):
 basis={}
 for x in A:
  while x:
   p=x.bit_length()-1
   if p in basis:x^=basis[p]
   else:basis[p]=x;break
 return len(basis)
def mats(u,v):return list(itertools.product(range(1<<u),repeat=v))
def identity(n):return tuple(1<<i for i in range(n))
def inv(A):
 n=len(A);return tuple(next(x for x in range(1<<n) if apply(A,x)==1<<i) for i in range(n))
def module_orbits():
 report=[];indecomposable=[];total=0
 for u in range(3):
  for v in range(3):
   glU=[(P,inv(P)) for P in mats(u,u) if rank(P)==u]
   glV=[(Q,inv(Q)) for Q in mats(v,v) if rank(Q)==v]
   seen=set();rank_orbits={}
   for a in mats(u,v):
    for b in mats(v,u):
     if any(mul(b,a)):continue
     total+=1
     orbit=min((mul(Pi,mul(a,Q)),mul(Qi,mul(b,P))) for P,Pi in glU for Q,Qi in glV)
     if orbit in seen:continue
     seen.add(orbit);r,s,t=rank(a),rank(b),rank(mul(a,b));key=(r,s,t)
     need(key not in rank_orbits,'ranks distinguish all enumerated orbits');rank_orbits[key]=orbit
     mult=(u-r-s,v-r-s+t,r-t,s-t,t)
     need(min(mult)>=0,'nonnegative multiplicities')
     split=False
     for E in mats(u,u):
      if mul(E,E)!=E:continue
      for F in mats(v,v):
       if mul(F,F)!=F or mul(E,a)!=mul(a,F) or mul(F,b)!=mul(b,E):continue
       if (E,F) not in ((tuple(0 for _ in range(u)),tuple(0 for _ in range(v))),(identity(u),identity(v))):split=True
     actual_indecomp=(u+v>0 and not split)
     need(actual_indecomp==(sum(mult)==1),'idempotent test versus decomposition')
     if actual_indecomp:indecomposable.append({'u':u,'v':v,'r':r,'s':s,'t':t})
   report.append({'u':u,'v':v,'orbits':len(seen)})
 need(len(indecomposable)==5,'five small indecomposables')
 return {'universe':'All F2 maps a:V->U,b:U->V with ba=0, dimU,dimV<=2; change-of-basis orbits and commuting idempotents.','representations':total,'orbits':report,'indecomposables':indecomposable}

def hp(rows,vertices=None):
 if vertices is None:vertices=list(range(len(rows)))
 return sum(all(rows[p[i]]>>p[i+1]&1 for i in range(len(p)-1)) for p in itertools.permutations(vertices))
def polynomial(rows,cut):
 n=len(rows);K=[0]*(n+1);O=K.copy();I=K.copy();J=K.copy()
 for l in range(1,n+1):
  for p in itertools.permutations(range(n),l):
   if not all(rows[p[i]]>>p[i+1]&1 for i in range(l-1)):continue
   w=hp(rows,[v for v in range(n) if v not in p]);K[l]+=w
   a=bool(cut>>p[0]&1);b=not(cut>>p[-1]&1)
   if a:O[l]+=w
   if b:I[l]+=w
   if a and b:J[l]+=w
 return K,O,I,J

def g(n):
 if n<2:return 1
 q=n//2;return g(q)+(q+(q&1) if not(n&1) else 0)
def child(core,cut,a,b):
 left=list(map(int,bin(a)[2:])) if a else []
 right=[1-int(c) for c in bin(b)[2:][::-1]] if b else []
 l=len(left);n=len(core);N=l+n+len(right);r=[0]*(N+1)
 for i in range(N):
  for j in range(i+1,N):
   ok=(core[i-l]>>(j-l)&1) if l<=i<j<l+n else True
   u,v=(i,j) if ok else (j,i);r[u]|=1<<v
 for i,bit in enumerate(left+[(cut>>v)&1 for v in range(n)]+right):r[N if bit else i]|=1<<(i if bit else N)
 return r
def dp(r):
 n=len(r);full=(1<<n)-1;table=[[0]*n for _ in range(full+1)]
 for i in range(n):table[1<<i][i]=1
 for S in range(1,full+1):
  for last,c in enumerate(table[S]):
   todo=r[last]&~S&full
   while c and todo:
    b=todo&-todo;todo-=b;table[S|b][b.bit_length()-1]+=c
 return sum(table[-1])
def main():
 result={'modules':module_orbits()};r=(6,12,24,17,3)
 p=polynomial(r,3);q=polynomial(r,5)
 need(p[:3]==q[:3],'proper observers unchanged')
 diff=[b-a for a,b in zip(p[-1],q[-1])]
 need(diff==[0,0,-2,-2,3,1],'factored central polynomial')
 tests=0
 for a,b in itertools.product(range(5),repeat=2):
  base=41+15*(g(a)+g(b)-2)+120*a*b+63*a+42*b+9*(a&1)+6*(b&1)
  for cut,offset in ((3,0),(5,2)):
   need(dp(child(r,cut,a,b))==base+offset,'literal two-flank child');tests+=1
 result['central_switch']={'rows':r,'edge_cut':3,'diagonal_cut':5,'edge_polynomials':p,'diagonal_polynomials':q,'difference':diff,'literal_child_checks':tests}
 (HERE/'pentad_checks.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
if __name__=='__main__':main()
