"""Independent exact audit of root's cycle-expansion determinant inequality.
Uses this lane's independently written tournament canonicalizer and HP DP only.
"""
from verify import adj,canon,hp as nonempty_hp,cycles,check

def hp(a):
 return nonempty_hp(a) if a else 1
from itertools import permutations
import json

def det(M):
 a=[list(r) for r in M];n=len(a);sign=1;den=1
 if not n:return 1
 for k in range(n-1):
  if not a[k][k]:
   r=next((r for r in range(k+1,n) if a[r][k]),None)
   if r is None:return 0
   a[k],a[r]=a[r],a[k];sign=-sign
  pivot=a[k][k]
  for i in range(k+1,n):
   for j in range(k+1,n):
    num=a[i][j]*pivot-a[i][k]*a[k][j]
    check(num%den==0,'Bareiss integral division')
    a[i][j]=num//den
   a[i][k]=0
  den=pivot
 return sign*a[-1][-1]
def det_perm(M):
 n=len(M);total=0
 for p in permutations(range(n)):
  x=(-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))
  for i in range(n):x*=M[i][p[i]]
  total+=x
 return total

def all_cycles(a):
 # Count all literal cycles modulo cyclic rotation; no odd-length filter.
 out=[]
 def rec(p,mask):
  first=p[0];last=p[-1]
  if len(p)>=3 and a[last]>>first&1:out.append((tuple(p),mask))
  opts=a[last]&~mask&~((1<<(first+1))-1)
  while opts:
   bit=opts&-opts;opts-=bit;rec(p+[bit.bit_length()-1],mask|bit)
 for i in range(len(a)):rec([i],1<<i)
 return out

def uv(a):
 n=len(a);A=[[a[i]>>j&1 for j in range(n)] for i in range(n)]
 B=[[A[i][j]+int(i==j) for j in range(n)] for i in range(n)]
 return det(B),(-1)**n*det(A)
def induce(a,mask):
 v=[i for i in range(len(a)) if mask>>i&1]
 return tuple(sum(((a[i]>>j)&1)<<k for k,j in enumerate(v)) for i in v)
def extend(a,S):
 n=len(a);b=list(a)+[S]
 for i in range(n):
  if not S>>i&1:b[i]|=1<<n
 return tuple(b)
def paths(a):
 out=[]
 def rec(p,s):
  if len(p)>=2:out.append((tuple(p),s))
  opts=a[p[-1]]&~s
  while opts:
   bit=opts&-opts;opts-=bit;rec(p+[bit.bit_length()-1],s|bit)
 for i in range(len(a)):rec([i],1<<i)
 return out

def main():
 reps={0};rows=[];checked=0;equality=[];insertions=0
 # Empty sidecar is deliberately NOT tested against the full inequality.
 check(uv(())==(1,1),'empty determinant sidecar')
 for n in range(1,9):
  nxt=set()
  for code in reps:
   for S in range(1<<(n-1)):nxt.add(canon(extend(adj(n-1,code),S)))
  reps=nxt;eq=0;cyclic=0;mindiff=None;packing_two=0;refined_eq=0;refined_min=None
  for code in sorted(reps):
   a=adj(n,code);u,v=uv(a);H=hp(a);cs=all_cycles(a);cyc=bool(cs)
   odd=sum(len(p)%2 for p,mask in cs);even=len(cs)-odd
   triangles=[mask for p,mask in cs if len(p)==3]
   has_disjoint=any(not(x&y) for i,x in enumerate(triangles) for y in triangles[i+1:])
   ham=sum(len(p)==n for p,mask in cs);proper=len(cs)-ham
   check(H-u>=len(cs),'weak gap counts all directed cycles')
   excess=H-u-abs(v)-proper
   check(excess>=0,'proper-cycle refined determinant inequality')
   check((excess==0)==(not has_disjoint),'refined equality iff triangle packing at most one')
   if has_disjoint:
    check(excess>=2,'disjoint triangles force extra two')
    packing_two+=1
    refined_min=min(refined_min,excess) if refined_min is not None else excess
   else:
    check((u,v,H)==(1+odd-even,-ham,1+2*odd),'packing one exact expansions')
    refined_eq+=1
   check(u>abs(v),'strict determinant domination')
   check(H>=u+abs(v),'global nonempty H inequality')
   check((H==u+abs(v))==((not cyc)or(n==3)),'equality classification')
   check((H==u)==(not cyc),'weak inequality equality classification')
   eq+=H==u+abs(v);cyclic+=cyc;checked+=1
   if cyc:mindiff=min(mindiff,H-u-abs(v)) if mindiff is not None else H-u-abs(v)
   if n<=5:
    A=[[a[i]>>j&1 for j in range(n)]for i in range(n)]
    B=[[A[i][j]+int(i==j)for j in range(n)]for i in range(n)]
    check((det_perm(B),(-1)**n*det_perm(A))==(u,v),'independent determinant permutation check')
    pp=paths(a);full=(1<<n)-1;comps={mask:(hp(induce(a,full^mask)),*uv(induce(a,full^mask)))for _,mask in pp}
    for S in range(1<<n):
     selected=[(p,comps[mask])for p,mask in pp if S>>p[0]&1 and not S>>p[-1]&1]
     E=[r for p,r in selected if len(p)%2==0];O=[r for p,r in selected if len(p)%2]
     b=extend(a,S);uu,vv=uv(b);hh=hp(b)
     check(uu==u+sum(r[1]for r in E)-sum(r[1]for r in O),'u child cycle expansion')
     check(vv==-sum(r[2]for p,r in selected),'v child cover expansion')
     check(hh==H+2*sum(r[0]for r in E),'H child OCF expansion')
     for sign in(-1,1):
      rhs=H-u+sum(2*h-q-sign*w for h,q,w in E)+sum(q-sign*w for h,q,w in O)
      check(hh-uu+sign*vv==rhs,'signed gap recurrence')
      check(all(q-sign*w>=0 for h,q,w in E+O),'nonnegative determinant kernels incl empty')
     insertions+=1
  rows.append({'n':n,'classes':len(reps),'cyclic_classes':cyclic,'equality_classes':eq,'minimum_cyclic_gap':mindiff,'packing_at_least_two_classes':packing_two,'refined_equality_classes':refined_eq,'minimum_packing_two_refined_excess':refined_min})
 print(json.dumps({'nonempty_classes_checked':checked,'literal_parent_cut_checks':insertions,'rows':rows,'empty_boundary':{'H':1,'u':1,'v':1,'full_inequality':False}},indent=2,sort_keys=True))
if __name__=='__main__':main()
