"""Independent literal-graph controls for clone matching and blow-up sparsity."""
from graphs import strong,clone,counts,edges
from itertools import combinations,product
from math import factorial
from pathlib import Path
import json

def need(b,s):
 if not b:raise RuntimeError(s)
def matching(E):
 best=0
 def rec(k,O,I):
  nonlocal best
  best=max(best,len(O))
  for j in range(k,len(E)):
   o,i=E[j]
   if o not in O and i not in I:rec(j+1,O|{o},I|{i})
 rec(0,set(),set());return best
def stirling(m,j):
 D=[[0]*(m+1) for _ in range(m+1)];D[0][0]=1
 for a in range(1,m+1):
  for b in range(1,a+1):D[a][b]=D[a-1][b-1]+b*D[a-1][b]
 return D[m][j] if j<=m else 0
def extract(seq):
 out=[]
 for m,x in enumerate(seq,1):
  residual=x-sum(a*factorial(j)*stirling(m,j) for j,a in enumerate(out,1))
  need(residual%factorial(m)==0,'nonintegral run coefficient');out.append(residual//factorial(m))
 return out
def control(r,a):
 nu=matching(edges(r,a));need(nu>=1,'strong root has no triangle')
 values=[counts(clone(r,a,m)) for m in range(1,nu+3)]
 H=extract([h for h,c in values]);C=extract([c for h,c in values])
 need(H[nu]>0 and all(x==0 for x in H[nu+1:]),'Hamilton path matching degree')
 need(C[nu-1]>0 and all(x==0 for x in C[nu:]),'Hamilton cycle matching degree')
 need(all(x>=0 for x in H+C),'negative run coefficient')
 for m,(h,c) in enumerate(values,1):
  need(c>=factorial(nu-1)*factorial(nu)*stirling(m,nu),'cycle packing lower bound')
 if nu==1:
  need(H[1]>=values[0][1],'binary amplitude cycle lower bound')
  need(all(c==values[0][1] for h,c in values),'star cycle preservation')
 return dict(nu=nu,H_coefficients=H,C_coefficients=C,counts=values)
def blowup(S,m):
 origin=[i for i,k in enumerate(m) for _ in range(k)];out=[0]*len(origin)
 for a in range(len(origin)):
  for b in range(a+1,len(origin)):
   f=origin[a]==origin[b] or S[origin[a]]>>origin[b]&1
   out[a if f else b]|=1<<(b if f else a)
 return out,origin

def main():
 tested=0;hist={};hostiles={}
 for n in range(3,6):
  P=list(combinations(range(n),2))
  for bits in range(1<<len(P)):
   r=[0]*n
   for k,(u,v) in enumerate(P):r[u if bits>>k&1 else v]|=1<<(v if bits>>k&1 else u)
   if not strong(r):continue
   for a in range(n):
    z=control(r,a);tested+=1;hist[z['nu']]=hist.get(z['nu'],0)+1
    h,c=z['counts'][0]
    if c<=2 and z['nu']==1 and h in (3,5,9,11,13):
     B=2*h-z['counts'][1][0];hostiles.setdefault((h,B),dict(rows=r,root=a,**z))
 R7=[sum(1<<((i+j)%7) for j in (1,2,3)) for i in range(7)];rank3=control(R7,0);need(rank3['nu']==3,'rank3 control')
 clock=[]
 for m in range(1,11):
  h,c=counts(clone([2,4,1],0,m));need(h==2**m+1 and c==1,'induced clock');clock.append([m,h,c])
 # Repeated replacements can occur at different current vertices and inside old clone modules.
 sequences=[[(0,2),(0,3),(3,2)],[(2,3),(0,2),(4,2)],[(1,2),(3,2),(0,2),(4,2)]]
 checked=0
 for seed in ([2,4,1],[6,12,8,1]):
  for steps in sequences:
   r=list(seed);origin=list(range(len(seed)))
   for a,m in steps:
    a%=len(r);r=clone(r,a,m);origin=origin[:a]+origin[a+1:]+[origin[a]]*m
    # All original-seed modules remain transitive; topologically sort each one.
    order=[];mult=[]
    for v in range(len(seed)):
     block=[i for i,x in enumerate(origin) if x==v];mult.append(len(block))
     snapshot=tuple(block);block.sort(key=lambda i:-sum(bool(r[i]>>j&1) for j in snapshot))
     need(all(r[block[i]]>>block[j]&1 for i in range(len(block)) for j in range(i+1,len(block))),'module not transitive')
     order+=block
    expected,_=blowup(seed,mult);actual=[sum(1<<j for j,v in enumerate(order)if r[u]>>v&1)for u in order]
    need(actual==expected,'changed-root clone escaped original blowup');checked+=1
    h,_=counts(r);need(all(h>=2**m+1 for m in mult),'multiplicity lower bound')
 out=dict(status='PASS',complete_strong_rooted_controls=tested,matching_histogram=hist,rank3_rows=R7,rank3=rank3,clock_controls=clock,changed_root_clone_controls=checked,small_binary_offsets=[dict(H=h,B=B,**v)for (h,B),v in sorted(hostiles.items())])
 data=json.dumps(out,sort_keys=True,indent=2)+'\n';(Path(__file__).parent/'VALIDATION.json').write_text(data);print(data,end='')
if __name__=='__main__':main()
