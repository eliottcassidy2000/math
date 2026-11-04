from itertools import permutations
import json,pathlib
P=pathlib.Path(__file__).parent
r=[48,5,9,3,14,30]
def count(r,S):return sum(all(r[a]>>b&1 for a,b in zip(t,t[1:])) for t in permutations(S))
def strong(r):
 for a in range(len(r)):
  s={a}
  while True:
   z=s|{j for i in s for j in range(len(r)) if r[i]>>j&1}
   if z==s:break
   s=z
  if len(s)!=len(r):return False
 return True
cases=[]
for edge in [None]+[(b,a) for b in [4,5] for a in [1,2,3]]:
 rr=r.copy()
 if edge:
  b,a=edge;rr[a]^=1<<b;rr[b]^=1<<a
 h=count(rr,range(6));q=count(rr,range(1,6));cases.append({'arc_reversed':edge,'rows':rr,'strong':strong(rr),'H':h,'Hminusroot':q,'delta':h-q})
 if not strong(rr) or (edge is None and (h,q)!=(33,3)) or (edge and (h,q)!=((37,5) if edge[0]==4 else(43,9))):raise RuntimeError('literal mismatch')
print(json.dumps(cases,indent=2));(P/'smoothing_literal_receipt.json').write_text(json.dumps(cases,indent=2)+'\n')
