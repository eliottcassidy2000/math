import itertools as it,json,math,random
from pathlib import Path
D=Path(__file__).resolve().parent;out={}
def need(b,s):
 if not b:raise RuntimeError(s)
def roots(p):
 z=[]
 for a,b in zip(p,p[1:]):
  r=math.isqrt(a+b);need(r*r==a+b,'edge');z.append(r)
 return z
# Independently validate every supplied vertex and edge, including both twin placements.
c=json.loads((D/'certificate.json').read_text());t=json.loads((D/'twins.json').read_text())
for key,p in c['Hamilton_paths'].items():need(sorted(p)==list(range(1,int(key)+1)),'vertex permutation');roots(p)
p,q=t['path'],t['second'];need(sorted(p)==sorted(q)==list(range(1,37)),'twins permutation')
need(roots(p)==roots(q),'same roots');need([y-x for x,y in zip(p,q)]==[6*(-1)**i for i in range(36)],'explicit twin shift')
out['certified_path_sizes']=sorted(map(int,c['Hamilton_paths']));out['twin_shift']=6
# All ordered words of arbitrary edge sums, not merely positive square sums.
count=0
for m in range(6):
 for sums in it.product((-2,0,3,7),repeat=m):
  constants=[0]
  for a in sums:constants.append(a-constants[-1])
  E=constants[::2];O=constants[1::2]
  for L,U in ((-5,5),(1,12),(3,7)):
   actual=[]
   for x in range(L,U+1):
    a=[x]
    for s in sums:a.append(s-a[-1])
    if len(a)==len(set(a)) and all(L<=v<=U for v in a):actual.append(x)
   expected=[]
   if len(set(E))==len(E) and len(set(O))==len(O):
    lo=max([L-e for e in E]+[o-U for o in O]);hi=min([U-e for e in E]+[o-L for o in O])
    expected=[x for x in range(lo,hi+1) if all(2*x!=o-e for e in E for o in O)]
   need(actual==expected,'general interval placement');count+=1
out['arbitrary_sum_interval_controls']=count
# Include odd n, absent from author finite controls; arbitrary even-position value sets.
count=0
for n in range(1,13):
 U=set(range(1,n+1))
 for S in it.combinations(U,(n+1)//2):
  A=set(S)
  for t0 in range(1,n+2):
   translated=[x+t0 for x in A]+[x-t0 for x in U-A]
   actual=sorted(translated)==list(range(1,n+1))
   predicted=n%(2*t0)==0 and A=={x for x in U if (x-1)//t0%2==0}
   need(actual==predicted,'shift classification');count+=1
out['shift_controls_including_odd_orders']=count
# Verify all anchor edges and component lower bound using union-find, not author's DFS.
for r in c['optimal_anchor_augmentations']:
 P=r['path'];A=set(r['anchors']);need(len(P)==len(set(P)) and A<=set(P),'anchors');roots(P)
 parent={a:a for a in A}
 def find(x):
  while parent[x]!=x:x=parent[x]
  return x
 for a,b in it.combinations(A,2):
  k=math.isqrt(a+b)
  if k*k==a+b:parent[find(a)]=find(b)
 components=len({find(a) for a in A});need(len(set(P)-A)==components-1,'optimal cardinality')
out['optimal_connector_sizes']=[2,3]
# Regression for the missing distinctness condition, communicated to author.
sol=[z for z in range(1,1000) if math.isqrt(z+2)**2==z+2 and math.isqrt(z+7)**2==z+7]
need(sol==[2],'degenerate connector regression');out['degenerate_connector_a2_b7']={'positive_solutions':sol,'simple_connector':False}
# Independent exhaustive z scan, not factor enumeration. Every solution has
# 2*sqrt(z+min(a,b)) < abs(a-b), so z < abs(a-b)^2 suffices.
height=[]
occupied={4,2,5,3,7,13,44,527,8}
for a in (7,527):
 for b in (3,4,13,44,8):
  candidates=[z for z in range(1,(a-b)**2+1) if z not in occupied and math.isqrt(z+a)**2==z+a and math.isqrt(z+b)**2==z+b]
  height.append({'a':a,'b':b,'connectors':candidates})
need(height==c['minimum_count_height_rows'],'all finite height rows')
need(min(z for row in height for z in row['connectors'])==317,'height minimum')
out['minimum_count_height']=317
out['height_rows']=height
print(json.dumps(out,sort_keys=True,indent=2))
