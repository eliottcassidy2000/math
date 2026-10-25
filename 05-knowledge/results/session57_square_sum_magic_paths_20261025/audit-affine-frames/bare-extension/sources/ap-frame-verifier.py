"""Exact AP-only frame counts. Standard library; all guards survive python -O."""
import itertools as it,json,pathlib,collections
P=pathlib.Path(__file__).parent

def require(x,msg):
 if not x: raise RuntimeError(msg)

class Field:
 def __init__(self,n,f):
  self.n=n;self.f=f;self.q=3**n
  self.d=[tuple((i//3**k)%3 for k in range(n)) for i in range(self.q)]
 def enc(self,t):return sum(a%3*3**k for k,a in enumerate(t))
 def add(self,a,b):return self.enc(x+y for x,y in zip(self.d[a],self.d[b]))
 def neg(self,a):return self.enc(-x for x in self.d[a])
 def sub(self,a,b):return self.add(a,self.neg(b))
 def mul(self,a,b):
  z=[0]*(2*self.n-1)
  for i,x in enumerate(self.d[a]):
   for j,y in enumerate(self.d[b]):z[i+j]+=x*y
  for k in range(len(z)-1,self.n-1,-1):
   for j in range(self.n):z[k-self.n+j]-=z[k]*self.f[j]
  return self.enc(z[:self.n])

def key(B,D):return tuple(sorted((tuple(sorted(B)),tuple(sorted(D)))))

def check(q):
 if q==81:
  r=next(r for r in json.loads((P/'sources/results.json').read_text())['fields'] if r['q']==81)
  w=r['spectral_hp_groups'][0]
 else:r=json.loads((P/f'sources/witness-{q}.json').read_text());w=r
 F=Field(r['n'],r['modulus']);add=F.add;sub=F.sub
 A=set(w['vertices']);C=w['center'];U={sub(a,C) for a in A};Q={F.mul(i,i) for i in range(q)}
 require(len(A)==9 and A<=Q,'square plane witness')
 require(all(add(u,v) in U for u in U for v in U),'linear plane')
 dirs=sorted({tuple(sorted((0,u,F.neg(u)))) for u in U-{0}})
 require(len(dirs)==4,'directions')
 frames=set();bad_by_direction=[];all_frames=set()
 for L,M in it.combinations(dirs,2):
  bad=[];local=set()
  for g in range(q):
   B={add(g,l) for l in L};D={add(sub(C,g),m) for m in M}
   require({add(b,d) for b in B for d in D}==A,'cross sums')
   k=key(B,D);require(k not in all_frames,'parameter duplicate');all_frames.add(k)
   collision=bool(B&D)
   require(collision==(add(g,g) in A),'bad gauge predicate')
   require(len(B&D)==int(collision),'intersection singleton')
   if collision:bad.append(g)
   else:frames.add(k);local.add(k)
  require(len(bad)==9 and len(local)==q-9,'per-direction count');bad_by_direction.append(bad)
 require(len(frames)==6*(q-9),'frame count')
 # Changing representative C by h=l_h+m_h reindexes g by l_h.
 L,M=dirs[:2]
 for h in U:
  lh,mh=next((l,m) for l in L for m in M if add(l,m)==h)
  for g in range(q):
   gp=add(g,lh);Cp=add(C,h)
   require({add(g,l) for l in L}=={add(gp,l) for l in L},'C representative B')
   require({add(sub(C,g),m) for m in M}=={add(sub(Cp,gp),m) for m in M},'C representative D')
 # All 432 parametrized magic arrays, directly checking eight line sums.
 lines=((0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6))
 offsets=((1,0),(-1,-1),(0,1),(-1,1),(0,0),(1,-1),(0,-1),(1,1),(-1,0))
 layouts=set()
 scale=lambda s,x:0 if s==0 else x if s==1 else F.neg(x)
 for c in A:
  for u in U-{0}:
   for v in U-{0,u,F.neg(u)}:
    m=tuple(add(c,add(scale(s,u),scale(t,v))) for s,t in offsets)
    require(set(m)==A and all(add(add(m[i],m[j]),m[k])==0 for i,j,k in lines),'magic layout')
    layouts.add(m)
 require(len(layouts)==432,'432 layouts')
 # Explicit bare-K33 exception: one horizontal line and nonlinear transversal.
 u=next(x for x in dirs[0] if x);v=next(x for x in dirs[1] if x)
 B0=set(dirs[0]);D0={C,add(C,v),add(C,add(F.neg(v),u))}
 require(add(add(*tuple(D0)[:2]),tuple(D0)[2])!=0,'noncollinear transversal')
 g=next(g for g in range(q) if not ({add(g,b) for b in B0}&{sub(d,g) for d in D0}))
 B={add(g,b) for b in B0};D={sub(d,g) for d in D0}
 require(len({add(b,d) for b in B for d in D})==9 and {add(b,d) for b in B for d in D}==A,'bare witness')
 require(key(B,D) not in frames,'bare witness excluded')
 bare_witness={'B':sorted(B),'D':sorted(D)}
 brute=None
 if q==81:
  # Independent all ambient affine lines, without assuming directions belong to U.
  table=[[add(a,b) for b in range(q)] for a in range(q)]
  ambient=sorted({tuple(sorted((a,b,F.neg(table[a][b])))) for a,b in it.combinations(range(q),2)})
  require(len(ambient)==1080,'ambient line census')
  found=set()
  for B,D in it.combinations(ambient,2):
   if set(B)&set(D):continue
   if {table[b][d] for b in B for d in D}==A:found.add(key(B,D))
  require(found==frames,'independent ambient brute census')
  fibers=collections.Counter()
  for left,right in frames:
   for B,D in ((left,right),(right,left)):
    for b in it.permutations(B):
     for d in it.permutations(D):
      e=tuple(table[b[i]][d[j]] for i,j in ((2,2),(1,1),(0,0),(0,1),(2,0),(1,2),(1,0),(0,2),(2,1)))
      require(e in layouts,'ordered frame maps to magic layout');fibers[e]+=1
  require(len(fibers)==432 and set(fibers.values())=={q-9},'layout gauge fibers')
  brute={'ambient_lines':len(ambient),'unordered_line_pairs_examined':len(ambient)*(len(ambient)-1)//2,'matching_disjoint_frames':len(found),'fully_ordered_frames_checked':sum(fibers.values()),'gauge_fibers_per_layout':q-9}
 return {'q':q,'modulus_low_first':r['modulus'],'square_plane':sorted(A),'C':C,'directions':dirs,'bad_gauges_by_direction_pair':bad_by_direction,'unordered_disjoint_AP_frames':len(frames),'ordered_and_internally_ordered_AP_frames':72*len(frames),'magic_layouts':len(layouts),'cross_skeleton_HP_occurrences':36*len(frames),'bare_non_AP_witness':bare_witness,'independent_brute':brute}

rows=[]
for q in (81,729,2187):
 rows.append(check(q));print('PASS',q,rows[-1]['unordered_disjoint_AP_frames'],flush=True)
# K33 Hamilton paths, orientations discarded, with fixed six vertex labels.
paths={min(p,p[::-1]) for p in it.permutations(range(6)) if all((a<3)!=(b<3) for a,b in zip(p,p[1:]))}
require(len(paths)==36,'Hamilton paths per skeleton')
result={'status':'PASS','rows':rows,'unoriented_HP_per_cross_skeleton':len(paths),'scope':'AP-only frames with fixed nine-point cross-sum plane; paths are occurrence counts, not deduplicated ambient paths'}
(P/'results.json').write_text(json.dumps(result,indent=2)+'\n')
