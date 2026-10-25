import itertools,json,math,hashlib
from fractions import Fraction
from pathlib import Path
OUT=Path(__file__).parent

def require(ok,msg):
 if not ok: raise RuntimeError(msg)

def elimination(a):
 a=[[Fraction(v) for v in row] for row in a];r=0
 for c in range(len(a[0])):
  p=next((k for k in range(r,len(a)) if a[k][c]),None)
  if p is None:continue
  a[p],a[r]=a[r],a[p];z=a[r][c];a[r]=[v/z for v in a[r]]
  for k in range(r+1,len(a)):
   z=a[k][c]
   if z:a[k]=[v-z*w for v,w in zip(a[k],a[r])]
  r+=1
 return r

def census(n):
 E=list(itertools.combinations(range(n),2));ix={e:k for k,e in enumerate(E)};P=list(itertools.permutations(range(n)));N=1<<len(E);full=N-1
 maps=[]
 for p in P:
  targets=[1<<ix[tuple(sorted((p[a],p[b])))] for a,b in E]
  flips=sum(targets[k] for k,(a,b) in enumerate(E) if p[a]>p[b])
  maps.append((targets,flips))
 def image(x,t):
  z=0
  while x:
   bit=x&-x;z|=t[bit.bit_length()-1];x-=bit
  return z
 graphs={};rows=[];ng=0
 for x in range(N):
  if x in graphs:continue
  ng+=1;signs={};ambiguous=False;sm=0;aut=0
  for t,flip in maps:
   y=image(x,t);s=1-2*((y&flip).bit_count()%2);sm+=s
   if y in signs and signs[y]!=s:ambiguous=True
   signs[y]=s
   if y==x:aut+=1
  graphs.update(dict.fromkeys(signs,x))
  if not ambiguous and x.bit_count()%2==0:rows.append((x,sm,aut))
 ts={};hist=[];reps=[]
 for x in range(N):
  if x in ts:continue
  orbit={image(x,t)^flip for t,flip in maps};orbit|={y^full for y in orbit};ts.update(dict.fromkeys(orbit,x));reps.append(x)
  h={}
  for t,flip in maps:
   y=image(x,t)^flip;q=min(graphs[y],graphs[y^full]);h[q]=h.get(q,0)+1
  hist.append(h)
 cols=sorted(set().union(*hist));M=[[h.get(c,0) for c in cols] for h in hist];r=elimination(M)
 require(len(rows)==len(reps),'Walsh dimension versus direct class count')
 require(r==sum(s!=0 for x,s,a in rows),'diagonal versus direct rank')
 require(len({tuple(v) for v in M})==len(M),'profile class separation')
 return {'n':n,'graph_classes':ng,'converse_classes':len(reps),'rank':r,'distinct_profiles':len(M),'zero_modes':[{'mask':x,'sum':s,'automorphisms':a,'edges':[E[k] for k in range(len(E)) if x>>k&1]} for x,s,a in rows if s==0],'rows':rows,'matrix':M,'reps':reps,'cols':cols}

def exact_integral():
 # 2 d(1-d)(1-2d)^3, independently coefficient expanded.
 poly={1:2,2:-14,3:36,4:-40,5:16}
 v=sum(Fraction(c,k+1) for k,c in poly.items());require(v==0,'connected-six integral');return str(v)

def component_check():
 E=[(0,1),(0,2),(0,3),(4,5),(4,6),(4,7)];total=0
 for p in itertools.permutations(range(8)):
  total+=(-1)**sum(p[a]>p[b] for a,b in E)
 require(total==0,'two three-leaf stars coefficient')
 return total

def main():
 data=[census(n) for n in range(1,7)]
 require([(x['converse_classes'],x['rank']) for x in data]==[(1,1),(1,1),(2,2),(3,3),(10,10),(34,32)],'expected table')
 require([z['mask'] for z in data[-1]['zero_modes']]==[255,2015],'exact two zero modes')
 out={'census':data,'connected_integral':exact_integral(),'two_star_order_sum':component_check()}
 (OUT/'independent_results.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({'table':[{k:d[k] for k in ['n','graph_classes','converse_classes','rank','distinct_profiles','zero_modes']} for d in data],'integral':out['connected_integral'],'two_star_order_sum':out['two_star_order_sum']},indent=2))
if __name__=='__main__':main()
