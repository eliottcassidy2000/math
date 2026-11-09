#!/usr/bin/env python3
"""Exact finite controls; all-order tower formulas are proved separately."""
import itertools,json
from fractions import Fraction

def require(test,msg):
 if not test: raise RuntimeError(msg)
def matrix(n,mask):
 A=[[0]*n for _ in range(n)]
 for b,(i,j) in enumerate(itertools.combinations(range(n),2)):
  A[i][j]=(mask>>b)&1;A[j][i]=1-A[i][j]
 return A
def circ(n,S):return [[int((j-i)%n in S) for j in range(n)] for i in range(n)]
def join(A,B):
 a,b=len(A),len(B)
 return [A[i]+[1]*b for i in range(a)]+[[0]*a+B[i] for i in range(b)]
def ham(A):
 n=len(A);d=[[0]*n for _ in range(1<<n)]
 for i in range(n):d[1<<i][i]=1
 for m in range(1,1<<n):
  for v in range(n):
   z=d[m][v]
   if z:
    for w in range(n):
     if not m>>w&1 and A[v][w]:d[m|1<<w][w]+=z
 return sum(d[-1])
def autos(A):
 n=len(A);d=[sum(x) for x in A]
 return [p for p in itertools.permutations(range(n)) if all(d[v]==d[p[v]] for v in range(n)) and all(A[i][j]==A[p[i]][p[j]] for i in range(n) for j in range(i))]
def cycles(A):
 n=len(A);out=[]
 # Enumerate simple paths with their minimum vertex first, close odd cycles.
 for start in range(n):
  def dfs(path,used):
   if len(path)>=3 and len(path)%2 and A[path[-1]][start]:out.append(tuple(path))
   for w in range(start+1,n):
    if not used>>w&1 and A[path[-1]][w]:dfs(path+[w],used|1<<w)
  dfs([start],1<<start)
 return sorted(out)
def cancycle(c):
 j=c.index(min(c));return c[j:]+c[:j]
def packings(A):
 cs=cycles(A);out={()}
 def visit(i,mask,p):
  for j in range(i,len(cs)):
   cm=sum(1<<v for v in cs[j])
   if not cm&mask:
    for col in [0,1]:
     z=tuple(sorted(p+((cs[j],col),)));out.add(z);visit(j+1,mask|cm,z)
 visit(0,0,())
 return out
def stats(A,G=None):
 if G is None:G=autos(A)
 P=packings(A);seen=set();hist={}
 for p in sorted(P):
  if p in seen:continue
  orb={tuple(sorted((cancycle(tuple(g[v] for v in c)),b) for c,b in p)) for g in G}
  require(orb<=P,'action not closed');seen|=orb
  require(len(G)%len(orb)==0,'nonintegral stabilizer')
  s=len(G)//len(orb);hist[s]=hist.get(s,0)+1
 require(seen==P,'orbit census incomplete')
 h=ham(A);require(h==len(P),'OCF fails')
 mass=sum((Fraction(count,s) for s,count in hist.items()),Fraction())
 require(mass==Fraction(h,len(G)),'mass identity')
 require(mass.denominator==1,'path action divisibility')
 require(sum(hist.values())>=mass,'negative defect')
 require((sum(hist.values())==mass)==(len(G)==1),'strictness boundary')
 return {'n':len(A),'H':h,'aut_or_subgroup':len(G),'path_orbits':int(mass),'coarse_packing_orbits':sum(hist.values()),'stabilizer_histogram':dict(sorted(hist.items()))}
def product_group(gs):
 out=[()];n=0
 for G in gs:
  out=[p+tuple(n+x for x in q) for p in out for q in G];n+=len(G[0])
 return out

def main():
 tested=0
 for n in range(1,6):
  for m in range(1<<(n*(n-1)//2)):
   stats(matrix(n,m));tested+=1
 C=circ(3,{1});R=circ(5,{1,2});T=R;Gs=[autos(R)];tower=[]
 for k in range(3):
  G=product_group(Gs);z=stats(T,G)
  require(z['H']==15*3**k and z['path_orbits']==3 and z['coarse_packing_orbits']==7*3**k,'tower formula')
  if k<=1:require(set(G)==set(autos(T)),'full group control')
  tower.append(z);T=join(T,C);Gs.append(autos(C))
 c=stats(C)
 require(c['stabilizer_histogram']=={3:3},'C3 hostile')
 require(tower[0]['stabilizer_histogram']=={1:2,5:5},'R5 histogram')
 # Every quotient-map fibre has cardinality three by factor enumeration.
 base=range(7)
 for k in range(1,6):
  images={p:0 for p in itertools.product(base,*([range(3)]*(k-1)))}
  for p in itertools.product(base,*([range(3)]*k)):images[p[:-1]]+=1
  require(set(images.values())=={3},'tower restriction fibre')
 print(json.dumps({'status':'PASS','universe':'all labelled tournaments n=1..5; explicit tower k=0,1,2; Cartesian restriction controls k=1..5','labelled_tournaments':tested,'C3':c,'tower':tower,'no_all_order_coverage_claim':True},sort_keys=True,indent=2))
if __name__=='__main__':main()
