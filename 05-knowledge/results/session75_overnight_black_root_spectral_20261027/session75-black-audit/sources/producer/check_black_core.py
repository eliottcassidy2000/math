#!/usr/bin/env python3
"""Session 75: black core geodesics. No producer imports, stdlib only.
Bit row i contains precisely its named outneighbors. R reverses path arcs;
K complements NONpath arcs. Every rejection uses explicit exceptions (-O safe).
"""
import itertools,json,collections

def rows(n,code):
 a=[0]*n
 for i in range(n):
  for j in range(i+1,n):
   if code&1:a[i]|=1<<j
   else:a[j]|=1<<i
   code>>=1
 return tuple(a)
def arc(a,i,j):return (a[i]>>j)&1
def paths(a,end=None):
 n=len(a); full=(1<<n)-1
 def rec(p,used):
  if used==full:
   if end is None or p[-1]==end:yield tuple(p)
   return
  z=a[p[-1]] & ~used & full
  if end is not None and used.bit_count()<n-1:z &= ~(1<<end)
  while z:
   b=z&-z;z-=b; v=b.bit_length()-1
   yield from rec(p+[v],used|b)
 for v in range(n):
  if end is None or n==1 or v!=end:yield from rec([v],1<<v)
def blue(a,p):
 n=len(a)
 return all(arc(a,p[i],p[j])==arc(a,p[n-1-j],p[n-1-i]) for i in range(n) for j in range(i+1,n))
def rev(a,p):
 b=list(a)
 for i,j in zip(p,p[1:]):b[i]^=1<<j;b[j]^=1<<i
 return tuple(b)
def strong(a,vs=None):
 vs=tuple(range(len(a))) if vs is None else tuple(vs)
 if not vs:return False
 mask=sum(1<<v for v in vs)
 for backwards in (False,True):
  seen=1<<vs[0]
  while True:
   old=seen
   for v in vs:
    if seen>>v&1:seen|= ((mask ^ (1<<v) ^ a[v]) if backwards else a[v])&mask
   if old==seen:break
  if seen!=mask:return False
 return True
def target(a):
 n=len(a)
 return any(a[v].bit_count() in (0,n-1) and strong(a,[u for u in range(n) if u!=v]) for v in range(n))
def nearreg(n):
 lo=(n-2)//2; hi=n//2
 for code in range(1<<(n*(n-1)//2)):
  a=rows(n,code)
  if all(lo<=r.bit_count()<=hi for r in a):yield a
def op(a):
 full=(1<<len(a))-1
 return tuple(full ^ (1<<i) ^ r for i,r in enumerate(a))
def construct(a):
 n=len(a); k=(n-1)//2; ds=[r.bit_count() for r in a]
 p=next((p for p in paths(a) if not blue(a,p) and ds[p[-1]]>=ds[p[0]]),None)
 if p is None:raise ValueError(('no first black descent',a))
 dual=(n%2==0 and ds[p[-1]]!=n//2)
 if dual:a=op(a);p=p[::-1]
 initial=a; v=p[-1];first=p[0];route=[p];a=rev(a,p)
 if k>1:
  q=next((q for q in paths(a,v) if q[0]!=first),None)
  if q is None:raise ValueError(('no split',initial,p))
  if blue(a,q):raise ValueError(('split blue',initial,p,q))
  route.append(q);a=rev(a,q)
 for j in range(2,k):
  q=next(paths(a,v))
  if blue(a,q):raise ValueError(('later blue',initial,route,q))
  route.append(q);a=rev(a,q)
 if not target(a):raise ValueError(('badtarget',initial,route,a))
 return {'dual':dual,'v':v,'paths':[list(p) for p in route],'terminal':a}
def circulant(n,S=None):
 S=set(range(1,(n+1)//2)) if S is None else set(S)
 return tuple(sum(1<<j for j in range(n) if (j-i)%n in S) for i in range(n))

def require(ok,why):
 if not ok:raise RuntimeError(why)
def cut_values(a):
 n=len(a);full=(1<<n)-1
 return {s:sum((a[v]&(full^s)).bit_count() for v in range(n) if s>>v&1) for s in range(1,full)}
def regulars(n):
 es=[(u,v) for u in range(n) for v in range(u+1,n)]; k=(n-1)//2
 a=[0]*n;d=[0]*n;left=[n-1]*n
 def go(i):
  if i==len(es):yield tuple(a);return
  u,v=es[i];left[u]-=1;left[v]-=1
  for x,y in ((u,v),(v,u)):
   d[x]+=1;a[x]|=1<<y
   if all(d[z]<=k<=d[z]+left[z] for z in (u,v)):yield from go(i+1)
   d[x]-=1;a[x]^=1<<y
  left[u]+=1;left[v]+=1
 yield from go(0)
def failed_anti(a,p):
 n=len(a)
 for i in range(n):
  for j in range(i+1,n):
   r=(p[n-1-j],p[n-1-i])
   if arc(a,p[i],p[j])!=arc(a,*r):return [p[i],p[j],*r]
 return None
def replay(a,certificate):
 # R route was optionally constructed on a converse. Recover an R route
 # on original a, and then a literal K route via alternating converses.
 pp=[tuple(p[::-1] if certificate['dual'] else p) for p in certificate['paths']]
 original=a; k=len(pp); n=len(a); moves=[]
 for j,p in enumerate(pp):
  require(sorted(p)==list(range(n)),'not permutation')
  require(all(arc(a,u,v) for u,v in zip(p,p[1:])),'not directed R path')
  require(failed_anti(a,p) is not None,'blue R path')
  b=rev(a,p); ca=cut_values(a);cb=cut_values(b)
  require(all(cb[s]-ca[s]==((s>>p[-1]&1)-(s>>p[0]&1)) for s in ca),'cut current')
  source=op(a) if j%2 else a;q=p[::-1] if j%2 else p
  require(all(arc(source,u,v) for u,v in zip(q,q[1:])),'not directed K path')
  # independent literal NONpath complement uses unordered edge sets
  pe={frozenset((u,v)) for u,v in zip(q,q[1:])}
  actual=tuple(sum(1<<v for v in range(n) if u!=v and bool(arc(source,u,v))==(frozenset((u,v)) in pe)) for u in range(n))
  expected=op(b) if (j+1)%2 else b
  require(actual==expected,'R to literal K conversion')
  moves.append({'R_path':p,'K_path':q,'K_source_rows':source,'K_target_rows':actual,'failed_anti_named_arcs':failed_anti(source,q)})
  a=b
 require(target(a),'not singleton plus strong core')
 require(k==min(min(r.bit_count(),n-1-r.bit_count()) for r in original),'endpoint lower bound not tight')
 require(len({p[0] if not certificate['dual'] else p[-1] for p in pp[:2]})==min(2,k),'lost distinct starts')
 return {'initial_rows':original,'dual_for_construction':certificate['dual'],'named_port':certificate['v'],'distance':k,'moves':moves,'terminal_R_rows':a}
def join(blocks):
 n=sum(len(b) for b in blocks);a=[0]*n;offset=0
 for b in blocks:
  future=((1<<n)-1)^((1<<(offset+len(b)))-1)
  for i,r in enumerate(b):a[offset+i]=(r<<offset)|future
  offset+=len(b)
 return tuple(a)
def seed_tests():
 seeds={4:[next(nearreg(4))],6:[(48,33,3,7,14,28),(36,17,10,35,13,22),(36,33,18,7,11,28),(48,9,3,37,14,22),(36,17,10,19,37,14)],7:[(112,97,67,7,14,28,56),(112,97,19,7,74,28,44),(112,73,35,21,70,26,44)]}
 return [replay(a,construct(a)) for ss in seeds.values() for a in ss]
def main():
 import hashlib,random,sys
 report={'status':'PASS','scope':'exact balanced black core distances; all-color cut theorem separate','complete_balanced_counts':{},'case_receipt_hashes':{},'seeds':seed_tests()}
 for n,universe in ((4,nearreg(4)),(6,nearreg(6)),(7,regulars(7))):
  count=0;h=hashlib.sha256()
  for a in universe:
   r=replay(a,construct(a));h.update(json.dumps(r,sort_keys=True).encode());count+=1
  require(count=={4:24,6:2640,7:2640}[n],'universe cardinality')
  report['complete_balanced_counts'][n]=count;report['case_receipt_hashes'][n]=h.hexdigest()
 examples=[]
 for n in (7,9,11,13):examples.append((f'cyclic_{n}',circulant(n)))
 for n in (7,11):examples.append((f'Paley_{n}',circulant(n,{x*x%n for x in range(1,n)})))
 rng=random.Random(75027)
 for n in (8,9,10,11,12,13):
  N=n if n%2 else n+1;a=circulant(N)
  for sample in range(6):
   for z in range(30):
    x,y,z=rng.sample(range(N),3)
    if arc(a,x,y) and arc(a,y,z) and arc(a,z,x):
     b=list(a)
     for u,v in ((x,y),(y,z),(z,x)):b[u]^=1<<v;b[v]^=1<<u
     a=tuple(b)
   b=a if n%2 else tuple(r&((1<<n)-1) for r in a[:n])
   examples.append((f'balanced_random_{n}_{sample}',b))
 report['examples']=[{'name':name,**replay(a,construct(a))} for name,a in examples]
 # Pure-blue boundary controls enumerate EVERY Hamilton word.
 c3=circulant(3);r5=circulant(5);hostiles=[('C3',c3),('R5',r5),('C3_to_C3',join([c3,c3])),('C3_to_1_to_C3',join([c3,(0,),c3]))]
 report['pure_blue_controls']=[]
 for name,a in hostiles:
  pp=list(paths(a));require(pp and all(blue(a,p) for p in pp),'pure blue hostile changed')
  report['pure_blue_controls'].append({'name':name,'rows':a,'hamilton_paths':len(pp),'black_paths':0})
 a=(966,1004,1016,1009,995,193,384,768,544,96)
 pp=((2,3,4,6,7,8,9,5,0,1),(3,0,2,5,9,8,7,6,4,1))
 require(min(cut_values(a).values())==1,'cut-hostile lambda')
 r=replay(a,{'dual':False,'v':1,'paths':pp})
 report['nonregular_target_separation']={'lambda_nonstrong':1,'minimum_semidegree_core_lower_bound':2,**r}
 # Named-start guard hostile: a first black word, followed by two BLUE
 # words with the same endpoints, reaches source + R5 + sink, not one core.
 a=(112,37,73,67,14,28,50)
 same=((6,1,5,4,2,3,0),(6,4,5,3,1,2,0),(6,5,4,3,2,1,0));colors=[]
 for p in same:
  require(all(arc(a,u,v) for u,v in zip(p,p[1:])), 'same-start hostile directedness')
  colors.append('blue' if blue(a,p) else 'black');a=rev(a,p)
 require(colors==['black','blue','blue'] and not target(a),'same-start hostile failure')
 report['same_start_hostile']={'initial_rows':(112,37,73,67,14,28,50),'R_paths':same,'colors':colors,'terminal_rows':a,'terminal_shape':'1 -> R5 -> 1','one_core_target':False}
 report['not_claimed']=['all-order black connectivity','all-order black distance=lambda','fixed named port attachment for every regular port','forest or transitive target distance']
 print(json.dumps(report,sort_keys=True,indent=2))
if __name__=='__main__':main()
