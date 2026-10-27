#!/usr/bin/env python3
"""Exact independent verifier. Standard library; no removable assertions."""
import hashlib,itertools,json,math,pathlib,random
from fractions import Fraction
from functools import lru_cache
ROOT=pathlib.Path(__file__).resolve().parent

def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def prime_paley(n):
 need(n%4==3 and all(n%d for d in range(2,math.isqrt(n)+1)),"prime Paley input")
 sq={pow(t,2,n) for t in range(1,n)}
 return [[0 if x==y else 1 if (y-x)%n in sq else -1 for y in range(n)] for x in range(n)]
def gram(k):return [[sum(a*b for a,b in zip(r,s)) for s in k] for r in k]
def flat_check(k):
 n=len(k);need(all(k[i][j]==-k[j][i] and (k[i][j] in (-1,1) if i!=j else k[i][j]==0) for i in range(n) for j in range(n)),"skew tournament")
 need(gram(k)==[[n-1 if i==j else -1 for j in range(n)] for i in range(n)],"flat Gram")
def vertices(mask):
 while mask:
  bit=mask&-mask;mask^=bit;yield bit.bit_length()-1

def partitions(items, count):
 """Independent restricted-growth set partitions into exactly count blocks."""
 blocks=[[] for _ in range(count)]
 def rec(i,used):
  if len(items)-i<count-used:return
  if i==len(items):
   if used==count:yield tuple(tuple(b) for b in blocks)
   return
  for j in range(min(count,used+1)):
   blocks[j].append(items[i]);yield from rec(i+1,max(used,j+1));blocks[j].pop()
 yield from rec(0,0)

def matrix_from_parse(tree,n):
 k=[[0]*n for _ in range(n)]
 def rec(node):
  if 'vertex' in node:return [node['vertex']]
  blocks=[rec(c) for c in node['children']];kind=node['kind'];q=len(blocks);sign=node['orientation']
  need(sign in (-1,1) and ((kind=='join' and q==2) or (kind=='cycle' and q==3)),"parse node")
  for i in range(q):
   for j in range(i+1,q):
    a=sign*(1 if kind=='join' or j-i==1 else -1)
    for u in blocks[i]:
     for v in blocks[j]:k[u][v]=a;k[v][u]=-a
  return sum(blocks,[])
 labels=rec(tree);need(sorted(labels)==list(range(n)),"parse names exactly once")
 return k

def check_bellman(case):
 n=case['n'];k=prime_paley(n);bounds=case['subset_lower_bounds'];N=1<<n
 need(len(bounds)==N and bounds[0]==0 and all(type(x)==int and x>=0 for x in bounds),"bound vector")
 wins=[sum(1<<v for v in range(n) if k[u][v]>0) for u in range(n)]
 @lru_cache(None)
 def directed_cost(a,b):
  # Cost if every baseline vertex in a points to every vertex in b.
  return sum((b & ~wins[u]).bit_count() for u in vertices(a))
 counts=[0,0]
 for mask in range(1,N):
  vs=tuple(vertices(mask))
  if len(vs)==1:need(bounds[mask]==0,"singleton lower bound")
  for q in (2,3):
   for blocks in partitions(vs,q):
    masks=tuple(sum(1<<v for v in block) for block in blocks)
    internal=sum(bounds[m] for m in masks)
    if q==2:
     a,b=masks;cost=min(directed_cost(a,b),directed_cost(b,a))
    else:
     a,b,c=masks
     cost=min(directed_cost(a,b)+directed_cost(b,c)+directed_cost(c,a),directed_cost(b,a)+directed_cost(c,b)+directed_cost(a,c))
    need(bounds[mask]<=internal+cost,f"Bellman lower certificate n={n} subset={mask}");counts[q-2]+=1
 witness=matrix_from_parse(case['witness'],n)
 bad=[(u,v) for u in range(n) for v in range(u+1,n) if witness[u][v]!=k[u][v]]
 need(len(bad)==case['minimum']==bounds[-1],"attaining witness")
 need(counts==case['root_partition_counts'],"partition universe")
 return dict(n=n,minimum_arc_flips=len(bad),partition_counts=counts,witness_exceptions=bad)

def module(k,sub,whole):
 return all(len({k[v][u] for u in vertices(sub)})==1 for v in vertices(whole^sub))
def grammar_member(k):
 @lru_cache(None)
 def go(s):
  if s.bit_count()<=3:return True
  m=(s-1)&s
  while m:
   if 1<m.bit_count()<s.bit_count() and module(k,m,s):
    representative=m&-m
    if go(m) and go((s^m)|representative):return True
   m=(m-1)&s
  return False
 return go((1<<len(k))-1)
def ball7():
 k=prime_paley(7);pairs=list(itertools.combinations(range(7),2));cnt=0
 for number in range(3):
  for choice in itertools.combinations(pairs,number):
   h=[r[:] for r in k]
   for u,v in choice:h[u][v]*=-1;h[v][u]*=-1
   need(not grammar_member(h),"Paley7 radius-two hostile");cnt+=1
 need(cnt==232,"ball universe")
 full=(1<<7)-1;proper=[s for s in range(1,full) if 1<s.bit_count()<7]
 need(len(proper)==119 and not any(module(k,s,full) for s in proper),"Paley7 prime module control")
 return dict(radius_two_masks=cnt,proper_nontrivial_modules_tested=len(proper))

# Independent finite-field arithmetic, with x^3=x+1 in characteristic three.
def digits(t):return (t%3,(t//3)%3,t//9)
def fadd(a,b):return sum(((x+y)%3)*3**i for i,(x,y) in enumerate(zip(digits(a),digits(b))))
def fneg(a):return sum((-x%3)*3**i for i,x in enumerate(digits(a)))
def fmul(a,b):
 c=[0]*5
 for i,x in enumerate(digits(a)):
  for j,y in enumerate(digits(b)):c[i+j]=(c[i+j]+x*y)%3
 for j in (4,3):c[j-3]=(c[j-3]+c[j])%3;c[j-2]=(c[j-2]+c[j])%3
 return sum(c[i]*3**i for i in range(3))
def fpow(a,t):
 result=1
 while t:
  if t&1:result=fmul(result,a)
  a=fmul(a,a);t//=2
 return result
def field27():
 need(all((x*x*x-x-1)%3 for x in range(3)),"irreducible cubic has no root")
 need(all(fmul(x,fpow(x,25))==1 for x in range(1,27)),"field inverses")
 chars=[]
 for x in range(27):
  z=fpow(x,13) if x else 0;need(z in (0,1,2),"field quadratic character");chars.append(0 if z==0 else 1 if z==1 else -1)
 k=[[chars[fadd(y,fneg(x))] for y in range(27)] for x in range(27)]
 flat_check(k);return k

def rank(matrix):
 a=[[Fraction(x) for x in row] for row in matrix];m=len(a);n=len(a[0]) if m else 0;r=0
 for col in range(n):
  pivot=next((j for j in range(r,m) if a[j][col]),None)
  if pivot is None:continue
  a[r],a[pivot]=a[pivot],a[r];v=a[r][col];a[r]=[x/v for x in a[r]]
  for j in range(r+1,m):
   if a[j][col]:
    v=a[j][col];a[j]=[x-v*y for x,y in zip(a[j],a[r])]
  r+=1
  if r==m:break
 return r

def tower_parse(level,counter=None):
 if counter is None:counter=[0]
 if not level:
  v=counter[0];counter[0]+=1;return {'vertex':v}
 return {'kind':'cycle','orientation':1,'children':[tower_parse(level-1,counter) for _ in range(3)]}
def tower_control():
 k=matrix_from_parse(tower_parse(3),27);need(all(sum(row)==0 for row in k),"regular tower")
 g=gram(k);nullities={str(lam):27-rank([[g[i][j]-(lam if i==j else 0) for j in range(27)] for i in range(27)]) for lam in (3,27,243)}
 need(nullities=={'3':18,'27':6,'243':2},"tower exact squared spectrum")
 need(rank(k)==26,"tower rank")
 examples=[]
 for level in (1,3,5,7):
  n=3**level;bracket=(n-1)-2*level*3**((level-1)//2);edge=Fraction(n*bracket,4)
  clock=max(0,math.ceil(edge/Fraction(n-1)))
  examples.append(dict(level=level,n=n,arc_lower=max(0,math.ceil(edge)),complement_lower=clock))
 need(examples[1]['arc_lower']==54 and examples[1]['complement_lower']==3,"order27 bound")
 need(examples[2]['arc_lower']==9234 and examples[2]['complement_lower']==39,"order243 bound")
 return dict(n=27,squared_singular_nullities=nullities,examples=examples)

def tree_labels(t):
 return {t['vertex']} if 'vertex' in t else set().union(*(tree_labels(c) for c in t['children']))
def restrict(t,keep):
 if 'vertex' in t:return t if t['vertex'] in keep else None
 # A restricted cycle with two surviving children can have either orientation.
 indexes=[i for i,c in enumerate(t['children']) if tree_labels(c)&keep]
 children=[restrict(t['children'][i],keep) for i in indexes]
 if not children:return None
 if len(children)==1:return children[0]
 sign=t['orientation']
 if len(children)==2:
  if t['kind']=='cycle' and indexes==[0,2]:sign*=-1
  return {'kind':'join','orientation':sign,'children':children}
 return {'kind':t['kind'],'orientation':sign,'children':children}
def comb(n,cyclic=False):
 if cyclic:need(n%2==1,"odd cyclic comb")
 tree={'vertex':0};v=1
 while v<n:
  q=2 if cyclic else 1
  tree={'kind':'cycle' if cyclic else 'join','orientation':1,'children':[tree]+[{'vertex':u} for u in range(v,v+q)]};v+=q
 return tree

def balanced_control(tree,b):
 n=len(tree_labels(tree));k=matrix_from_parse(tree,n);cut_sizes=[]
 def rec(t,depth=0):
  labels=tree_labels(t);size=len(labels)
  if size==1:return 0,depth
  cur=t
  while True:
   large=[c for c in cur['children'] if (b+1)*len(tree_labels(c))>b*size]
   if not large:break
   cur=large[0]
  child=max(cur['children'],key=lambda c:len(tree_labels(c)));left=tree_labels(child);right=labels-left;m=len(left)
  need((b+1)*m>size and (b+1)*m<=b*size,"balanced module inequalities")
  need(all(len({k[v][u] for u in left})==1 for v in right),"actual module cut")
  a,da=rec(restrict(t,left),depth+1);c,dc=rec(restrict(t,right),depth+1)
  rad=m*(size-m);ceilroot=math.isqrt(rad)+(math.isqrt(rad)**2<rad);cut_sizes.append([size,m])
  return a+c+2*ceilroot,max(da,dc)
 budget,depth=rec(tree);need(len(cut_sizes)==n-1,"binary certificate leaves")
 return dict(n=n,b=b,cut_count=len(cut_sizes),max_balanced_depth=depth,integer_nuclear_upper=budget)

def native_controls():
 total=0
 for n in range(2,7):
  pairs=list(itertools.combinations(range(n),2));free=[p for p in pairs if p[1]!=p[0]+1];path={p for p in pairs if p[1]==p[0]+1}
  for bits in range(1<<len(free)):
   k=[[0]*n for _ in range(n)]
   for u,v in pairs:k[u][v]=1;k[v][u]=-1
   for j,(u,v) in enumerate(free):
    if bits>>j&1:k[u][v]=-1;k[v][u]=1
   sign=1 if bits%2 else -1;baseline=[[sign*(1 if u<v else -1) if u!=v else 0 for v in range(n)] for u in range(n)]
   d={p for p in pairs if k[p[0]][p[1]]!=baseline[p[0]][p[1]]}
   new=[row[:] for row in k]
   for u,v in free:new[u][v]*=-1;new[v][u]*=-1
   newd={p for p in pairs if new[p[0]][p[1]]!=-baseline[p[0]][p[1]]}
   need(newd==d^path and abs(len(newd)-len(d))<=n-1,"native XOR tariff")
   for u in range(n):
    for v in range(n):
     l=(1 if v==u+1 else -1 if u==v+1 else 0)
     need(new[u][v]==-k[u][v]+2*l,"native matrix identity")
   total+=1
 need(total==1099,"native marked universe")
 return dict(all_fixed_path_tournaments_n2_to_n6=total)

def miscellaneous():
 # Constants and sufficient threshold: purely integer comparisons.
 need(3**5<2**8 and 2700<2704,"rational C3 enclosure")
 need(169*512*81<=36*511**2,"n512 threshold")
 t=511;n=2*t+1;raw_lower=sum(2*math.isqrt(4*j-1) for j in range(1,t+1));true_upper=Fraction(13*n*((n-1).bit_length()),12)
 need(raw_lower>true_upper,"comb defeats naive node triangle budget")
 # Clone C3[T2,T2,T2] loses the three internal pairs on quotient-only reduction.
 weights=[2,2,2];cross_frob=sum(weights[i]*weights[j] for i in range(3) for j in range(3) if i!=j)
 need(cross_frob==24 and 6*5-cross_frob==6,"clone internal Frobenius loss")
 # A rank-two block can contain quadratically many pair exceptions.
 a,b=5,7;block=[[0]*(a+b) for _ in range(a+b)]
 for i in range(a):
  for j in range(a,a+b):block[i][j]=2;block[j][i]=-2
 need(rank(block)==2 and sum(x*x for row in block for x in row)==8*a*b,"block-rank currency")
 # Weighted dual identity, with positive unequal rational vertex costs.
 k=prime_paley(7);base=matrix_from_parse(comb(7,True),7);w=[Fraction(i+1,7) for i in range(7)]
 m=sum(w[i]*w[j] for i in range(7) for j in range(i+1,7) if k[i][j]!=base[i][j]);inner=sum(w[i]*w[j]*k[i][j]*base[i][j] for i in range(7) for j in range(7))
 need(inner==sum(w)**2-sum(x*x for x in w)-4*m,"weighted dual inner identity")
 return dict(comb_n=n,naive_triangle_sum_lower=raw_lower,true_nuclear_bound_upper=str(true_upper),clone_frobenius_loss=6,rank_two_block_arc_count=a*b,weighted_pair_debt=str(m))

def main():
 pins=json.loads((ROOT/'SOURCE_PINS.json').read_text())
 for pin in pins:need(hashlib.sha256((ROOT/pin['snapshot']).read_bytes()).hexdigest()==pin['sha256'],"pinned source hash")
 certificate=json.loads((ROOT/'certificate.json').read_text());distance=[check_bellman(c) for c in certificate['cases']]
 need([x['minimum_arc_flips'] for x in distance]==[0,3,11],"finite distance results")
 prime_orders=(3,7,11,19,31,43,59)
 for p in prime_orders:flat_check(prime_paley(p))
 field27()
 result=dict(status='PASS',source_pins=len(pins),exact_distances=distance,paley7_hostile=ball7(),flat_prime_orders=list(prime_orders),flat_field_order=27,tower=tower_control(),balanced_cut_controls=[balanced_control(comb(32),2),balanced_control(comb(31,True),3),balanced_control(tower_parse(3),3)]+[balanced_control(c['witness'],3) for c in certificate['cases']],native=native_controls(),miscellaneous=miscellaneous())
 print(json.dumps(result,sort_keys=True,indent=2))
if __name__=='__main__':main()
