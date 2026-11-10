#!/usr/bin/env python3
"""Independent literal DFS, orbit-completion and equality checks. No producer imports."""
import itertools,json,sys,math

def require(x,msg):
 if not x: raise RuntimeError(msg)

def hp(a):
 n=len(a)
 def walk(last,seen):
  if seen==(1<<n)-1:return 1
  return sum(walk(v,seen|(1<<v)) for v in range(n) if not seen>>v&1 and a[last][v])
 return sum(walk(v,1<<v) for v in range(n)) if n else 1

def invariant_tournaments(profile):
 blocks=[]; n=0
 for d in profile:blocks.append(list(range(n,n+d)));n+=d
 perm=list(range(n))
 for b in blocks:
  for i,v in enumerate(b):perm[v]=b[(i+1)%len(b)]
 unseen=set(itertools.combinations(range(n),2));orbits=[]
 while unseen:
  x,y=min(unseen);orb=[];u,v=x,y
  while (u,v) not in orb:
   require((v,u) not in orb,'pair reversal');orb.append((u,v));u,v=perm[u],perm[v]
  unseen.difference_update(tuple(sorted(pair)) for pair in orb);orbits.append(orb)
 for bits in itertools.product((0,1),repeat=len(orbits)):
  a=[[False]*n for _ in range(n)]
  for bit,orb in zip(bits,orbits):
   for u,v in orb:a[u if not bit else v][v if not bit else u]=True
  yield a,blocks

def submatrix(a,b):return [[a[u][v] for v in b] for u in b]
def ordinal(a,blocks):
 t=len(blocks);q=[[False]*t for _ in range(t)]
 for i,j in itertools.combinations(range(t),2):
  vals={a[u][v] for u in blocks[i] for v in blocks[j]}
  if len(vals)!=1:return False
  q[i][j]=vals.pop();q[j][i]=not q[i][j]
 return not any(q[i][j] and q[j][k] and q[k][i] for i,j,k in itertools.permutations(range(t),3))

def main():
 minima={1:1};rows=[]
 for d in (3,5,7,9):
  vals=[hp(a) for a,b in invariant_tournaments((d,))];minima[d]=min(vals)
  rows.append({'profile':[d],'cases':len(vals),'minimum':min(vals),'histogram':{str(h):vals.count(h) for h in sorted(set(vals))}})
 expected=176 if '--inject-failure' in sys.argv else 175
 require(minima[3]==3 and minima[5]==15 and minima[7]==expected,'prime minima')
 # Independently identify all six minimizing connection sets as affine Q7 images.
 actual={frozenset(j for j in range(1,7) if a[0][j]) for a,b in invariant_tournaments((7,)) if hp(a)==175}
 expectedsets={frozenset(u*j%7 for j in (1,2,3)) for u in range(1,7)}
 require(actual==expectedsets,'Q7 affine minimizers')
 # The inherited P7 affine group has no element of order21.
 orders=[]
 for a in (1,2,4):
  for b in range(7):
   perm=[(a*x+b)%7 for x in range(7)]; cur=list(range(7));order=0
   while True:
    order+=1;cur=[perm[x] for x in cur]
    if cur==list(range(7)):break
   orders.append(order)
 require(set(orders)=={1,3,7},'P7 cyclic subgroup orders')
 total=0
 for profile in ((3,3),(3,1,1),(5,1),(3,3,1),(3,5),(1,1,1)):
  cases=equal=0
  for a,blocks in invariant_tournaments(profile):
   h=hp(a); hs=[hp(submatrix(a,b)) for b in blocks];bound=math.prod(hs);floor=math.prod(minima[d] for d in profile)
   require(h>=bound>=floor,'product lower bound')
   require((h==bound)==ordinal(a,blocks),'product equality iff ordinal')
   require((h==floor)==(ordinal(a,blocks) and all(c==minima[d] for c,d in zip(hs,profile))),'minimum equality')
   require(h%math.lcm(*profile)==0,'actual subgroup divisibility')
   cases+=1;equal+=h==floor
  rows.append({'profile':profile,'cases':cases,'equality_cases':equal});total+=cases
 print(json.dumps({'status':'PASS','profile_cases':total,'rows':rows,'scope':'DFS through order9; inherited vertex-transitive classification not rerun'},indent=2))
if __name__=='__main__':main()
