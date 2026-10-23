#!/usr/bin/env python3
import itertools,json,math
def require(condition, message):
 if not condition: raise RuntimeError(message)
rows=[]
for n in range(1,7):
 pairs=list(itertools.combinations(range(n),2)); total=assoc=tri_free=0
 for mask in range(1<<len(pairs)):
  win=[[i for j in range(n)] for i in range(n)]
  for b,(i,j) in enumerate(pairs): win[i][j]=win[j][i]=j if mask>>b&1 else i
  a=all(win[win[i][j]][k]==win[i][win[j][k]] for i,j,k in itertools.product(range(n),repeat=3))
  t=all(not (win[i][j]==i and win[j][k]==j and win[k][i]==k) and not(win[i][j]==j and win[j][k]==k and win[k][i]==i) for i,j,k in itertools.combinations(range(n),3))
  require(a==t, 'explicit check at original line 13')
  assoc+=a;tri_free+=t;total+=1
 require(assoc==math.factorial(n), 'explicit check at original line 15')
 rows.append(dict(n=n,tournaments=total,associative=assoc,triangle_free=tri_free))
# normalized Z2 cocycle: pentagon on all 16 four-tuples.
w=lambda a,b,c:(-1)**(a*b*c)
for a,b,c,d in itertools.product(range(2),repeat=4):
 require(w(b,c,d)*w(a,b^c,d)*w(a,b,c)==w(a^b,c,d)*w(a,b,c^d), 'explicit check at original line 20')
# With r(a,b)=i^(ab), both scalar hexagons for associator w.
r=lambda a,b:1j**(a*b)
for a,b,c in itertools.product(range(2),repeat=3):
 require(w(a,b,c)*r(a,b^c)*w(b,c,a)==r(a,b)*w(b,a,c)*r(a,c), 'explicit check at original line 24')
 require(r(a^b,c)==w(a,b,c)*r(b,c)/w(a,c,b)*r(a,c)*w(c,a,b), 'explicit check at original line 25')
payload={'winner_product':rows,'Z2_pentagon_cases':16,'Z2_hexagon_cases':16}
# Independent quotient check: span of all associator differences versus SCCs.
quot=[]
for n in range(1,7):
 pairs=list(itertools.combinations(range(n),2)); counts={}
 for mask in range(1<<len(pairs)):
  win=[[i for j in range(n)] for i in range(n)]; reach=[[i==j for j in range(n)] for i in range(n)]
  for b,(i,j) in enumerate(pairs):
   wnr=j if mask>>b&1 else i; lsr=i+j-wnr
   win[i][j]=win[j][i]=wnr;reach[wnr][lsr]=True
  parent=list(range(n))
  def root(i):
   while parent[i]!=i: i=parent[i]
   return i
  for i,j,k in itertools.combinations(range(n),3):
   if (win[i][j]==i and win[j][k]==j and win[k][i]==k) or (win[i][j]==j and win[j][k]==k and win[k][i]==i):
    parent[root(j)]=root(i);parent[root(k)]=root(i)
  for k in range(n):
   for i in range(n):
    for j in range(n):reach[i][j]=reach[i][j] or (reach[i][k] and reach[k][j])
  require(all((root(i)==root(j)) == (reach[i][j] and reach[j][i]) for i in range(n) for j in range(n)), 'explicit check at original line 46')
  m=len({root(i) for i in range(n)});counts[m]=counts.get(m,0)+1
 quot.append({'n':n,'SCC_count_distribution':counts})
payload['associator_span_equals_SCC_kernel']=quot
print(json.dumps(payload,indent=2,sort_keys=True))
