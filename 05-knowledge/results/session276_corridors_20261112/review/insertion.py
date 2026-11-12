"""Independent core/clock proof of the H13-to-H27 insertion obstruction.
No producer automaton or transfer imports. Core lengths use proved Busch bound.
"""
import contextlib,io,itertools,json
with contextlib.redirect_stdout(io.StringIO()):
 from languages import graph
 from independent import subset_count

def strong(G):
 n=len(G)
 for i in range(n):
  reach={i};todo=[i]
  while todo:
   u=todo.pop()
   for v in range(n):
    if G[u][v] and v not in reach:reach.add(v);todo.append(v)
  if len(reach)!=n:return False
 return True

def one(w):
 n=len(w);G=[[False]*(n+1) for _ in range(n+1)]
 for i,b in enumerate(w):
  for j in range(i+1,n):G[i][j]=True
  G[n][i]=bool(b);G[i][n]=not b
 return G

def sub(p,w):
 j=0
 for a in w:
  if j<len(p) and a==p[j]:j+=1
 return j==len(p)
# B(27)=6, so one-apex spine<=6 and same-apex spine<=5.
a={1:[()]};c13=[];c27=[]
for n in range(2,7):
 for mid in itertools.product(range(2),repeat=n-2):
  w=(1,)+mid+(0,);h=subset_count(one(w))
  if h<=27:a.setdefault(h,[]).append(w)
for n in range(1,6):
 for w in itertools.product(range(4),repeat=n):
  G=graph(w)
  if not strong(G):continue
  h=subset_count(G)
  if h==13:c13.append(w)
  if h==27:c27.append(w)
assert set(c13)=={tuple(map(int,x)) for x in ['120','212','230','302','310','321']}
# For a three-letter subsequence, every neutral reservoir may be saturated at3.
maxwords=[(0,)*3+c+(3,)*3 for c in c27]
for d in a:
 if 27%d or 27//d not in a:continue
 for x in a[d]:
  for y in a[27//d]:maxwords.append((0,)*3+x+(1,)*3+tuple(1+2*b for b in y)+(3,)*3)
assert not any(sub(p,w) for p in c13 for w in maxwords)
# Positive repair uses literal adjacency matrices only.
G=graph((0,1,2,0));mid=[r[:] for r in G];mid[1][2],mid[2][1]=mid[2][1],mid[1][2]
upper=[r[:] for r in mid];upper[0][2],upper[2][0]=upper[2][0],upper[0][2]
vals=[subset_count(X) for X in (graph((1,2,0)),graph((1,1,2,0)),G,mid,upper)]
assert vals==[13,25,13,15,27]
p=[2,0,1,3,4,5];ren=[[upper[i][j] for j in p] for i in p];assert ren==graph((2,0,1,0))
print(json.dumps({'status':'PASS','independent_H13_strong_cores':len(c13),'independent_H27_strong_cores':len(c27),'maximal_clock_templates':len(maxwords),'forbidden_subsequence_pairs_checked':len(maxwords)*len(c13),'positive_repair_counts':vals,'scope':'all word lengths by finite-core theorem and pattern-length clock saturation'},indent=2))
