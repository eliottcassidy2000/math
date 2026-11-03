#!/usr/bin/env python3
"""Separate algorithms, no import from the producer."""
import itertools,json,math
from collections import Counter
C=Counter()
def check(v,k):
 C[k]+=1
 if not v:raise ValueError(k)
def graph(n,mask):
 A=[set() for _ in range(n)]
 for i,(a,b) in enumerate(itertools.combinations(range(n),2)):
  if mask>>i&1:A[a].add(b);A[b].add(a)
 return A
def stats(A):
 unused=set(range(len(A)));components=[]
 while unused:
  stack=[unused.pop()];V=set(stack)
  while stack:
   for w in A[stack.pop()]:
    if w in unused:unused.remove(w);V.add(w);stack.append(w)
  components.append((len(V),sum(len(A[v]) for v in V)//2))
 return all(e<v for v,e in components),all(e<=v for v,e in components)
def alpha(A):
 best=0
 for bits in range(1<<len(A)):
  if bits.bit_count()<=best:continue
  if all(not any(bits>>w&1 for w in A[v]) for v in range(len(A)) if bits>>v&1):best=bits.bit_count()
 return best
def colorable(A,k):
 order=sorted(range(len(A)),key=lambda v:-len(A[v]));colors={}
 def dfs(i):
  if i==len(order):return True
  v=order[i];used={colors[w] for w in A[v] if w in colors}
  for c in range(k):
   if c not in used:
    colors[v]=c
    if dfs(i+1):return True
  colors.pop(v,None);return False
 return dfs(0)
counts=[]
for n in range(1,7):
 F=P=0
 for mask in range(1<<(n*(n-1)//2)):
  A=graph(n,mask);f,p=stats(A)
  if f:F+=1;check(colorable(A,2),'forest_2_coloring')
  if p:P+=1;check(colorable(A,3),'pseudoforest_3_coloring')
 counts.append({'n':n,'forests':F,'pseudoforests':P})

# Distinct computation of the three-color K5 universe, inspecting degrees
# in a triangle, rather than the producer's combinations membership test.
pairs=list(itertools.combinations(range(5),2));valid=0;linear=[]
for mask in range(1024):
 A=graph(5,mask);f,p=stats(A)
 if f and max(map(len,A))<=2:linear.append(mask)
check(len(linear)==206,'linear_forest_count')
for a in linear:
 for b in linear:
  A=graph(5,a|b);check(colorable(A,4),'two_path_forest_coloring')
for code in range(3**10):
 d=code;R=[set() for _ in range(5)];B=[set() for _ in range(5)]
 for a,b in pairs:
  digit=d%3;d//=3
  if digit<2:
   E=(R,B)[digit];E[a].add(b);E[b].add(a)
 if any(any(A[u]&A[v] for u in range(5) for v in A[u]) for A in (R,B)):continue
 valid+=1;rf,rp=stats(R);bf,bp=stats(B);H=[R[i]|B[i] for i in range(5)]
 if rp and bp:check(colorable(H,5),'separate_pseudoforest_5_coloring')
 if (rf and bp) or (bf and rp):check(colorable(H,4),'mixed_4_coloring')
check(valid==29197,'K5_avoiding_census')

# Actual Hamilton word DFS. No subset dynamic program, and no arc-transitive
# reduction: the specified arc is tallied on every actual full word.
q=11;QR={i*i%q for i in range(1,q)};endpoint=[[0]*q for _ in range(q)];inc=0
A=[sum(1<<v for v in range(q) if (v-u)%q in QR) for u in range(q)]
def visit(start,last,unused,has_arc):
 global inc
 if not unused:
  endpoint[start][last]+=1;inc+=int(has_arc);return
 candidates=A[last]&unused
 while candidates:
  bit=candidates&-candidates;candidates-=bit;v=bit.bit_length()-1
  visit(start,v,unused^bit,has_arc or (last==0 and v==1))
for s in range(q):visit(s,s,((1<<q)-1)^(1<<s),False)
check(sum(map(sum,endpoint))==95095,'Paley_word_enumeration')
check(inc==17290,'Paley_direct_arc_enumeration')
for u in range(q):
 for v in range(q):check(endpoint[u][v]==(0 if u==v else 628 if (v-u)%q in QR else 1101),'Paley_full_endpoint_matrix')

# Independent matrix powers for the braid coloring hostile.
def multiply(A,B):return [[sum(A[i][k]*B[k][j] for k in range(2))%3 for j in range(2)] for i in range(2)]
M=[[2,2],[1,0]];M3=multiply(multiply(M,M),M)
check(M3==[[1,0],[0,1]],'Fox_matrix_cube')
def fixed(A):return sum(all(sum(A[i][j]*x[j] for j in range(2))%3==x[i] for i in range(2)) for x in itertools.product(range(3),repeat=2))
check((fixed(M),fixed(M3))==(3,9),'Fox_matrix_fixed_points')

# Closed-component correction: build the actual graph with m disjoint circles.
def pairings(vertices):
 if not vertices:
  yield []
  return
 a=vertices[0]
 for i in range(1,len(vertices)):
  for rest in pairings(vertices[1:i]+vertices[i+1:]):
   yield [(a,vertices[i])]+rest

def count_components(A):
 seen=set();total=0
 for v in range(len(A)):
  if v in seen:continue
  total+=1;stack=[v];seen.add(v)
  while stack:
   u=stack.pop()
   for w in A[u]:
    if w not in seen:seen.add(w);stack.append(w)
 return total

for r in range(1,4):
 P=[(2*i,2*i+1) for i in range(r)]
 for Q in pairings(list(range(2*r))):
  A=[set() for _ in range(2*r)]
  for u,v in P+Q:A[u].add(v);A[v].add(u)
  bare=count_components(A)
  for m in range(4):
   B=[set(x) for x in A]+[set() for _ in range(3*m)]
   for j in range(m):
    a=2*r+3*j
    for u,v in [(a,a+1),(a+1,a+2),(a+2,a)]:B[u].add(v);B[v].add(u)
   check(count_components(B)==m+bare,'closed_components_add_to_matching_loops')
   if m:check(count_components(B)!=bare,'extra_circle_hostile')

print(json.dumps({'status':'PASS','graph_counts':counts,'K5_avoiding':valid,'two_linear_forest_pairs':len(linear)**2,'paley11':{'H':sum(map(sum,endpoint)),'arc_0_1':inc,'endpoint_row_0':endpoint[0]},'fox_matrix_cube':M3,'checks':dict(sorted(C.items()))},indent=2,sort_keys=True))
