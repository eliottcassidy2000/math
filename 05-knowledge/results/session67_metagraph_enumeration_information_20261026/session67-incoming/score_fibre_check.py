"""Independent labelled census n<=5 plus small constructed tree-fibre controls."""
from itertools import combinations, permutations, product
from collections import defaultdict
from pathlib import Path
import json, math
ROOT=Path(__file__).resolve().parent

def c3key(M):
 return tuple(int(all(sum(M[i][j] for j in t if j!=i)==1 for i in t)) for t in combinations(range(len(M)),3))
def scorekey(M): return tuple(map(sum,M))
def modules(M):
 n=len(M); ans=[]
 for mask in range(1,1<<n):
  S=frozenset(i for i in range(n) if mask>>i&1)
  if all(len({M[v][u] for u in S})==1 for v in range(n) if v not in S):ans.append(S)
 return ans

def tree(M):
 mods=modules(M); strong=[S for S in mods if not any(S&T and S-T and T-S for T in mods)]
 def build(S):
  if len(S)==1:return {'vertices':sorted(S),'children':[],'kind':'leaf','b':0,'N':1}
  proper=[T for T in strong if T<S]; children=sorted([T for T in proper if not any(T<U<S for U in strong)],key=lambda T:min(T)); kids=[build(T) for T in children];k=len(kids)
  Q=[[M[min(children[i])][min(children[j])] if i!=j else 0 for j in range(k)] for i in range(k)]
  linear=sorted(map(sum,Q))==list(range(k)); weights=list(map(len,children));bal=not linear and all(sum(weights[j]*(2*Q[i][j]-1) for j in range(k) if j!=i)==0 for i in range(k))
  return {'vertices':sorted(S),'children':kids,'Q':Q,'kind':'linear' if linear else 'prime','b':sum(t['b'] for t in kids)+int(bal),'N':math.prod(t['N'] for t in kids)*(math.factorial(k) if linear else 2)}
 return build(frozenset(range(len(M))))

def tree_fibre(M,tr):
 n=len(M)
 def generate(t):
  if t['kind']=='leaf':return [dict()]
  children=t['children'];k=len(children); Q=t['Q']; outs=[]
  if t['kind']=='linear':
   orders=list(permutations(range(k)));qs=[[[int(order.index(i)<order.index(j)) if i!=j else 0 for j in range(k)] for i in range(k)] for order in orders]
  else:qs=[Q,list(map(list,zip(*Q)))]
  for parts in product(*(generate(c) for c in children)):
   base={ij:x for p in parts for ij,x in p.items()}
   for q in qs:
    arcs=dict(base)
    for i in range(k):
     for j in range(k):
      if i!=j:
       for a in children[i]['vertices']:
        for b in children[j]['vertices']:arcs[a,b]=q[i][j]
    outs.append(arcs)
  return outs
 return [[[arcs.get((i,j),0) for j in range(n)] for i in range(n)] for arcs in generate(tr)]

def matrix(n,mask):
 M=[[0]*n for _ in range(n)]
 for z,(i,j) in enumerate(combinations(range(n),2)):M[i][j]=mask>>z&1;M[j][i]=1-M[i][j]
 return M
summary=[]
for n in range(1,6):
 groups=defaultdict(list);cgroups=defaultdict(int); rows=[]
 for mask in range(1<<(n*(n-1)//2)):
  M=matrix(n,mask);h=c3key(M);s=scorekey(M);tr=tree(M);groups[h,s].append(mask);cgroups[h]+=1;rows.append((h,s,tr['b'],tr['N']))
 bad=sum(len(groups[h,s])!=2**b or cgroups[h]!=N for h,s,b,N in rows)
 summary.append({'n':n,'labelled_tournaments':len(rows),'C3_score_fibres':len(groups),'mismatches':bad,'max_b':max(r[2] for r in rows)})
 if bad:raise RuntimeError((n,bad))

def sub(Q,children):
 sizes=list(map(len,children));off=[sum(sizes[:i]) for i in range(len(sizes))];n=sum(sizes);M=[[0]*n for _ in range(n)]
 for a,A in enumerate(children):
  for b,B in enumerate(children):
   for i in range(len(A)):
    for j in range(len(B)):M[off[a]+i][off[b]+j]=A[i][j] if a==b else Q[a][b]
 return M

def TT(n):return [[int(i<j) for j in range(n)] for i in range(n)]
C=[[0,1,0],[0,0,1],[1,0,0]];A=sub(TT(2),[TT(1),C]);cases={'weighted_balanced_222':sub(C,[TT(2)]*3),'weighted_unbalanced_123':sub(C,[TT(1),TT(2),TT(3)]),'nested_balanced':sub(C,[C]*3),'aut_counterexample_T':sub(C,[A]*3)}
controls={}
for name,M in cases.items():
 tr=tree(M);fibre=tree_fibre(M,tr);same=sum(scorekey(N)==scorekey(M) for N in fibre);same_c3=all(c3key(N)==c3key(M) for N in fibre)
 controls[name]={'n':len(M),'b':tr['b'],'predicted_score_fibre':2**tr['b'],'enumerated_score_fibre':same,'predicted_C3_fibre':tr['N'],'enumerated_C3_fibre':len(fibre),'all_same_C3':same_c3}
 if same!=2**tr['b'] or len(fibre)!=tr['N'] or not same_c3:raise RuntimeError(name)
receipt={'status':'pass','scope':'All labelled tournaments n=1..5; four explicit generated modular-tree fibres, not all order6/9/12 tournaments','census':summary,'controls':controls}
print(json.dumps(receipt,indent=2));(ROOT/'score-fibre-check.json').write_text(json.dumps(receipt,indent=2)+'\n')
