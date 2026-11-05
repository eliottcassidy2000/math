"""Exact physical adjoint, midpoint inverse, and independent graph controls."""
from pathlib import Path
from itertools import product, permutations
from math import isqrt
from functools import lru_cache
import argparse, hashlib, heapq, json

T = (1,0,2,1)
B = (3,1,3,1)
Q = ((0,0,0,1),(0,1,0,-1),(0,0,1,-1),(1,-1,-1,1))
P = ((1,0,0,0),(0,-1,0,0),(0,0,1,0),(0,0,0,0))
RAW = {
 'O':((1,0,0,0),(0,0,0,0),(0,0,2,0),(0,0,0,1)),
 'Z':((1,-1,1,0),(0,0,0,1),(0,0,2,-1),(0,0,0,1)),
 'J':((1,-1,1,0),(1,-1,0,0),(1,0,1,0),(1,0,0,0)),
}
SIGMA = {'O':'Z','Z':'O','J':'J'}
CHANGE = ((1,0,2,1),(0,0,1,0),(1,1,2,2),(0,0,1,1))
QP = ((1,1,2,2),(1,1,1,1),(2,1,4,3),(2,1,3,3))

def need(ok, msg):
 if not ok: raise RuntimeError(msg)

def transpose(A): return tuple(zip(*A))
def dot(x,y): return sum(a*b for a,b in zip(x,y))
def mv(A,v): return tuple(dot(row,v) for row in A)
def mm(A,C): return tuple(tuple(dot(row,col) for col in transpose(C)) for row in A)
def rowmv(v,A): return tuple(dot(v,col) for col in transpose(A))
def pairing(v,w): return dot(v,mv(Q,w))
def dagger(word): return ''.join(SIGMA[c] for c in reversed(word))
def F(s): return s[0]-s[1]+s[2]
def G(s): return s[0]+2*s[2]

def step(s,c):
 D,E,U,d=s
 if c=='O': return D,0,2*U,d
 if c=='Z': return D-E+U,d,2*U-d,d
 return D-E+U,D-E,D+U,D

def state(word, seed=T):
 for c in word: seed=step(seed,c)
 return seed

def source_bank(cap, seed=T):
 """Topological minimum-J relaxation. Every predecessor has smaller U."""
 if seed[2]>cap: return {}
 out={seed:(0,'')}; heap=[(seed[2],seed)]; done=set()
 while heap:
  _,s=heapq.heappop(heap)
  need(s not in done,'one topological extraction')
  done.add(s);cost,word=out[s]
  for c in 'OZJ':
   q=step(s,c)
   need(s[2]<q[2]<=2*s[2],'physical cover step bounds')
   if q[2]>cap: continue
   candidate=(cost+(c=='J'),word+c)
   if q not in out:
    out[q]=candidate;heapq.heappush(heap,(q[2],q))
   elif candidate[0]<out[q][0]:
    need(q not in done,'minimum-cost update before extraction')
    out[q]=candidate
 return out

def construct_pair(vword,wword,balanced=False):
 if balanced: return vword+dagger(wword)+'OJ'
 return wword+dagger(vword)+'J'

def inverse(n, balanced=False, budget=None):
 if balanced and n==3:
  return dict(status='ACCEPT',target=3,boundary='singleton OJ gives triangle',rows=[2,4,1],J_cost=0)
 if n<3 or n%2==0: return dict(status='REJECT',target=n,scope='selected C1 scalar class')
 cap=isqrt(4*(n-1));left=source_bank(cap);right=source_bank(cap,B) if balanced else left
 ls=sorted(left,key=lambda s:(s[2],s));rs=sorted(right,key=lambda s:(s[2],s));tests=0
 for i,v in enumerate(ls):
  cv,va=left[v]
  if budget is not None and cv>budget: continue
  row=rowmv(v,Q)
  for j,w in enumerate(rs):
   if not balanced and j<i: continue
   prod=v[2]*w[2]
   if prod<n:continue
   if prod>2*(n-1):break
   cw,wa=right[w]
   if budget is not None and cv+cw>budget:continue
   tests+=1
   if dot(row,w)==n:
    word=construct_pair(va,wa,balanced)
    need(state(word)[0]==n,'literal source-word Hamilton readout')
    return dict(status='ACCEPT',target=n,balanced=balanced,budget=budget,cover_cap=cap,
                left=v,right=w,left_word=va,right_word=wa,J_cost=cv+cw,
                graph_word_from_transitive_edge=word,graph_order=len(word)+2,
                tested_pairs=tests,left_inventory=len(left),right_inventory=len(right))
 return dict(status='REJECT',target=n,balanced=balanced,budget=budget,cover_cap=cap,
             tested_pairs=tests,left_inventory=len(left),right_inventory=len(right),
             scope='Complete inverse for the selected C1 class only, not the unrestricted tournament spectrum.')

def graph(word):
 rows=[2,0];root=0
 for c in word:
  n=len(rows);rows.append(0)
  for v in range(n):
   if c!='O' and v==root: rows[n]|=1<<v
   else: rows[v]|=1<<n
  if c=='J':root=n
 return rows,root

@lru_cache(None)
def graph_counts(word):
 rows,root=graph(word);n=len(rows);full=(1<<n)-1
 end=[[0]*n for _ in range(full+1)];fixed=[[0]*n for _ in range(full+1)]
 for i in range(n):end[1<<i][i]=1
 fixed[1][0]=1
 for mask in range(1,full+1):
  for j in range(n):
   if not mask>>j&1:continue
   rest=mask^(1<<j)
   if rest:
    end[mask][j]=sum(end[rest][i] for i in range(n) if rest>>i&1 and rows[i]>>j&1)
    if mask&1:fixed[mask][j]=sum(fixed[rest][i] for i in range(n) if rest>>i&1 and rows[i]>>j&1)
 H=sum(end[full]);cycles=sum(fixed[full][j]for j in range(1,n)if rows[j]&1)
 return H,cycles

def determinant(A):
 total=0
 for p in permutations(range(4)):
  sign=(-1)**sum(p[i]>p[j] for i in range(4)for j in range(i+1,4))
  value=sign
  for i in range(4):value*=A[i][p[i]]
  total+=value
 return total

def verify():
 need(determinant(Q)==-1,'unimodular pairing')
 for c in RAW:
  need(mm(transpose(RAW[c]),Q)==mm(Q,RAW[SIGMA[c]]),'physical adjoint '+c)
 need(mm(Q,RAW['J'])==P,'one-entry factor')
 need(rowmv(T,Q)==(1,-1,1,0)and rowmv(B,Q)==(1,0,2,0),'boundary readouts')
 need(mm(mm(transpose(CHANGE),Q),CHANGE)==QP,'positive coordinate substitution')
 u=(1,1,2,2);d=(0,0,1,1)
 upper=tuple(tuple(u[i]*u[j]-QP[i][j]for j in range(4))for i in range(4))
 lower=tuple(tuple(2*QP[i][j]-u[i]*u[j]-2*d[i]*d[j]for j in range(4))for i in range(4))
 need(min(x for row in upper for x in row)>=0 and min(x for row in lower for x in row)>=0,'positive product bounds')
 midpoint_checks=0
 for k in range(8):
  for letters in product('OZJ',repeat=k):
   word=''.join(letters);s=state(word)
   need(state(dagger(dagger(word)))==s,'word involution')
   for seed,readout in [(T,F),(B,G)]:
    n=readout(s);cap=isqrt(4*(n-1));found=False
    for cut in range(k+1):
     v=state(word[:cut]);w=state(dagger(word[cut:]),seed)
     need(pairing(v,w)==n,'exact midpoint transfer')
     need(word[:cut].count('J')+dagger(word[cut:]).count('J')==word.count('J'),'root-change budget')
     need(v[2]*w[2]<=2*(n-1)and n<=v[2]*w[2],'physical product interval')
     if max(v[2],w[2])<=cap:found=True
    need(found,'two-sided square-root midpoint')
    midpoint_checks+=1
 small=[''.join(w)for k in range(3)for w in product('OZJ',repeat=k)]
 literal=0
 for alpha in small:
  for beta in small:
   word=construct_pair(alpha,beta)
   n=pairing(state(alpha),state(beta))
   need(graph_counts(word)==(n,1),'literal full C1 pair')
   word=construct_pair(alpha,beta,True)
   n=pairing(state(alpha),state(beta,B))
   need(graph_counts(word)==(n,1),'literal balanced pair')
   literal+=2
 N=5001;cap=isqrt(4*(N-1));big=source_bank(N-1);left=source_bank(cap);right=source_bank(cap,B)
 direct_full={F(s)for s in big if F(s)<=N}
 direct_bal={3}|{G(s)for s in big if G(s)<=N}
 direct_budget=[{3}|{G(s)for s,(c,w)in big.items()if c<=r and G(s)<=N}for r in range(4)]
 pair_full=set();pair_bal={3};pair_budget=[{3}for _ in range(4)];tests=0
 ss=list(left)
 for i,v in enumerate(ss):
  row=rowmv(v,Q)
  for w in ss[i:]:
   if v[2]*w[2]>2*(N-1):continue
   n=dot(row,w);tests+=1
   if n<=N:pair_full.add(n)
  for w,(cw,_)in right.items():
   if v[2]*w[2]>2*(N-1):continue
   n=dot(row,w);tests+=1
   if n<=N:
    pair_bal.add(n);cost=left[v][0]+cw
    for r in range(cost,4):pair_budget[r].add(n)
 need(pair_full==direct_full,'full C1 complete inverse image')
 need(pair_bal==direct_bal,'balanced complete inverse image')
 need(pair_budget==direct_budget,'bounded-root-change exact inverse images')
 samples=[]
 for n in [3,5,7,11,13,19,67,129,541,2647,3187,4999,5001]:
  rec=inverse(n)
  need((rec['status']=='ACCEPT')==(n in direct_full),'individual C1 inverse')
  samples.append(rec)
 for n in [3,5,7,11,19,3187,4999,5001]:
  rec=inverse(n,True,3)
  need((rec['status']=='ACCEPT')==(n in direct_budget[3]),'individual four-block inverse')
  samples.append(rec)
 here=Path(__file__).resolve().parent
 return dict(status='PASS',determinant=-1,adjoint_letters=3,positive_matrix=QP,
             product_upper_certificate=upper,product_lower_certificate=lower,
             midpoint_word_checks=midpoint_checks,literal_pair_graphs=literal,
             finite_comparison=dict(target_cap=N,square_root_cover_cap=cap,forward_profiles=len(big),
                 ordinary_sqrt_profiles=len(left),triangle_sqrt_profiles=len(right),pair_tests=tests,
                 C1_support=len(pair_full),balanced_support=len(pair_bal),
                 bounded_root_support=[len(s)for s in pair_budget],
                 C1_gaps=[n for n in range(3,N+1,2)if n not in pair_full]),
             exact_target_examples=samples,
             dependency_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest()for p in sorted((here/'dependencies').glob('*.md'))},
             scope='All-order physical adjoint and complete square-root-height inverses; universal odd coverage remains open.')

if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--target',type=int);parser.add_argument('--balanced',action='store_true');parser.add_argument('--budget',type=int)
 args=parser.parse_args()
 if args.budget is not None:need(args.budget>=0,'nonnegative root budget')
 result=verify()if args.target is None else inverse(args.target,args.balanced or args.budget is not None,args.budget)
 print(json.dumps(result,indent=2,sort_keys=True))
