from pathlib import Path
from functools import lru_cache
from itertools import product
import hashlib,json,math
HERE=Path(__file__).resolve().parent;SRC=HERE.parent/'even_profile_linearization'
def need(c,m):
 if not c:raise RuntimeError(m)
def arcrows(Q,a,c):
 m=len(Q);R=list(Q)+[1<<(m+1),0]
 for u in range(m):
  for bits,p in[(a,m),(c,m+1)]:
   if bits>>u&1:R[u]|=1<<p
   else:R[p]|=1<<u
 return R

def run(Q,a,c,caps):
 R=arcrows(Q,a,c);m=len(Q)
 @lru_cache(None)
 def f(t,last):
  if not t[last]:return 0
  rem=list(t);rem[last]-=1;rem=tuple(rem)
  if not any(rem):return 1
  return sum(f(rem,u) for u in range(m+2) if R[u]>>last&1)
 return {e:sum(f(e+(1,1),u) for u in range(m+2)) for e in product(*(range(v+1) for v in caps))}
def physical(Q,sizes,a,c):
 off=[sum(sizes[:i]) for i in range(len(Q))];n=sum(sizes);T=[0]*n
 for i,k in enumerate(sizes):
  for u in range(k):T[off[i]+u]=sum(1<<(off[i]+j) for j in range(u+1,k))|sum(((1<<sizes[j])-1)<<off[j] for j in range(len(Q)) if Q[i]>>j&1)
 aa=sum(((1<<k)-1)<<off[i] for i,k in enumerate(sizes) if a>>i&1);cc=sum(((1<<k)-1)<<off[i] for i,k in enumerate(sizes) if c>>i&1);R=arcrows(T,aa,cc);N=len(R)
 @lru_cache(None)
 def f(mask,first):
  rest=mask^(1<<first)
  return 1 if not rest else sum(f(rest,v) for v in range(N) if rest>>v&1 and R[first]>>v&1)
 return sum(f((1<<N)-1,v) for v in range(N))
def main():
 src=json.loads((SRC/'normal.json').read_text());e=(0,1,14,4,8,0);types=[(0,1,b,c,8,d) for b,c,d in product([6,14],[4,12],[0,8])];formal=physicaln=queries=0;records=[]
 for rr in src['records']:
  Q=rr['Q'];m=len(Q);qc=[]
  for r in rr['full_block_query_controls']:
   a,c=r['a'],r['c'];K=run(Q,a,c,[5]*m);need([[*v,k] for v,k in sorted(K.items()) if k]==r['complete_integer_run_coefficients'],'complete integer run-word kernel');D=[]
   s=physical(Q,[1]*m,a,c);physicaln+=1;need(s==r['actual_singleton_response'],'independent physical singleton count')
   for i in range(m):
    t=[1]*m;t[i]=2;d=K[tuple(t)];D.append(d);v=physical(Q,t,a,c);physicaln+=1;need(v==r['actual_pair_responses'][i] and v-s==2*d,'exact physical derivative')
    for k in [3,5]:t[i]=k;need(K[tuple(t)]%2==0,'both higher corrections even')
   baseline=sum(v*math.prod(e[k] for k in t) for t,v in K.items());need(baseline==r['TT4_baseline_functional'],'baseline full integer functional')
   if m<=2:need(physical(Q,[4]*m,a,c)%16==baseline%16,'independent actual TT4 baseline');physicaln+=int(m<=2)
   for ids in product(range(8),repeat=m):
    gs=[types[i] for i in ids];v=sum(w*math.prod(gs[i][t[i]] for i in range(m)) for t,w in K.items())%16;pred=(baseline+8*sum(D[i] for i in range(m) if gs[i][2]==6))%16;need(v==pred,'complete formal profile states');formal+=1
   qc.append([a,c,baseline,D]);queries+=1
  records.append({'Q':Q,'query_count':len(qc),'complete_summary_sha256':hashlib.sha256(json.dumps(qc,separators=(',',':')).encode()).hexdigest()})
 h=src['higher_order_return_control'];Q=h['Q'];a,c=h['query'];swaps=[]
 for k,expected in h['coefficients']:
  caps=[k]+[1]*6;K=run(Q,a,c,caps);need(K[tuple(caps)]==expected,'nonzero seven-child control exact');R=arcrows(Q,a,c);left=caps+[1,1];words=[]
  def dfs(w):
   if not any(left):words.append(tuple(w));return
   for v in range(9):
    if left[v] and (not w or R[w[-1]]>>v&1):left[v]-=1;dfs(w+[v]);left[v]+=1
  dfs([]);S=set(words);need(len(S)==expected,'literal valid run words complete')
  def swap(w):
   p=[j for j,v in enumerate(w) if v==0];i,j,l=p[:3];return w[:i+1]+w[j+1:l]+(0,)+w[i+1:j]+w[l:]
  need(all(swap(w)!=w and swap(w) in S and swap(swap(w))==w for w in words),'explicit complete free return-swap orbits');swaps.append({'multiplicity':k,'word_count':expected,'word_set_sha256':hashlib.sha256(repr(sorted(S)).encode()).hexdigest()})
 need(formal==src['formal_response_checks'] and physicaln==src['literal_physical_counts'],'all source control totals');need((SRC/'normal.json').read_bytes()==(SRC/'optimized.json').read_bytes(),'paired source modes');pins={f:hashlib.sha256((SRC/f).read_bytes()).hexdigest() for f in ['PROOF.md','verify.py','normal.json','optimized.json']}
 print(json.dumps({'status':'PASS','formal_states':formal,'independent_physical_counts':physicaln,'queries':queries,'records':records,'explicit_return_involution_controls':swaps,'source_pins':pins,'audit_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'All-order Hamilton-only one-bit reduction. Not ordered covers and not universal noncollision.'},sort_keys=True,indent=2))
if __name__=='__main__':main()
