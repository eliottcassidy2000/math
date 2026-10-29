"""Literal sparse ordered path-cover audit; no determinant computations."""
from check_avoiding_dag import *
from itertools import permutations,combinations

def check(A,C):
 n=len(A);E,order,s=dag(A,C);need(len(order)==n,'unsafe input')
 f={v:i for i,v in enumerate(order)};ci={v:i for i,v in enumerate(C)};pred={s[v]:v for v in s};receipt=[]
 for k in range(1,n+1):
  ports=list(permutations(range(n),k));canon={}
  for a in ports:
   S=sorted(a,key=ci.get);nexts={S[i-1]:S[i] for i in range(k)}
   t=tuple(pred[nexts[v]] for v in a);canon[a]=t
   ends=sorted(t,key=ci.get);prev={ends[i]:ends[i-1] for i in range(k)}
   need(tuple(s[prev[v]] for v in t)==a,'inverse canonical matching failed')
  need(len(set(canon.values()))==len(ports),'not bijective')
  entries=collections.Counter()
  for p in permutations(range(n)):
   for cuts in combinations(range(1,n),k-1):
    cutset=set(cuts)
    if any(not A[p[i-1]][p[i]] for i in range(1,n) if i not in cutset):continue
    bounds=(0,)+cuts+(n,)
    a=tuple(p[bounds[i]] for i in range(k));t=tuple(p[bounds[i+1]-1] for i in range(k))
    entries[a,t]+=1
  for (a,t),v in entries.items():
   if t==canon[a]:need(v==1,'canonical entry not one')
   else:need(sum(f[u] for u in a)<sum(f[s[u]] for u in t),'potential failed')
  need(all(entries[a,t]==1 for a,t in canon.items()),'canonical missing')
  receipt.append({'k':k,'dimension':len(ports),'nonzero_entries':len(entries),'ordered_covers':sum(entries.values())})
 return receipt

def main():
 here=Path(__file__).resolve().parent
 data=json.loads((here.parent/'singularity/pinned/tournament_results.json').read_text());out=[]
 for ns,layer in data['layers'].items():
  if int(ns)>6:continue
  for rec in layer['records']:
   w=rec['word']
   if 'be' in w or 'ee' in w:continue
   A,C,_=build(w);out.append({'word':w,'order':len(A),'checks':check(A,C)})
 result={'status':'FINITE-EXACT independent literal enumeration','universe':'Every avoiding grammar word through order6, every k from1 to n','tournaments':len(out),'results':out}
 (here/'all_k_sparse_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
 print('PASS',len(out),'tournaments;',sum(len(x['checks']) for x in out),'all-k observers; literal permutations and cuts, sparse entries')
if __name__=='__main__':main()
