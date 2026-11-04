import json,itertools,collections
from pathlib import Path
D=Path(__file__).parent

def subs(m):
 s=m
 while True:
  yield s
  if not s:break
  s=(s-1)&m

def readmodels(models,vals):
 it=iter(map(int,vals.split()));out=[]
 for m in models:
  n=next(it);L=next(it);masks=[s for s in range(1<<n) if s.bit_count()<=2]
  if len(masks)!=L or n!=len(m['rows']):raise RuntimeError('dimension mismatch')
  values=[[next(it) for _ in masks] for _ in masks];pos={a:i for i,a in enumerate(masks)}
  coeff=[];coarse=collections.defaultdict(int)
  which={v:i for i,b in enumerate(m['blocks']) for v in b}
  bm={a:sum(1<<which[v] for v in range(n) if a>>v&1) for a in masks}
  # Two physical names in one child map to the same idempotent bit.
  bm={a:sum(1<<i for i,b in enumerate(m['blocks']) if any(a>>v&1 for v in b)) for a in masks}
  for a in masks:
   row=[]
   for c in masks:
    v=sum((-1)**((a.bit_count()+c.bit_count()-s.bit_count()-t.bit_count()))*values[pos[s]][pos[t]] for s in subs(a) for t in subs(c))
    row.append(v);coarse[bm[a],bm[c]]+=v
   coeff.append(row)
  bs=[a for a in range(128) if a.bit_count()<=2]
  cb=[[coarse[a,c] for c in bs] for a in bs]
  cv=[[sum(coarse[s,t] for s in subs(a) for t in subs(c)) for c in bs] for a in bs]
  out.append(dict(**m,named_masks=masks,named_values=values,named_coefficients=coeff,block_masks=bs,block_coefficients=cb,block_values=cv))
 try:next(it);raise RuntimeError('trailing integers')
 except StopIteration:pass
 return out

def flatten(x):return sum(x,[])
def bitvec(x):
 if any(v%4 for v in x):raise RuntimeError('not binary contrast')
 return sum(((v%8)//4)<<i for i,v in enumerate(x))
def solve(base,cols):
 basis={}
 for i,c in enumerate(cols):
  x=c;w=1<<i
  while x:
   j=x.bit_length()-1
   if j not in basis:basis[j]=(x,w);break
   x^=basis[j][0];w^=basis[j][1]
 x=base;w=0
 while x:
  j=x.bit_length()-1
  if j not in basis:return len(basis),None
  x^=basis[j][0];w^=basis[j][1]
 return len(basis),w

def summarize(M):
 byQ=collections.defaultdict(list)
 for m in M:byQ[tuple(m['Q'])].append(m)
 out=[]
 for Q,ms in byQ.items():
  if len(ms)!=72:continue
  d={}
  for m in ms:d[tuple(m['sizes']),m['side']]=flatten(m['block_coefficients'])
  both=[]
  for side in range(2):
   base=d[(1,)*7,side];k2=[];k3=[];k22=[]
   for i in range(7):
    ss=[1]*7;ss[i]=2;v2=d[tuple(ss),side];ss[i]=3;v3=d[tuple(ss),side]
    if any((a-b)%2 for a,b in zip(v2,base)):raise RuntimeError('K2 divisibility')
    if any((a-3*b+2*c)%12 for a,b,c in zip(v3,v2,base)):raise RuntimeError('K3 even divisibility')
    k2.append([(a-b)//2 for a,b in zip(v2,base)]);k3.append([(a-3*b+2*c)//6 for a,b,c in zip(v3,v2,base)])
   for i,j in itertools.combinations(range(7),2):
    ss=[1]*7;ss[i]=ss[j]=2;v=d[tuple(ss),side]
    num=[a-b-2*c-2*e for a,b,c,e in zip(v,base,k2[i],k2[j])]
    if any(x%4 for x in num):raise RuntimeError('K22 divisibility')
    k22.append([x//4 for x in num])
   both.append(dict(K0=base,K2=k2,K3=k3,K22=k22))
  K={key:[[a-b for a,b in zip(x,y)] for x,y in zip(both[0][key],both[1][key])] if key!='K0' else [a-b for a,b in zip(both[0][key],both[1][key])] for key in both[0]}
  base=[(K['K0'][k]+2*sum(z[k] for z in K['K2'])+4*sum(z[k] for z in K['K22']))%8 for k in range(841)]
  cols=[bitvec([4*x for x in v]) for v in K['K2']]+[bitvec([2*x for x in v]) for v in K['K3']]
  rank,sol=solve(bitvec(base),cols)
  oo=dict(Q=Q,full_run_coefficients=both,rank=rank,blind_solution=sol,blind_type_count=0 if sol is None else 2**(14-rank),all_pairs_contrast=base)
  if sol is not None:oo['one_blind_type']=['PECO'[((sol>>i)&1)+2*((sol>>(i+7))&1)] for i in range(7)]
  out.append(oo)
 return out
if __name__=='__main__':
 M=readmodels(json.load(open(D/'models.json')),(D/'named_values.txt').read_text());json.dump(M,open(D/'named_banks.json','w'),separators=(',',':'))
 S=summarize(M);json.dump(S,open(D/'affine_banks.json','w'),separators=(',',':'))
 print(json.dumps([{k:v for k,v in s.items() if k not in ['full_run_coefficients','all_pairs_contrast']} for s in S],indent=2))
