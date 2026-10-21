import json,pathlib,itertools,hashlib
P=pathlib.Path(__file__).resolve().parent;checks=0
def need(x):
 global checks
 checks+=1
 if not x:raise RuntimeError('independent check failed')
def sc(n,D,V):
 reach={x:{x}|{y for y in V if (y-x)%n in D} for x in V}
 for k in V:
  for x in V:
   if k in reach[x]:reach[x]|=reach[k]
 return all(len(reach[x])==len(V) for x in V)
def cut(n,D):
 for k in range(n-1):
  for S in itertools.combinations(range(n),k):
   if not sc(n,D,set(range(n))-set(S)):return k
 raise RuntimeError('no cut')
r=json.loads((P.parent/'session22-primecenters/receipt_normal.json').read_text());o=json.loads((P.parent/'session22-primecenters/receipt.json').read_text());nv=lambda x:{k:v for k,v in x.items() if k not in ('start_utc','end_utc','wall_seconds')};need(nv(r)==nv(o));need(r['checks']==17159)
for row in r['prime_connectivity_controls']:
 n=row['p'];D=set(row['connection_set']);need(cut(n,D)==(n-1)//2==row['minimum_cut_size'])
need(cut(9,{1,3,4,7})==3)
Y={0,3,6};D={1,3,4,7};S={1,4,7};X={2,5,8}
need(all((y-x)%9 in D for x in X for y in Y));need({(y+d)%9 for y in Y for d in D|{0}}==Y|S)
# Independent sumset representation: rotate residue bitmasks and union.
paircounts={}
for n in (3,5,7):
 full=(1<<n)-1;count=0
 for a in range(1,full+1):
  for b in range(1,full+1):
   s=0
   for i in range(n):
    if a>>i&1:s|=((b<<i)|(b>>(n-i)))&full
   need(s.bit_count()>=min(n,a.bit_count()+b.bit_count()-1));count+=1
 paircounts[str(n)]=count
need(paircounts==r['prime_subset_pair_counts'])
expected=set(range(11,48,2))|{12,18,20,24,28};expected-={37,41,43,47};need(sorted(expected)==r['unknown_support'] and len(expected)==20 and max(expected)==45)
out={'checks':checks,'producer_checks_per_mode':r['checks'],'nonvolatile_equal':True,'prime_additive_pair_counts':paircounts,'prime_connectivity_orientations':len(r['prime_connectivity_controls']),'composite9_minimum_cut':3,'degree_bound':45,'unknown_coefficients':20,'safe_recurrence_onset':58,'scope':'tiny controls and arithmetic only, not large center census'}
(P/'independent_results.json').write_text(json.dumps(out,indent=2)+'\n');print(out)
