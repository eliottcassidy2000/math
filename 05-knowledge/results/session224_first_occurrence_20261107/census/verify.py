import csv,json,itertools,collections,math,pathlib,hashlib,sys
ROOT=pathlib.Path(sys.argv[1]) if len(sys.argv)>1 else pathlib.Path(__file__).resolve().parents[2]/'actual_hamilton_multiplicities_20261021'/'session24-multiplicities'
OUT=pathlib.Path(__file__).parent
checks=0
def need(v):
 global checks
 checks+=1
 if not v:raise ValueError('check '+str(checks))
def dual(a):return tuple(x^((1<<len(a))-1)^(1<<i) for i,x in enumerate(a))
def canon(a):
 n=len(a);cells=[[i for i,x in enumerate(a) if x.bit_count()==d] for d in range(n)];best=None;aut=0
 for pp in itertools.product(*(itertools.permutations(c) for c in cells if c)):
  p=sum(pp,());code=sum(1<<(8*i+j) for i in range(n) for j in range(n) if a[p[i]]>>p[j]&1)
  if best is None or code<best:best=code;aut=1
  elif code==best:aut+=1
 return best,aut
def paths(a):
 n=len(a);d=[[0]*n for _ in range(1<<n)]
 for i in range(n):d[1<<i][i]=1
 for s in range(1,1<<n):
  for i in range(n):
   if not d[s][i]:continue
   t=a[i]&~s
   while t:
    b=t&-t;t-=b;d[s|b][b.bit_length()-1]+=d[s][i]
 return sum(d[-1])
rows=[]
# Independent labelled brute force for all n<=5, without imported atlases.
for n in range(1,6):
 classes={};edges=list(itertools.combinations(range(n),2))
 for mask in range(1<<len(edges)):
  a=[0]*n
  for k,(i,j) in enumerate(edges):
   if mask>>k&1:a[i]|=1<<j
   else:a[j]|=1<<i
  code,aut=canon(a)
  if code not in classes:classes[code]=(a,aut)
 for code,(a,aut) in classes.items():rows.append(dict(n=n,code=code,H=paths(a),aut=aut,self_converse=int(canon(dual(a))[0]==code),a=a))
# Existing complete atlas, independently checked canonical/HP/aut channels.
src=json.load(open(ROOT/'census.json'))
for r in src['records']:
 if r['n'] not in (6,7):continue
 a=r['adjacency_outmasks']
 for b in ([a] if r['self_converse'] else [a,dual(a)]):
  code,aut=canon(b);need(aut==r['aut']);need(paths(b)==r['H']);need((canon(dual(b))[0]==code)==r['self_converse'])
  rows.append(dict(n=len(b),code=code,H=r['H'],aut=aut,self_converse=int(r['self_converse']),a=list(b)))
for r in csv.DictReader(open(ROOT/'n8_classes.csv')):
 r={k:int(v) for k,v in r.items()};a=[(r['code']>>(8*i))&255 for i in range(8)]
 code,aut=canon(a);need(code==r['code']);need(aut==r['aut']);need(paths(a)==r['H']);need((canon(dual(a))[0]==code)==bool(r['self_converse']))
 rows.append(dict(n=8,code=code,H=r['H'],aut=aut,self_converse=r['self_converse'],a=a))
expected=[1,1,2,4,12,56,456,6880];stats={};birth={}
for n in range(1,9):
 rr=[r for r in rows if r['n']==n];need(len(rr)==expected[n-1]);need(len({r['code'] for r in rr})==len(rr));need(sum(math.factorial(n)//r['aut'] for r in rr)==2**(n*(n-1)//2));need(sum(r['H']//r['aut'] for r in rr)==2**((n-1)*(n-2)//2));need(all(r['H']%2 and r['H']%r['aut']==0 for r in rr))
 codes={r['code'] for r in rr};need(all(canon(dual(r['a']))[0] in codes for r in rr))
 for h in sorted({r['H'] for r in rr}):
  if h in birth:continue
  fiber=[r for r in rr if r['H']==h];sc=sum(r['self_converse'] for r in fiber);need((len(fiber)-sc)%2==0)
  birth[h]=dict(H=h,first_n=n,ordinary=len(fiber),self_converse=sc,converse_pairs=(len(fiber)-sc)//2,converse_orbits=(len(fiber)+sc)//2)
 stats[n]=dict(classes=len(rr),SC=sum(r['self_converse'] for r in rr),new_H=sum(r['first_n']==n for r in birth.values()))
with open(OUT/'birth_fibers.csv','w') as f:
 w=csv.DictWriter(f,fieldnames=list(next(iter(birth.values())).keys()));w.writeheader();w.writerows(birth[h] for h in sorted(birth))
summary=dict(status='PASS',checks=checks,scope='All tournaments of orders 1 through 8; no assertion of global realizability outside observed values',by_order=stats,birth_values=len(birth),multi_orbit_birth_values=sum(r['converse_orbits']>1 for r in birth.values()),source_sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in ['census.json','n8_classes.csv','n7_representatives.txt']})
(OUT/'receipt.json').write_text(json.dumps(summary,indent=2)+'\n')
(OUT/'small_witnesses.json').write_text(json.dumps([r for r in rows if r['n']<=6 and r['H'] in [9,15,33] and birth[r['H']]['first_n']==r['n']],indent=2)+'\n')
print(json.dumps(summary,indent=2))
