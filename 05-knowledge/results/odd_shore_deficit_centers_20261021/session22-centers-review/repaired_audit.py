import json, pathlib, copy, hashlib, math, time
P=pathlib.Path(__file__).resolve().parent
source=(P.parent/'session22-centers/check.py').read_text()
ns={'__file__':str(P.parent/'session22-centers/check.py')}
exec(compile(source.split('start=datetime.datetime')[0],'<producer definitions>','exec'),ns)
R=json.loads((P.parent/'session22-centers/receipt_normal.json').read_text()); O=json.loads((P.parent/'session22-centers/receipt.json').read_text())
volatile={'start_utc','end_utc','wall_seconds'}
assert {k:v for k,v in R.items() if k not in volatile}=={k:v for k,v in O.items() if k not in volatile}
def verify(r):
 A=[[int(x) for x in row] for row in r['source_rows']];n=len(A);C=r['cycle'];v=r['terminal']
 assert r['direction'] in ('exit','entry')
 assert all(len(x)==n for x in A)
 assert all(x in (0,1) for row in A for x in row)
 assert all(A[i][i]==0 for i in range(n))
 assert all(A[i][j]+A[j][i]==1 for i in range(n) for j in range(i))
 assert len(C)==len(set(C))==8 and set(C)<=set(range(n)) and n>=16 and v not in C and v in range(n)
 B=A if r['direction']=='exit' else list(zip(*A))
 assert all(B[C[i]][C[(i+1)%8]] for i in range(8))
 suffix=None; prefixes=set()
 for w in r['words']:
  assert len(w)==n and set(w)==set(range(n)) and w[8]==v
  assert all(B[x][y] for x,y in zip(w,w[1:]))
  assert tuple(w[:8]) in {tuple(C[i:]+C[:i]) for i in range(8)}
  if suffix is None:suffix=w[8:]
  assert suffix==w[8:]
  prefixes.add(tuple(w[:8]))
 assert len(prefixes)==len(r['words'])>=5
for r in R['prefix_receipts']:verify(r)
r=copy.deepcopy(R['prefix_receipts'][0]);n=len(r['source_rows'])
r['words']=[[(-1 if x==n-1 else x) for x in w] for w in r['words']]
r['terminal']=-1 if r['terminal']==n-1 else r['terminal']
negative_accepted=ns['audit_receipt'](r)
try:verify(r); strict_rejected=False
except AssertionError:strict_rejected=True
# Independent permutation composition and order check of support threshold.
def compose(p,q):return [p[q[i]] for i in range(len(p))]
def order(p):
 seen=set();o=1
 for i in range(len(p)):
  if i in seen:continue
  x=i;l=0
  while x not in seen:seen.add(x);l+=1;x=p[x]
  o=math.lcm(o,l)
 return o
boundary=[]
for r in (2,4,8,16):
 for n in (2*r-1,2*r,2*r+1):
  u=list(range(n));u[:r]=list(range(1,r))+[0];rho=list(reversed(range(n)))
  p=compose(u,compose(rho,compose(u,rho)))
  got=order(p);assert got==(2*r-1 if n==2*r-1 else r)
  boundary.append({'r':r,'n':n,'square_order':got})
# Boundary cut count is 2(t-1), not 2t.
assert 2*(4-1)==6 and 2*(5-1)==8
out={'source_sha256':hashlib.sha256(source.encode()).hexdigest(),'normal_checks':R['checks'],'optimized_checks':O['checks'],'nonvolatile_receipts_equal':True,'independently_valid_prefix_receipts':len(R['prefix_receipts']),'producer_accepts_negative_label_alias':negative_accepted,'strict_consumer_rejects_alias':strict_rejected,'half_support_boundaries':boundary,'entry_receipts':sum(x['direction']=='entry' for x in R['prefix_receipts']),'exit_receipts':sum(x['direction']=='exit' for x in R['prefix_receipts'])}
(P/'repaired_audit_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
