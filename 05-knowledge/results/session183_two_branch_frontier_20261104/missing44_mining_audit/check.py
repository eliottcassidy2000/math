import json,hashlib,math
from pathlib import Path
D=Path(__file__).parent;B=D.parent/'missing44';S=B/'source.json'
def ck(x,msg):
 if not x:raise RuntimeError(msg)
ck(hashlib.sha256(S.read_bytes()).hexdigest()=='fd41d36a6bd6256c37d028a42babcdd269a3b4de8bbeebe680807e16678bcf7a','frozen Session182 source byte lineage')
s=json.loads(S.read_text());c=json.loads((B/'candidate44.json').read_text());v=json.loads((B/'verification.json').read_text())['candidates'][0]
ck(s['names']==c['names'],'physical name equality')
ns=c['names'];sn=c['suffix_names'];n=len(sn);N=1<<n
r=[sum(1<<j for j,b in enumerate(sn) if c['rows'][ns.index(a)]>>ns.index(b)&1) for a in sn]
def forward(rows):
 d=[[0]*n for _ in range(N)]
 for i in range(n):d[1<<i][i]=1
 for mask in range(1,N):
  for u in range(n):
   q=d[mask][u]
   if not q:continue
   avail=rows[u]&~mask
   while avail:
    bit=avail&-avail;avail-=bit;w=bit.bit_length()-1
    d[mask|bit][w]+=q
 return d
end=forward(r);start=forward([sum(1<<j for j in range(n) if r[j]>>i&1)for i in range(n)])
h=[sum(x)for x in end];h[0]=1
F=[0]*n;J=[[0]*n for _ in range(n)]
for mask in range(1,N):
 co=(N-1)^mask
 for i in range(n):
  q=end[mask][i]
  if not q:continue
  F[i]+=q*h[co]
  for j in range(n):J[i][j]+=q*start[co][j]
ck({str(sn[i]):F[i]for i in range(n)}==v['all_port_responses'],'all endpoint fields')
ck(all(J[sn.index(z['end'])][sn.index(z['start'])]==z['value']for z in v['all_joint']),'all joint fields')
W=[sn.index(x)for x in c['in_ports']];f=sum(F[i]for i in W);j=sum(J[i][k]for i in W for k in W);h0=h[-1];H=h0-j+8*f
ck((h0,f,j,H)==(40085,633897,499420,4611841),'target')
pairs={frozenset((ns[i],ns[k]))for i in range(len(ns))for k in range(i) if (s['rows'][i]>>k&1)!=(c['rows'][i]>>k&1)}
ck(pairs=={frozenset(p)for p in c['all_changed_pairs']},'exact change list')
ck(all(not(set(p)&set(c['prefix']))for p in pairs),'prefix untouched')
ck(h0-j-3*(s['h']-s['J'])==-503924 and f-3*s['F']==62991,'arbitrary prefix coefficients')
ck(math.gcd(62991,4)==1 and math.factorial(7)<62992,'isolation bounds')
rpt={'status':'PASS','source_byte_equal':True,'names':len(ns),'changed_pairs':len(pairs),'full_joint_coefficients':n*n,'all_endpoint_fields':n,'h':h0,'F':f,'J':j,'H':H,'arbitrary_prefix_difference':[-503924,62991],'isolation':'p >= (n+1)a and p <= 8a force n <= 7; congruence a=1 mod62991 with a <=5040 forces a=1; p=8 then TT3','native_scope':'not audited here','source_sha256':hashlib.sha256(S.read_bytes()).hexdigest()}
(D/'receipt.json').write_text(json.dumps(rpt,indent=2)+'\n');print(json.dumps(rpt))
