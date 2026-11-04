"""Independent parity profiles, Hamilton cycles and actual ambient checks."""
from pathlib import Path
from itertools import permutations,combinations
from functools import lru_cache
import json,hashlib,subprocess,sys
HERE=Path(__file__).resolve().parent;S=HERE.parent/'mod8_transfer'
def need(x,msg):
 if not x:raise RuntimeError(msg)
@lru_cache(None)
def Hs(R):
 n=len(R);dp=[{} for _ in range(1<<n)]
 for u in range(n):dp[1<<u][u]=1
 for used in range(1,1<<n):
  for u,z in dp[used].items():
   rem=R[u]&~used
   while rem:
    bit=rem&-rem;rem-=bit;v=bit.bit_length()-1;dest=dp[used|bit];dest[v]=dest.get(v,0)+z
 return (1,)+tuple(sum(dp[M].values()) for M in range(1,1<<n))
def H(R):return Hs(tuple(R))[-1]
def prof(R):
 R=tuple(R);n=len(R);h=Hs(R);p=[[0]*(n+1) for _ in h];p[0][0]=1
 for M in range(1,len(h)):
  A=M
  while A:
   for k in range(1,n+1):p[M][k]+=h[A]*p[M^A][k-1]
   A=(A-1)&M
 return p[-1]
def cycles(R):
 n=len(R)
 if n<3:return 0
 return sum(all(R[u]>>v&1 for u,v in zip((0,)+p,p+(0,))) for p in permutations(range(1,n)))
def root(Q,r,A):
 outside=[v for v in range(len(Q)) if v!=r];R=list(A)+[0]*len(outside);n=len(A)
 for i,v in enumerate(outside):
  for j,w in enumerate(outside):
   if Q[v]>>w&1:R[n+i]|=1<<(n+j)
  if Q[v]>>r&1:R[n+i]|=(1<<n)-1
  else:
   for u in range(n):R[u]|=1<<(n+i)
 return tuple(R)
def rows(n,code):
 R=[0]*n
 # Match the record convention: unset bit points from lower to higher name.
 for k,(u,v) in enumerate(combinations(range(n),2)):
  if code>>k&1:R[v]|=1<<u
  else:R[u]|=1<<v
 return tuple(R)
def main():
 dependencies=[S/'PROOF.md',S/'verify.py',HERE.parent/'covered_response'/'verify.py',HERE.parent/'response'/'verify.py'];pins={str(p.relative_to(HERE.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in dependencies}
 record=json.loads((S/'receipt.json').read_text());profiles=literalcycles=0
 for rec in record['profile_records']:
  if rec['order']>5:continue
  R=tuple(rec['rows']);p=prof(R);need(p==rec['q'],'independent full profile');z=cycles(R) if len(R)%2 else 0;need(z==rec['hamilton_cycles_if_odd'],'literal rotational Hamilton cycles');profiles+=1;literalcycles+=len(R)%2
  q=p+[0,0,0];N=len(R);need((q[2]-2*q[1]-(1<<N)+4+4*z)%8==0,'q2 including odd correction')
  s=[[0]*4 for _ in range(N+1)];s[0][0]=1
  for n in range(1,N+1):
   for k in range(1,4):s[n][k]=k*s[n-1][k]+s[n-1][k-1]
  need((q[3]-6*s[N][3]-(2*q[1]-2-4*z if N%2 else 0))%8==0,'q3 including odd correction')
  need(all(x%8==0 for x in q[4:]),'higher profile residues')
 controls=[]
 for rec in record['cyclic_port_core_controls']:
  X,Y=(tuple(rec[k]['rows']) for k in ('port_core_T','port_core_U'));px,py=prof(X),prof(Y);need(px==rec['port_core_T']['q'] and py==rec['port_core_U']['q'],'complete actual port profiles')
  N=len(X);z=cycles(X)-cycles(Y) if N%2 else 0;d=px[1]-py[1];actual=H(rec['ambient_T_rows'])-H(rec['ambient_U_rows']);need(actual==rec['ambient_difference'],'actual ambient core counts')
  need(actual==3*d+px[2]-py[2],'cyclic run coefficients')
  if N%2:need(z%2==1 and d%8==4 and actual%8==0 and actual!=0,'odd obstruction not exact collision')
  else:need(d%8==actual%8==4,'even residue four persists')
  controls.append({'order':N,'core_gap':d,'ambient_gap':actual,'odd_cycle_gap':z})
 pairs=[(next(rows(4,c) for c in range(64) if H(rows(4,c))==3),rows(4,0)),((8,1,3,6),rows(4,0))];ambients=0
 for n,code,r,d,factor,actual in record['ambient_records']:
  Q=rows(n,code);B,A=next((B,A) for B,A in pairs if H(B)-H(A)==d)
  need(H(root(Q,r,(2,0)))==factor and factor%2==1,'literal odd unit H(Q[TT2])')
  gap=H(root(Q,r,B))-H(root(Q,r,A));need(gap==actual and (gap-factor*d)%8==0,'literal even ambient replacement');ambients+=1
 normal=subprocess.check_output([sys.executable,str(S/'verify.py')]);optimized=subprocess.check_output([sys.executable,'-O',str(S/'verify.py')]);need(normal==optimized and json.loads(normal)==record,'producer saved/normal/optimized agreement')
 need(pins=={str(p.relative_to(HERE.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in dependencies},'source and helper freeze')
 out={'status':'PASS','independent_profiles_through5':profiles,'literal_odd_order_cycle_counts':literalcycles,'independent_even_ambient_substitutions':ambients,'physical_port_controls':controls,'producer_receipt_sha256':hashlib.sha256(normal).hexdigest(),'producer_normal_optimized_equal':True,'source_sha256':pins,'audit_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 (HERE/'receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,sort_keys=True))
if __name__=='__main__':main()
