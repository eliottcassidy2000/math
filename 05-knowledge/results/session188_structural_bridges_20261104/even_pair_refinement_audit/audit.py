"""Independent endpoint identities, literal quotient finite differences and physical witness."""
from pathlib import Path
from itertools import combinations
import json,hashlib,subprocess,tempfile,sys
HERE=Path(__file__).resolve().parent;BASE=HERE.parent;SOURCE=BASE/'even_pair_refinement'
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def rows(n,code):
 T=[0]*n
 for k,(i,j) in enumerate(combinations(range(n),2)):
  T[i if code>>k&1 else j]|=1<<(j if code>>k&1 else i)
 return tuple(T)
def op(T):return tuple(((1<<len(T))-1)^r^(1<<i) for i,r in enumerate(T))
def ends(T):
 n=len(T);d=[[0]*n for _ in range(1<<n)]
 for i in range(n):d[1<<i][i]=1
 for S in range(1,1<<n):
  for i,z in enumerate(d[S]):
   if not z:continue
   avail=T[i]&~S
   while avail:
    b=avail&-avail;avail-=b;d[S|b][b.bit_length()-1]+=z
 return d
def masks(D):return [sum((v&1)<<i for i,v in enumerate(row)) for row in D]
def attach(T,a,c):
 n=len(T);R=list(T)+[1<<(n+1),0]
 for i in range(n):
  for v,m in ((n,a),(n+1,c)):
   if m>>i&1:R[i]|=1<<v
   else:R[v]|=1<<i
 return R
def hp(T):return sum(ends(T)[-1])
def features(T):
 n=len(T);E,S=masks(ends(T)),masks(ends(op(T)));b=0;h=[0]*n
 for M in range(1,1<<n):
  v=E[M]^S[M];k=M.bit_count()
  if k>=2 and n-k>=2 and (n-k)%2==0:b^=v
  if k%2==0:
   for j in range(n):
    if not M>>j&1:h[j]^=v
 return b,h

def coefficient(bank,m):
 total=0;s=m
 while True:
  total+=(-1)**((m^s).bit_count())*bank[s]
  if not s:return total
  s=(s-1)&m

def main():
 paths=[SOURCE/x for x in ('PROOF.md','verify.py','run_dp.cpp','literal_hostile.cpp')]+[BASE/d/'PROOF.md' for d in ('odd_pair_quotient','mod8_transfer')];pins={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
 source=json.loads((SOURCE/'normal.json').read_text());run=json.loads((SOURCE/'run_dp.json').read_text());literal_controls=0;quotient_queries=0;shared=0
 for row in run['small_trace_controls']:
  n,code=row[:2];T=rows(n,code);b,h=features(T)
  need(row[2:2+n]==[4*(1^((b>>j)&1)) for j in range(n)],'independent integer subset detector versus physical trace')
  need(row[2+n:]==[0 if i==j else 4*((h[j]>>i)&1) for j in range(n) for i in range(n)],'independent physical shared cubics')
  # Recover delta_Q itself from complete actual integer attachment banks.
  N=1<<n;bank=[]
  for c in range(N):
   for a in range(N):
    delta=hp(attach(T,a,c))-hp(attach(op(T),a,c));need(delta%2==0,'integer half contrast');bank.append(delta//2);quotient_queries+=2
  for j in range(n):
   trace=0
   for i in range(n):
    if i==j:continue
    hh=coefficient(bank,(1<<i)|(1<<j)|(1<<(n+j)))&1;need(hh==(h[j]>>i&1),'independent shared cubic cycle identity');shared+=1
    trace+=coefficient(bank,(1<<i)|(1<<(n+j)))
    if T[j]>>i&1:trace+=hh
   need(trace%2==1,'actual E1 anchor identity')
  literal_controls+=1
 Q=tuple(source['hostile_rows']);b,h=features(Q);need(b==763 and sum(x.bit_count() for x in h)==44 and h[0]>>3&1,'hostile exact independent detector')
 # Direct integer endpoint sums over every eligible ordered disjoint support.
 E,S=ends(Q),ends(op(Q));N=(1<<len(Q))-1;direct=0;pairs=0
 for A in range(1,N+1):
  if A.bit_count()<2:continue
  ha=sum(E[A]);B=N^A
  while B:
   if B.bit_count()>=2:
    direct^=sum(((ha*(E[B][j]+S[B][j]))&1)<<j for j in range(len(Q)));pairs+=1
   B=(B-1)&(N^A)
 need(direct==b,'hostile independent support-weight correction')
 with tempfile.TemporaryDirectory() as d:
  binary=Path(d)/'literal';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(HERE/'literal.cpp'),'-o',str(binary)],check=True)
  physical=json.loads(subprocess.check_output([str(binary)]))
 need(physical['values_mod8']==[3,7,3,7],'independent forward physical witness and pair-bit controls')
 normal=subprocess.check_output([sys.executable,str(SOURCE/'verify.py')]);optimized=subprocess.check_output([sys.executable,'-O',str(SOURCE/'verify.py')]);need(normal==optimized==(SOURCE/'normal.json').read_bytes(),'fresh producer modes replay')
 need(pins=={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},'source freeze')
 print(json.dumps({'status':'PASS','actual_quotient_attachment_values':quotient_queries,'complete_quotient_trace_controls':literal_controls,'independent_shared_cubic_coefficients':shared,'hostile_disjoint_support_pairs':pairs,'hostile_b_mask':b,'hostile_shared_cubic_columns':h,'physical_forward_dp':physical,'producer_modes_equal':True,'source_sha256':pins,'audit_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (HERE/'audit.py',HERE/'literal.cpp')}},indent=2,sort_keys=True))
if __name__=='__main__':main()
