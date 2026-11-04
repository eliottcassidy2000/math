from pathlib import Path
from itertools import permutations
import json,hashlib,subprocess,tempfile
HERE=Path(__file__).resolve().parent;SRC=HERE.parent/'tt3_circulation'
def need(c,m):
 if not c:raise RuntimeError(m)
def aug(T,a,c):
 n=len(T);R=list(T)+[1<<(n+1),0]
 for v in range(n):
  for bits,p in[(a,n),(c,n+1)]:
   if bits>>v&1:R[v]|=1<<p
   else:R[p]|=1<<v
 return tuple(R)
def flip(T):return [r^((7^(1<<i)) if i<3 else 0) for i,r in enumerate(T)]
def add(*ps):
 out={}
 for P in ps:
  for k,v in P.items():out[k]=out.get(k,0)+v
 return {k:v for k,v in out.items() if v}
def scale(c,P):return {k:c*v for k,v in P.items() if c*v}
def mul(A,B):
 out={}
 for (m,t),v in A.items():
  for (s,u),w in B.items():
   if not m&s:out[m|s,t+u]=out.get((m|s,t+u),0)+v*w
 return {k:v for k,v in out.items() if v}
def matmul(A,B):return [[add(*(mul(A[i][k],B[k][j]) for k in range(3))) for j in range(3)] for i in range(3)]
def formal():
 records=[]
 for T,pairs in [([6,4,0],[(0,1),(1,2),(0,2)]),([4,1,2],[(0,2)])]:
  Z=[[{(1<<i,0):1} if i==j else{} for j in range(3)] for i in range(3)];A=[[{(1<<i,0):int(bool(T[i]>>j&1)),(1<<i,1):1} for j in range(3)] for i in range(3)];AZ=matmul(A,Z);AAZ=matmul(A,AZ);need(all(not x for row in matmul(A,AAZ) for x in row),'nilpotent matrix tail');M=[[add(Z[i][j],AZ[i][j],AAZ[i][j]) for j in range(3)] for i in range(3)];C=[add(*(M[i][j] for i in range(3))) for j in range(3)];R=[add(*row) for row in M];H=add(*C)
  for p,q in pairs:
   S=add(mul(C[p],C[q]),scale(-1,mul(R[p],R[q])));D=add(mul(C[p],R[q]),scale(-1,mul(C[q],R[p])));E=add(M[p][q],scale(-1,M[q][p]));HE=mul(H,E);need(not mul(H,S) and not mul(H,D) and not mul(H,HE),'all omitted terms vanish identically');records.append({'rows':T,'p':p,'q':q,'coefficients':[{str(t):v for (m,t),v in P.items() if m==7} for P in [S,D,E,HE]]})
 return records

def extform(R,ao,i,j):
 n=len(R);G=[[{} for _ in R] for _ in R]
 for k in range(1,n+1):
  for w in permutations(range(n),k):
   if all(R[u]>>v&1 for u,v in zip(w,w[1:])):
    key=(sum(1<<v for v in w),0);D=G[w[0]][w[-1]];D[key]=D.get(key,0)+1
 E=[add(*(G[u][v] for u in range(n))) for v in range(n)];S=[add(*row) for row in G];O=[v for v in range(n) if ao>>v&1];I=[v for v in range(n) if not ao>>v&1];one={(0,0):1};alpha=add(one,*(E[v] for v in I));beta=add(one,*(S[v] for v in O));tau=add(*(G[u][v] for u in O for v in I));ell=[add(*(G[u][v] for v in I)) for u in range(n)];rho=[add(*(G[v][u] for v in O)) for u in range(n)]
 KS=add(mul(alpha,add(mul(ell[i],S[j]),mul(ell[j],S[i]))),scale(-1,mul(beta,add(mul(E[i],rho[j]),mul(E[j],rho[i])))));KD=mul(mul(alpha,beta),add(G[i][j],scale(-1,G[j][i])));B0=add(mul(E[i],S[j]),scale(-1,mul(E[j],S[i])));B1=add(mul(alpha,add(mul(rho[i],S[j]),scale(-1,mul(rho[j],S[i])))),mul(beta,add(mul(E[i],ell[j]),scale(-1,mul(E[j],ell[i])))));z=mul(add(one,scale(2,tau)),B0);P=[add(scale(-2,KS),B1),add(scale(2,KS),B1),add(scale(-2,KD),z,B1),add(scale(2,KD),scale(-1,z),B1)];return [p.get(((1<<n)-1,0),0) for p in P]
def H(T):
 n=len(T);D=[[0]*n for _ in range(1<<n)]
 for M in range(1,1<<n):
  for v in range(n):
   if M>>v&1:
    S=M^(1<<v);D[M][v]=1 if not S else sum(D[S][u] for u in range(n) if S>>u&1 and T[u]>>v&1)
 return sum(D[-1])
def induced(T,vs):return [sum(1<<j for j,v in enumerate(vs) if T[u]>>v&1) for u in vs]
def main():
 src=json.loads((SRC/'normal.json').read_text());ft=formal();need(ft==src['formal_effective_path_table'],'independent Neumann matrix table including tau');expected={};records=[]
 def take(T,U,qs):
  for q in qs:
   for side,R in enumerate([T,U]):
    g=aug(R,q['a'],q['c']);v=q['values'][side];need(g not in expected or expected[g]==v,'consistent duplicate physical values');expected[g]=v
  return sum(q['sign']*(q['values'][0]-q['values'][1]) for q in qs)
 for r in src['universal_controls']:
  need([x&7 for x in r['T'][:3]]==[6,4,0] and [x&7 for x in r['V'][:3]]==[4,1,2],'physical TT3 and virtual cycle orientation');need(all((r['T'][k]^r['V'][k])&~7==0 for k in range(3)) and r['T'][3:]==r['V'][3:],'only internal root differences')
  ds=[take(r['T'],flip(r['T']),qs) for qs in r['tt_quartics']];cy=take(r['V'],flip(r['V']),r['cycle_quartic']);circ=take(r['T'],flip(r['T']),r['circulation']);formula=extform(r['R'],r['Aout'],r['i'],r['j']);need(ds+[cy]==formula and ds[0]+ds[1]-ds[2]==cy==circ,'all arbitrary-exterior exact quartics and circulation');records.append(formula)
 bounds=[]
 for r in src['positive_witness_controls']:
  local=take(r['T'],flip(r['T']),r['local']);quot=take(r['T'],r['U'],r['quotient']);cy=take(r['V'],flip(r['V']),r['cycle']);E=r['E'];ao=r['Aout'];R=[ao<<1]+[(v<<1)|(0 if ao>>i&1 else 1) for i,v in enumerate(E)];bound=2*H(r['blocks'][1])*H(r['blocks'][2])*H(R);need(local==quot==cy==r['coefficient'] and local>=bound>0,'all positive grouped witnesses and independent lower bounds');bounds.append([local,bound])
 amb=[]
 for r in src['strong_ambient_controls']:
  local=take(r['T'],flip(r['T']),r['local']);quot=take(r['T'],r['U'],r['quotient']);cy=take(r['V'],flip(r['V']),r['cycle']);vs=[0]+list(range(6,len(r['T'])));bound=2*H(induced(r['T'],[3,4]))*H(induced(r['T'],[5]))*H(induced(r['T'],vs));need(local==quot==cy==r['coefficient'] and local>=bound>0,'strong ambient circulation and enlarged exterior bound');amb.append([local,bound])
 bank=src['complete_named_bank'];T,U=bank['T'],bank['U'];lookup={}
 for q in bank['coordinates']:
  take(T,U,[dict(q,sign=1)]);lookup[q['a'],q['c']]=q['values']
 def sub(a):
  x=a
  while True:
   yield x
   if not x:break
   x=(x-1)&a
 coeff={};hist={d:0 for d in range(5)}
 for q in bank['coordinates']:
  a,c=q['a'],q['c'];vv=[sum((-1)**((a^x).bit_count()+(c^y).bit_count())*lookup[x,y][s] for x in sub(a) for y in sub(c)) for s in [0,1]];need(vv==q['coefficients'],'every841 complete named coefficient');coeff[a,c]=vv;hist[a.bit_count()+c.bit_count()]+=vv[0]!=vv[1]
 for q in bank['full_mask_controls']:
  take(T,U,[dict(q,sign=1)]);a,c=q['a'],q['c'];need([sum(v[s] for (x,y),v in coeff.items() if x&a==x and y&c==y) for s in [0,1]]==q['values'],'full named polynomial reconstruction')
 need({str(k):v for k,v in hist.items()}==bank['nonzero_contrast_by_degree'],'complete cubic/quartic histograms');graphs=list(expected);need(len(graphs)==src['distinct_actual_graphs']==43434,'every retained actual graph independently reconstructed');lines=[]
 for g in graphs:lines.extend([str(len(g)),' '.join(map(str,g))])
 with tempfile.TemporaryDirectory() as td:
  exe=Path(td)/'physical';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(HERE/'physical.cpp'),'-o',str(exe)],check=True);actual=list(map(int,subprocess.check_output([str(exe)],input=('\n'.join(lines)+'\n').encode()).split()))
 need(actual==[expected[g] for g in graphs],'all43434 physical values by independent predecessor DP');need((SRC/'normal.json').read_bytes()==(SRC/'optimized.json').read_bytes(),'paired producer receipts');pins={f:hashlib.sha256((SRC/f).read_bytes()).hexdigest() for f in ['PROOF.md','verify.py','positive.py','kernel.py','count.cpp','normal.json','optimized.json']}
 print(json.dumps({'status':'PASS','source_pins':pins,'formal_table':ft,'universal_formula_controls':len(records),'universal_formulas_sha256':hashlib.sha256(repr(records).encode()).hexdigest(),'positive_witness_bounds':bounds,'strong_ambient_bounds':amb,'complete_named_coordinates':len(coeff),'nonzero_named_difference_histogram':hist,'distinct_actual_graphs':len(graphs),'maximum_physical_order':max(map(len,graphs)),'full_physical_value_vector_sha256':hashlib.sha256(repr(actual).encode()).hexdigest(),'audit_pins':{f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['audit.py','physical.cpp']}},sort_keys=True,indent=2))
if __name__=='__main__':main()
