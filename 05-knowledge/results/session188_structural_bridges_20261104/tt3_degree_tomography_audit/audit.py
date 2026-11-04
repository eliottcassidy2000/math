from pathlib import Path
from itertools import combinations,permutations,product
from math import comb
import hashlib,json,subprocess,tempfile
ROOT=Path(__file__).resolve().parent;SRC=ROOT.parent/'tt3_degree_tomography'
def need(x,m):
 if not x:raise RuntimeError(m)
def rows(n,code):
 R=[0]*n
 for k,(i,j)in enumerate(combinations(range(n),2)):
  if code>>k&1:R[i]|=1<<j
  else:R[j]|=1<<i
 return R
def op(R):return [((1<<len(R))-1)^r^(1<<i)for i,r in enumerate(R)]
def tt(n):return [sum(1<<j for j in range(i+1,n))for i in range(n)]
def tq(n,kmax):return [0]+[sum((-1)**(k-j)*comb(k,j)*j**n for j in range(k+1))for k in range(1,kmax+1)]
def blow(Q,B):
 off=[];n=0
 for b in B:off.append(n);n+=len(b)
 R=[0]*n
 for i,b in enumerate(B):
  for u in range(len(b)):
   for j,c in enumerate(B):
    for v in range(len(c)):
     if(i==j and b[u]>>v&1)or(i!=j and Q[i]>>j&1):R[off[i]+u]|=1<<(off[j]+v)
 return R,off
def aug(R,a,c):
 n=len(R);out=list(R)+[1<<(n+1),0]
 for i in range(n):
  for v,M in [(n,a),(n+1,c)]:
   if M>>i&1:out[i]|=1<<v
   else:out[v]|=1<<i
 return out
def outer(R,a,c):
 V,o=blow([2,4,1],[R,[0],[0]]);return aug(V,a|1<<o[2],c|1<<o[2])
def mul(A,B):
 N=len(A);out=[0]*N
 for a,x in enumerate(A):
  for b,y in enumerate(B):
   if not a&b:out[a|b]+=x*y
 return out
def inputs(R,out):
 n=len(R);N=1<<n;S=[[0]*N for _ in R];E=[[0]*N for _ in R];tau=[0]*N;ell=[[0]*N for _ in R];rho=[[0]*N for _ in R]
 for k in range(1,n+1):
  for w in permutations(range(n),k):
   if not all(R[u]>>v&1 for u,v in zip(w,w[1:])):continue
   M=sum(1<<u for u in w);u,v=w[0],w[-1];S[u][M]+=1;E[v][M]+=1
   if out>>u&1 and not out>>v&1:tau[M]+=1
   if not out>>v&1:ell[u][M]+=1
   if out>>u&1:rho[v][M]+=1
 al=[int(M==0)+sum(E[v][M]for v in range(n)if not out>>v&1)for M in range(N)];be=[int(M==0)+sum(S[v][M]for v in range(n)if out>>v&1)for M in range(N)]
 return S,E,al,be,tau,ell,rho
def setup(E,ao,mo):
 Q=[2,4,1]+[v<<3 for v in E]
 for i in range(len(E)):
  if ao>>i&1:Q[0]|=1<<(i+3)
  else:Q[i+3]|=1
  if mo>>i&1:Q[1]|=1<<(i+3);Q[2]|=1<<(i+3)
  else:Q[i+3]|=6
 return Q
def h_literal(R):return sum(all(R[u]>>v&1 for u,v in zip(w,w[1:]))for w in permutations(range(len(R))))
def main():
 src=json.loads((SRC/'normal.json').read_text());need((SRC/'normal.json').read_bytes()==(SRC/'optimized.json').read_bytes(),'paired source receipts');need(src['proof_sha256']==hashlib.sha256((SRC/'PROOF.md').read_bytes()).hexdigest()and src['verifier_sha256']==hashlib.sha256((SRC/'verify.py').read_bytes()).hexdigest(),'frozen source pins')
 graphs=[];where={}
 def add(R):
  key=tuple(R)
  if key not in where:where[key]=len(graphs);graphs.append(R)
  return where[key]
 tom=[];cache={}
 for n,co,out,ac,r,i,expected in src['tomography']:
  R=rows(n,co);A=rows(3,ac)
  T=[v|sum(1<<(3+j)for j in range(n)if out>>j&1)for v in A]+[(v<<3)|(0 if out>>j&1 else 7)for j,v in enumerate(R)]
  key=n,co,out
  if key not in cache:cache[key]=inputs(R,out)
  S,E,al,be,tau,ell,rho=cache[key];d=sum(A[v]>>r&1 for v in range(3));as_=mul(al,S[i]);be_=mul(be,E[i]);rhs=d*as_[-1]+(2-d)*be_[-1]+2*mul(as_,tau)[-1]+2*mul(be_,tau)[-1]+2*mul(mul(al,be),[x+y for x,y in zip(ell[i],rho[i])])[-1]
  need(rhs==expected,'independent exterior literal-path formula');p,q=[v for v in range(3)if v!=r];roles=[('a',p),('a',q),('c',r),('c',3+i)];ids=[]
  for bits in range(16):
   a=sum(1<<v for j,(t,v)in enumerate(roles)if t=='a'and bits>>j&1);c=sum(1<<v for j,(t,v)in enumerate(roles)if t=='c'and bits>>j&1);ids.append(add(aug(T,a,c)))
  tom.append((ids,expected))
 need(len(tom)==5040,'complete tomography domain')
 physical=src['physical'];fullids=[[add(aug(R,a,c))for a,c in product(range(32),repeat=2)]for R in [physical['T'],physical['U']]]
 families=[]
 for side_size in [1,2]:
  B=[tt(3),tt(side_size),tt(side_size)];T,off=blow([2,4,1],B);n=len(T);X=((1<<side_size)-1)<<off[1];Y=((1<<side_size)-1)<<off[2];qs={(0,0)}
  for r in range(3):
   p,q=[i for i in range(3)if i!=r]
   for Z in [X,Y]:
    for bits in range(16):qs.add(((1<<p if bits&1 else 0)|(1<<q if bits&2 else 0),(1<<r if bits&4 else 0)|(Z if bits&8 else 0)))
  blocks=[1,2,4,X,Y];anchor=0;out_blocks=[C for C in blocks if C!=1 and T[0]&C==C];d=len(out_blocks)
  for C in blocks[1:]:qs.add((C|1 if C in out_blocks else C,1))
  if d%2:qs.add((1,1))
  if(4+d)%2:qs.add((0,1))
  qs.add((1,0));qs=sorted(qs);modules=[3,6,7]+([X,Y]if side_size==2 else[])
  for outerflag in [False,True]:
   ext=outer if outerflag else aug;sourceids=[add(ext(T,a,c))for a,c in qs]
   for bits in range(1<<len(modules)):
    V=op(T)
    for j,M in enumerate(modules):
     if bits>>j&1:
      for u,v in combinations(range(n),2):
       if M>>u&1 and M>>v&1:V[u]^=1<<v;V[v]^=1<<u
    targetids=[add(ext(V,a,c))for a,c in qs];families.append({'side_size':side_size,'ambient':outerflag,'family_bits':bits,'source_ids':sourceids,'target_ids':targetids,'query_count':len(qs),'queries':qs})
 with tempfile.TemporaryDirectory()as td:
  ex=Path(td)/'count';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(ROOT/'count.cpp'),'-o',str(ex)],check=True);wire=str(len(graphs))+'\n'+''.join(str(len(R))+' '+' '.join(map(str,R))+'\n'for R in graphs);values=list(map(int,subprocess.check_output([str(ex)],input=wire.encode()).split()));need(len(values)==len(graphs),'all literal graphs')
  for ids,w in tom:need(sum((-1)**(4-b.bit_count())*values[k]for b,k in enumerate(ids))==w,'every literal physical quartic')
  need([values[i]for i in fullids[0]]==physical['source_values']and[values[i]for i in fullids[1]]==physical['target_values'],'complete actual response banks')
  half=[a for a in range(32)if a.bit_count()<=2];coeff=[]
  def subsets(M):
   S=M
   while True:
    yield S
    if not S:break
    S=(S-1)&M
  for a,c in product(half,repeat=2):
   z=[sum((-1)**(a.bit_count()+c.bit_count()-x.bit_count()-y.bit_count())*values[fullids[side][32*x+y]]for x in subsets(a)for y in subsets(c))for side in range(2)];coeff.append([a,c,*z])
  need(coeff==physical['complete_named_bank'],'all256 named coefficients including zeros')
  for a,c in product(range(32),repeat=2):
   need([sum(v[side]for x,y,*v in coeff if x&a==x and y&c==y)for side in range(2)]==[values[fullids[side][32*a+c]]for side in range(2)],'complete source/target full-mask reconstruction')
  for f in families:
   sv=[values[i]for i in f.pop('source_ids')];tv=[values[i]for i in f.pop('target_ids')];need(any(a!=b for a,b in zip(sv,tv)),'fixed root-family battery');f['source_values']=sv;f['target_values']=tv
  configs=[(n,co,ao,mo)for n in range(4)for co in range(1<<(n*(n-1)//2))for ao in range(1<<n)for mo in range(1<<n)];wire=str(len(configs))+'\n'+''.join(str(n+3)+' '+' '.join(map(str,setup(rows(n,co),ao,mo)))+'\n'for n,co,ao,mo in configs);ex=Path(td)/'run';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(ROOT/'run_covers.cpp'),'-o',str(ex)],check=True);raw=[list(map(int,s.split()))for s in subprocess.check_output([str(ex)],input=wire.encode()).decode().splitlines()]
 cursor=0;lookup={}
 for config in configs:
  n,co,ao,mo=config;count=(n+3)**2;lookup[config]=raw[cursor:cursor+count];cursor+=count
 need(cursor==len(raw),'complete direct ordered-cover run tables')
 for n,co,ao,mo,qi,kx,ky,lo in src['equal_profile']:
  q=tq(qi+1,n+3)if qi<8 else[0,3,6,6]+[0]*n
  vx=sum((a-b)*q[i]*q[j]for i,j,a,b,c,d in lookup[n,co,ao,mo]);vy=sum((c-d)*q[i]*q[j]for i,j,a,b,c,d in lookup[n,co,ao,mo]);E=rows(n,co);EM=[sum(1<<(i+1)for i in range(n)if mo>>i&1)]+[(r<<1)|(0 if mo>>i&1 else 1)for i,r in enumerate(E)]
  need([vx,vy,2*q[1]**2*h_literal(EM)]==[kx,ky,lo]and vx-vy>=lo>0,'independent cover-run skew and lower bound')
 qprofiles=[tq(n,n)for n in range(1,21)]+[[0,3,6,6]]+[[0]+[comb(n,k)for k in range(n+1)]for n in range(1,16)]
 for qi,c,m,w in src['binomial']:
  q=qprofiles[qi];val=lambda i:q[i]if 0<=i<len(q)else 0;p={(0,0):1}
  for _ in range(m):
   z={}
   for(i,j),a in p.items():z[i+1,j]=z.get((i+1,j),0)+a;z[i,j+1]=z.get((i,j+1),0)+a
   p=z
  result=sum(a*(2*val(i+c+1)*val(j+c+1)-val(i+c+2)*val(j+c)-val(i+c)*val(j+c+2))for(i,j),a in p.items());need(result==w and result>=0,'independent binomial polynomial contraction')
 need(len(src['equal_profile'])==4941 and len(src['binomial'])==3456,'full numerical scopes')
 return {'status':'PASS','source_proof':src['proof_sha256'],'source_verifier':src['verifier_sha256'],'auditor_pins':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest()for f in ['audit.py','count.cpp','run_covers.cpp']},'literal_quartic_coefficients':len(tom),'distinct_actual_graphs':len(graphs),'literal_full_response_values':2048,'independent_cover_run_configurations':len(configs),'equal_profile_controls':4941,'binomial_controls':3456,'root_family_controls':families,'scope':'Analytical identities, equal-log-concave positivity and source-relative root-family integration; no arbitrary-prime completeness.'}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
