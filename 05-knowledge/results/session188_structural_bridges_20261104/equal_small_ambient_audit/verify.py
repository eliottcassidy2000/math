"""Independent exact polarization and canonical-partition cover audit."""
from pathlib import Path
from itertools import product
from fractions import Fraction
from math import factorial,comb
import json,hashlib
ROOT=Path(__file__).resolve().parent;SOURCE=ROOT.parent/'equal_small_ambient'
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def plus(*ps):
 out={}
 for p in ps:
  for e,c in p.items():out[e]=out.get(e,0)+c
 return {e:c for e,c in out.items() if c}
def times(*ps):
 out={(0,0,0):1}
 for p in ps:
  new={}
  for e,c in out.items():
   for f,d in p.items():
    g=tuple(a+b for a,b in zip(e,f));new[g]=new.get(g,0)+c*d
  out={e:c for e,c in new.items() if c}
 return out
def scale(p,c):return {e:c*v for e,v in p.items()}
ONE={(0,0,0):1};X={(1,0,0):1};Y={(0,1,0):1};Z={(0,0,1):1};P=times(X,Y,Z)
E=plus(X,Y,Z,times(X,Y),times(X,Z),times(Y,Z),scale(P,3))
K=plus(times(X,X,Y,Y),times(X,X,Z,Z),times(Y,Y,Z,Z),times(X,X,Y),times(X,X,Z),times(X,Y,Y),times(X,Z,Z),times(Y,Y,Z),times(Y,Z,Z),scale(times(P,plus(X,Y,Z)),-1),scale(P,-6))
def ttq(n):return [0]+[sum((-1)**(k-j)*comb(k,j)*j**n for j in range(k+1)) for k in range(1,n+1)]
def trilinear(num,h,qs):
 out=0
 for k in range(min(map(len,qs))):
  for ex,c in num.items():
   t=c*comb(k+h,h)
   for i,q in enumerate(qs):t*=q[k+ex[i]] if k+ex[i]<len(q) else 0
   out+=t
 return Fraction(h*out,3)
def poly(small,big,h):
 A=ttq(small)+[0]*(big-small);B=ttq(big);num=K
 for _ in range(h-1):num=times(num,E)
 out=[Fraction(0)]*4
 for choices in product(range(2),repeat=3):out[sum(choices)]+=trilinear(num,h,[A if b else B for b in choices])
 need(all(x.denominator==1 for x in out),'polarized coefficient not integral');return [int(x) for x in out],A,B,num
def profile(T):
 n=len(T);N=1<<n;D=[[0]*n for _ in range(N)];H=[0]*N;H[0]=1
 for S in range(1,N):
  rem=S
  while rem:
   bit=rem&-rem;rem-=bit;u=bit.bit_length()-1;rest=S^bit
   if not rest:D[S][u]=1;continue
   nxt=T[u]&rest
   while nxt:
    b=nxt&-nxt;nxt-=b;D[S][u]+=D[rest][b.bit_length()-1]
  H[S]=sum(D[S])
 U=[[0]*(n+1) for _ in range(N)];U[0][0]=1
 for S in range(1,N):
  least=S&-S;rest=S^least;A=rest
  while True:
   block=A|least;tail=S^block
   for k in range(1,n+1):U[S][k]+=H[block]*U[tail][k-1]
   if not A:break
   A=(A-1)&rest
 return [factorial(k)*v for k,v in enumerate(U[-1])]
def op(T):return [r^((1<<len(T))-1)^(1<<i) for i,r in enumerate(T)]
def blow(Q,blocks):
 offs=[sum(map(len,blocks[:i])) for i in range(len(blocks))];out=[0]*sum(map(len,blocks))
 for i,B in enumerate(blocks):
  for u,r in enumerate(B):
   out[offs[i]+u]=r<<offs[i]
   for j,C in enumerate(blocks):
    if Q[i]>>j&1:out[offs[i]+u]|=((1<<len(C))-1)<<offs[j]
 return out
def extend(T,a,c):
 n=len(T);R=list(T)+[1<<(n+1),0]
 for u in range(n):
  if a>>u&1:R[u]|=1<<n
  else:R[n]|=1<<u
  if c>>u&1:R[u]|=1<<(n+1)
  else:R[n+1]|=1<<u
 return R
def sos(q):
 a=q+[0,0];n=len(q)-1;z=a[1]**3
 for j in range(1,n+2):z+=j*(a[j-1]-a[j+1])**2*(a[j-1]+3*a[j]+2*a[j+1])
 for j in range(1,n+1):z+=(a[j]-a[j+1])**2*((j+2)*a[j]+(2*j+1)*a[j+1])
 need(z%3==0 and z>0,'strict square certificate');return z//3

def main():
 import tempfile,subprocess
 data=json.loads((SOURCE/'normal.json').read_text());coordinates=[]
 for c in data['coordinate_controls']:
  T=c['child_rows'];q=profile(T);n=len(T);need(q==c['q'],'independent child ordered-cover profile');num=K;gaps=[]
  for h in range(1,3*n+3):
   num={e:v for e,v in num.items() if max(e)<=n};d=trilinear(num,h,[q]*3);need(d.denominator==1,'integer cover gap');gaps.append(int(d));num=times(num,E)
  P=sos(q);need(gaps==c['all_augmented_cover_gaps'] and gaps[0]==-P==-c['P'] and all(d<=0 for d in gaps),'complete coordinate domination and scalar strictness')
  coordinates.append({'name':c['name'],'profile':q,'gaps':gaps,'positive_scalar_certificate':P})
 jobs=[];meta=[]
 for r in data['literal_controls']['records']:
  A=r['child_rows'];Q=r['outer_rows'];n=len(A);M=blow([2,4,1],[A,A,A]);N=blow([4,1,2],[op(A),A,op(A)]);T=blow(Q,[M]+[[0] for _ in Q[1:]]);U=blow(Q,[N]+[[0] for _ in Q[1:]])
  need(T==r['source_old_rows'] and U==r['target_old_rows'],'fixed physical old labels and target interiors');outside=sum(1<<(3*n+j-1) for j in range(1,len(Q)) if Q[j]&1);a=((1<<n)-1)|outside;c=(((1<<n)-1)<<n)|outside;need((a,c)==(r['a'],r['c']),'physical exterior port lift')
  X,Y=extend(T,a,c),extend(U,a,c);augmented=((1<<(3*n))-1)|(1<<len(T))|(1<<(len(T)+1))
  need(all(Z[v]&augmented in(0,augmented) for Z in(X,Y) for v in range(3*n,len(T))),'common actual augmented module')
  jobs.extend([X,Y]);P=sos(profile(A));H=profile(Q)[1];need(r['lower_bound']==H*P,'exact quantitative bound');meta.append((r,H*P))
 lines=[str(len(jobs))]+[str(len(T))+' '+' '.join(map(str,T)) for T in jobs]
 with tempfile.TemporaryDirectory() as d:
  exe=Path(d)/'count';subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(ROOT/'count.cpp'),'-o',str(exe)],check=True);counts=list(map(int,subprocess.check_output([str(exe)],input=('\n'.join(lines)+'\n').encode()).split()))
 need(len(counts)==72,'all independent physical counts');physical=[]
 for i,(r,bound) in enumerate(meta):
  f,g=counts[2*i:2*i+2];need((f,g)==(r['F_source'],r['F_target']) and g-f==r['target_minus_source']>=bound>0,'literal strict ambient sign');physical.append({'child_rows':r['child_rows'],'outer_rows':r['outer_rows'],'physical_query':[r['a'],r['c']],'values':[f,g],'positive_lower_bound':bound,'augmented_order':len(jobs[2*i])})
 return {'status':'PASS','coordinates':coordinates,'independent_physical_counts':72,'physical_controls':physical,'source_pins':{f:hashlib.sha256((SOURCE/f).read_bytes()).hexdigest() for f in ['PROOF.md','verify.py','literal.cpp','normal.json','optimized.json']},'scope':'Actual common child profiles of order1..7, fixed cyclic quotient converse and arbitrary fixed ambient with the prescribed common-module port lift; not arbitrary profile or quotient reconstruction.'}
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
