"""Independent exact polarization and canonical-partition cover audit."""
from pathlib import Path
from itertools import product
from fractions import Fraction
from math import factorial,comb
import json,hashlib
ROOT=Path(__file__).resolve().parent;SOURCE=ROOT.parent/'profile_cone_obstruction'
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
def tt(n):return [((1<<n)-1)^((1<<(i+1))-1) for i in range(n)]
def build(sizes,reverse):
 offs=[sum(sizes[:i]) for i in range(3)];n=sum(sizes);T=[0]*n
 for i,m in enumerate(sizes):
  j=(i+(2 if reverse else 1))%3
  for u,r in enumerate(tt(m)):T[offs[i]+u]=(r<<offs[i])|(((1<<sizes[j])-1)<<offs[j])
 a=(1<<sizes[0])-1;c=((1<<sizes[1])-1)<<offs[1];R=T+[1<<(n+1),0]
 for u in range(n):
  if a>>u&1:R[u]|=1<<n
  else:R[n]|=1<<u
  if c>>u&1:R[u]|=1<<(n+1)
  else:R[n+1]|=1<<u
 return R
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
def main():
 producer=json.loads((SOURCE/'receipt.json').read_text());out=[]
 for rec in producer['records']:
  small,big,h,w=map(rec.get,['small_generator','large_generator','h','weight']);co,A,B,num=poly(small,big,h);need(co==rec['coefficient_vector'],'complete cubic polarization')
  q=[b+w*a for a,b in zip(A,B)];val=sum(x*w**j for j,x in enumerate(co));need(q==rec['artificial_profile'] and val==rec['positive_contrast']>0,'exact profile evaluation');need(q[1]>factorial(big),'unrealizability bound');controls=[]
  for original in rec['literal_polarized_controls']:
   sizes=original['sizes'];T,U=build(sizes,False),build(sizes,True);need(T==original['T_port'] and U==original['U_port'],'physical label construction');p,r=profile(T),profile(U)
   need(p==original['q_T_port'] and r==original['q_U_port'],'independent complete ordered-cover profiles');controls.append({'sizes':sizes,'h':h,'T_profile':p,'U_profile':r,'difference':p[h]-r[h]})
  need(sum(c['difference'] for c in controls)==co[2],'literal polarized coefficient sum')
  out.append({'small':small,'big':big,'h':h,'weight':w,'coefficient_vector':co,'profile':q,'value':val,'polarized_controls':controls})
 odd=out[0]['coefficient_vector'];need(odd[0]+1580*odd[1]+1580**2*odd[2]<0,'first even threshold lower endpoint');need(odd[1]>0 and odd[2]>0 and odd[3]==0,'monotone odd threshold')
 return {'status':'PASS','records':out,'independent_full_cover_profiles':12,'source_pins':{f:hashlib.sha256((SOURCE/f).read_bytes()).hexdigest() for f in ['PROOF.md','verify.py','receipt.json','receipt.optimized.json']},'primary_reference':{'url':'https://www.cip.ifi.lmu.de/~grinberg/algebra/redeiberge-long.pdf','version':'Detailed version, June 13, 2026','theorem':'1.39','printed_page':14},'scope':'Artificial positive same-parity Stirling-cone obstruction only, in both parities; no actual tournament or ambient-response counterexample.'}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
