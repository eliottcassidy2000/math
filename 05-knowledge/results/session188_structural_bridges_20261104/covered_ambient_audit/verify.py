"""Independent exact checks of the ordered-cover ambient reduction."""
from pathlib import Path
from itertools import permutations
import hashlib,json,math,sys
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'covered_quotient'))
from verify import cycle,extend,profile,need

def plus(*polys):
 out={}
 for p in polys:
  for e,c in p.items():out[e]=out.get(e,0)+c
 return {e:c for e,c in out.items() if c}
def scale(p,c):return {e:v*c for e,v in p.items() if v*c}
def times(p,q):
 out={}
 for e,c in p.items():
  for f,d in q.items():
   z=tuple(x+y for x,y in zip(e,f));out[z]=out.get(z,0)+c*d
 return {e:c for e,c in out.items() if c}
def mult(*p):
 q={(0,0,0):1}
 for a in p:q=times(q,a)
 return q
ONE={(0,0,0):1};X={(1,0,0):1};Y={(0,1,0):1};Z={(0,0,1):1};XYZ=mult(X,Y,Z);R=plus(ONE,scale(XYZ,-1))
N=[[X,mult(X,Y),XYZ],[XYZ,Y,mult(Y,Z)],[mult(X,Z),XYZ,Z]]
E=plus(X,Y,Z,mult(X,Y),mult(X,Z),mult(Y,Z),scale(XYZ,3))
K=plus(mult(X,X,Y,Y),mult(X,X,Z,Z),mult(Y,Y,Z,Z),mult(X,X,Y),mult(X,X,Z),mult(Y,Y,X),mult(Y,Y,Z),mult(Z,Z,X),mult(Z,Z,Y),scale(mult(XYZ,plus(X,Y,Z)),-1),scale(XYZ,-6))
def singles(N):
 def left(j):return plus(R,*[N[i][j] for i in range(3)])
 def right(omit):return plus(R,*[N[i][j] for i in range(3) if i!=omit for j in range(3)])
 return mult(left(0),right(0),left(1),right(1))
def perm(p,pi):
 out={}
 for e,c in p.items():
  f=[0]*3
  for i,j in enumerate(pi):f[j]=e[i]
  out[tuple(f)]=c
 return out

def contrast(q,h):
 num=K
 for _ in range(h-1):num=times(num,E)
 z=0
 for e,c in num.items():
  for k in range(len(q)):
   if max(e)+k>=len(q):continue
   z+=h*c*math.comb(k+h,h)*math.prod(q[t+k] for t in e)
 need(z%3==0,'integer contrast');return z//3

def main():
 D=plus(singles(N),scale(singles(list(map(list,zip(*N)))),-1))
 need(plus(D,perm(D,(1,0,2)))=={},'exact separated-port exchange')
 need(plus(*[perm(D,p) for p in permutations(range(3))])=={},'exact symmetrized separated ports')
 cases=[]
 for A in ([0],[2,0],[6,4,0],[2,4,1]):
  T,_=cycle([A,A,A]);U,_=cycle([A,A,A],True);m=len(A);mask=(1<<m)-1
  T=extend(T,mask,mask<<m);U=extend(U,mask,mask<<m);q=profile(A)
  p=profile(T);r=profile(U);actual=[a-b for a,b in zip(p,r)]
  predicted=[0]+[contrast(q,h) for h in range(1,len(T)+1)]
  need(actual==predicted,'literal whole ordered profile')
  cases.append({'child':A,'q':q,'T_minus_U_profile':actual})
 artificial=[0,9,0,1,1,1,0];need(contrast(artificial,2)==146,'artificial sign obstruction')
 root=Path(__file__).resolve().parent.parent
 print(json.dumps({'status':'PASS','literal_ordered_profiles':cases,'artificial_h2':146,'independent_separated_numerator_terms':len(D),'source_sha256':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [root/'covered_ambient'/'PROOF.md',root/'covered_ambient'/'verify.py',root/'covered_quotient'/'verify.py',Path(__file__).resolve()]}},indent=2,sort_keys=True))
if __name__=='__main__':main()
