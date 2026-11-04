"""Exact polynomial and literal controls for same-parity profile-cone obstruction."""
from pathlib import Path
import importlib.util,itertools,json,math,hashlib
ROOT=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('covered',ROOT.parent/'covered_response'/'verify.py');v=importlib.util.module_from_spec(s);s.loader.exec_module(v)
E=v.add(v.X,v.Y,v.Z,v.prod(v.X,v.Y),v.prod(v.X,v.Z),v.prod(v.Y,v.Z),v.scale(v.P,3))
def convolve(a,b):
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return c
def polynomial(small,big,h):
 A=v.profile(v.tt(small))+[0]*(big-small);B=v.profile(v.tt(big));power=v.ONE
 for _ in range(h-1):power=v.mul(power,E)
 num=v.mul(v.N,power);out=[0]*4
 for k in range(big+1):
  factor=math.comb(k+h,h)
  for ex,c in num.items():
   term=[c*factor*h]
   for r in ex:
    j=k+r;term=convolve(term,[B[j],A[j]] if j<=big else [0,0])
   for j,x in enumerate(term):out[j]+=x
 v.need(all(x%3==0 for x in out),'integral functional polynomial coefficients');return [x//3 for x in out],A,B,num
def build(sizes,reverse=False):
 off=[sum(sizes[:i]) for i in range(3)];rows=[0]*sum(sizes)
 for i,n in enumerate(sizes):
  j=(i+(2 if reverse else 1))%3
  for u,r in enumerate(v.tt(n)):rows[off[i]+u]=(r<<off[i])|(((1<<sizes[j])-1)<<off[j])
 return rows,off
def literal(small,big,h):
 controls=[]
 for position in range(3):
  sizes=[small]*3;sizes[position]=big;T,off=build(sizes);U,_=build(sizes,True);a=(1<<sizes[0])-1;c=((1<<sizes[1])-1)<<off[1];X=v.extend(T,a,c);Y=v.extend(U,a,c);p=v.profile(X);q=v.profile(Y);controls.append({'sizes':sizes,'query':[a,c],'T_port':X,'U_port':Y,'q_T_port':p,'q_U_port':q,'h':h,'contrast':p[h]-q[h]})
 return controls
records=[]
for small,big,h,w,expected in [(1,5,2,1582,[-2910041940,1746120,60,0]),(2,6,4,8192,[-223898369213472,41383812864,834144,-192])]:
 coeff,A,B,num=polynomial(small,big,h);v.need(coeff==expected,'formal one-variable polynomial identity');q=[b+w*a for a,b in zip(A,B)];value=sum(c*w**j for j,c in enumerate(coeff));v.need(value>0 and h*v.functional(num,h+1,q)==3*value,'exact positive cone violation');v.need(q[1]>math.factorial(big),'artificial profile is not a tournament profile');controls=literal(small,big,h);v.need(sum(c['contrast'] for c in controls)==coeff[2],'literal polarized mixed coefficient');records.append({'small_generator':small,'large_generator':big,'h':h,'weight':w,'coefficient_vector':coeff,'artificial_profile':q,'positive_contrast':value,'literal_polarized_controls':controls})
v.need(sum(c*1580**j for j,c in enumerate(records[0]['coefficient_vector']))<0,'first even crossing boundary')
print(json.dumps({'status':'PASS','scope':'Positive same-parity Stirling cones alone do not imply higher ordered-cover domination. Artificial cone profiles, not tournament or full-response counterexamples.','records':records,'source_sha256':{str(p.relative_to(ROOT.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'PROOF.md',ROOT/'verify.py',ROOT.parent/'covered_ambient'/'PROOF.md',ROOT.parent/'covered_response'/'verify.py']}},sort_keys=True,indent=2))
