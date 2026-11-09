import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from pathlib import Path
import json,math
from core import compose,states
from verify import literal,need
w=json.loads(Path(__file__).with_name('LOOP_WITNESS.json').read_text());core=states(w['adj'])[w['mask']]
need(core==literal(w['adj'])[w['mask']],'literal core');h,q,A,G,T,S,E=core
X=list(map(int,w['X']));Y=list(map(int,w['Y']));d=len(X)
need(len(Y)==d and X[-1]==0 and Y[0]==0,'boundary');a,g,p=A,G,Y[-1];gain=0
for x,y in zip(X,Y):
 need(x==g%2,'parity')
 aa=2*a+h*y-q*x;gg=(g+q*y-h*x)//2
 need(-h<aa<q and -h<gg<q,'bounds')
 wt=a*y+h*p*y-x*gg-h*x*(gg%2);need(wt%2==0,'edge parity');need((aa*gg-a*g-q*wt)%h==0,'cocycle congruence')
 gain+=wt;a,g,p=aa,gg,y
need((a,g,p)==(A,G,Y[-1]),'not a closed walk');need(gain==190,'gain');need(gain%(2*h//math.gcd(h,q))==0,'quantization')
one=[(1,2,b,b,0,b,b) for b in (0,1)]
def fold(ps):
 s=ps[0]
 for z in ps[1:]:s=compose(s,z)
 return s
left=fold([one[0],core]+[one[b] for b in Y]+[one[0]])
right=fold([one[0]]+[one[b] for b in X]+[core,one[0]])
need(left[:4]+left[5:]==right[:4]+right[5:],'base interface');need(right[4]-left[4]==-190,'base difference')
B=2**d;R=int(w['Y'],2);V=int(w['Y'][::-1],2)
checks=[]
for k in range(11):
 vals=[];iface=None
 for j in range(k+1):
  ps=[one[1],one[0]]+[one[b] for b in X]*j+[core]+[one[b] for b in Y]*(k-j)+[one[0]];z=fold(ps)
  expected=(95,1276*2**(d*k+3),2822*2**(d*k+1)+190*R*(B**k-1)//(B-1),1583+5104*V*(B**k-1)//(B-1),95,0)
  need(z[:4]+z[5:]==expected,'explicit interface');vals.append(2*z[2]+h-z[4])
 need(all(v==vals[0]+190*j for j,v in enumerate(vals)),'progression')
 checks.append(dict(k=k,vertices=d*k+12,count=len(vals),first_count=vals[0],last_count=vals[-1]))
print(json.dumps(dict(status='PASS',length=d,core=core,cycle_gain=gain,step=190,checks=checks),sort_keys=True,indent=2))
