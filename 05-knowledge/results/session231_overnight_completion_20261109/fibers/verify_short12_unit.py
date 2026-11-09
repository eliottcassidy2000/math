import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from pathlib import Path
import json,math
from core import compose,states
from verify import literal,need
w=json.loads(Path(__file__).with_name('SHORT12_UNIT_LOOP_WITNESS.json').read_text());core=states(w['adj'])[w['mask']];need(core==literal(w['adj'])[w['mask']],'literal core')
h,q,A,G,T,S,E=core;X=list(map(int,w['X']));Y=list(map(int,w['Y']));d=len(X);need(len(Y)==d,'length');u=X[-1];v=Y[0]
a,g,p=A,G,Y[-1];gain=0
for x,y in zip(X,Y):
 need(x==g%2,'parity');aa=2*a+h*y-q*x;gg=(g+q*y-h*x)//2
 need(-h<aa<q and -h<gg<q,'bounds');wt=a*y+h*p*y-x*gg-h*x*(gg%2);gain+=wt;a,g,p=aa,gg,y
need((a,g,p)==(A,G,Y[-1]),'not closed');need(gain==w['cycle_gain'] and abs(gain)==2,'gain')
one=[(1,2,b,b,0,b,b) for b in (0,1)]
def fold(ps):
 s=ps[0]
 for z in ps[1:]:s=compose(s,z)
 return s
left=fold([one[u],core]+[one[b] for b in Y]+[one[v]])
right=fold([one[u]]+[one[b] for b in X]+[core,one[v]])
need(left[:4]+left[5:]==right[:4]+right[5:],'base interface');need(right[4]-left[4]==-gain,'base difference')
checks=[]
for k in range(9):
 vals=[];iface=None
 for j in range(k+1):
  ps=([one[1]] if not u else [])+[one[u]]+[one[b] for b in X]*j+[core]+[one[b] for b in Y]*(k-j)+[one[v]];z=fold(ps)
  sig=z[:4]+z[5:]
  if iface is None:iface=sig
  need(sig==iface,'interface');vals.append(2*z[2]+h-z[4])
 need(all(c==vals[0]+gain*j for j,c in enumerate(vals)),'consecutive odds')
 checks.append(dict(k=k,vertices=d*k+len(w['adj'])+4+(not u),count=len(vals),first_mod_10000=vals[0]%10000,last_mod_10000=vals[-1]%10000,first_bitlength=vals[0].bit_length()))
print(json.dumps(dict(status='PASS',length=d,core=core,first_y=v,last_x=u,cycle_gain=gain,step=2,checks=checks),sort_keys=True,indent=2))
