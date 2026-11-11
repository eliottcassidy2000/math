"""Literal completion of N(35,r,2)=5+140*2**r for every r>=1."""
from pathlib import Path
import argparse,json
from head_lift_probe import inflate,ends,converse
HERE=Path(__file__).resolve().parent
H5=(2,4,9,3)
H57=(2,64,67,39,15,23,57)
H113=(34,48,11,67,13,28,55)
SEED=tuple(inflate((2,40,3,5,15,29),4,2))
SELECTED=(2,3)

def ordinal(A,B):
 n=len(A);new=((1<<len(B))-1)<<n
 return [row|new for row in A]+[row<<n for row in B]

def tail(B,O,t):
 if t<1:raise ValueError('The tail must be nonempty.')
 n=len(B);z=n+t;T=list(B)+[0]*(t+1)
 for i in range(n):
  for j in range(n,z):T[i]|=1<<j
  if i in O:T[z]|=1<<i
  else:T[i]|=1<<z
 for j in range(n,z):
  for k in range(j+1,z+1):T[j]|=1<<k
 return T

def seed_metrics(B,O):
 n=len(B);M=(1<<n)-1;e=ends(B);s=ends(converse(B));h=[sum(x) for x in e];h[0]=1
 U=sum(sum(s[m][i] for i in O)*h[M^m] for m in range(1,M+1))
 W=sum(sum(e[m][i] for i in O)*sum(s[M^m][j] for j in O) for m in range(1,M))
 return h[M],U,W

def realize(r):
 if r<1:raise ValueError('r must be positive')
 inner=H57 if r==1 else H113 if r==2 else tail(SEED,SELECTED,r-2)
 return ordinal(H5,inner)

def check():
 H,U,W=seed_metrics(SEED,SELECTED)
 assert H-W==1 and U==112
 records=[]
 for r in range(1,7):
  T=realize(r);actual=sum(ends(T)[-1]);expected=5+140*2**r
  assert actual==expected
  records.append({'r':r,'order':len(T),'H':actual,'rows':T})
 receipt={'status':'FINITE-EXACT controls; all-r proof in PROOF.md',
  'seed':{'rows':SEED,'selected':SELECTED,'H_U_W':[H,U,W]},
  'H5_rows':H5,'H57_rows':H57,'H113_rows':H113,'literal_controls':records}
 (HERE/'realize35_s2_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
 print(json.dumps({'status':'pass','seed_H_U_W':[H,U,W],'H':[r['H'] for r in records]}))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--r',type=int);a=p.parse_args()
 if a.r is None:check()
 else:
  T=realize(a.r);print(json.dumps({'r':a.r,'H':5+140*2**a.r,'order':len(T),'rows':T},indent=2))
