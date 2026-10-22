import json,numpy as np,pathlib,random,subprocess
from experiment import rot,flip,bits
from test_lollipop import lollipop
ROOT=pathlib.Path(__file__).resolve().parent
p=subprocess.Popen([str(ROOT/'count')],stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True)
S=rot(7);rng=random.Random(30500)
while True:
 p.stdin.write(f'7 {bits(S)}\n');p.stdin.flush()
 if int(p.stdout.readline())==171:break
 S=flip(S,rng,10)
p.stdin.close();p.wait()
out=[]
for k in [1,3,5,9,15,31,63]:
 T=np.kron(S,np.ones((k,k),dtype=int))+np.kron(np.eye(7,dtype=int),rot(k));n=len(T)-1
 for block in [0,1]:
  A=np.delete(np.delete(T,k*block,0),k*block,1);L=float(lollipop(A));r=dict(k=k,n=n,block=block,L=L,nL=n*L,n2L=n*n*L);out.append(r);print(r,flush=True)
json.dump(out,open(ROOT/'lollipop_scaling.json','w'),indent=2)
