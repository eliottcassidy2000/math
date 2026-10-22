import json,random,subprocess,pathlib
import numpy as np
from experiment import frombits,bits,rot,flip
ROOT=pathlib.Path(__file__).resolve().parent

def charpoly(S):
 n=len(S);P=np.eye(n,dtype=object);tr=[n]
 for k in range(1,n+1):P=P@S;tr.append(sum(P[i,i] for i in range(n)))
 c=[1]
 for k in range(1,n+1):
  num=-sum(c[k-i]*tr[i] for i in range(1,k+1));assert num%k==0;c.append(int(num//k))
 return tuple(c)

def main():
 proc=subprocess.Popen([str(ROOT/'count')],stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True)
 def count(S):proc.stdin.write(f'{len(S)} {bits(S)}\n');proc.stdin.flush();return int(proc.stdout.readline())
 out=[]
 groups={}
 for n,rs in json.load(open(ROOT/'small_classes.json')).items():
  for r in rs:groups.setdefault((int(n),charpoly(frombits(int(n),r['bits']))),[]).append(r)
 print('small_groups',groups,flush=True)
 for (n,cp),group in groups.items():
  for a in group:
   for b in group:
    if int(a['H'])<int(b['H']):out.append(dict(n=n,charpoly=list(cp),first=a,second=b))
 for n in [8,10,12]:
  seen={};S=rot(n+1);rng=random.Random(30400+n)
  for trial in range(2000):
   S=flip(S,rng,10);T=S[1:,1:];b=bits(T);cp=charpoly(T)
   if cp in seen:
    prior=seen[cp]
    if b==prior['bits']:continue
    H=count(T)
    if H!=int(prior['H']):
     row=dict(n=n,charpoly=list(cp),first=prior,second=dict(bits=b,H=str(H),trial=trial));out.append(row);print('FOUND',json.dumps(row),flush=True);break
   else:seen[cp]=dict(bits=b,H=str(count(T)),trial=trial)
  else:print('NO_COLLISION',n,len(seen),flush=True)
 json.dump(out,open(ROOT/'isospectral_witnesses.json','w'),indent=2)
 proc.stdin.close();proc.wait()
if __name__=='__main__':main()
