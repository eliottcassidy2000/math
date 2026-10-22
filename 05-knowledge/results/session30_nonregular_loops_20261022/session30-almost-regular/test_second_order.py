"""Independent H targets versus explicitly producer-engine predictions."""
import numpy as np,json,pathlib,random,subprocess,sys,time
from experiment import pred,rot,flip,bits,frombits
ROOT=pathlib.Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'session30-decorated-engine'))
from predict import predict

def main():
 cases=[]
 for n in [12,16,20]:
  rng=random.Random(30600+n);S=rot(n+1)
  for j in range(3):S=flip(S,rng,1000);cases.append((f'fresh-chain-{n}-{j}',S[1:,1:]))
 for n in [6,10,18]:
  q=n+1;sq={i*i%q for i in range(1,q)};S=np.array([[0 if i==j else 1 if (j-i)%q in sq else -1 for j in range(n)] for i in range(n)]);cases.append((f'retrospective-paley-{n}',S))
 for r in [json.loads(x) for x in open(ROOT/'lollipop_heldout.jsonl')][:2]:cases.append(('retrospective-'+r['name'],frombits(r['n'],r['bits'])))
 for n in [10,16,20]:cases.append((f'retrospective-rotation-{n}',rot(n+1)[1:,1:]))
 process=subprocess.Popen([str(ROOT/'count')],stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True)
 with open(ROOT/'second_order.jsonl','w') as f:
  for name,S in cases:
   n=len(S);b=bits(S);process.stdin.write(f'{n} {b}\n');process.stdin.flush();H=int(process.stdout.readline());target=pred(S,H);t=time.time();p=predict(S);assert abs(p['rho']-target['rho'])<1e-11
   r=dict(name=name,n=n,bits=b,H=str(H),exact=target,prediction=p,residual_first=target['y']-p['y1'],residual_second=target['y']-p['y2'],seconds=time.time()-t);f.write(json.dumps(r)+'\n');f.flush();print(name,'y',target['y'],'y1',p['y1'],'y2',p['y2'],'nres1',n*r['residual_first'],'n2res2',n*n*r['residual_second'],flush=True)
 process.stdin.close();process.wait()
if __name__=='__main__':main()
