import numpy as np,json,pathlib,random,subprocess
from experiment import pred,rot,flip,bits,frombits
ROOT=pathlib.Path(__file__).resolve().parent

def lollipop(S):
 n=len(S);M=S/n;Q=M@M;G=-Q@np.linalg.inv(np.eye(n)+Q);B=Q@np.linalg.inv(np.eye(n)-Q)
 return -2*np.diag(B-G)@G@np.ones(n)

def main():
 rows=[json.loads(x) for x in open(ROOT/'data.jsonl')]
 for r in rows:
  S=frombits(r['n'],r['bits']);L=float(lollipop(S));r.update(lollipop=L,y_lollipop=r['y_tree']-r['n']*L,residual_lollipop=r['residual']+r['n']*L)
 with open(ROOT/'lollipop_retrospective.jsonl','w') as f:
  for r in rows:f.write(json.dumps(r)+'\n')
 proc=subprocess.Popen([str(ROOT/'count')],stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True)
 def count(S):proc.stdin.write(f'{len(S)} {bits(S)}\n');proc.stdin.flush();return int(proc.stdout.readline())
 S=rot(7);rng=random.Random(30500)
 while count(S)!=171:S=flip(S,rng,10)
 T=np.kron(S,np.ones((3,3),dtype=int))+np.kron(np.eye(7,dtype=int),rot(3));assert np.all(T.sum(axis=1)==0)
 with open(ROOT/'lollipop_heldout.jsonl','w') as f:
  for v in range(0,21,3):
   A=np.delete(np.delete(T,v,0),v,1);H=count(A);r=dict(name=f'nonvt7-lex3-delete-{v}',n=20,bits=bits(A),H=str(H),**pred(A,H));L=float(lollipop(A));r.update(lollipop=L,y_lollipop=r['y_tree']-20*L,residual_lollipop=r['residual']+20*L);f.write(json.dumps(r)+'\n');f.flush();print(r['name'],'y',r['y'],'path',r['y_tree'],'new',r['y_lollipop'],'L',L,flush=True)
 proc.stdin.close();proc.wait()
if __name__=='__main__':main()
