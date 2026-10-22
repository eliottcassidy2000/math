import numpy as np, itertools, math, random, json, subprocess, pathlib
ROOT=pathlib.Path(__file__).resolve().parent

def bits(S):return ''.join('1' if S[i,j]>0 else '0' for i in range(len(S)) for j in range(i+1,len(S)))
def frombits(n,b):
 S=np.zeros((n,n),dtype=int)
 for v,(i,j) in zip(b,itertools.combinations(range(n),2)):S[i,j]=1 if v=='1' else -1;S[j,i]=-S[i,j]
 return S

def rot(n,signs=None):
 if signs is None:signs=[1]*(n//2)
 S=np.zeros((n,n),dtype=int)
 for i in range(n):
  for d,s in enumerate(signs,1):S[i,(i+d)%n]=s;S[(i+d)%n,i]=-s
 return S

def flip(S,rng,count=1000):
 S=S.copy(); done=0
 while done<count:
  a,b,c=rng.sample(range(len(S)),3)
  if S[a,b]==S[b,c]==S[c,a]:
   for i,j in [(a,b),(b,c),(c,a)]:S[i,j]*=-1;S[j,i]*=-1
   done+=1
 return S

def pred(S,H):
 n=len(S);M=S/n;M2=M@M;I=np.eye(n);G=-M2@np.linalg.inv(I+M2);B=M2@np.linalg.inv(I-M2);O=M@np.linalg.inv(I-M2);d=np.diag(B-G)
 f8=np.sum(2*np.diag(G)*np.diag(B)-1.5*np.diag(G)**2);theta=np.sum(G*(B**2-O**2-2/3*G**2));dumb=-d@G@d;X1=f8+theta+dumb
 xs=np.linalg.eigvalsh(1j*M)**2
 D1=np.sum(2*xs/(1-xs**2));D2=np.sum(2*xs/(1-xs)**2+2*xs/(1+xs)**2);corr=(D2+D1**2-D1)/2
 rho=math.exp(.5*(np.linalg.slogdet(I-M2)[1]-np.linalg.slogdet(I+M2)[1]));kappa=H*2**(n-1)/math.factorial(n);delta=M@np.ones(n);tree=n*delta@np.linalg.solve(I+M2,delta);yr=-n*X1-corr;y=n*(1-kappa/rho)
 return dict(rho=rho,kappa=kappa,y=y,y_regular=yr,tree=tree,y_tree=yr+tree,residual=y-yr-tree,n_residual=n*(y-yr-tree),n2_residual=n*n*(y-yr-tree),X1=X1,f8=f8,theta=theta,dumbbell=dumb,D1=D1,D2=D2,score_imbalance_sq=int(np.sum(S.sum(axis=1)**2)))

def brute(S):return sum(all(S[p[i],p[i+1]]==1 for i in range(len(S)-1)) for p in itertools.permutations(range(len(S))))

def run():
 cases=[]
 for n in [2,4,6]:
  for m in range(1<<(n*(n-1)//2)):
   b=format(m,'0%db'%(n*(n-1)//2));S=frombits(n,b)
   if np.all(abs(S.sum(axis=1))==1):cases.append((f'exhaustive-labelled-{n}-{m}',S))
 for n in [8,10,12,14,16,18,20]:
  R=rot(n+1);cases.append((f'rot-delete-{n}',R[1:,1:]));rng=random.Random(30001+n);C=rot(n+1,[rng.choice([-1,1]) for _ in range((n+1)//2)]);cases.append((f'circulant-delete-{n}',C[1:,1:]));Q=flip(R,random.Random(30002+n));cases.append((f'chain-delete-{n}',Q[1:,1:]))
 for q in [7,11,19]:
  sq={i*i%q for i in range(1,q)};S=np.array([[0 if i==j else (1 if (j-i)%q in sq else -1) for j in range(q)] for i in range(q)]);cases.append((f'paley-delete-{q-1}',S[1:,1:]))
 for n in [4,6,8,10,12]:cases.append((f'hostile-transitive-{n}',np.triu(np.ones((n,n),dtype=int),1)-np.tril(np.ones((n,n),dtype=int),-1)))
 for n in [5,7,9,11]:cases.append((f'control-regular-{n}',rot(n)))
 process=subprocess.Popen([str(ROOT/'count')],stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True)
 counts={};out=[]
 with open(ROOT/'data.jsonl','w') as f:
  for name,S in cases:
   n=len(S);b=bits(S);process.stdin.write(f'{n} {b}\n');process.stdin.flush();H=int(process.stdout.readline());assert H%2==1
   if n<=6 and (n not in counts or name.startswith('paley')):assert H==brute(S);counts[n]=1
   row=dict(name=name,n=n,bits=b,H=str(H),**pred(S,H));out.append(row);f.write(json.dumps(row)+'\n');f.flush()
   if not name.startswith('exhaustive'):print(name,H,'y',round(row['y'],9),'pred',round(row['y_tree'],9),'nres',round(row['n_residual'],9),flush=True)
 process.stdin.close();process.wait()
 with open(ROOT/'summary.json','w') as f:
  summary={}
  for n in [2,4,6]:
   group=[r for r in out if r['name'].startswith(f'exhaustive-labelled-{n}-')];res=np.array([r['residual'] for r in group]);summary[str(n)]=dict(count=len(group),H_distinct=sorted(set(int(r['H']) for r in group)),residual_min=float(min(res)),residual_max=float(max(res)),residual_mean=float(np.mean(res)),residual_sd=float(np.std(res)))
  json.dump(summary,f,indent=2);print(summary)
if __name__=='__main__':run()
