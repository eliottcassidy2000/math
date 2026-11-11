from module_tools import *
import json

def convolution(a,b,n):
 out=[0]*(1<<n)
 for mask in range(1<<n):
  s=mask
  while True:
   out[mask]+=a[s]*b[mask^s]
   if not s:break
   s=(s-1)&mask
 return out

def tensor(T,O):
 n=len(T);e=m.endpoints(T);rev=tuple(sum(1<<i for i in range(n) if T[i]>>j&1) for j in range(n));s=m.endpoints(rev)
 H=[sum(z) for z in e];E=[sum(z[i] for i in range(n) if O>>i&1) for z in e];S=[sum(z[i] for i in range(n) if O>>i&1) for z in s]
 power=[0]*(1<<n);power[0]=1;K=convolution(E,S,n);ret=[]
 for r in range(n+1):
  ret.append((power[-1],convolution(S,power,n)[-1],convolution(E,power,n)[-1],convolution(K,power,n)[-1]))
  power=convolution(H,power,n)
 return tuple(ret)

seeds=json.loads(Path(__file__).with_name('even_births.json').read_text())['seeds'];D={}
for T,O,q in seeds:D.setdefault(tensor(tuple(T),O),[]).append((T,O,q))
hits=[{'tensor':k,'seeds':v} for k,v in D.items() if len(v)>1]
for h in hits:print(h,flush=True)
Path(__file__).with_name('marked_trades.json').write_text(json.dumps(hits,indent=2))
