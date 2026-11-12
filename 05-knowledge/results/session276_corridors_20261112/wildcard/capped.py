"""Exact capped positive automaton for two-apex Hamiltonian path counts."""
from collections import deque,Counter
import json,sys,time
K=36
OUT=tuple(i for a in range(2) for r in range(3) for t in range(3) for z in range(2) if r!=2 and t!=2 and not(t==0 and a==1) for i in [a*18+r*6+t*2+z])
TRANS=[]
for l in range(4):
 a=l&1;b=l>>1;tr=[]
 for arr in range(2):
  x,y=(a,b) if not arr else (b,a)
  for r in range(3):
   for t in range(3):
    for z in range(2):
     dst=[arr*18+(1 if x==0 else 2)*6+t*2+z]
     if t!=0 or x==1:dst.append(arr*18+r*6+(1 if y==0 else 2)*2+z)
     if z!=0 or y==1:dst.append(arr*18+r*6+t*2+1)
     tr.append(tuple(dst))
 TRANS.append(tr)
def step(v,l,cap):
 o=[0]*36
 for c,ds in zip(v,TRANS[l]):
  if c:
   for j in ds:o[j]=min(cap,o[j]+c)
 return tuple(o)
def total(v):return sum(v[i] for i in OUT)
def init():
 v=[0]*36;v[0]=v[18]=1;return tuple(v)
def count(w):
 v=init()
 for l in w:v=step(v,l,10**100)
 return total(v)
def build(N,maxstates=2000000):
 vs=[init()];ix={vs[0]:0};ts=[];words=[()];hs=[1];i=0
 while i<len(vs):
  row=[]
  for l in range(4):
   z=step(vs[i],l,N+1)
   if total(z)>N:row.append(-1);continue
   if z not in ix:
    if len(vs)>=maxstates:raise RuntimeError(('state budget',N,len(vs)))
    ix[z]=len(vs);vs.append(z);words.append(words[i]+(l,));hs.append(total(z))
   row.append(ix[z])
  ts.append(tuple(row));i+=1
  if i%10000==0:print('progress',N,i,len(vs),flush=True)
 return vs,ts,words,hs
# Moore refinement keeps the entire count output, not just H=N.
def minimize(ts,hs,dead_output=-1):
 tr=ts+[(len(ts),)*4];out=hs+[dead_output]
 tr=[tuple(len(ts) if x<0 else x for x in r) for r in tr]
 part=[{h:i for i,h in enumerate(sorted(set(out)))}[h] for h in out]
 while True:
  sig=[(out[i],*(part[j] for j in r)) for i,r in enumerate(tr)]
  dic={};nxt=[]
  for s in sig:
   if s not in dic:dic[s]=len(dic)
   nxt.append(dic[s])
  if len(set(nxt))==len(set(part)):break
  part=nxt
 rep={}
 for i,p in enumerate(nxt):rep.setdefault(p,i)
 return [tuple(nxt[j] for j in tr[rep[p]]) for p in range(len(rep))],[out[rep[p]] for p in range(len(rep))],nxt
if __name__=='__main__':
 N=int(sys.argv[1]) if len(sys.argv)>1 else 21
 t=time.time();vs,ts,ws,hs=build(N);mt,mh,mp=minimize(ts,hs)
 print(json.dumps({'N':N,'raw_states':len(vs),'moore_states':len(mt),'max_shortest_word':max(map(len,ws)),'values':sorted(set(hs)),'minimal_witnesses':{h:''.join(map(str,ws[hs.index(h)])) for h in sorted(set(hs))},'seconds':round(time.time()-t,3)},indent=2))
 dest=sys.argv[2] if len(sys.argv)>2 else '/tmp/capped'+str(N)+'.json'
 with open(dest,'w') as f:json.dump({'N':N,'transitions':mt,'outputs':mh,'initial':mp[0],'raw_states':len(vs),'max_shortest_word':max(map(len,ws))},f)
