
from pathlib import Path
from runpy import run_path
from sys import argv
base=run_path(str(Path(__file__).with_name('gewirtz-transport.py')))
globals().update({k:v for k,v in base.items() if not k.startswith('__')})
out=argv[1] if len(argv)>1 else 'rainbow-data.bin'
s=0; prev=0; enumerator=Counter()
for i in range(1<<21):
    g=i^(i>>1)
    if i:s^=line_bits[(g^prev).bit_length()-1]
    if s==0:enumerator[g.bit_count()]+=1
    prev=g
assert enumerator=={0:1,6:168,8:210,10:1008,12:280,14:360,16:21}
print('line code weight enumerator',dict(sorted(enumerator.items())),flush=True)
points=[sum(1<<j for j,h in enumerate(orbit) if h>>i&1) for i in range(21)]
Pi={s:i for i,s in enumerate(points)}
Pbits=[[j for j,h in enumerate(orbit) if h>>i&1] for i in range(21)]
def pts(p):
    return bytes(Pi[sum(1<<p[j] for j in v)] for v in Pbits)
Candidates={}; Free={}
for p in G1:
    if any(i==j for i,j in enumerate(p)):continue
    t=pts(p)
    k=sum(i==j for i,j in enumerate(t))
    if k==1:Candidates[p]=t
    if k==0:Free[p]=t
assert len(Candidates)==7560
Gen=G+[pl_perm(frob),du]
def inv(p):return bytes(p.index(i) for i in range(56))
# Use only point-preserving conjugation to preserve named source points.
tt=lambda p:p+bytes(range(56,256))
gen=[(tt(g),tt(inv(g))) for g in Gen[:-1]]
cs=[];rest=set(Candidates); reps=[]
while rest:
    p=min(rest);reps.append(p);Q=[p];rest.remove(p)
    seen={p}
    for u in Q:
      for g,h in gen:
        v=bytes(g[u[h[i]]] for i in range(56))
        assert v in Candidates
        if v not in seen:seen.add(v);rest.remove(v);Q.append(v)
    cs.append(len(Q))
assert sorted(cs)==[2520,5040]
print('class sizes',cs,flush=True)
position={p:i for i,p in enumerate(Candidates)}
with open(out,'wb') as fp:
  import struct
  fp.write(struct.pack('<II',len(Candidates),len(reps)))
  for p in reps:fp.write(struct.pack('<I',position[p]))
  for p,q in Candidates.items():fp.write(p+q)

rest=set(Free);reps=[];sizes=[]
while rest:
    p=min(rest);reps.append(p);Q=[p];rest.remove(p);seen={p}
    for u in Q:
      for g,h in gen:
        v=bytes(g[u[h[i]]] for i in range(56))
        assert v in Free
        if v not in seen:seen.add(v);rest.remove(v);Q.append(v)
    sizes.append(len(Q))
assert sizes==[2880]*4
print('zero-color class sizes',sizes,flush=True)
with open(out,'ab') as fp:
    fp.write(struct.pack('<I',len(reps)))
    for p in reps:fp.write(p+Free[p])

    fp.write(struct.pack('<21I',*line_bits))
