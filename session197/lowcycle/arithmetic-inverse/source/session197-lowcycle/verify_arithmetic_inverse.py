#!/usr/bin/env python3
"""Exact arithmetic inverse controls and finite positive coverage; -O safe."""
from pathlib import Path
from functools import lru_cache
from math import gcd
import itertools,json,hashlib,subprocess,tempfile,importlib.util
P=Path(__file__).resolve().parent
GAPS=[7,13,21,41,43,55,67,77,83,163,173,231,233,247,251,277,299,421,449,503,533,569,583,637,821,1067,1339,1861,2647]
HARD=[5099741,5477461,8231473,9551177]
def need(v,s):
 if not v:raise RuntimeError(s)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
@lru_cache(None)
def gadget(L):
 if L==1:return 1,1,1
 ell=(L+1)//2;g,Q,B=gadget(ell)
 return (g+ell-ell%2,Q+B,2*B)if L%2 else(g,Q,2*B)
def bits(L):return [1-((L-1)>>i&1)for i in range((L-1).bit_length())]
def transfer(x,L):
 h,e,R=x;g,Q,B=gadget(L)
 return g*h-(L%2)*e+L*R,h-(L==1)*e,Q*h+B*R
def build(rows,entry,L):
 rows=list(rows);b=len(rows);entry=set(entry);rows.append(sum(1<<u for u in entry))
 for u in range(b):
  if u not in entry:rows[u]|=1<<b
 for bit in bits(L):
  v=len(rows);rows.append(0)
  for u in range(v):
   if u==b and not bit:rows[v]|=1<<u
   else:rows[u]|=1<<v
 return rows,b

def tables(rows):
 n=len(rows);N=1<<n;end=[[0]*n for _ in range(N)];start=[[0]*n for _ in range(N)]
 for u in range(n):end[1<<u][u]=start[1<<u][u]=1
 for S in range(1,N):
  for u in range(n):
   if not(S>>u&1):continue
   T=S^(1<<u)
   if not T:continue
   end[S][u]=sum(end[T][v]for v in range(n)if T>>v&1 and rows[v]>>u&1)
   start[S][u]=sum(start[T][v]for v in range(n)if T>>v&1 and rows[u]>>v&1)
 H=[sum(x)for x in end];H[0]=1
 return end,start,H

def profile(rows,r):
 end,start,H=tables(rows);F=(1<<len(rows))-1
 R=sum(start[S][r]*H[F^S]for S in range(1,F+1))
 return H[F],end[F][r],R,start[F][r]
def literal(rows):
 n=len(rows);full=(1<<n)-1;h=0;c=0;starts=[0]*n;ends=[0]*n
 def go(first,last,S):
  nonlocal h,c
  if S==full:
   h+=1;starts[first]+=1;ends[last]+=1
   if first==0 and rows[last]&1:c+=1
   return
  q=rows[last]&~S
  while q:
   bit=q&-q;q-=bit;go(first,bit.bit_length()-1,S|bit)
 for u in range(n):go(u,u,1<<u)
 return h,c,starts,ends

def prefixes(B):
 out=[(1,1,1)];seen=set(out)
 for x in out:
  for L in range(1,B//x[2]+1):
   y=transfer(x,L)
   if y[2]<=B and y not in seen:seen.add(y);out.append(y)
 return out

def inverse(N,x):
 h,e,R=x;out=[]
 if transfer(x,1)[0]==N:out.append(1)
 for t in (0,1):
  low=max(2,(N+h+R-1)//(h+R));high=(N-h+t*e)//R
  lo=max(1,(low-t+1)//2);hi=(high-t)//2
  if lo>hi:continue
  d=gcd(2*R,h);rhs=N+t*e-R*t
  if rhs%d:continue
  mod=h//d;res=0 if mod==1 else rhs//d*pow(2*R//d,-1,mod)%mod
  for k in range(lo+(res-lo)%mod,hi+1,mod):
   L=2*k+t
   if transfer(x,L)[0]==N:out.append(L)
 return sorted(out)

def module(name,path):
 s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def main():
 constant_tests=100000
 for L in range(1,constant_tests+1):
  g,Q,B=gadget(L);need(g%2==1 and 1<=g<=L and Q>=g and B>=L,'gadget inequalities')
  need(2*g-Q>=(3 if L>=3 and L%2 else 1),'gadget deficit')
  if L>=2:need(B<=2*L-2,'dyadic bound')
 physical=0
 for n in range(1,4):
  pairs=list(itertools.combinations(range(n),2))
  for mask in range(1<<len(pairs)):
   rows=[0]*n
   for j,(u,v)in enumerate(pairs):rows[u if mask>>j&1 else v]|=1<<(v if mask>>j&1 else u)
   for r in range(n):
    x=profile(rows,r)[:3]
    for L in range(1,17):
     out,b=build(rows,[r],L);y=profile(out,b)[:3]
     need(y==transfer(x,L),'physical one-seam transfer');need(y[0]<=y[2]<=2*(y[0]-y[1]),'physical cone');physical+=1
 bank=prefixes(50);inverse_tests=0
 for x in bank:
  expected={}
  for L in range(1,2001//x[2]+1):
   N=transfer(x,L)[0]
   if N<=2001:expected.setdefault(N,[]).append(L)
  for N in range(1,2002,2):need(inverse(N,x)==expected.get(N,[]),'exact congruence inverse');inverse_tests+=1
 need(transfer((1,1,1),3)==(5,1,7),'physical inverse hostile')
 need((3,2,4)[2]>2*((3,2,4)[0]-(3,2,4)[1]),'hostile cone')
 c2=[]
 for r in json.loads((P/'arithmetic-c2-witnesses.json').read_text()):
  rows=r['source_rows'];a,b=r['pair'];end,start,H=tables(rows);F=(1<<len(rows))-1
  lit=literal(rows);need(lit[1]==1 and start[F][a]==start[F][b]==1,'pair source')
  R=sum((start[S][a]+start[S][b])*H[F^S]for S in range(1,F+1))
  K=sum(end[S][a]*start[F^S][b]+end[S][b]*start[F^S][a]for S in range(1,F))
  need([H[F],end[F][a]+end[F][b],R,K]==r['response'],'pair response')
  rows,root=build(rows,[a,b],r['pair_L'])
  for L in r['ordinary_Ls']:rows,root=build(rows,[root],L)
  p=profile(rows,root);h,c,s,e=literal(rows)
  need(h==p[0]==r['target'] and c==2 and s[root]==p[3]==2 and e[root]==p[1],'literal C2 witness')
  c2.append(dict(target=h,order=len(rows),rows=rows,root=root,profile=list(p),cycles=c))
 need(sorted(r['target']for r in c2)==sorted(set(GAPS)-{7,21,67}),'all26 repaired')
 D=P.parent/'session197/lowcycle-two-port';compiler=module('arithmetic_compiler',D/'compiler.py');core=module('arithmetic_twoport',D/'twoport_core.py')
 hard=[]
 for r in json.loads((P/'arithmetic-hard-witnesses.json').read_text()):
  rows=[0];root=0;x=(1,1,1)
  for L in r['recipe']:
   rows,root=build(rows,[root],L);x=transfer(x,L)
  need(rows==r['rows']and root==r['root']and x[0]==r['target'],'hard physical recipe')
  comp=compiler.compile_one(rows,{root});state=core.brute([2,0],[0,1]);C=0
  for step in comp['steps'][2:]:
   mask=step['mask'];C=(0,state[1]+state[2],state[3]+state[5],state[2]+state[5])[mask]
   state=core.step(state,2,mask,step['drop'])
  need(core.H(state,2)==r['target']and C==1,'independent uncapped physical interface')
  need(r['recipe'][-1]in inverse(r['target'],tuple(r['prefix_state'])),'hard terminal inverse')
  hard.append(dict(target=r['target'],order=len(rows),profile=list(x),cycles=C))
 coverage=[]
 with tempfile.TemporaryDirectory(prefix='arithmetic-verify-',dir=P)as tmp:
  tmp=Path(tmp);exe=tmp/'coverage';subprocess.run(['g++','-O3','-std=c++17',str(P/'arithmetic_coverage.cpp'),'-o',str(exe)],check=True)
  for B in (500,1000):
   image=tmp/f'image-{B}.bin';proc=subprocess.run([str(exe),'10000001',str(B),str(image)],check=True,text=True,capture_output=True);v=json.loads(proc.stdout)
   need(v['missing_odd']==sorted(GAPS+(HARD if B==500 else[])),'complete positive prefix image')
   need(v['prefix_states']==(2717 if B==500 else 7497),'complete prefix count')
   v['packed_image_sha256']=sha(image);coverage.append(v)
 deps={str(q.relative_to(P.parent)):sha(q)for q in [P/'ROOTED.md',P/'C1_COMPLETE.md',P/'unique_start.py',D/'PROOF.md',D/'compiler.py',D/'twoport_core.py']}
 out=dict(status='PASS',gadget_constant_controls=constant_tests,physical_cone_controls=physical,inverse_controls=inverse_tests,prefix_bank50=len(bank),coverage=coverage,C2_witnesses=c2,large_C1_witnesses=hard,dependencies=deps)
 (P/'arithmetic-inverse-certificate.json').write_text(json.dumps(out,sort_keys=True,indent=2)+'\n')
 print(json.dumps(out,sort_keys=True,indent=2))
if __name__=='__main__':main()
