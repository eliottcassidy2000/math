from pathlib import Path
from itertools import combinations
from functools import lru_cache
from math import factorial
import hashlib,json,subprocess,tempfile
ROOT=Path(__file__).resolve().parent;SOURCE=ROOT.parent/'monotone_root_response'
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def profile(T):
 n=len(T);N=1<<n;D=[[0]*n for _ in range(N)];H=[0]*N;H[0]=1
 for S in range(1,N):
  rem=S
  while rem:
   bit=rem&-rem;rem-=bit;u=bit.bit_length()-1;rest=S^bit
   if not rest:D[S][u]=1;continue
   nxt=T[u]&rest
   while nxt:
    b=nxt&-nxt;nxt-=b;D[S][u]+=D[rest][b.bit_length()-1]
  H[S]=sum(D[S])
 U=[[0]*(n+1) for _ in range(N)];U[0][0]=1
 for S in range(1,N):
  least=S&-S;rest=S^least;A=rest
  while True:
   block=A|least;tail=S^block
   for k in range(1,n+1):U[S][k]+=H[block]*U[tail][k-1]
   if not A:break
   A=(A-1)&rest
 return [factorial(k)*v for k,v in enumerate(U[-1])]
def op(T):return [r^((1<<len(T))-1)^(1<<i) for i,r in enumerate(T)]
def blow(Q,blocks):
 offs=[sum(map(len,blocks[:i])) for i in range(len(blocks))];out=[0]*sum(map(len,blocks))
 for i,B in enumerate(blocks):
  for u,r in enumerate(B):
   out[offs[i]+u]=r<<offs[i]
   for j,C in enumerate(blocks):
    if Q[i]>>j&1:out[offs[i]+u]|=((1<<len(C))-1)<<offs[j]
 return out
def extend(T,a,c):
 n=len(T);R=list(T)+[1<<(n+1),0]
 for u in range(n):
  if a>>u&1:R[u]|=1<<n
  else:R[n]|=1<<u
  if c>>u&1:R[u]|=1<<(n+1)
  else:R[n+1]|=1<<u
 return R
def sos(q):
 a=q+[0,0];n=len(q)-1;z=a[1]**3
 for j in range(1,n+2):z+=j*(a[j-1]-a[j+1])**2*(a[j-1]+3*a[j]+2*a[j+1])
 for j in range(1,n+1):z+=(a[j]-a[j+1])**2*((j+2)*a[j]+(2*j+1)*a[j+1])
 need(z%3==0 and z>0,'strict square certificate');return z//3

def hamilton(T):
 n=len(T);d=[[0]*n for _ in range(1<<n)]
 for S in range(1,1<<n):
  for v in range(n):
   if not S>>v&1:continue
   R=S^(1<<v)
   d[S][v]=1 if not R else sum(d[R][w] for w in range(n) if R>>w&1 and T[v]>>w&1)
 return sum(d[-1])
def modules(T):
 out=[]
 for M in range(1,(1<<len(T))-1):
  verts=[v for v in range(len(T)) if M>>v&1]
  if len(verts)<2:continue
  if all(all(bool(T[z]>>v&1)==bool(T[z]>>verts[0]&1) for v in verts) for z in range(len(T)) if z not in verts):out.append(M)
 return out
def xor_family(T,family):
 n=len(T);R=[0]*n
 for i,j in combinations(range(n),2):
  direction=(T[i]>>j)&1
  for M in family:direction^=bool(M>>i&1 and M>>j&1)
  R[i if direction else j]|=1<<(j if direction else i)
 return R
def lift(T,ambient):
 Q=[2,4,1] if ambient else [0];A=blow(Q,[T]+[[0] for _ in Q[1:]]);ext=1<<10 if ambient else 0;return A,7|ext,56|ext

def main():
 data=json.loads((SOURCE/'normal.json').read_text());jobs=[];expected=[];records=[];sources=[];strictcounts={};basis=[]
 lookup={(r['child'],r['cyclic_exterior'],tuple(r['family'])):r for r in data['records']};source_lookup={(r['child'],r['cyclic_exterior']):r for r in data['sources']};need(len(lookup)==1040,'unique complete family records')
 for label,child in [('TT3',[6,4,0]),('C3',[2,4,1])]:
  T=blow([2,4,1],[child]*3);mods=modules(T);need(len(mods)==(9 if label=='TT3' else 3),'independent original modules');basis.append([label,mods]);q=profile(child);P=sos(q);H=hamilton(T);qcache={}
  for ambient in [False,True]:
   A,a,c=lift(T,ambient);source=source_lookup[label,ambient];need(A==source['old_rows'] and [a,c]==source['query'],'fixed source graph and physical lift');sources.append([label,ambient,source['value']]);jobs.append(extend(A,a,c));expected.append(source['value']);strict=0
   for bits in range(1<<len(mods)):
    family=[511]+[M for i,M in enumerate(mods) if bits>>i&1];target=xor_family(T,family);r=lookup[label,ambient,tuple(family)];targetprofiles=[]
    for b in range(3):
     C=tuple((target[3*b+i]>>(3*b))&7 for i in range(3))
     if C not in qcache:qcache[C]=profile(C)
     targetprofiles.append(qcache[C]);need(all(x>=y for x,y in zip(qcache[C],q)),'canonical-partition target profile domination')
    HU=hamilton(target);need(HU>=H and (HU>H)==r['scalar_strict'],'independent scalar status');strict+=HU>H
    B,a,c=lift(target,ambient);need(B==r['target_old_rows'] and [a,c]==r['query'] and targetprofiles==r['target_child_profiles'],'edge-XOR target and child profiles');G=extend(B,a,c);full=511|(1<<len(B))|(1<<(len(B)+1));need(all(G[v]&full in(0,full) for v in range(9,len(B))),'actual augmented core module');bound=(3 if ambient else 1)*P;need(r['minimum_gap']==bound and r['source_value']==source['value'] and r['target_value']-source['value']==r['gap']>=bound>0,'fixed-query strict direction and lower bound');jobs.append(G);expected.append(r['target_value']);records.append([label,ambient,bits,HU,source['value'],r['target_value'],bound])
   strictcounts[f'{label}_{ambient}']=strict
 lines=[str(len(jobs))]+[str(len(T))+' '+' '.join(map(str,T)) for T in jobs]
 with tempfile.TemporaryDirectory() as td:
  exe=Path(td)/'count';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(ROOT/'count.cpp'),'-o',str(exe)],check=True);counts=list(map(int,subprocess.check_output([str(exe)],input=('\n'.join(lines)+'\n').encode()).split()))
 need(counts==expected and len(counts)==1044,'all independent backward-DP literal counts');need(strictcounts['TT3_False']==296 and strictcounts['TT3_True']==296 and strictcounts['C3_False']==strictcounts['C3_True']==0,'complete scalar-strict census')
 return {'status':'PASS','independent_original_module_lists':basis,'physical_counts':len(counts),'family_contexts':len(records),'scalar_strict_counts':strictcounts,'record_columns':['child','cyclic_exterior','child_module_bitmask','target_core_H','source_F','target_F','minimum_gap'],'records':records,'source_pins':{f:hashlib.sha256((SOURCE/f).read_bytes()).hexdigest() for f in ['PROOF.md','verify.py','count.cpp','normal.json','optimized.json']},'scope':'Fixed source equal-profile cyclic root and all root-containing original-module XOR families; arbitrary ambient conditional on the separately certified baseline sign. Actual controls at common child order3.'}
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
