#!/usr/bin/env python3
"""Independent exact all-cover and physical audit of the pair-child ambient lift."""
from pathlib import Path
from itertools import combinations,product
from functools import lru_cache
import hashlib,json,subprocess,tempfile
ROOT=Path(__file__).resolve().parent;SRC=ROOT.parent/'ordered_pair_ambient';BASE=ROOT.parent
PINS={'PROOF.md':'40bce2a4ee3f302781949af2e9c2ad4b93f40810d43340ccf43db3bfcdb470ff','verify.py':'37ca8c79a72e56340fd3b92d6f94f5bb1db16ae458a7e43585b1a11f66339b91','count.cpp':'4f34e4aede4e9a14564f61f47fd609ae1324ee009dece48f36e4abc1eb227a5f','normal.json':'d3c716273d77425c33cd061b4d0e2c180c963d9679b79158b90cc697f78bf061','optimized.json':'d3c716273d77425c33cd061b4d0e2c180c963d9679b79158b90cc697f78bf061'}
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def op(T):return [((1<<len(T))-1)^(1<<i)^r for i,r in enumerate(T)]
def reverse(T,mm):
 R=list(T)
 for M in mm:
  for i in range(len(T)):
   if M>>i&1:R[i]^=M^(1<<i)
 return R
def blow(Q,B):
 off=[sum(map(len,B[:i])) for i in range(len(B))];R=[0]*sum(map(len,B))
 for i,C in enumerate(B):
  for v,r in enumerate(C):
   R[off[i]+v]=r<<off[i]
   for j,D in enumerate(B):
    if Q[i]>>j&1:R[off[i]+v]|=((1<<len(D))-1)<<off[j]
 return R,off
def augment(T,a,c):
 n=len(T);R=list(T)+[1<<(n+1),0]
 for i in range(n):
  if a>>i&1:R[i]|=1<<n
  else:R[n]|=1<<i
  if c>>i&1:R[i]|=1<<(n+1)
  else:R[n+1]|=1<<i
 return tuple(R)
def endpoints(T):
 n=len(T);D=[[0]*n for _ in range(1<<n)]
 for i in range(n):D[1<<i][i]=1
 for M in range(1,1<<n):
  for i in range(n):
   if M>>i&1 and M!=(1<<i):D[M][i]=sum(D[M^(1<<i)][j] for j in range(n) if T[j]>>i&1)
 return D
@lru_cache(None)
def profile(T):
 n=len(T);N=1<<n;H=list(map(sum,endpoints(T)));q=[[0]*N for _ in range(n+1)];q[0][0]=1
 for h in range(1,n+1):
  for M in range(1,N):
   part=M
   while part:q[h][M]+=H[part]*q[h-1][M^part];part=(part-1)&M
 return [q[h][-1] for h in range(n+1)]
def modules(T):return [M for M in range(1,1<<len(T)) if M.bit_count()>1 and all(T[i]&M in(0,M) for i in range(len(T)) if not M>>i&1)]
def subsets(M):
 S=M
 while True:
  yield S
  if not S:break
  S=(S-1)&M

def hp(T):return sum(endpoints(T)[-1])
def induced(T,names):return [sum(1<<j for j,w in enumerate(names) if T[v]>>w&1) for v in names]
def skew(T,X,Y):
 R=[r>>2 for r in T[2:]];X>>=2;Y>>=2;n=len(R);M=X|Y;need(all(R[v]&M in(0,M) for v in range(n) if not M>>v&1) and all(R[v]&Y==Y for v in range(n) if X>>v&1),'global ordered module')
 E=endpoints(R);S=endpoints(op(R));full=(1<<n)-1;xy=yx=0
 for A in range(1,full):
  xy+=sum(E[A][v] for v in range(n) if X>>v&1)*sum(S[full^A][v] for v in range(n) if Y>>v&1)
  yx+=sum(E[A][v] for v in range(n) if Y>>v&1)*sum(S[full^A][v] for v in range(n) if X>>v&1)
 qx=profile(tuple(induced(R,[v for v in range(n) if X>>v&1])));qy=profile(tuple(induced(R,[v for v in range(n) if Y>>v&1])));reps=[(M&-M).bit_length()-1]+[v for v in range(n) if not M>>v&1];lower=hp(induced(R,reps))*(qx[2]*qy[1]+qx[1]*qy[2]+qx[2]*qy[2]);need(yx-xy>=lower>0,'exact skew and positive lower bound');return yx-xy,lower

def main():
 for f,h in PINS.items():need(hashlib.sha256((SRC/f).read_bytes()).hexdigest()==h,'frozen source '+f)
 d=json.loads((SRC/'normal.json').read_text());need(d['status']=='PASS','source status')
 for f,h in d['source_sha256'].items():need(hashlib.sha256((BASE/f).read_bytes()).hexdigest()==h,'bound source/dependency '+f)
 exact={};cover=[]
 def retain(T,a,c,v):
  G=augment(T,a,c);need(G not in exact or exact[G]==v,'exact physical consistency');exact[G]=v
 for r in d['all_cover_restoration_controls']:
  T,U,W=[r[k] for k in('source_rows','quotient_rows','local_rows')];need(W==reverse(T,[3]),'local pair target');coeff=[[0]*(len(T)+3) for _ in range(3)]
  for q in r['corners']:
   for j,(key,V) in enumerate(zip(('source','quotient','local'),(T,U,W))):
    p=profile(augment(V,q['a'],q['c']));need(p==q[key],'every complete ordered-cover coordinate')
    for h,z in enumerate(p):coeff[j][h]+=q['sign']*z
  delta=[a-b for a,b in zip(coeff[0],coeff[1])];need(delta==r['all_cover_coefficients']==[a-b for a,b in zip(coeff[0],coeff[2])] and delta[1]>0,'full-h exact restoration');cover.append(coeff)
 for r in d['physical_ambient_quartics']:
  T,U,W=[r[k] for k in('source_rows','quotient_rows','local_rows')];need(W==reverse(T,[3]),'ambient local pair target');gap,lower=skew(T,r['X'],r['Y']);need(gap==r['skew']==r['actual_coefficient'] and lower==r['lower_bound'],'full exterior endpoint skew');K=r['outer'];n=len(T)-len(K)+1;ext=sum(1<<(n+j-1) for j in range(1,len(K)) if K[j]&1);need(ext==r['outside_bits'],'fixed root incidences')
  coremask=((1<<n)-1)|(1<<len(T))|(1<<(len(T)+1));tot=[0]*3
  for bits,q in enumerate(r['corners']):
   a=ext|(2 if bits&1 else 0)|(r['X'] if bits&2 else 0);c=ext|(1 if bits&4 else 0)|(r['Y'] if bits&8 else 0);sg=(-1)**bits.bit_count();need((a,c,sg)==(q['a'],q['c'],q['sign']),'saturated physical finite difference')
   for j,(key,V) in enumerate(zip(('source','quotient','local'),(T,U,W))):
    G=augment(V,a,c);need(all(G[v]&coremask in(0,coremask) for v in range(n,len(T))),'same augmented ambient module');retain(V,a,c,q[key]);tot[j]+=sg*q[key]
  need(tot[0]-tot[1]==tot[0]-tot[2]==gap,'exact ambient master equality')
 branches=dict(strict=0,pair_selected=0,pair_unselected=0);families=0
 for b in d['root_family_contexts']:
  T=b['core_rows'];n=len(T);K=b['outer'];source,_=blow(K,[T]+[[0]]*(len(K)-1));need(source==b['source_rows'],'actual root substitution');mods=[M for M in modules(T) if M!=(1<<n)-1];need(mods==b['original_proper_modules'],'complete original modules');selections=[[M for j,M in enumerate(mods) if bits>>j&1] for bits in range(1<<len(mods))];need(selections==[r['selected_modules'] for r in b['records']],'complete original family cube')
  ext=sum(1<<(n+j-1) for j in range(1,len(K)) if K[j]&1);need(ext==b['outside_bits'],'unchanged physical exterior incidence');qs=list(map(tuple,b['battery']));at={q:i for i,q in enumerate(qs)};base=b['baseline_values'];X=b['quartet'][2][0];Y=b['quartet'][8][1];gap,_=skew(source,X,Y);need(gap==b['expected_selected_pair_coefficient'],'expected ambient selected coefficient')
  # Independently recover child sets from the maximal proper modules.
  maximal=[M for M in mods if not any(M!=N and M&N==M for N in mods)];maximal.sort(key=lambda M:(M&-M).bit_length());parts=[2]+maximal[1:];out=sum(bool(T[0]&M) for M in parts);anchor=[(M|(1 if T[0]&M else 0),1) for M in parts]
  if out%2:anchor.append((1,1))
  if (len(parts)+out)%2:anchor.append((0,1))
  anchor.append((1,0));quart=[((2 if k&1 else 0)|(X if k&2 else 0),(1 if k&4 else 0)|(Y if k&8 else 0),(-1)**k.bit_count()) for k in range(16)];expected=list(dict.fromkeys([(a,c) for a,c,s in quart]+anchor));need(anchor==list(map(tuple,b['anchor'])) and quart==list(map(tuple,b['quartet'])) and expected==qs and len(qs)<=len(maximal)+18,'physical anchor/battery and odd-prime bound')
  for (a,c),v in zip(qs,base):retain(source,a|ext,c|ext,v)
  H=hp(T)
  for r,sel in zip(b['records'],selections):
   U=reverse(op(T),sel);V,_=blow(K,[U]+[[0]]*(len(K)-1));need(U==r['target_core_rows'] and V==r['target_rows'],'physical original XOR target');h=hp(U);need(h==r['core_H'] and H==r['source_core_H'] and h>=H,'exact core scalar comparison');vals=r['values'];need(vals!=base and len(vals)==len(qs),'complete integer response vector separation');need(r['pair_selected']==(3 in sel),'selected-pair parity')
   for (a,c),v in zip(qs,vals):retain(V,a|ext,c|ext,v)
   if h>H:branch='strict';need(vals[at[(0,0)]]>base[at[(0,0)]],'lifted internal-zero strictness')
   elif 3 in sel:branch='pair_selected';z=sum(s*(base[at[(a,c)]]-vals[at[(a,c)]]) for a,c,s in quart);need(z==r['coefficient']==gap>0,'equal-profile selected-pair quartic')
   else:branch='pair_unselected';res=[(base[at[q]]-vals[at[q]])%4 for q in anchor];need(res==r['anchor_residues'] and 2 in res,'fixed-ambient exposed anchor')
   need(branch==r['branch'],'branch label');branches[branch]+=1;families+=1
 need(families==d['family_count']==128 and branches==d['branches']==dict(strict=16,pair_selected=56,pair_unselected=56),'all family branches');need(len(exact)==d['distinct_actual_augmented_graphs']==3142,'complete literal graph census')
 with tempfile.TemporaryDirectory(prefix='ordered-pair-audit-') as td:
  exe=Path(td)/'count';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(ROOT/'count.cpp'),'-o',str(exe)],check=True);G=sorted(exact);wire=str(len(G))+'\n'+'\n'.join(str(len(T))+' '+' '.join(map(str,T)) for T in G)+'\n';vals=list(map(int,subprocess.check_output([str(exe)],input=wire.encode()).split()));need(vals==[exact[T] for T in G],'every literal value by independent predecessor DP')
 # Excluded |X union Y|=2 boundary: exact zero even in a strong outer context.
 core,_=blow([2,4,1],[[2,0],[0],[0]]);quot,_=blow([4,1,2],[[2,0],[0],[0]]);boundary_old=[blow([2,4,1],[V,[0],[0]])[0] for V in(core,quot,reverse(core,[3]))];boundary=[];co=[0,0,0]
 for bits in range(16):
  a=32|(2 if bits&1 else 0)|(4 if bits&2 else 0);c=32|(1 if bits&4 else 0)|(8 if bits&8 else 0);z=[hp(augment(V,a,c)) for V in boundary_old];boundary.append([a,c,*z])
  for j,v in enumerate(z):co[j]+=(-1)**bits.bit_count()*v
 need(co[0]==co[1]==co[2],'two-singleton-side strictness exclusion')
 return dict(status='PASS',excluded_two_singleton_boundary=boundary,source_sha256=PINS,dependencies=d['source_sha256'],complete_cover_coefficients=cover,actual_graphs=len(exact),families=families,branches=branches,complete_exact_graph_values=[[list(T),v] for T,v in zip(G,vals)])
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
