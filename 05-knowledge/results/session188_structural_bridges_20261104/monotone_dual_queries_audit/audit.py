#!/usr/bin/env python3
"""Independent predecessor-DP audit of the full physical dual-query packet."""
from pathlib import Path
from itertools import combinations,product
from functools import lru_cache
import hashlib,json,subprocess,sys,tempfile
ROOT=Path(__file__).resolve().parent;SRC=ROOT.parent/'monotone_dual_queries';BASE=ROOT.parent
PINS={'PROOF.md':'f617b74aa6c9545e945d6c5f6d8258049a9dd68f3127c5252542b2592f746d5f','verify.py':'32c6d95c2824386e5add6655c5859e902f93093c46ade399b90997e932bd9559','literal.cpp':'ab6bc0b847283f2031b3d46c1b9690731f1497c02159152476fae8fbba27b8ec','normal.json':'0152c7cc566cef65704c281b69cae9caf5e986953ee0400c106482ca945ba4a6','optimized.json':'0152c7cc566cef65704c281b69cae9caf5e986953ee0400c106482ca945ba4a6'}
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

def main():
 for name,h in PINS.items():need(hashlib.sha256((SRC/name).read_bytes()).hexdigest()==h,'source pin '+name)
 data=json.loads((SRC/'normal.json').read_text());need(data['status']=='PASS','source status')
 for name,h in data['dependency_pins'].items():need(hashlib.sha256((BASE/name).read_bytes()).hexdigest()==h,'dependency pin '+name)
 records=data['family_controls'];Q=[4,13,24,17,3];need(modules(Q)==[31],'prime named quotient');E=endpoints(Q)[-1];S=endpoints(op(Q))[-1];need((E[2]+S[2])%2==1,'quotient endpoint anchor');need({tuple(E),tuple(S)}=={tuple(data['quotient_start_counts']),tuple(data['quotient_end_counts'])},'complete endpoint counts')
 exact={};number=0
 def retain(T,a,c,v):
  nonlocal number
  G=augment(T,a,c);need(G not in exact or exact[G]==v,'same physical graph consistent count');exact[G]=v;number+=1
 dual=[]
 for r in data['all_coordinate_duality']:
  T,U=r['source_old_rows'],r['baseline_old_rows'];qs=r['core_queries'];p=[profile(augment(A,a,c)) for A in(T,U) for a,c in qs];need(p==r['complete_augmented_profiles'],'all exact augmented profile coordinates');need(p[0]==p[3] and p[1]==p[2],'augmented profile duality');dual.append(p)
 count=strict=0;family_sizes=[];negative_indices=[]
 for r in records:
  B=r['child_rows'];T,off=blow(Q,B);U0,_=blow(op(Q),B);R=r['outer_rows'];source,_=blow(R,[T]+[[0]]*(len(R)-1));baseline,_=blow(R,[U0]+[[0]]*(len(R)-1));need(source==r['source_old_rows'] and baseline==r['baseline_old_rows'],'fixed ambient source/baseline')
  n=len(T);M=(1<<n)-1;I=((1<<len(B[2]))-1)<<off[2];ext=sum(1<<(n+j-1) for j in range(1,len(R)) if R[j]&1);qs=[(I|ext,ext),(M|ext,(M^I)|ext)];need(qs==list(map(tuple,r['physical_queries'])),'only core bits complemented; exterior unchanged')
  vals=r['baseline_values'];d=[a-b for a,b in vals];need(d[0]==-d[1] and d[0]%4==2,'opposite exact nonzero ambient contrasts');chosen=0 if d[0]<0 else 1;need(chosen==r['fixed_negative_query_index'],'baseline-sign selection');negative_indices.append(chosen)
  for (a,c),(s,b) in zip(qs,vals):retain(source,a,c,s);retain(baseline,a,c,b)
  mm=[modules(b) for b in B];choices=[[tuple(m for k,m in enumerate(z) if bits>>k&1) for bits in range(1<<len(z))] for z in mm];expected=list(product(*choices));need(expected==[tuple(map(tuple,f['selected_original_child_modules'])) for f in r['families']],'all original-module families, not a sample');family_sizes.append(len(expected))
  for f,selected in zip(r['families'],expected):
   new=[op(reverse(b,sel)) for b,sel in zip(B,selected)];need(new==f['target_child_rows'],'original-module XOR target interiors');qa=[profile(tuple(b)) for b in B];qb=[profile(tuple(b)) for b in new];need(all(all(y>=x for x,y in zip(a,b)) for a,b in zip(qa,qb)),'complete child-profile monotonicity');strictflag=qa!=qb;need(strictflag==f['profile_strict'],'profile strictness')
   U,_=blow(op(Q),new);target,_=blow(R,[U]+[[0]]*(len(R)-1));actual=f['physical_values'];need(len(actual)==2,'both full actual queries')
   coremask=((1<<n)-1)|(1<<len(source))|(1<<(len(source)+1))
   for (a,c),(s,t),(bs,bt) in zip(qs,actual,vals):
    retain(source,a,c,s);retain(target,a,c,t);need(s==bs and t>=bt,'all target query monotonicity')
    for old in(source,baseline,target):
     G=augment(old,a,c);need(all(G[v]&coremask in(0,coremask) for v in range(n,len(source))),'same physical augmented core module')
   need(actual[chosen][0]-actual[chosen][1]<=d[chosen]<0,'one fixed member separates every permitted target');count+=1;strict+=strictflag
 need(count==data['families']==208 and strict==data['profile_strict_families']==76,'complete family census')
 bank=data['complete_named_bank'];T,U,V=[bank[k] for k in('source_old_rows','baseline_old_rows','strict_crossing_target_old_rows')];TT3=[6,4,0];need(T==blow(Q,[TT3]+[[0]]*4)[0] and U==blow(op(Q),[TT3]+[[0]]*4)[0],'named bank source/baseline');cross=op(reverse(TT3,[3,6]));need(V==blow(op(Q),[cross]+[[0]]*4)[0] and profile(tuple(cross))[1]>1,'literal crossing target')
 masks=[M for M in range(128) if M.bit_count()<=2];entries=bank['records'];need([(a,c) for a,c,*_ in entries]==list(product(masks,repeat=2)) and len(entries)==841,'complete named bidegree(2,2) including zeros');table={(a,c):vs for a,c,*vs in entries};hist=[0]*5
 for a,c,*vs in entries:
  for A,value in zip((T,U,V),vs[:3]):retain(A,a,c,value)
  coeff=[sum((-1)**(a.bit_count()+c.bit_count()-x.bit_count()-y.bit_count())*table[x,y][j] for x in subsets(a) for y in subsets(c)) for j in range(3)];need(coeff==vs[3:],'full exact Boolean coefficient vector');hist[a.bit_count()+c.bit_count()]+=coeff[0]!=coeff[2]
 for a,c,*vs in bank['full_mask_reconstructions']:
  calc=[sum(z[3+j] for (x,y),z in table.items() if x&~a==0 and y&~c==0) for j in range(3)];need(calc==vs,'complete polynomial arbitrary-mask reconstruction')
  for A,value in zip((T,U,V),vs):retain(A,a,c,value)
 need(hist==[bank['source_target_nonzero_histogram'][str(i)] for i in range(5)] and hist[3]==290 and hist[4]==266,'cubic quartic retained')
 # Source counted repeated source values in the two pair banks: +881 repeats.
 need(number+881==data['literal_counts']==4380,'retained physical count semantics')
 with tempfile.TemporaryDirectory(prefix='dual-query-audit-') as td:
  exe=Path(td)/'count';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(ROOT/'count.cpp'),'-o',str(exe)],check=True);graphs=sorted(exact);wire=str(len(graphs))+'\n'+'\n'.join(str(len(G))+' '+' '.join(map(str,G)) for G in graphs)+'\n';values=list(map(int,subprocess.check_output([str(exe)],input=wire.encode()).split()));need(values==[exact[G] for G in graphs],'all literal values by independent predecessor DP')
  args=[sys.executable]+(['-O'] if sys.flags.optimize else [])+[str(SRC/'verify.py')];replay=subprocess.check_output(args,cwd=td);need(replay==(SRC/'normal.json').read_bytes(),'fresh source full replay')
 return dict(status='PASS',source_sha256=PINS,dependency_sha256=data['dependency_pins'],all_augmented_profiles=dual,complete_family_sizes=family_sizes,families=count,strict_families=strict,baseline_negative_query_indices=negative_indices,full_named_coordinates=841,contrast_histogram=hist,actual_counts_including_source_repeats=number+881,distinct_actual_graphs=len(graphs),complete_exact_graph_values=[[list(G),v] for G,v in zip(graphs,values)],source_replay_matches=True)
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
