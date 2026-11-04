"""Complete symbolic run certificate and literal named controls; no sampled proof step."""
from pathlib import Path
from itertools import combinations,product
import hashlib,json,subprocess,tempfile
ROOT=Path(__file__).resolve().parent
Q=[238,4,240,118,195,146,34,74]
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def op(R):return [((1<<len(R))-1)^r^(1<<i)for i,r in enumerate(R)]
def tt(n):return [sum(1<<j for j in range(i+1,n))for i in range(n)]
C3=[2,4,1]
def blow(Q,blocks):
 off=[];n=0
 for B in blocks:off.append(n);n+=len(B)
 R=[]
 for i,B in enumerate(blocks):
  ext=sum(((1<<len(C))-1)<<off[j]for j,C in enumerate(blocks)if Q[i]>>j&1)
  R += [(r<<off[i])|ext for r in B]
 return R,off
def expand(mask,blocks):
 off=0;out=0
 for i,B in enumerate(blocks):
  if mask>>i&1:out|=((1<<len(B))-1)<<off
  off+=len(B)
 return out
def aug(R,a,c):
 n=len(R);out=list(R)+[1<<(n+1),0]
 for i in range(n):
  for v,m in [(n,a),(n+1,c)]:
   if m>>i&1:out[i]|=1<<v
   else:out[v]|=1<<i
 return out
def physical_query(core,blocks,a,c,K=None):
 aa,cc=expand(a,blocks),expand(c,blocks)
 if K is None:return aug(core,aa,cc)
 R,off=blow(K,[core]+[tt(1)]*(len(K)-1));incoming=sum(1<<off[i]for i in range(1,len(K))if K[i]&1)
 return aug(R,aa|incoming,cc|incoming)
def poly(rows,h):
 out={}
 for row in rows:
  s=row[2];d=[(s>>(2*i))&3 for i in range(8)];cur={():row[3+2*h]-row[4+2*h]}
  for i,r in enumerate(d):
   if r==1:continue
   factor={():2,(2*i,):4}if r==2 else{(2*i,):4,(2*i+1,):2};nxt={}
   for m,z in cur.items():
    for b,w in factor.items():key=tuple(sorted(m+b));nxt[key]=nxt.get(key,0)+z*w
   cur=nxt
  for m,z in cur.items():out[m]=(out.get(m,0)+z)%8
 return [[list(m),z]for m,z in sorted(out.items())if z]
def kernel_value(rows,h,profiles,side):
 total=0
 for row in rows:
  d=[(row[2]>>(2*i))&3 for i in range(8)];z=row[3+2*h+side]
  for q,k in zip(profiles,d):z*=q[k]if k<len(q)else 0
  total+=z
 return total%8
def blockq(B):
 if B==C3:return [0,3,6,6]
 n=len(B);S=[[0]*(n+1)for _ in range(n+1)];S[0][0]=1
 import math
 for i in range(1,n+1):
  for k in range(1,i+1):S[i][k]=S[i-1][k-1]+k*S[i-1][k]
 return [0]+[math.factorial(k)*S[n][k]for k in range(1,n+1)]
def reverse_inside(R,M):
 return [r^((M^(1<<i))if M>>i&1 else 0)for i,r in enumerate(R)]
def main():
 with tempfile.TemporaryDirectory()as td:
  ex={}
  for f in ['covers','physical']:
   exe=Path(td)/f;subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(ROOT/(f+'.cpp')),'-o',str(exe)],check=True);ex[f]=exe
  kr=[list(map(int,line.split()))for line in subprocess.check_output([str(ex['covers'])]).decode().splitlines()]
  rows={q:[r for r in kr if r[:2]==list(q)]for q in [(1,0),(255,254)]};states=[]
  for s in range(65536):
   d=[(s>>(2*i))&3 for i in range(8)]
   if min(d)>0 and sum(x>=2 for x in d)<=2:states.append(s)
  need(len(states)==129 and all([r[2]for r in rs]==states for rs in rows.values()),'complete relevant run domain')
  symbolic=[]
  for q,rs in rows.items():
   for h in range(3):
    p=poly(rs,h);need(p==([[[],4]]if h==0 else[]),'universal symbolic profile contrast');symbolic.append({'query':q,'cover_slots':h+1,'nonzero_residue_coefficients':p})
  graphs=[];meta=[]
  def add(R,cov,label):graphs.append((R,cov));meta.append(label);return len(graphs)-1
  small=[];B=[tt(1)]*8;small.append(('singletons',B))
  for i in range(8):
   for name,A in [('pair',tt(2)),('TT3',tt(3)),('C3',C3)]:
    C=list(B);C[i]=A;small.append((f'{name}_{i}',C))
  for i,j in combinations(range(8),2):
   C=list(B);C[i]=C[j]=tt(2);small.append((f'pairs_{i}_{j}',C))
  smallrec=[]
  for name,blocks in small:
   T,_=blow(Q,blocks);U,_=blow(op(Q),blocks);profiles=list(map(blockq,blocks));ids=[]
   for a,c in rows:
    for side,R in enumerate([T,U]):ids.append(add(physical_query(R,blocks,a,c),1,[name,a,c,side]))
   smallrec.append({'name':name,'blocks':blocks,'profiles':profiles,'ids':ids})
  large=[]
  for name,changes,K in [('all_pairs',{},None),('C3_child0',{0:C3},None),('TT3_child3',{3:tt(3)},None),('TT4_child4',{4:tt(4)},None),('mixed_odd',{1:C3,7:tt(3)},None),('C3_ambient',{},C3),('ordered_ambient',{},tt(2))]:
   blocks=[changes.get(i,tt(2))for i in range(8)];T,_=blow(Q,blocks);U,_=blow(op(Q),blocks);ids=[]
   for a,c in rows:
    for side,R in enumerate([T,U]):ids.append(add(physical_query(R,blocks,a,c,K),0,[name,a,c,side]))
   large.append({'name':name,'blocks':blocks,'ambient':K,'ids':ids})
  blocks=[tt(2)]*8;T,_=blow(Q,blocks);U,_=blow(op(Q),blocks);batteries=[]
  for name,roles in [('shared_cubic',[('a',1),('a',4),('c',4)]),('quartic',[('a',1),('a',4),('c',0),('c',6)])]:
   ids=[];qs=[]
   for bits in range(1<<len(roles)):
    a=sum(1<<v for k,(kind,v)in enumerate(roles)if kind=='a'and bits>>k&1);c=sum(1<<v for k,(kind,v)in enumerate(roles)if kind=='c'and bits>>k&1);qs.append([a,c]);ids.append([add(aug(R,a,c),0,[name,bits,side])for side,R in enumerate([T,U])])
   batteries.append({'name':name,'roles':roles,'queries':qs,'ids':ids})
  family=[];blocks=[tt(3)]+[tt(2)]*7;T,_=blow(Q,blocks);baseline,_=blow(op(Q),blocks);sourceids=[add(physical_query(R,blocks,a,c),0,['family_baseline',a,c,side])for a,c in rows for side,R in enumerate([T,baseline])]
  for bits in range(8):
   V=op(T)
   for j,M in enumerate([3,6,7]):
    if bits>>j&1:V=reverse_inside(V,M)
   ids=[add(physical_query(V,blocks,a,c),0,['root_family',bits,a,c])for a,c in rows]
   family.append({'original_modules':[3,6,7],'selected_bits':bits,'target_rows':V,'ids':ids})
  data=str(len(graphs))+'\n'+''.join(f'{len(R)} {cov}\n'+ ' '.join(map(str,R))+'\n'for R,cov in graphs)
  vals=[list(map(int,line.split()))for line in subprocess.check_output([str(ex['physical'])],input=data.encode()).decode().splitlines()];need(len(vals)==len(graphs),'physical row count')
 for rec in smallrec:
  actual=[vals[i]for i in rec.pop('ids')];expected=[]
  for q,rs in rows.items():
   for side in [0,1]:expected.append([kernel_value(rs,h,rec['profiles'],side)for h in range(3)])
  need([[z%8 for z in v]for v in actual]==expected,'literal all-cover model versus skeleton');rec['values']=actual
 for rec in large:
  actual=[vals[i][0]for i in rec.pop('ids')];need((actual[0]-actual[1])%8==(actual[2]-actual[3])%8==4,'large physical/ambient separator');need(actual[0]-actual[1]==actual[3]-actual[2],'exact dual baseline');rec['values']=actual
 for rec in batteries:
  actual=[[vals[i][0]for i in pair]for pair in rec.pop('ids')];degree=len(rec['roles']);coef=[sum((-1)**(degree-b.bit_count())*v[side]for b,v in enumerate(actual))for side in range(2)];rec['values']=actual;rec['coefficients']=coef;rec['contrast']=coef[0]-coef[1]
 need(all(rec['contrast']!=0 for rec in batteries),'retained nonzero cubic and quartic controls')
 base=[vals[i][0]for i in sourceids];neg=0 if base[0]<base[1]else 1;need(base[2*neg]<base[2*neg+1],'negative dual choice')
 for rec in family:
  actual=[vals[i][0]for i in rec.pop('ids')];need(base[2*neg]<actual[neg],'root family monotone separator');rec['values']=actual
 deps=['five_profile_types/PROOF.md','mod8_transfer/PROOF.md','monotone_dual_queries/PROOF.md','module_decomposition/PROOF.md','even_endpoint_obstruction/PROOF.md']
 return {'status':'PASS','quotient':Q,'kernel_columns':['a','c','base4_run_count_code','q1_source','q1_target','q2_source','q2_target','q3_source','q3_target'],'full_kernels':kr,'symbolic_identity':symbolic,'small_all_cover_models':smallrec,'large_physical_controls':large,'full_named_batteries':batteries,'root_family_baseline_values':base,'root_family_negative_query_index':neg,'root_family_controls':family,'physical_graph_occurrences':len(graphs),'distinct_physical_graphs':len({(tuple(R),cov)for R,cov in graphs}),'full_physical_values':vals,'source_pins':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest()for f in ['PROOF.md','verify.py','covers.cpp','physical.cpp']},'dependency_pins':{f:hashlib.sha256((ROOT.parent/f).read_bytes()).hexdigest()for f in deps},'scope':'One explicit prime-eight quotient; arbitrary nontrivial children and fixed ambient; no arbitrary-prime theorem or response collision.'}
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
