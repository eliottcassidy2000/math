#!/usr/bin/env python3
"""Source-bound independent physical and algebra audit. No producer code imports."""
from pathlib import Path
from itertools import combinations
import hashlib,json,subprocess,sys,tempfile
import independent as I
ROOT=Path(__file__).resolve().parent;SRC=ROOT.parent/'ordered_cyclic_module';BASE=ROOT.parent
PINS={
 'PROOF.md':'d528ad7a713b16bb88971466bf1d7ba1f4e2260a1efe7abeafb790ddbcea01cf',
 'README.md':'c2b754e1f05b7112a74c8f434567af71384acb843292a50d67d90013f2eb80c5',
 'verify.py':'217687f959b620d0393ee9c5f09fe71e62a9c8026b73d503a9e61c40c38a72ff',
 'algebra.py':'34b15db6cd33f4e560c536435784c999b68c26a1358676e823021dd07fbed2c6',
 'kernel.py':'6d34008f36ce92ab7c01c3655ff9c2232065246e8ab2de4136aef201bfa75061',
 'count.cpp':'4f34e4aede4e9a14564f61f47fd609ae1324ee009dece48f36e4abc1eb227a5f',
 'normal.json':'4429f0957d66c0737fdde2e7ab3672bd64beea2d5c74382e97c06d099ab029d3',
 'optimized.json':'4429f0957d66c0737fdde2e7ab3672bd64beea2d5c74382e97c06d099ab029d3'}
need=I.need

def flip(T,M):return [r^(M^(1<<i) if M>>i&1 else 0) for i,r in enumerate(T)]
def augment(T,a,c):
 n=len(T);R=list(T)+[1<<(n+1),0]
 for v in range(n):
  if a>>v&1:R[v]|=1<<n
  else:R[n]|=1<<v
  if c>>v&1:R[v]|=1<<(n+1)
  else:R[n+1]|=1<<v
 return tuple(R)
def blow(Q,B):
 off=[sum(map(len,B[:i])) for i in range(len(B))];T=[0]*sum(map(len,B))
 for i,C in enumerate(B):
  for u,r in enumerate(C):
   T[off[i]+u]=r<<off[i]
   for j,D in enumerate(B):
    if Q[i]>>j&1:T[off[i]+u]|=((1<<len(D))-1)<<off[j]
 return T,off
def cyclic(r):
 n=2*r+1;return [sum(1<<j for j in range(n) if 1<=(j-i)%n<=r) for i in range(n)]
def ufamily(r):
 n=2*r+1;T=[0]*n
 for u,v in combinations(range(n),2):
  if v<r:a,b=v,u
  elif u>=r:a,b=u,v
  elif v-r<=u:a,b=u,v
  else:a,b=v,u
  T[a]|=1<<b
 return T
def wfamily(r):
 n=2*r+1;T=[0]*n
 for u,v in combinations(range(n),2):
  if u==0:a,b=(0,v) if v%2 else(v,0)
  else:a,b=u,v
  T[a]|=1<<b
 return T
def ismodule(T,M):return all(T[v]&M in(0,M) for v in range(len(T)) if not M>>v&1)
def tournament(T):
 n=len(T);need(all(0<=r<1<<n and not r>>i&1 for i,r in enumerate(T)),'physical names in range')
 need(all(bool(T[i]>>j&1)!=bool(T[j]>>i&1) for i,j in combinations(range(n),2)),'literal tournament')

def main():
 for name,h in PINS.items():need(hashlib.sha256((SRC/name).read_bytes()).hexdigest()==h,'frozen source '+name)
 data=json.loads((SRC/'normal.json').read_text());need(data['status']=='PASS','producer status')
 for name,h in data['source_sha256'].items():need(hashlib.sha256((BASE/name).read_bytes()).hexdigest()==h,'dependency source '+name)
 # Fresh full independent algebra/physical controls, including non-VT balanced marks.
 ours=I.main();raw=(json.dumps(ours,sort_keys=True,indent=2)+'\n').encode()
 need(raw==(ROOT/'independent.normal.json').read_bytes()==(ROOT/'independent.optimized.json').read_bytes(),'full independent paired replay')
 graph_values={}
 def retain(T,a,c,value):
  tournament(T);G=augment(T,a,c)
  need(G not in graph_values or graph_values[G]==value,'identical graph same exact value')
  graph_values[G]=value
 # Match producer's entire arbitrary-exterior census by physical incidences, not hashes.
 controls={(tuple(z['E']),((1<<len(z['E']))-1)^z['IA'],((1<<len(z['E']))-1)^z['IM']):z for z in ours['records'] if z['tag']=='arbitrary_exterior'}
 need(len(controls)==549==len(data['small_arbitrary_exterior_controls']),'complete small census')
 for c in data['small_arbitrary_exterior_controls']:
  z=controls[(tuple(c['E']),c['Aout'],c['Mout'])];need(c['coefficient']==z['actual_coefficient'] and c['lower_bound']==z['strict_lower_bound'],'all small interface values')
 for T,b in data['balanced_profiles']:need(I.marked(T,0)==b,'independent ordered-cover profiles')
 need(len(data['literal_positive_controls'])==21,'literal control census')
 for c in data['literal_positive_controls']:
  E=c['E'];allE=(1<<len(E))-1;T,U,ox,oy=I.build(E,allE^c['Aout'],allE^c['Mout'],c['X'],c['Y']);need(T==c['T'] and U==c['local_target'],'physical local graph construction')
  blocks=[[2,4,1],c['X'],c['Y']]+[[0]]*len(E);root,off=blow(I.op(c['Q']),blocks);need(root==c['quotient_target'],'fixed physical quotient converse')
  need(c['physical_marks']==[0,ox,1,oy],'physical marked names')
  pred,lower,bx,by=I.positive_value(E,allE^c['Aout'],allE^c['Mout'],c['X'],c['Y'],0,0)
  need(pred==c['positive_formula']['coefficient']==sum(c['kernel_terms']) and lower==c['positive_formula']['lower_bound'],'independent complete positive functional')
  for key,target in [('local_queries',U),('quotient_queries',root)]:
   need(len(c[key])==16,'complete physical quartic finite difference');gap=0
   for bits,q in enumerate(c[key]):
    a=(1 if bits&1 else 0)|(1<<ox if bits&2 else 0);cc=(2 if bits&4 else 0)|(1<<oy if bits&8 else 0);sign=(-1)**bits.bit_count()
    need((a,cc,sign)==(q['a'],q['c'],q['sign']),'physical quartic mask/sign');retain(T,a,cc,q['source']);retain(target,a,cc,q['target']);gap+=sign*(q['source']-q['target'])
   need(gap==pred>0,'exact local/quotient quartic')
 funcs={'T':cyclic,'U':ufamily,'W':wfamily};structures=data['critical_witness_structures'];need(len(structures)==45,'critical batches')
 for z in structures:
  r=z['r'];Q=funcs[z['family']](r);need(len(Q)==z['order']==2*r+1,'critical order')
  expected=[(0,r,r+1)] if z['family']=='T' else([(i,r+i,r+i+1) for i in range(r)] if z['family']=='U' else[(0,2*s+1,2*s+2) for s in range(r)])
  need(expected==list(map(tuple,z['witnesses'])),'complete named witnesses')
  for a,b,c in expected:
   need(Q[a]>>b&1 and Q[b]>>c&1 and Q[c]>>a&1,'cyclic witness');need(all(bool(Q[t]>>b&1)==bool(Q[t]>>c&1) for t in range(len(Q)) if t not in(a,b,c)),'deleted-root module')
 family_count=0
 for rec in data['root_families']:
  Q=rec['Q'];m=len(Q);T,off=blow(Q,[[2,4,1]]*m);need(T==rec['T'],'physical all-C3 inflation');n=len(T);mods=[M for M in range(1,(1<<n)-1) if M.bit_count()>1 and ismodule(T,M)];blocks=[7<<(3*t) for t in range(m)];need(mods==blocks,'all original proper modules')
  z,u,v=rec['witness'];p=3*z;q=p+1;i=3*u;j=3*v
  quart=[((1<<p if bits&1 else 0)|(1<<i if bits&2 else 0),(1<<q if bits&4 else 0)|(1<<j if bits&8 else 0),(-1)**bits.bit_count()) for bits in range(16)]
  need(quart==list(map(tuple,rec['quartic'])),'root physical quartic')
  compressed=[1<<(p+1),1<<(p+2)]+[blocks[t] for t in range(m) if t!=z];d=sum(bool(T[p]&B) for B in compressed);anchor=[(B|((1<<p) if T[p]&B else 0),1<<p) for B in compressed]
  if d%2:anchor.append((1<<p,1<<p))
  if (len(compressed)+d)%2:anchor.append((0,1<<p))
  anchor.append((1<<p,0));battery=list(dict.fromkeys([(a,c) for a,c,s in quart]+anchor));need(anchor==list(map(tuple,rec['anchor'])) and battery==list(map(tuple,rec['battery'])) and len(battery)<=m+20,'exact named anchor/battery')
  base=rec['baseline_values'];need(len(base)==len(battery),'complete baseline');at={q:i for i,q in enumerate(battery)}
  for (a,c),v0 in zip(battery,base):retain(T,a,c,v0)
  need([f['selected_children_bits'] for f in rec['families']]==list(range(1<<m)),'complete original root family cube')
  for f in rec['families']:
   bits=f['selected_children_bits'];target=I.op(T)
   for t,B in enumerate(blocks):
    if bits>>t&1:target=flip(target,B)
   need(target==f['target_rows'],'original-module XOR semantics');vals=f['values'];need(len(vals)==len(battery),'complete target values')
   for (a,c),v0 in zip(battery,vals):retain(target,a,c,v0)
   need(vals!=base and vals[at[(0,0)]]==base[at[(0,0)]],'full exact vectors separate scalar-equal root family')
   if bits>>z&1:
    gap=sum(s*(base[at[(a,c)]]-vals[at[(a,c)]]) for a,c,s in quart);need(gap==f['quartic_gap']>0,'restored triangle positive quartic')
   else:
    residues=[(base[at[q0]]-vals[at[q0]])%4 for q0 in anchor];need(residues==f['anchor_residues'] and 2 in residues,'unrestored triangle literal anchor')
   family_count+=1
 need(family_count==104 and len(graph_values)==data['distinct_literal_augmented_tournaments']==3466,'full family and graph counts')
 with tempfile.TemporaryDirectory(prefix='ordered-cyclic-audit-') as td:
  exe=Path(td)/'independent';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(ROOT/'replay.cpp'),'-o',str(exe)],check=True)
  graphs=sorted(graph_values);wire=str(len(graphs))+'\n'+'\n'.join(str(len(G))+' '+' '.join(map(str,G)) for G in graphs)+'\n';actual=list(map(int,subprocess.check_output([str(exe)],input=wire.encode()).split()));need(actual==[graph_values[G] for G in graphs],'all retained exact physical values by independent predecessor DP')
  args=[sys.executable]+(['-O'] if sys.flags.optimize else [])+[str(SRC/'verify.py')];replay=subprocess.check_output(args,cwd=td);need(replay==(SRC/'normal.json').read_bytes(),'fresh source replay matches complete retained data')
 return dict(status='PASS',source_sha256=PINS,dependency_sha256=data['source_sha256'],independent_full_controls_sha256=hashlib.sha256(raw).hexdigest(),independent_physical_values=ours['literal_physical_values'],independent_extra_balanced_controls=48,independent_excluded_transitive_control=1,source_exterior_controls=549,source_local_quotient_controls=21,critical_batches=45,root_families=family_count,retained_exact_graph_count=len(graph_values),retained_exact_graph_values=[[list(G),v] for G,v in zip(graphs,actual)],source_replay_matches=True)
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
