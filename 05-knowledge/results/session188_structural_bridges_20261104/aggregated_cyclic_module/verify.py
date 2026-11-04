#!/usr/bin/env python3
"""Aggregate-mark theorem: exact profile functional and literal named responses."""
from pathlib import Path
from itertools import combinations
import hashlib,json,subprocess,tempfile
from run_algebra import symbolic_identity,positive
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent

def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def rows(n,code):
 T=[0]*n
 for k,(i,j) in enumerate(combinations(range(n),2)):
  a,b=(i,j) if code>>k&1 else(j,i);T[a]|=1<<b
 return T
def op(T):return [((1<<len(T))-1)^(1<<i)^r for i,r in enumerate(T)]
def flip(T,M):return [r^(M^(1<<i) if M>>i&1 else 0) for i,r in enumerate(T)]
def module(T,M):return all(T[v]&M in(0,M) for v in range(len(T)) if not M>>v&1)
def covers(T):
 n=len(T);N=1<<n;d=[[0]*n for _ in range(N)]
 for i in range(n):d[1<<i][i]=1
 for M in range(1,N):
  for v in range(n):
   if not M>>v&1:continue
   rest=M^(1<<v)
   if rest:d[M][v]=sum(d[rest][u] for u in range(n) if T[u]>>v&1)
 H=list(map(sum,d));q=[[0]*N for _ in range(n+1)];q[0][0]=1
 for h in range(1,n+1):
  for M in range(1,N):
   S=M
   while S:q[h][M]+=H[S]*q[h-1][M^S];S=(S-1)&M
 return [q[h][-1] for h in range(n+1)]
def blow(Q,blocks):
 off=[sum(map(len,blocks[:i])) for i in range(len(blocks))];T=[0]*sum(map(len,blocks))
 for i,B in enumerate(blocks):
  for v,r in enumerate(B):
   T[off[i]+v]=r<<off[i]
   for j,C in enumerate(blocks):
    if Q[i]>>j&1:T[off[i]+v]|=((1<<len(C))-1)<<off[j]
 return T,off
def augment(T,a,c):
 n=len(T);U=list(T)+[1<<(n+1),0]
 for v in range(n):
  if a>>v&1:U[v]|=1<<n
  else:U[n]|=1<<v
  if c>>v&1:U[v]|=1<<(n+1)
  else:U[n+1]|=1<<v
 return tuple(U)
def setup(E,ao,mo):
 Q=[2,4,1]+[r<<3 for r in E]
 for j in range(len(E)):
  if ao>>j&1:Q[0]|=1<<(j+3)
  else:Q[j+3]|=1
  if mo>>j&1:Q[1]|=1<<(j+3);Q[2]|=1<<(j+3)
  else:Q[j+3]|=6
 return Q
def family(kind,r=2):
 n=2*r+1
 if kind=='T':return [sum(1<<((i+d)%n) for d in range(1,r+1)) for i in range(n)]
 if kind=='W':return [sum(1<<v for v in range(n) if(i and v and i<v)or(i==0 and v and v%2)or(v==0 and i and i%2==0)) for i in range(n)]
 Q=[0]*n
 for i,j in combinations(range(n),2):
  a,b=(j,i) if j<r else((i,j) if i>=r or j-r<=i else(j,i));Q[a]|=1<<b
 return Q
def battery(T,block_masks,z,u,v):
 Z=block_masks[z];p=(Z&-Z).bit_length()-1;q=next(t for t in range(len(T)) if Z>>t&1 and T[p]>>t&1);X=block_masks[u];Y=block_masks[v]
 quartet=[((1<<p if bits&1 else 0)|(X if bits&2 else 0),(1<<q if bits&4 else 0)|(Y if bits&8 else 0),(-1)**bits.bit_count()) for bits in range(16)]
 parts=[1<<t for t in range(len(T)) if Z>>t&1 and t!=p]+[M for i,M in enumerate(block_masks) if i!=z];d=sum(bool(T[p]&M) for M in parts);anchor=[(M|((1<<p) if T[p]&M else 0),1<<p) for M in parts]
 if d%2:anchor.append((1<<p,1<<p))
 if (len(parts)+d)%2:anchor.append((0,1<<p))
 anchor.append((1<<p,0));qs=list(dict.fromkeys([(a,c) for a,c,s in quartet]+anchor));need(len(qs)<=len(block_masks)+20,'m+20 root battery')
 return quartet,anchor,qs

def main():
 identity=symbolic_identity();C3=[2,4,1];one=[0];TT2=[2,0];TT3=[6,4,0];graphs=[];gid={}
 def value(T,a,c):
  G=augment(T,a,c)
  if G not in gid:gid[G]=len(graphs);graphs.append(G)
  return gid[G]
 def quartic_record(T,U,X,Y,p=0,q=1):
  return [{'a':(1<<p if k&1 else 0)|(X if k&2 else 0),'c':(1<<q if k&4 else 0)|(Y if k&8 else 0),'sign':(-1)**k.bit_count(),'source_id':value(T,(1<<p if k&1 else 0)|(X if k&2 else 0),(1<<q if k&4 else 0)|(Y if k&8 else 0)),'target_id':value(U,(1<<p if k&1 else 0)|(X if k&2 else 0),(1<<q if k&4 else 0)|(Y if k&8 else 0))} for k in range(16)]
 small=[rows(n,code) for n in range(1,4) for code in range(1<<(n*(n-1)//2))];specs=[([],0,0,X,Y) for X in small for Y in small]
 for n in range(3):
  for code in range(1<<(n*(n-1)//2)):
   for ao in range(1<<n):
    for mo in range(1<<n):specs.append((rows(n,code),ao,mo,TT2,one))
 for code in range(2):
  for ao in range(4):
   for mo in range(4):specs.append((rows(2,code),ao,mo,TT3,TT2))
 physical=[]
 for E,ao,mo,X,Y in specs:
  Q=setup(E,ao,mo);T,off=blow(Q,[C3,X,Y]+[one]*len(E));UL=flip(T,7);UQ,_=blow(op(Q),[C3,X,Y]+[one]*len(E));bx,by=covers(X),covers(Y);f=positive(E,ao,mo,bx,by);XM=((1<<len(X))-1)<<off[1];YM=((1<<len(Y))-1)<<off[2]
  physical.append(dict(E=E,Aout=ao,Mout=mo,X=X,Y=Y,qX=bx,qY=by,T=T,local_target=UL,quotient_target=UQ,positive_formula=f,local_queries=quartic_record(T,UL,XM,YM),quotient_queries=quartic_record(T,UQ,XM,YM)))
 # Complete original root-containing families with non-VT marked witness children.
 roots=[]
 for label,Q,z,u,v in [('C3',C3,0,1,2),('T5',family('T'),0,2,3),('U5',family('U'),0,2,3),('W5',family('W'),0,1,2)]:
  B=[one for _ in Q];B[z]=C3;B[u]=TT3;B[v]=TT2;T,off=blow(Q,B);n=len(T);masks=[((1<<len(b))-1)<<o for b,o in zip(B,off)];proper=[M for M in range(1,(1<<n)-1) if M.bit_count()>1 and module(T,M)];quart,anchor,qs=battery(T,masks,z,u,v);rec=dict(label=label,Q=Q,blocks=B,block_masks=masks,witness=[z,u,v],T=T,proper_original_modules=proper,quartic=quart,anchor=anchor,battery=qs,baseline_ids=[value(T,a,c) for a,c in qs],families=[])
  for selection in range(1<<len(proper)):
   U=op(T)
   for k,M in enumerate(proper):
    if selection>>k&1:U=flip(U,M)
   child_profiles=[]
   for old,o in zip(B,off):
    b=[(U[o+i]>>o)&((1<<len(old))-1) for i in range(len(old))];qold,qnew=covers(old),covers(b);need(all(y>=x for x,y in zip(qold,qnew)),'hereditary child q monotonicity');child_profiles.append(qnew)
   rec['families'].append(dict(selection=selection,target_rows=U,target_child_profiles=child_profiles,value_ids=[value(U,a,c) for a,c in qs]))
  roots.append(rec)
 # All sixteen quotient/whole-child orientation states: complete exact finite battery vectors.
 Q=C3;B=[C3,TT2,TT3];T,off=blow(Q,B);masks=[((1<<len(b))-1)<<o for b,o in zip(B,off)];states=[];union=[]
 for bits in range(16):
  q0=op(Q) if bits&1 else Q;b0=[op(b) if bits>>(i+1)&1 else b for i,b in enumerate(B)];S,_=blow(q0,b0);u,v=(2,1) if bits&1 else(1,2);_,_,qs=battery(S,masks,0,u,v);union.extend(qs);states.append(dict(bits=bits,T=S))
 union=list(dict.fromkeys(union))
 # Add physical isolated-child probes for both quotient orientations and every core attachment.
 for S in states:
  for M in masks:
   names=[i for i in range(len(T)) if M>>i&1];p=names[0];extA=sum(1<<j for j in range(len(T)) if not M>>j&1 and S['T'][j]>>p&1)
   for a in range(1<<len(names)):
    for c in range(1<<len(names)):union.append((extA|sum(1<<v for j,v in enumerate(names) if a>>j&1),extA|sum(1<<v for j,v in enumerate(names) if c>>j&1)))
 union=list(dict.fromkeys(union))
 for S in states:S['value_ids']=[value(S['T'],a,c) for a,c in union]
 cube=dict(blocks=B,queries=union,states=states)
 # Full named bank with an unbalanced TT3 witness and a scalar-strict crossing target.
 Q=C3;B=[C3,TT3,one];T,off=blow(Q,B);U,_=blow(op(Q),B);V,_=blow(op(Q),[C3,C3,one]);n=len(T);masks=[M for M in range(1<<n) if M.bit_count()<=2];bank=[]
 for a in masks:
  for c in masks:bank.append(dict(a=a,c=c,value_ids=[value(S,a,c) for S in(T,U,V)]))
 extra=[]
 for k in range(40):a=(19*k+7)&127;c=(43*k+31)&127;extra.append(dict(a=a,c=c,value_ids=[value(S,a,c) for S in(T,U,V)]))
 with tempfile.TemporaryDirectory(prefix='aggregate-cyclic-') as td:
  exe=Path(td)/'count';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(ROOT/'count.cpp'),'-o',str(exe)],check=True);wire=str(len(graphs))+'\n'+'\n'.join(str(len(G))+' '+' '.join(map(str,G)) for G in graphs)+'\n';values=list(map(int,subprocess.check_output([str(exe)],input=wire.encode()).split()));need(len(values)==len(graphs),'all exact physical counts')
 for r in physical:
  for key in('local_queries','quotient_queries'):
   gap=0
   for q in r[key]:q['source']=values[q.pop('source_id')];q['target']=values[q.pop('target_id')];gap+=q['sign']*(q['source']-q['target'])
   need(gap==r['positive_formula']['coefficient']>0,'aggregate exact quartic and alias')
 strict=equal=0
 for r in roots:
  base=r['baseline_values']=[values[i] for i in r.pop('baseline_ids')];at={q:i for i,q in enumerate(r['battery'])};z=r['witness'][0];Z=r['block_masks'][z];zi=r['proper_original_modules'].index(Z)
  for f in r['families']:
   vals=f['values']=[values[i] for i in f.pop('value_ids')];need(vals!=base,'root-family full vector separation');hgap=vals[at[(0,0)]]-base[at[(0,0)]];need(hgap>=0,'scalar monotonicity');f['scalar_increase']=hgap
   if hgap:strict+=1;continue
   equal+=1;need(f['target_child_profiles']==list(map(covers,r['blocks'])),'scalar-equality complete profile restoration')
   if f['selection']>>zi&1:
    gap=sum(s*(base[at[(a,c)]]-vals[at[(a,c)]]) for a,c,s in r['quartic']);need(gap>0,'restored hub positive grouped quartic');f['quartic_gap']=gap
   else:
    residues=[(base[at[q]]-vals[at[q]])%4 for q in r['anchor']];need(2 in residues,'unrestored exposed triangle');f['anchor_residues']=residues
 for S in states:S['values']=[values[i] for i in S.pop('value_ids')]
 need(all(A['values']!=B0['values'] for A,B0 in combinations(states,2)),'all sixteen exact whole-child cube vectors distinct')
 table={}
 for rec in bank:
  rec['values']=[values[i] for i in rec.pop('value_ids')];table[rec['a'],rec['c']]=rec['values']
 def subsets(M):
  S=M
  while True:
   yield S
   if not S:break
   S=(S-1)&M
 for rec in bank:
  a,c=rec['a'],rec['c'];rec['coefficients']=[sum((-1)**((a.bit_count()+c.bit_count()-x.bit_count()-y.bit_count()))*table[x,y][j] for x in subsets(a) for y in subsets(c)) for j in range(3)]
 for rec in extra:
  rec['values']=[values[i] for i in rec.pop('value_ids')];calc=[sum(z['coefficients'][j] for z in bank if z['a']&~rec['a']==0 and z['c']&~rec['c']==0) for j in range(3)];need(calc==rec['values'],'complete named polynomial reconstruction')
 alias_terms=[[i,j,next(z['coefficients'][0]-z['coefficients'][1] for z in bank if z['a']==(1|(1<<i)) and z['c']==(2|(1<<j)))] for i in range(3,6) for j in [6]]
 alias_gap=sum(z[2] for z in alias_terms);matching=next(r for r in physical if r['T']==T and r['quotient_target']==U);need(alias_gap==matching['positive_formula']['coefficient']>0,'full named alias sum equals actual grouped quartic')
 selected=[127,7,24,48];target=list(T)
 for M in selected:need(module(T,M),'strict crossing bank uses original modules');target=flip(target,M)
 need(target==V,'literal root/triangle/crossing XOR bank target')
 counts={degree:sum(z['coefficients'][0]!=z['coefficients'][2] for z in bank if z['a'].bit_count()+z['c'].bit_count()==degree) for degree in range(5)}
 paths=[ROOT/n for n in('PROOF.md','verify.py','run_algebra.py','kernel.py','count.cpp')]+[BASE/'ordered_cyclic_module'/'PROOF.md',BASE/'module_decomposition'/'PROOF.md',BASE/'returns_next'/'UNIFORM.md',BASE/'returns_next'/'PROOF.md']
 return dict(status='PASS',identity=identity,physical_controls=physical,root_families=roots,root_family_count=strict+equal,scalar_strict_families=strict,scalar_equal_families=equal,whole_child_cube=cube,complete_named_bank=dict(source=T,baseline=U,strict_target=V,strict_target_original_selection=selected,aggregate_alias_terms=alias_terms,aggregate_alias_gap=alias_gap,coordinates=bank,full_mask_controls=extra,contrast_nonzeros_by_degree=counts),distinct_actual_augmented_graphs=len(graphs),source_sha256={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
