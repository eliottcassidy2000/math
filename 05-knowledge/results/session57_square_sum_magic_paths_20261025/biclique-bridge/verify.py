#!/usr/bin/env python3
"""Exact decorated-biclique certificate; standard library; safe under -O."""
import json,sys,hashlib,itertools
from collections import Counter
from pathlib import Path
from math import isqrt
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'sources'))
from field import Field,magic,valid

def require(ok,msg='certificate mismatch'):
 if not ok: raise RuntimeError(msg)
def sq(n):return n>=0 and isqrt(n)**2==n

def field_case(F,m):
 require(valid(F,m),'invalid input magic square')
 c=m[4];x=F.sub(m[0],c);y=F.sub(m[2],c)
 S=[0,x,y,F.add(x,y)]
 D=[F.sub(t,s) for s in S for t in S]; ss={F.add(s,t) for s in S for t in S}
 require(len(S)==len(set(S))==4 and len(set(D))==len(ss)==9)
 require(ss=={F.add(F.add(x,y),d) for d in D})
 gauges=[];forbidden=[]
 for g in range(F.q):
  L=[F.sub(F.sub(c,s),g) for s in S];R=[F.add(s,g) for s in S]
  E=[F.add(l,r) for l in L for r in R]
  require(set(E)==set(m) and set(E)<=F.squares)
  require(sorted(Counter(E).values())==[1,1,1,1,2,2,2,2,4])
  bad=F.add(g,g) in {F.sub(c,z) for z in ss}
  require(bad==(len(set(L+R))!=8))
  (forbidden if bad else gauges).append(g)
 require(len(forbidden)==9 and len(gauges)==F.q-9)
 g=gauges[0];L=[F.sub(F.sub(c,s),g) for s in S];R=[F.add(s,g) for s in S]
 return dict(q=F.q,p=F.p,n=F.n,modulus=F.modulus,matrix=m,c=c,x=x,y=y,S=S,good_gauge_count=len(gauges),forbidden_gauges=forbidden,example_gauge=g,L=L,R=R,edge_labels=[[F.add(l,r) for r in R] for l in L])

def main():
 data=json.loads((ROOT/'sources/results.json').read_text());cases=[]
 for q in (29,59,81):
  rec=next(r for r in data['fields'] if r['q']==q);F=Field(rec['p'],rec['n'],rec['modulus'])
  for group in rec['spectral_hp_groups']: cases.append(field_case(F,group['matrix']))
 r=json.loads((ROOT/'sources/witness-2187.json').read_text());cases.append(field_case(Field(r['p'],r['n'],r['modulus']),r['matrix']))
 # Integer near-miss: exact magic and reflected decoration, but nonsquare labels.
 c,x,y=16,12,1;S=[0,x,y,x+y];L=[4*c-4*s-1 for s in S];R=[4*s+1 for s in S]
 E=[l+r for l in L for r in R]
 require(min(L+R)>0 and len(set(L+R))==8 and {z%4 for z in L}=={3} and {z%4 for z in R}=={1})
 require(len(set(E))==9 and not all(map(sq,E)))
 near=dict(c=c,x=x,y=y,L=L,R=R,labels=sorted(set(E)),nonsquares=sorted({z for z in E if not sq(z)}))
 # Bare square-sum K4,4, but thirteen labels: decoration cannot be recovered.
 L0=[0,1,2,6];R0=[3,15,19,20];Q={i*i%59 for i in range(59)};E0=[(l+r)%59 for l in L0 for r in R0]
 require(len(set(L0+R0))==8 and set(E0)<=Q and len(set(E0))==13)
 # Symbolic linear algebra represented by integer triples (C,x,y).
 add=lambda a,b:tuple(u+v for u,v in zip(a,b));neg=lambda a:tuple(-u for u in a)
 C=(1,0,0);X=(0,1,0);Y=(0,0,1);O=(0,0,0);s=[O,X,Y,add(X,Y)]
 edge=Counter(add(C,add(t,neg(v))) for v in s for t in s)
 M=[(1,1,0),(1,-1,-1),(1,0,1),(1,-1,1),C,(1,1,-1),(1,0,-1),(1,1,1),(1,-1,0)]
 require(set(edge)==set(M) and sorted(edge.values())==[1,1,1,1,2,2,2,2,4])
 lines=[(0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6)]
 require(all(add(add(M[a],M[b]),M[d])==(3,0,0) for a,b,d in lines))
 # Three Hamiltonian paths cover the entire K4,4, preserving every edge label.
 paths=[];uncovered={(i,j) for i in range(4) for j in range(4)}
 for a in itertools.permutations(range(4)):
  for b in itertools.permutations(range(4)):
   path=[v for ij in zip(a,tuple(4+j for j in b)) for v in ij]
   edges={tuple((u,v-4) if u<4 else (v,u-4)) for u,v in zip(path,path[1:])}
   paths.append((path,edges))
 chosen=[]
 while uncovered:
  path,edges=max(paths,key=lambda t:len(t[1]&uncovered));chosen.append(path);uncovered-=edges
 require(len(chosen)==3)
 return dict(field_witnesses=cases,integer_near_miss=near,bare_K44_F59=dict(L=L0,R=R0,labels=sorted(set(E0)),distinct_label_count=13),symbolic_grid=M,symbolic_edge_multiplicities=[dict(label=list(k),multiplicity=v) for k,v in sorted(edge.items())],covering_Hamiltonian_paths=chosen,hamiltonian_path_count_up_to_reversal=len(paths))

if __name__=='__main__':
 out=json.dumps(main(),indent=2,sort_keys=True)+'\n';target=ROOT/'certificate.json'
 if '--check' in sys.argv:require(target.read_text()==out,'stored certificate differs')
 else:target.write_text(out)
 print('PASS: symbolic identities, finite-field gauges, square labels, near-miss, bare K4,4, and path cover')
