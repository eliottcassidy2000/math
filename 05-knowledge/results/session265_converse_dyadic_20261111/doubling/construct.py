#!/usr/bin/env python3
"""Exact new Fano/converse parity switch controls; no external packages."""
from itertools import combinations, permutations
from collections import Counter
from math import gcd
from pathlib import Path
import json
LINES=[(a-1,b-1,(a^b)-1) for a in range(1,8) for b in range(a+1,8) if b<(a^b)]
def fano(mask):
 T=[0]*7
 for k,(a,b,c) in enumerate(LINES):
  if mask>>k&1:b,c=c,b
  for u,v in [(a,b),(b,c),(c,a)]:T[u]|=1<<v
 return tuple(T)
def rebase(T,r):
 p=[r]+[i for i in range(len(T)) if i!=r]
 return tuple(sum(((T[p[i]]>>p[j])&1)<<j for j in range(len(T))) for i in range(len(T)))
def converse(T):return tuple(sum(((T[j]>>i)&1)<<j for j in range(len(T))) for i in range(len(T)))
def hp(T,root=None):
 n=len(T);dp=[[0]*n for _ in range(1<<n)];H=[0]*(1<<n);H[0]=1 if root is None else 0
 for v in range(n):
  if root is None or root==v:dp[1<<v][v]=1
 for S in range(1,1<<n):
  H[S]=sum(dp[S])
  for v,c in enumerate(dp[S]):
   if c:
    rem=T[v]&~S
    while rem:
     b=rem&-rem;rem-=b;dp[S|b][b.bit_length()-1]+=c
 return H
def rooted(T):
 H=hp(T);R=hp(T,0);full=(1<<len(T))-1;P=[0]*(len(T)+1)
 for S,c in enumerate(R):P[S.bit_count()]+=c*H[full^S]
 return P
def profile(T):
 H=hp(T);full=len(H)-1;R=full^1;U=rooted(T);V=rooted(converse(T));Q=k=0;S=R
 while True:
  v=H[S]*H[R^S];Q+=v;k+=(-1)**S.bit_count()*v
  if S==0:break
  S=(S-1)&R
 return dict(H=H[full],h=H[R],Q=Q,k=k,U=sum(U),V=sum(V),d=sum((-1)**(i+1)*x for i,x in enumerate(U)),e=sum((-1)**(i+1)*x for i,x in enumerate(V)))
def matrix(p):
 h=p['h'];return [[h,p['U']-h,p['d']-h,p['H']-h],[0,p['Q'],0,p['V']-h],[0,0,p['k'],p['e']-h],[0,0,0,h]]
def mm(A,B):return [[sum(A[i][k]*B[k][j] for k in range(4)) for j in range(4)] for i in range(4)]
def fuse(A,B):
 n=len(A);m=len(B);T=[0]*(n+m-1)
 for X,p in [(A,list(range(n))),(B,[0]+list(range(n,n+m-1)))]:
  for i in range(len(X)):
   for j in range(len(X)):
    if X[i]>>j&1:T[p[i]]|=1<<p[j]
 for i in range(1,n):
  for j in range(n,n+m-1):T[i]|=1<<j
 return tuple(T)
def literal(T):
 n=len(T);W=[[0]*n for _ in T];E=[[0]*n for _ in T];H=0
 for p in permutations(range(n)):
  if all(T[u]>>v&1 for u,v in zip(p,p[1:])):
   H+=1;E[p[0]][p[-1]]+=1
   for u,v in zip(p,p[1:]):W[u][v]+=1;W[v][u]+=1
 return H,W,E
def walsh(vals):return [sum((-1)**((s&m).bit_count())*v for m,v in enumerate(vals)) for s in range(len(vals))]
def main():
 A=rebase(fano(0),2);B=converse(A);pa,pb=profile(A),profile(B)
 assert pa==dict(H=171,h=43,Q=424,k=0,U=260,V=257,d=24,e=25)
 assert pb==dict(H=171,h=43,Q=424,k=0,U=257,V=260,d=25,e=24)
 la,lb=literal(A),literal(B);assert la[:2]==lb[:2];assert la[0]==hp(A)[-1]
 M,N=matrix(pa),matrix(pb);C=[[a-b for a,b in zip(r,s)] for r,s in zip(mm(M,N),mm(N,M))]
 assert C==[[0,1143,43,1330],[0,0,0,1143],[0,0,0,43],[0,0,0,0]]
 AB,BA=fuse(A,B),fuse(B,A);tri=(2,4,1)
 counts={name:hp(T)[-1] for name,T in [('AB',AB),('BA',BA),('C3AB',fuse(tri,AB)),('C3BA',fuse(tri,BA))]}
 assert counts['AB']-counts['BA']==1330
 assert counts['C3AB']-counts['C3BA']==3616
 assert gcd(150,3616)==2
 common={name:hp(T)[-1] for name,T in [
  ('CABC',fuse(fuse(tri,AB),tri)),('CBAC',fuse(fuse(tri,BA),tri)),
  ('ABCC',fuse(fuse(AB,tri),tri)),('BACC',fuse(fuse(BA,tri),tri))]}
 assert common==dict(CABC=1191647,CBAC=1185745,ABCC=1087199,BACC=1074439)
 J=matrix(profile(tri))
 expected={'CABC':mm(mm(mm(J,M),N),J),'CBAC':mm(mm(mm(J,N),M),J),'ABCC':mm(mm(mm(M,N),J),J),'BACC':mm(mm(mm(N,M),J),J)}
 assert all(v[0][0]+v[0][3]==common[k] for k,v in expected.items())
 assert sorted([hp(fuse((2,0),(0,1)))[-1],hp(fuse((0,1),(2,0)))[-1]])==[1,3]
 for T in [A,B,tri]:
  for i in range(len(T)):
   seen=1<<i
   while True:
    nxt=seen
    for j in range(len(T)):
     if seen>>j&1:nxt|=T[j]
    if nxt==seen:break
    seen=nxt
   assert seen==(1<<len(T))-1
 assert gcd(common['CABC']-common['CBAC'],common['ABCC']-common['BACC'])==2
 # Finite full seven-line orientation cube. This checks the new Fourier bookkeeping,
 # not a new claim to the inherited octonion classification.
 vals=[];diff=[];w01=[]
 for m in range(128):
  T=fano(m);p=profile(T);h,w,e=literal(T);assert h==p['H'];vals.append(h);diff.append(p['U']-p['V']);w01.append(w[0][1])
 Hf,Wf,Df=map(walsh,[vals,w01,diff]);assert all(Hf[s]==Wf[s]==0 for s in range(128) if s.bit_count()%2)
 assert all(Df[s]==0 for s in range(128) if s.bit_count()%2==0)
 assert Counter(vals)=={171:112,189:16}
 # Context proof controls: all labelled tournaments of orders1..3 on both sides.
 bank=[]
 for n in range(1,4):
  ed=list(combinations(range(n),2))
  for m in range(1<<len(ed)):
   T=[0]*n
   for i,(a,b) in enumerate(ed):
    if m>>i&1:a,b=b,a
    T[a]|=1<<b
   bank.append(profile(tuple(T)))
 gaps=[]
 for L in bank:
  for R in bank:
   gap=mm(mm(matrix(L),C),matrix(R))[0][3]
   formula=1143*(L['h']*R['V']+L['U']*R['h'])+43*(L['h']*R['e']+L['d']*R['h'])-1042*L['h']*R['h']
   assert gap==formula and gap>=1158*L['h']*R['h'];gaps.append(gap)
 out=dict(status='FINITE-EXACT controls of newly constructed switch',A=list(A),B=list(B),profile_A=pa,profile_B=pb,profile_C=profile(tri),rooted_A=rooted(A),rooted_B=rooted(B),weighted_graph_A=la[1],endpoint_matrix_A=la[2],commutator=C,literal_counts=counts,common_inventory_counts=common,common_inventory_gap_gcd=2,independent_matrix_check_common_inventory=True,all_input_modules_strong=True,small_nonstrong_boundary_counts=[1,3],gcd_with_150=gcd(150,3616),context_pairs=len(gaps),fano_cube_counts=dict(Counter(vals)),walsh_H_numerators={s:v for s,v in enumerate(Hf) if v},walsh_W01_support=[s for s,v in enumerate(Wf) if v],walsh_UminusV_support=[s for s,v in enumerate(Df) if v],walsh_denominator=128)
 path=Path(__file__).with_name('receipt.json');path.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k not in ['weighted_graph_A','endpoint_matrix_A']},indent=2))
if __name__=='__main__':main()
