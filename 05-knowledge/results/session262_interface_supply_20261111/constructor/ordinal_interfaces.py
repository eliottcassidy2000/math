#!/usr/bin/env python3
"""Literal arbitrary-depth ordinal source compiler; no universal spectrum claim."""
import json
from pathlib import Path

def ordinal(A,OA,B,OB):
 n=len(A);m=len(B);return tuple([x+(((1<<m)-1)<<n) for x in A]+[x<<n for x in B]),OA|(OB<<n)
def profile(T,O):
 n=len(T);N=1<<n;end=[[0]*n for _ in range(N)]
 for v in range(n):end[1<<v][v]=1
 for mask in range(1,N):
  for v in range(n):
   z=end[mask][v]
   if z:
    avail=T[v]&~mask
    while avail:
     b=avail&-avail;avail-=b;end[mask|b][b.bit_length()-1]+=z
 h=[sum(a) for a in end];h[0]=1
 E=[sum(z for i,z in enumerate(a) if O>>i&1) for a in end]
 rev=tuple(sum(1<<j for j in range(n) if T[j]>>i&1) for i in range(n))
 start=[[0]*n for _ in range(N)]
 for v in range(n):start[1<<v][v]=1
 for mask in range(1,N):
  for v in range(n):
   z=start[mask][v]
   if z:
    avail=rev[v]&~mask
    while avail:
     b=avail&-avail;avail-=b;start[mask|b][b.bit_length()-1]+=z
 S=[sum(z for i,z in enumerate(a) if O>>i&1) for a in start];full=N-1
 return dict(H=h[full],Q=sum(h[a]*h[full^a] for a in range(N)),S=S[-1],E=E[-1],U=sum(h[a]*S[full^a] for a in range(N)),V=sum(E[a]*h[full^a] for a in range(N)),W=sum(E[a]*S[full^a] for a in range(N)))
def combine(a,b):
 return dict(H=a['H']*b['H'],Q=a['Q']*b['Q'],S=a['S']*b['H'],E=a['H']*b['E'],U=a['U']*b['Q']+a['H']*b['U'],V=a['Q']*b['V']+a['V']*b['H'],W=a['W']*b['H']+a['H']*b['W']+a['U']*b['V']+a['E']*b['S'])
def compile_word(atoms):
 if not atoms:raise ValueError('nonempty word required')
 T,O=atoms[0];p=profile(T,O)
 for B,OB in atoms[1:]:
  p=combine(p,profile(B,OB));T,O=ordinal(T,O,B,OB)
 return T,O,p

def balanced_trade_source(k,choices=None,padding=0):
 if k<0 or padding<0:raise ValueError('nonnegative k and padding required')
 choices=[0]*(2*k) if choices is None else list(choices)
 if len(choices)!=2*k or any(x not in (0,1) for x in choices):raise ValueError('need 2k binary choices')
 bank=((2,4,9,19,39,15),(2,4,57,3,11,27))
 atoms=[((0,),0)]*padding+[(bank[x],0) for x in choices[:k]]+[((0,),1)]+[(bank[x],0) for x in choices[k:]]+[((0,),0)]*padding
 return compile_word(atoms)
def selected_tail(B,O,t):
 if t<1:raise ValueError('tail length must be positive')
 n=len(B);out=list(B)+[O]+[0]*t
 for i in range(n):
  out[i]|=((1<<t)-1)<<(n+1)
  if not O>>i&1:out[i]|=1<<n
 for j in range(t):out[n+1+j]=(1<<n)|sum(1<<(n+1+a) for a in range(j+1,t))
 return tuple(out)

def main():
 atoms={'unmarked':((0,),0),'marked':((0,),1),'triangle':((2,4,1),0),'unit32':((2,4,9,19,7),9),'offset9':((2,28,25,17,1),1)}
 words=[['marked'],['unit32'],['offset9'],['unit32','unmarked','unmarked'],['triangle','unit32'],['unit32','triangle'],['marked','unmarked','marked'],['triangle','marked','triangle'],['unit32','marked'],['marked','unit32'],['offset9','unmarked'],['unit32','unit32']]
 receipt=[]
 for word in words:
  T,O,p=compile_word([atoms[k] for k in word]);actual=profile(T,O)
  if p!=actual:raise RuntimeError((word,p,actual))
  receipt.append(dict(word=word,rows=T,marks=O,profile=p,D=p['H']-p['W']))
 families=[]
 for k in (0,1):
  T,O,p=balanced_trade_source(k,choices=[0,1]*k)
  if p!=profile(T,O) or (p['H'],p['W'],p['U'],p['V'])!=(289**k,0,3808**k,3808**k):raise RuntimeError('balanced family')
  G=selected_tail(T,O,1);target=289**k+2*3808**k
  if profile(G,0)['H']!=target:raise RuntimeError('tail family')
  families.append(dict(k=k,source_order=len(T),source_profile=p,tail_order=len(G),tail_H=target))
 out={'balanced_initial_cases':families,'status':'PASS','scope':'Twelve literal controls of proved arbitrary-depth ordinal interface law; not a complete source census','controls':receipt}
 Path(__file__).with_name('ordinal_receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
if __name__=='__main__':
 import argparse
 parser=argparse.ArgumentParser();parser.add_argument('--k',type=int);parser.add_argument('--tail',type=int,default=1);parser.add_argument('--padding',type=int,default=0);args=parser.parse_args()
 if args.k is None:main()
 else:
  T,O,p=balanced_trade_source(args.k,padding=args.padding)
  print(json.dumps(dict(k=args.k,padding=args.padding,source_rows=T,marks=O,profile=p,tail_rows=selected_tail(T,O,args.tail),tail_H=p['H']-p['W']+2**args.tail*p['U'])))
