#!/usr/bin/env python3
"""Literal two-port constructions for every odd 23-smooth N except 7 and21."""
from pathlib import Path
import argparse,json
from interface import physical
ROOT=Path(__file__).resolve().parent
PRIMES=(3,5,7,11,13,17,19,23)

def factor_parts(N):
 if N<1:raise ValueError('N must be positive')
 rem=N;power={}
 for p in PRIMES:
  power[p]=0
  while rem%p==0:rem//=p;power[p]+=1
 if rem!=1:raise ValueError('N must be odd and have no prime factor above 23')
 if N in (7,21):raise ValueError('7 and 21 are excluded')
 parts=[];a=power[7];power[7]=0
 if a%2:
  if a>=3:parts.append(343);a-=3
  else:
   q=next((p for p in PRIMES if p not in (3,7) and power[p]),None)
   if q is not None:parts.append(7*q);power[q]-=1;a=0
   elif power[3]>=2:parts.append(63);power[3]-=2;a=0
   else:raise RuntimeError('unhandled factorization')
 parts.extend([49]*(a//2))
 for p in PRIMES:parts.extend([p]*power[p])
 prod=1
 for p in parts:prod*=p
 if prod!=N:raise RuntimeError('factorization mismatch')
 return parts

def compile_word(N):
 bank=json.loads((ROOT/'smooth_seeds.json').read_text());parts=factor_parts(N);word=[]
 for i,p in enumerate(parts):
  if i:word.extend([[0,0],[0,0]])
  word.extend(bank[str(p)]['word'])
 return word,parts

def compile_graph(N):
 word,parts=compile_word(N);rows=[2,0];ports=[0,1]
 for S,d in word:rows,ports=physical(rows,ports,S,d)
 return {'target':N,'factors':parts,'word':word,'rows':rows,'ports':ports}

if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('N',type=int);args=parser.parse_args()
 print(json.dumps(compile_graph(args.N),indent=2))
