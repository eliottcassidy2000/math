#!/usr/bin/env python3
"""Literal block profiles and the exact two-integer singleton-context formula."""
import json
from itertools import permutations
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'square_sum_spectrum_overnight_20261031b'))
from audit_inverse_witnesses import adjacency
from audit_feedback import subset_dp,require

def defect(n):
 d=0
 while n:
  a=n//2
  if n%2==0:d+=(a+1)//2
  n=a
 return d

def block_profiles(rows):
 n=len(rows);starts=[0]*n;ends=[0]*n;W=[[0]*n for _ in rows];H=0
 def path(p):return all(rows[u]>>v&1 for u,v in zip(p,p[1:]))
 for word in permutations(range(n)):
  if path(word):H+=1;starts[word[0]]+=1;ends[word[-1]]+=1
  for k in range(1,n):
   if path(word[:k]) and path(word[k:]):W[word[k-1]][word[k]]+=1
 E=2*H+sum(map(sum,W));weights=[[(W[u][v]+W[v][u])//2 for v in range(n)] for u in range(n)]
 require(all((W[u][v]+W[v][u])%2==0 for u in range(n) for v in range(n)),'integral cut weights')
 result=[]
 for mask in range(1<<n):
  s=sum(starts[v] for v in range(n) if mask>>v&1)
  t=sum(ends[v] for v in range(n) if not(mask>>v&1))
  p=sum(ends[u]+sum(W[u]) for u in range(n) if not(mask>>u&1))
  q=sum(starts[v]+sum(W[u][v] for u in range(n)) for v in range(n) if mask>>v&1)
  F=H+2*sum(weights[u][v] for u in range(n) for v in range(n)
            if mask>>u&1 and not(mask>>v&1) and rows[u]>>v&1)
  result.append({'H':H,'E':E,'p':p,'q':q,'s':s,'t':t,'F':F,'mask':mask})
 return result

def formula(b,x,y):
 H,E,p,q,s,t,F=(b[k] for k in ('H','E','p','q','s','t','F'))
 return F+(p+E*y)*x+q*y+(H-s)*(x%2)+(H-t)*(y%2)+2*H*(defect(x)+defect(y))

def main():
 D=[8,1,3,6];blocks=block_profiles(D)
 double=[r|240 for r in D]+[r<<4 for r in D];two=block_profiles(double)
 require(all(b['H']==5 and b['E']==32 for b in blocks),'single D profile')
 require(all(b['H']==25 and b['E']==1024 for b in two),'two adjacent D profile')
 checks=0
 for b in blocks:
  for x in range(4):
   for y in range(4):
    word=[('S',int(c)) for c in bin(x)[2:]] if x else []
    word+=[('D',b['mask'])]+[('S',1-(y>>j&1)) for j in range(y.bit_length())]
    require(subset_dp(adjacency(word))==formula(b,x,y),'independent full adjacency formula');checks+=1
 for name,records in [('D',blocks),('DD',two)]:
  (HERE/f'{name}_profiles.json').write_text(json.dumps(records,indent=2,sort_keys=True)+'\n')
  (HERE/f'{name}_profiles.tsv').write_text('\n'.join(' '.join(str(b[k]) for k in ('H','E','p','q','s','t','F','mask')) for b in records)+'\n')
 print(json.dumps({'status':'passed','literal_D_controls':checks,'D_profiles':16,'DD_profiles':256,
                   'D_recurrence':'D(0)=0; D(2a)=D(a)+ceil(a/2); D(2a+1)=D(a)',
                   'scope':'Every fixed block with arbitrary singleton contexts; a terminal trailing singleton1 does not change F.'},indent=2,sort_keys=True))

if __name__=='__main__':main()
