#!/usr/bin/env python3
"""Actual witnesses in the finite prefix-cap64 digital-form bank."""
import argparse,json
from pathlib import Path
from one_block_formula import defect
from audit_inverse_witnesses import adjacency
HERE=Path(__file__).resolve().parent
def f(n):return 2*defect(n)+1
def seed_word(x):return [('S',int(c)) for c in bin(x)[2:]] if x else []
def suffix_word(y):return [('S',1-(y>>j&1)) for j in range(y.bit_length())]
def block_word(name,mask):
 if name=='D':return [('D',mask)]
 k=0 if name=='DD' else int(name[2:-1])
 return [('D',mask&15)]+[('S',(mask>>(4+j))&1) for j in range(k)]+[('D',mask>>(4+k))]

def search(N,cap):
 profiles=[]
 for name in ['D','DD']+[f'DT{k}D' for k in range(1,6)]:
  for line in(HERE/f'{name}_profiles.tsv').read_text().splitlines():profiles.append((name,*map(int,line.split())))
 evaluations=0
 for x in range(cap+1):
  fx=f(x)
  for name,H,E,p,q,s,t,F,mask in profiles:
   F0=H*fx+p*x+(H-s)*(x%2)+F-H;Q=E*x+q;g=H-t
   if F0>N:continue
   if Q==0:
    if N%H:continue
    m=N//H
    candidates=[0] if m==1 else (2*a for a in range((m-1)//2,m) if f(a)+a+a%2==m)
   else:
    lo=max(0,(N-F0-g+Q+H-1)//(Q+H));hi=(N-F0)//Q
    candidates=range(lo,hi+1)
   for y in candidates:
    evaluations+=1
    if F0+Q*y+g*(y%2)+H*(f(y)-1)!=N:continue
    word=seed_word(x)+block_word(name,mask)+suffix_word(y);rows=adjacency(word)
    return {'target':N,'x':x,'y':y,'parent':name,'mask':mask,'word':word,'vertices':len(rows),'rows':rows,'evaluations':evaluations}
 # The declared bank also includes words with no non-singleton block.
 candidates=[0] if N==1 else (2*a for a in range((N-1)//2,N) if f(a)+a+a%2==N)
 for y in candidates:
  word=seed_word(y);rows=adjacency(word)
  return {'target':N,'x':0,'y':y,'parent':'singleton','mask':None,'word':word,'vertices':len(rows),'rows':rows,'evaluations':evaluations}
 return {'target':N,'x':None,'prefix_cap':cap,'evaluations':evaluations}

def main():
 p=argparse.ArgumentParser();p.add_argument('targets',nargs='*',type=int);p.add_argument('--cap',type=int,default=63);a=p.parse_args()
 if a.cap<0 or any(N<1 or N%2==0 for N in a.targets):raise ValueError('nonnegative cap and positive odd targets required')
 records=[search(N,a.cap) for N in a.targets or [154711,9999991,9999999]]
 folder=HERE/'formula_witnesses';folder.mkdir(exist_ok=True)
 for r in records:
  if r['x'] is not None:
   n=r['vertices'];(folder/f"H{r['target']}.txt").write_text('\n'.join([f"{n} {r['target']}"]+[''.join(str(row>>j&1) for j in range(n)) for row in r['rows']])+'\n')
 print(json.dumps({'status':'complete target search','prefix_cap':a.cap,'results':records},indent=2,sort_keys=True))

if __name__=='__main__':main()
