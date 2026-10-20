#!/usr/bin/env python3
import itertools,json
from pathlib import Path

def T(x,m): return ((x<<1)&((1<<m)-1)) | (1^(x>>(m-1)))
def R(x,m): return ((x<<1)&((1<<m)-1)) | (x>>(m-1))
def orbit(x,m,tw=False):
 y=x;out=[]
 while y not in out:out.append(y);y=(T if tw else R)(y,m)
 return out

def valid(x,m):
 d=len(orbit(x,m));return sum((x>>i)&1 for i in range(d))%2==0

def generate(m,d):
 reps=sorted({min(orbit(x,m,True)) for x in range(1<<m)})
 targets=sorted({min(orbit(x,m)) for x in range(1<<m) if valid(x,m)})
 mons=[sum(1<<i for i in c) for k in range(d+1) for c in itertools.combinations(range(m-1),k)]
 points=[]
 for x in reps:
  for i in range(m):
   z=x ^ (((x>>(m-1))&1)*((1<<m)-1));points.append(z)
   x=(x>>1) | ((1^(x&1))<<(m-1))
 basis=[sum(int((z&s)==s)<<i for i,z in enumerate(points)) for s in mons]
 can=[targets.index(min(orbit(x,m))) if valid(x,m) else -1 for x in range(1<<m)]
 return reps,targets,mons,points,basis,can

def rref(rows):
 rows=[x for x in rows if x];piv=[];k=0
 for c in range(max(rows,default=0).bit_length()):
  hit=next((j for j in range(k,len(rows)) if rows[j]>>c&1),None)
  if hit is None:continue
  rows[k],rows[hit]=rows[hit],rows[k]
  for j in range(len(rows)):
   if j!=k and rows[j]>>c&1:rows[j]^=rows[k]
  piv.append(c);k+=1
 return rows[:k],piv

def nullspace(rows,n):
 reduced,piv=rref(rows)
 return [sum((r>>f&1)<<p for r,p in zip(reduced,piv))|(1<<f) for f in range(n) if f not in piv]

if __name__=='__main__':
 m=8;reps,tg,mons,pts,basis,can=generate(m,2)
 text=f'{len(basis)} {len(reps)}\n'+''.join(f'{b&((1<<64)-1)} {b>>64}\n' for b in basis)+' '.join(map(str,can))+'\n'
 old=Path(__file__).with_name('input8.txt').read_text()
 if text.split()!=old.split():raise RuntimeError('input8 basis mismatch')
 parity=[sum(((b>>(m*j))&((1<<m)-1)).bit_count()%2<<k for k,b in enumerate(basis)) for j in range(len(reps))]
 print('basis input independently regenerated, parity rank',len(rref(parity)[0]))
 for j,rep in enumerate(reps):
  rows=parity+[sum((b>>(m*j+i)&1)<<k for k,b in enumerate(basis)) for i in range(m)]
  print(rep,'null dimension',len(nullspace(rows,len(mons))))
