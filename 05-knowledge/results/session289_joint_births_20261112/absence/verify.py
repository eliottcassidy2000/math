#!/usr/bin/env python3
"""Fresh exact controls: complete labelled tournaments through6, SCC monoid,
triangle-module lower floor, and99-ray7-adic obstruction. Not a uniform proof
of Chan–Pak; mathematical all-order statements are in PROOF.md.
"""
import json

def hp(rows):
 n=len(rows);d=[[0]*n for _ in range(1<<n)]
 for i in range(n):d[1<<i][i]=1
 for mask in range(1,1<<n):
  for i,x in enumerate(d[mask]):
   if not x:continue
   bits=rows[i]&~mask
   while bits:
    bit=bits&-bits;bits-=bit;j=bit.bit_length()-1
    d[mask|bit][j]+=x
 return sum(d[-1])
def strong(rows):
 n=len(rows)
 for s in range(n):
  seen=1<<s;new=seen
  while new:
   nxt=0
   while new:
    b=new&-new;new-=b;nxt|=rows[b.bit_length()-1]
   new=nxt&~seen;seen|=new
  if seen!=(1<<n)-1:return False
 return True
def floor(n):
 a,r=divmod(n,3)
 return 3*5**(a-1) if r==0 else 5**a if r==1 else 9*5**(a-1)
def rows_code(n,code):
 r=[0]*n;k=0
 for i in range(n):
  for j in range(i+1,n):
   if code>>k&1:r[i]|=1<<j
   else:r[j]|=1<<i
   k+=1
 return r
def run():
 alph={1};census=[]
 for n in range(1,7):
  vals=set();number=0
  for code in range(1<<(n*(n-1)//2)):
   rows=rows_code(n,code)
   if strong(rows):number+=1;vals.add(hp(rows))
  alph|=vals;census.append({'n':n,'strong_labelled':number,'strong_values':sorted(vals)})
 cutoff=floor(7);mon={1}
 for x in range(1,cutoff):
  if x in mon:
   for a in alph:
    if x*a<cutoff:mon.add(x*a)
 holes=[x for x in range(1,cutoff,2) if x not in mon]
 if holes!=[7,21]:raise RuntimeError(holes)
 triangle=[]
 for m in range(1,9):
  n=m+2;rows=[0]*n
  for i in range(m):
   rows[i]=sum(1<<j for j in range(i+1,m))|(1<<m)
  rows[m]=1<<(m+1);rows[m+1]=(1<<m)-1
  h=hp(rows)
  if h!=2**m+1 or not strong(rows):raise RuntimeError(m)
  triangle.append({'m':m,'H':h})
 for n in range(3,101):
  if floor(n)**3<3**n:raise RuntimeError(('floor',n))
 residues=[(99*2**k+23)%7 for k in range(1,13)]
 if set(residues)!={3,4,6}:raise RuntimeError(residues)
 primes=[n for n in range(2,100) if all(n%d for d in range(2,int(n**.5)+1))]
 if len(primes)!=25:raise RuntimeError(primes)
 return {'status':'PASS','fresh_scope':'all labelled tournaments through6; fixed examples and algebraic controls',
 'census':census,'strong_alphabet':sorted(alph),'frozen_cutoff':cutoff,'frozen_holes':holes,
 'triangle_twins':triangle,'ray_mod7':residues,'primes_below100':primes}
if __name__=='__main__':print(json.dumps(run(),indent=2,sort_keys=True))
