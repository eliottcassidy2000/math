"""Leading multiplicity classification on all 3,5,7-smooth H; finite formula controls."""
import json,math
from pathlib import Path
from functools import lru_cache
from fractions import Fraction as Q
ROOT=Path(__file__).parent
# (v7,v3,v5), total number of strong isomorphism classes, inherited.
ATOMS=[((0,1,0),1),((0,0,1),1),((1,0,1),2),((2,0,0),8),((1,2,0),2),((3,0,0),103)]
@lru_cache(None)
def dp(a,b,c):
 if (a,b,c)==(0,0,0):return 0,1
 best,count=-1,0
 for (u,v,w),mass in ATOMS:
  if a<u or b<v or c<w:continue
  r,k=dp(a-u,b-v,c-w)
  if r<0:continue
  r+=1
  if r>best:best,count=r,mass*k
  elif r==best:count+=mass*k
 return best,count

def formula(a,b,c):
 if a==1 and c==0:
  if b<2:return -1,Q(0)
  return b-1,Q(2,math.factorial(b-2))
 m=a//2;r=b+c+m
 if a%2==0:return r,Q(8**m,math.factorial(m)*math.factorial(b)*math.factorial(c))
 lead=Q(2*c*8**m,math.factorial(m))
 if m:lead+=Q(103*8**(m-1),math.factorial(m-1))
 return r,lead/(math.factorial(b)*math.factorial(c))

def main():
 checks=0
 for a in range(10):
  for b in range(10):
   for c in range(10):
    r,k=dp(a,b,c);s,L=formula(a,b,c)
    if r!=s or (Q(k,math.factorial(r)) if r>=0 else Q(0))!=L:raise RuntimeError((a,b,c,r,k,s,L))
    checks+=1
 examples=[]
 for a,b,c in [(0,2,0),(1,2,0),(1,3,0),(2,0,0),(3,0,0),(5,0,0),(2,0,2),(2,4,0),(3,0,1)]:
  r,L=formula(a,b,c);examples.append(dict(H=7**a*3**b*5**c,exponents_7_3_5=[a,b,c],degree=r,leading_coefficient=str(L),top_core_words=dp(a,b,c)[1]))
 out=dict(status='PASS',checks=checks,scope='Formula versus ordered-six-atom dynamic program on exponents 0..9 in each variable; no inherited inventory recount.',examples=examples)
 (ROOT/'receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
