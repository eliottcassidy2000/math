"""Exact weighted AM-GM certificates on the palindromic descent cone, orders13..16."""
from pathlib import Path
from itertools import combinations,combinations_with_replacement,permutations
import hashlib,json,math
ROOT=Path(__file__).resolve().parent

def need(ok,label):
 if not ok:raise RuntimeError(label)
def choose(n,k):return math.comb(n,k) if 0<=k<=n else 0
def generator(n,d):return [sum(choose(n-1-e,k-1-e) for e in {d,n-1-d}) for k in range(n+1)]
def mul(P,Q,n):
 out={}
 for a,x in P.items():
  for b,y in Q.items():
   e=tuple(u+v for u,v in zip(a,b))
   if max(e)<=n:out[e]=out.get(e,0)+x*y
 return {e:c for e,c in out.items() if c}
E={(1,0,0):1,(0,1,0):1,(0,0,1):1,(1,1,0):1,(1,0,1):1,(0,1,1):1,(1,1,1):3}
K={(2,2,0):1,(2,0,2):1,(0,2,2):1,(2,1,0):1,(2,0,1):1,(1,2,0):1,(1,0,2):1,(0,2,1):1,(0,1,2):1,(2,1,1):-1,(1,2,1):-1,(1,1,2):-1,(1,1,1):-6}
def mixed(num,h,profiles,n):
 return sum(c*math.comb(j+h,h)*profiles[0][j+e[0]]*profiles[1][j+e[1]]*profiles[2][j+e[2]] for e,c in num.items() for j in range(n-max(e)+1))
def coefficient(num,h,profiles,n,t):
 mult=6//math.prod(math.factorial(t.count(i)) for i in set(t));v=h*mult*mixed(num,h,[profiles[i] for i in t],n);need(v%3==0,'integer polarized coefficient');return v//3

SCALES={
 (13,2):[1,2,4,8,8,16,32],
 (14,2):[1,2,4,8,8,16,32],(14,3):[1,2,4,8,8,16,16],
 (15,2):[1,2,4,8,16,16,32,64],(15,3):[1,2,4,8,8,16,32,64],(15,4):[1,2,4,8,8,16,32,64],
 (16,2):[1,2,4,8,16,16,32,64],(16,3):[1,2,4,8,8,16,32,32],(16,4):[1,2,4,8,8,16,32,32]}

def main():
 records=[];total=coordinate_total=0;used=[]
 for n in range(13,17):
  profiles=[generator(n,d) for d in range((n+1)//2)];triples=list(combinations_with_replacement(range(len(profiles)),3));power={(0,0,0):1};coordinates=[]
  for h in range(1,3*n+3):
   num=mul(K,power,n)
   if h>=2:
    coeffs={t:coefficient(num,h,profiles,n,t) for t in triples};positive={t:c for t,c in coeffs.items() if c>0};record={'h':h,'coefficients':[[list(t),c] for t,c in coeffs.items()]};total+=len(coeffs);coordinate_total+=1
    if positive:
     need((n,h) in SCALES,'every positive-coefficient coordinate has a certificate');d=SCALES[n,h];need(len(d)==len(profiles) and all(isinstance(v,int) and v>0 for v in d),'positive integer scales');A=[-coeffs[i,i,i] for i in range(len(d))];need(all(a>=0 for a in A),'nonpositive pure cubes');usage=[0]*len(d)
     for t,c in positive.items():
      z=c*math.prod(d[i] for i in t)
      for i in t:usage[i]+=z
     capacity=[3*A[i]*d[i]**3 for i in range(len(d))];slack=[a-b for a,b in zip(capacity,usage)];need(min(slack)>=0,'all exact weighted AM-GM budgets');record['amgm']={'scales':d,'negative_cube_magnitudes':A,'positive_terms':[[list(t),c] for t,c in positive.items()],'usage':usage,'capacity':capacity,'slack':slack};used.append([n,h])
    else:record['all_coefficients_nonpositive']=True
    coordinates.append(record)
   power=mul(power,E,n)
  need(not mul(K,power,n),'vanishing all later support');records.append({'n':n,'generators':profiles,'coordinates':coordinates})
 need({tuple(x) for x in used}==set(SCALES),'exact certificate coverage')
 pins=[ROOT/'PROOF.md',ROOT/'verify.py']+[ROOT.parent/d/'PROOF.md' for d in('descent_profile_cone','covered_ambient','covered_positivity','equal_small_ambient')]
 print(json.dumps({'status':'PASS','scope':'Entire nonnegative palindromic descent cone at orders13..16; prior packet covers1..12. No inference above16.','coefficient_count':total,'higher_cover_coordinates':coordinate_total,'positive_coefficient_coordinates':used,'records':records,'source_sha256':{str(p.relative_to(ROOT.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in pins}},sort_keys=True,indent=2))
if __name__=='__main__':main()
