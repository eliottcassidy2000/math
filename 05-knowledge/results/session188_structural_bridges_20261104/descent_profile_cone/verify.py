"""Complete palindromic-descent-cone cubic certificate, orders one through twelve."""
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

def small_transform_controls():
 records=[]
 for n in range(1,6):
  masks=range(1<<n);full=(1<<n)-1;subwords={M:list(permutations(v for v in range(n) if M>>v&1)) for M in masks};seen=[]
  for code in range(1<<(n*(n-1)//2)):
   T=[0]*n
   for bit,(u,v) in enumerate(combinations(range(n),2)):
    a,b=(u,v) if code>>bit&1 else(v,u);T[a]|=1<<b
   c=[0]*n
   for p in subwords[full]:c[sum(not(T[u]>>v&1) for u,v in zip(p,p[1:]))]+=1
   need(c==c[::-1] and sum(c)==math.factorial(n),'literal palindrome and mass')
   H=[0]*(1<<n)
   for M in range(1,1<<n):H[M]=sum(all(T[u]>>v&1 for u,v in zip(p,p[1:])) for p in subwords[M])
   cover=[[0]*(1<<n) for _ in range(n+1)];cover[0][0]=1
   for k in range(1,n+1):
    for M in range(1,1<<n):
     S=M
     while S:cover[k][M]+=H[S]*cover[k-1][M^S];S=(S-1)&M
   q=[cover[k][full] for k in range(n+1)];transform=[sum(c[d]*choose(n-1-d,k-1-d) for d in range(n)) for k in range(n+1)]
   need(q==transform,'independent ordered-cover convolution versus word/cut transform')
   gs=[generator(n,d) for d in range((n+1)//2)];need(q==[sum(c[d]*gs[d][k] for d in range(len(gs))) for k in range(n+1)],'paired generators with central singleton')
   seen.append([code,c,q])
  records.append({'n':n,'labelled_tournaments':len(seen),'complete_descent_cover_records_sha256':hashlib.sha256(json.dumps(seen,separators=(',',':')).encode()).hexdigest()})
 return records

def main():
 records=[];total=0;coordinate_total=0
 for n in range(1,13):
  profiles=[generator(n,d) for d in range((n+1)//2)];triples=list(combinations_with_replacement(range(len(profiles)),3));power={(0,0,0):1};coordinates=[]
  for h in range(1,3*n+3):
   num=mul(K,power,n)
   if h>=2:
    coeffs=[[list(t),coefficient(num,h,profiles,n,t)] for t in triples];need(all(v<=0 for _,v in coeffs),'all complete cone coefficients nonpositive');coordinates.append({'h':h,'coefficients':coeffs});total+=len(coeffs);coordinate_total+=1
   power=mul(power,E,n)
  need(not mul(K,power,n),'all higher support vanishes');records.append({'n':n,'generators':profiles,'coordinates':coordinates})
 n=13;gs=[generator(n,d) for d in range(7)];positive=coefficient(mul(K,E,n),2,gs,n,(0,0,6));need(positive==1930071780,'coefficient-sign method boundary')
 pins=[ROOT/'PROOF.md',ROOT/'verify.py']+[ROOT.parent/d/'PROOF.md' for d in('covered_ambient','covered_positivity','equal_small_ambient')]
 print(json.dumps({'status':'PASS','scope':'Entire nonnegative palindromic descent cone through order12; actual equal-profile cyclic children are separated in every fixed ambient. No inference above12.','coefficient_count':total,'higher_cover_coordinates':coordinate_total,'records':records,'literal_transform_controls':small_transform_controls(),'order13_positive_mixed_coefficient':{'h':2,'generator_triple':[0,0,6],'coefficient':positive,'meaning':'Sufficient coefficientwise sign test fails; this is not a positive full cubic or an actual-profile counterexample.'},'source_sha256':{str(p.relative_to(ROOT.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in pins}},sort_keys=True,indent=2))
if __name__=='__main__':main()
