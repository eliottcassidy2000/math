from pathlib import Path
import json,random,csv,re,ast
P=Path(__file__).resolve().parent

def ck(b,s):
 if not b:raise RuntimeError(s)
def integrate(a):
 out=[];s=0
 for v in a:s+=v;out.append(s)
 return out
def twice(a):return integrate(integrate(a))
rng=random.Random(527)
N=100
# Generate actual truncated series by repeated multiplication by 1/(1-z).
pair=[0]*N;pair[12]=8;pair[13]=8
pair=integrate(twice(pair))
for n in range(N):
 ck(pair[n]==(8*(n-11)**2 if n>=11 else 0),'decomposable series')
checks=0
for _ in range(1000):
 cuts=sorted([0,162]+[rng.randrange(163) for _ in range(3)])
 counts=[cuts[i+1]-cuts[i] for i in range(4)]
 a=[0]*N;a[8]=22
 for m,c in zip(range(9,13),counts):a[m]=c
 W=sum(i*c for i,c in enumerate(a));one=twice(a)
 ck(1634<=W<=2120,'moment range')
 for n in range(11,N):
  ck(one[n]+pair[n]==8*n*n+8*n+1152-W,'independent convolution');checks+=1
 for n in range(11,N-1):ck(one[n+1]+pair[n+1]-one[n]-pair[n]==16*(n+1),'difference')
a=[0]*N;a[7]=14;a[8]=7;a[9]=2;seq=twice(a)
ck(seq[7]==14 and 23*7-149==12,'minimal89 onset')
for n in range(8,N):ck(seq[n]==23*n-149,'89 series')
# Exhaust all residue shifts, independently from explicit polynomial sequence.
s=[(8*n*n+8*n)%72 for n in range(144)]
periods=[d for d in range(1,73) if all(s[i]==s[i+d] for i in range(72))]
ck(periods==list(range(9,73,9)),'exact quadratic periods')
s=[(23*n-149)%72 for n in range(144)]
ck([d for d in range(1,73) if all(s[i]==s[i+d] for i in range(72))]==[72],'exact linear period')
print(json.dumps({'status':'PASS','random_admissible_distributions':1000,'coefficient_comparisons':checks,'quadratic_period':9,'linear_period':72,'exhaustive_tournament_census_replayed':False},indent=2))
