import importlib.util,itertools,math,json,random
from pathlib import Path
ROOT=Path(__file__).resolve().parent.parent
p=ROOT/'multi-exit/verify.py'
if not p.exists(): p=ROOT/'session215-multi-exit/verify.py'
s=importlib.util.spec_from_file_location('m',p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
rng=random.Random(215);cases=0
for k in range(1,7):
 for trial in range(10):
  A=[0]*k
  for a,b in itertools.combinations(range(k),2):
   i,j=(a,b) if rng.randrange(2) else (b,a);A[i]|=1<<j
  # Finite differences with scaling isolate highest clock from the EGF engine.
  vals=[m.egf(A,k,t) for t in range(1,k+2)]
  for root in range(1,k+1):vals=[vals[i+1]-root*vals[i] for i in range(len(vals)-1)]
  expected=math.factorial(k)**3*(k+1)
  if vals!=[expected]:raise RuntimeError((k,A,vals,expected))
  cases+=1
print(json.dumps({'status':'PASS','independent_random_leading_clock_cases':cases,'k_range':[1,6]},indent=2))
