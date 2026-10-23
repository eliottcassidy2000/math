import runpy,itertools,math,json
from pathlib import Path
f=runpy.run_path(str(Path(__file__).with_name('augmentation.py')));pf=f['pf'];n=7
free=[(i,j) for i in range(n) for j in range(i+2,n)];profiles={};bydet={};wit=None
for mask in range(1<<len(free)):
 M=[[0]*n for _ in range(n)]
 for i in range(n-1):M[i][i+1]=1;M[i+1][i]=-1
 for b,(i,j) in enumerate(free):M[i][j]=1 if mask>>b&1 else -1;M[j][i]=-M[i][j]
 K=[[M[i][j]-M[i][-1]+M[j][-1] for j in range(n-1)] for i in range(n-1)]
 p=abs(pf(tuple(sum(K,[]))));g4=math.gcd(*(pf(tuple(K[i][j] for i in sub for j in sub)) for sub in itertools.combinations(range(6),4)))
 key=(1,g4,p//g4);profiles[key]=profiles.get(key,0)+1
 if p in bydet and bydet[p][0]!=key and wit is None:wit={'abs_pf':p,'first':bydet[p],'second':(key,mask)}
 bydet.setdefault(p,(key,mask))
print('n7 full fixedHP profile census',profiles,'samePf_differentSmith_witness',wit)
