import importlib.util,json
from pathlib import Path
p=Path(__file__).resolve().parents[1]/'verify_ternary_natural_boundary.py'
spec=importlib.util.spec_from_file_location('author',p);a=importlib.util.module_from_spec(spec);spec.loader.exec_module(a)
counts,digits,*_=a.compile_counts(130)
matrices=[((1,-2,2),(2,-1,2),(2,-2,3)),((1,2,2),(2,1,2),(2,2,3)),((-1,2,2),(-2,1,2),(-2,2,3))]
level=[(3,4,5)];got=[];checks=0
for n in range(9):
 c=sum(z*z>4*x*y for x,y,z in level);got.append(c)
 if c!=counts[n]:raise RuntimeError(('Berggren mismatch',n,c,counts[n]))
 checks+=1
 if n<8: level=[tuple(sum(row[i]*t[i] for i in range(3)) for row in m) for t in level for m in matrices]
b=0
for n in range(129):
 if counts[n]!=2*b+digits[n]:raise RuntimeError(('digit',n))
 b=3*b+digits[n];checks+=1
for n in range(1,129):
 if counts[n]-3*counts[n-1]!=digits[n]-digits[n-1]:raise RuntimeError(('defect',n))
 checks+=1
print(json.dumps({'checks':checks,'direct_berggren_depth':8,'direct_counts':got,'status':'PASS'},indent=2))
