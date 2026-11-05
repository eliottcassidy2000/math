from pathlib import Path
import sys,csv,json,hashlib
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P/'source'));from spectrum import inverse
expect={}
for h,f,u,v in csv.reader((P/'inverse-normal.csv').open()):expect.setdefault(int(h),set()).add((f,int(u),int(v)))
solutions=0
for h in range(1,1862):
 got=set()
 for x in inverse(h):
  got.add(('R',0,0)if x['family']=='regular5'else(x['family'],x['u'],x['v']))
 if got!=expect.get(h,set()):raise RuntimeError(('complete inverse mismatch',h,got,expect.get(h)))
 solutions+=len(got)
print(json.dumps({'status':'PASS','all_integer_targets':1861,'exact_parameter_solutions':solutions,'independent_rectangle_sha256':hashlib.sha256((P/'inverse-normal.csv').read_bytes()).hexdigest()},indent=2,sort_keys=True))
