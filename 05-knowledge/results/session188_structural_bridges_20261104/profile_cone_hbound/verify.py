"""Exact H-bounded cone obstruction; no claim of tournament realizability."""
from pathlib import Path
import importlib.util,itertools,json,math,hashlib
ROOT=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('v',ROOT.parent/'covered_response'/'verify.py');v=importlib.util.module_from_spec(s);s.loader.exec_module(v)
E=v.add(v.X,v.Y,v.Z,v.prod(v.X,v.Y),v.prod(v.X,v.Z),v.prod(v.Y,v.Z),v.scale(v.P,3))
def conv(a,b):
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return c
def polynomial(small,n,h):
 a=v.profile(v.tt(small))+[0]*(n-small);b=v.profile(v.tt(n));p=v.ONE
 for _ in range(h-1):p={e:c for e,c in v.mul(p,E).items() if max(e)<=n}
 num=v.mul(v.N,p);out=[0]*4
 for e,c in num.items():
  for k in range(n-max(e)+1):
   term=[h*c*math.comb(k+h,h)]
   for r in e:term=conv(term,[b[k+r],a[k+r]])
   for j,x in enumerate(term):out[j]+=x
 v.need(all(x%3==0 for x in out),'integer coefficient');return [x//3 for x in out],a,b
co,a,b=polynomial(5,11,9)
v.need(co==[-2756654916959793711894994556905685760,72943657112400893970868442880,533607768637358688000,-10019483308800],'formal polynomial')
def f(w):return sum(c*w**j for j,c in enumerate(co))
w=35000000;q=[y+w*x for x,y in zip(a,b)];gap=f(w);bound=math.factorial(11)
v.need(gap==20457251690201969885400943894314240 and gap>0,'positive gap')
v.need(q[1]<=bound and w%2==0,'Hamilton bound and even cone coefficient')
v.need(q[0]==0 and all(x>0 for x in q[1:]),'positive full support')
v.need(all(q[k]%math.factorial(k)==0 for k in range(1,12)),'factorial divisibility')
v.need(q[-1]==bound and q[-2]==bound*10//2,'top cover coordinates')
# The derivative is strictly positive throughout [0,11!-1].
v.need(co[1]-3*abs(co[3])*(bound-1)**2>0 and co[2]>0,'monotonicity certificate')
lo=0;hi=(bound-1)//2
while lo<hi:
 mid=(lo+hi)//2
 if f(2*mid)>0:hi=mid
 else:lo=mid+1
threshold=2*lo;v.need(f(threshold)>0 and f(threshold-2)<=0,'first even crossing')
cover_bounds=[bound*math.comb(10,k-1) for k in range(1,12)]
violations=[{'k':k,'q_k':q[k],'universal_word_cut_bound':cover_bounds[k-1]} for k in range(1,12) if q[k]>cover_bounds[k-1]]
v.need([x['k'] for x in violations]==[2,3,4],'explicit nonrealizability obstruction')
# Every word contributes choose(10-d,k-1-d) cuts if it has d backward adjacencies.
descent=[]
for k in range(1,12):descent.append(q[k]-sum(descent[d]*math.comb(10-d,k-1-d) for d in range(k-1)))
v.need(descent==descent[::-1] and min(descent)<0,'palindromic but impossible descent distribution')
print(json.dumps({'status':'PASS','scope':'H-bounded artificial same-parity cone counterexample, not a tournament profile or full-response collision.','order':11,'small_generator':5,'h':9,'weight':w,'coefficient_vector':co,'q':q,'H_bound':bound,'positive_contrast':gap,'first_even_positive_weight_in_this_ray':threshold,'contrast_before_threshold':f(threshold-2),'contrast_at_threshold':f(threshold),'factorial_divisibility':True,'cover_bound_violations':violations,'descent_coefficients':descent,'source_sha256':{str(p.relative_to(ROOT.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'PROOF.md',ROOT/'verify.py',ROOT.parent/'covered_ambient'/'PROOF.md',ROOT.parent/'covered_response'/'verify.py']}},sort_keys=True,indent=2))
