from itertools import product,permutations
from pathlib import Path
import runpy,json
v=runpy.run_path(str(Path(__file__).with_name('verify.py')));a=[18,8,3,5,14];zs=[]
for c in (3,6):
 g=v['closure'](a,c)
 total=sum(all(g[i]>>j&1 for i,j in zip(p,p[1:])) for p in permutations(range(len(g))))
 v['need'](total==v['tables'](g)[0][-1],'permutation versus subset DP')
 zs.append(v['mul'](v['interface'](a,c),(1,2,0,0,0,0,0)))
rows=[]
for k in range(1,9):
 vals=[]
 for bits in product((0,1),repeat=k):
  z=zs[bits[0]]
  for bit in bits[1:]:z=v['mul'](z,zs[bit])
  vals.append(z[0]+2*z[2]-z[4])
 steps=[]
 for j in range(k):
  hp=13**j;Ap=0 if j==0 else 64*13**(j-1)*(16**j-1)//15
  steps.append(-13**(k-j-1)*(3*Ap+2*hp))
 for bits,total in zip(product((0,1),repeat=k),vals):v['need'](total==vals[0]+sum(b*s for b,s in zip(bits,steps)),'exact affine count')
 if k>=2:v['need'](all(s%(2*13**(k-2))==0 for s in steps),'scoped divisor')
 rows.append(dict(k=k,choices=len(vals),steps=steps))
from math import gcd
bs=[]
for cut in (17,33):
 z=v['interface']([32,17,3,7,13,30],cut)
 bs.append(v['mul'](z,(1,2,0,0,0,0,0)))
v['need'](tuple(y-x for x,y in zip(*bs))==(0,0,0,-6,20,-2,0),'h37 difference')
mixed=[]
for k in range(1,13):
 vals=[]
 for bits in product((0,1),repeat=k):
  z=zs[bits[0]]
  for bit in bits[1:]:z=v['mul'](z,bs[bit])
  vals.append(z[0]+2*z[2]-z[4])
 steps=[-2*37**(k-1)];hp,Ap=13,64
 for i in range(2,k+1):
  steps.append(37**(k-i)*(6*Ap-20*hp));Ap,hp=704*Ap+192*hp,37*hp
 for bits,total in zip(product((0,1),repeat=k),vals):v['need'](total==vals[0]+sum(b*s for b,s in zip(bits,steps)),'mixed affine count')
 v['need'](gcd(*steps)==2,'mixed all-length gcd')
 ordered=sorted(vals)
 v['need'](len(set(vals))==2**k,'mixed injectivity')
 spacing=min(b-a for a,b in zip(ordered,ordered[1:]))
 v['need'](spacing==(2 if k==1 else 50*37**(k-2)),'exact mixed minimum spacing')
 mixed.append(dict(k=k,steps=steps,choices=len(vals),gcd=2,minimum_spacing=spacing))
print(json.dumps(dict(status='PASS',mixed_chains=mixed,permutation_control='both seven-vertex closures',chains=rows),indent=2))
