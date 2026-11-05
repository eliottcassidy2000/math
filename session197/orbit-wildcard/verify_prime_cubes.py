from verify import *
from collections import Counter
hist={};total=0
for p in (5,7,11,13):
 c=Counter();F=set()
 for bits in range(1<<((p-1)//2)):
  ds=[i if bits>>(i-1)&1 else p-i for i in range(1,(p+1)//2)]
  R=cyclic(p,ds);h=h_all(R);other=h_backward(R)
  gate(h[-1]==other,'independent orbit H');c[other]+=1;total+=1
  if p<=7:F.add(sum(h[s]*h[((1<<p)-1)^s] for s in range(1,(1<<p)-1)))
 hist[str(p)]={str(k):v for k,v in sorted(c.items())}
 if p==5:gate(c=={15:4} and F=={90},'five orbit bank')
 if p==7:gate(c=={175:6,189:2} and F=={1638},'seven orbit bank')
 if p==11:gate(c=={92411:10,93027:10,93467:10,95095:2} and min(c)>p**3,'eleven orbit bank')
 if p==13:gate(c=={3669497:12,3683797:12,3703011:4,3704857:12,3707483:12,3711175:12} and min(c)>p**3,'thirteen orbit bank')
gate(5**16>17**9 and 5*17**9>18**9,'all-order integer bound')
gate(15**2>125 and 3*15+90>125 and 125%15!=0,'p5 module inequalities')
gate(175**2>343 and 3*175+1638>343 and all(343%h for h in (175,189)),'p7 module inequalities')
R=[0]*9
for a in range(9):
 for b in range(9):
  if a//3<b//3 or a//3==b//3 and (b-a)%3==1:R[a]|=1<<b
gate(h_all(R)[-1]==h_backward(R)==27,'p3 boundary')
p=[1,2,0,3,4,5,6,7,8]
gate(all(((R[i]>>j)&1)==((R[p[i]]>>p[j])&1) for i in range(9) for j in range(9)),'p3 nontrivial automorphism')
print(json.dumps({'status':'PASS','circulants_checked':total,'histograms':hist,'large_prime_base':{'5^16':5**16,'17^9':17**9,'5*17^9':5*17**9,'18^9':18**9},'p3_boundary_rows':R},sort_keys=True,indent=2))
