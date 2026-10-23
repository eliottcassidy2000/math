from pathlib import Path
import json
P=[2,3,5,7,11,13]; e=[5,7,3,4,1,2]
N=1
for p,a in zip(P,e):N*=p**a
assert N==39046339332000
f=dict(zip(P,e))
def v(n,p):
 a=0
 while n%p==0:n//=p;a+=1
 return a
assert [v(N,p) for p in P]==e
bad={}
for m in [3,9]:
 bad[m]=[(p,q,f[p],f[q]) for p in P for q in P if p<q and (p-q)%m==0 and (f[p]-f[q])%m]
 assert bad[m]
assert (5-2)%3==0 and (3-5)%3!=0
quot=lambda x:13 if x==4 else x
for number,period in [(N,5),(N//2,2)]:
 d={quot(p):quot(v(number,p)) for p in P};x=13;orbit=[]
 while x not in orbit:orbit.append(x);x=d[x]
 assert x==13 and len(orbit)==period
 print(number,orbit)
out={'status':'PASS exact integer checks','N':N,'valuations':f,'modular_descent_witnesses':bad,'set_quotient_cycles':{'N':[13,2,5,3,7],'N_over_2':[13,2]},'scope':'The identification4=13 is a set quotient, not integer equality or a ring-compatible dynamical law.'}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
