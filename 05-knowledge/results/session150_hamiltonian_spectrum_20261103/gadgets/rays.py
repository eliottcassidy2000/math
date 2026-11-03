from test import *
r={}
for n in range(2,6):
 for code in range(1<<(n*(n-1)//2)):
  a=graph(n,code)
  for v in range(n):
   h,s,q,f=profile(a,v); A=h-s+2*q-f; B=f-q
   if B>0:r.setdefault((A,B),(n,code,v,h,s,q,f))
for k in sorted(r,key=lambda z:(z[1],z[0])): print(k, r[k])
