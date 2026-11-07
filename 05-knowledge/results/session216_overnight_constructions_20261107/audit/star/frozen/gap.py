from itertools import permutations

def rows(n,bits):
 a=[0]*n;k=0
 for i in range(n):
  for j in range(i+1,n):
   if bits>>k&1:a[i]|=1<<j
   else:a[j]|=1<<i
   k+=1
 return a

def coeff(a,v):
 n=len(a); c=[0]*(n+1)
 for p in permutations([i for i in range(n) if i!=v]):
  pol=[1]
  for i in range(n):
   l=p[i-1] if i else None;r=p[i] if i<n-1 else None
   e=int(l is None or r is None or a[l]>>r&1)
   f=int((l is None or a[l]>>v&1) and (r is None or a[v]>>r&1))
   q=[0]*(len(pol)+1)
   for j,x in enumerate(pol):q[j]+=e*x;q[j+1]+=f*x
   pol=q
  for j,x in enumerate(pol):c[j]+=x
 return c
if __name__=='__main__':
 for n in range(3,6):
  seen={}
  for b in range(1<<(n*(n-1)//2)):
   a=rows(n,b)
   for v in range(n):
    c=coeff(a,v)
    if c[2] and not any(c[3:]):seen.setdefault((c[1],c[2]),(a,v,c))
  print(n,seen)
