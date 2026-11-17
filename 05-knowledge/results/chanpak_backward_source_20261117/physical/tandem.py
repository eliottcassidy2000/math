from collections import defaultdict
a=(174,124,120,240,65,80,1,118)
b=(174,124,120,240,65,16,33,118)

def exact(a,p,q):
 n=len(a);k=n+2;M=1<<n
 old=a
 arr=[[defaultdict(int) for i in range(k)] for j in range(M)]
 A=n;B=n+1
 a=list(old)
 for j in range(n):
  if not p>>j&1:a[j]|=1<<A
  if not q>>j&1:a[j]|=1<<B
 a.extend([p|1<<B,q])
 for j in range(n):arr[1<<j][j][0,0]=1
 arr[0][A][1,0]=1;arr[0][B][0,1]=1
 # Track increasing r+s then old mask; r,s<=9 since external barriers.
 for rs in range(20):
  for s in range(M):
   for j in range(k):
    for r in range(min(rs,9)+1):
     t=rs-r
     if t>9:continue
     m=arr[s][j].get((r,t),0)
     if not m:continue
     for i in range(k):
      if not (a[j]>>i&1):continue
      if i<n:
       if s>>i&1:continue
       arr[s|1<<i][i][r,t]+=m
      elif i==A and r<9:arr[s][i][r+1,t]+=m
      elif i==B and t<9:arr[s][i][r,t+1]+=m
 C=[]
 for s in range(M):
  z=defaultdict(int)
  for d in arr[s]:
   for key,v in d.items():z[key]+=v
  C.append(z)
 C[0][0,0]=1
 P=defaultdict(int)
 for s in range(M):
  for (r,t),u in C[s].items():
   for (v,w),x in C[M-1-s].items():
    if r+v and t+w:P[r+v,t+w]+=u*x
 return C[-1],P

def coeff(p,q):
 X=exact(a,p,q);Y=exact(b,p,q)
 return [{k:B.get(k,0)-A.get(k,0) for k in (A.keys()|B.keys()) if min(k)>0 and B.get(k,0)!=A.get(k,0)} for A,B in zip(X,Y)]
if __name__=='__main__':
 for p,q in [(12,31),(43,75)]:
  print('PAIR',p,q)
  for d in coeff(p,q):
   for r in range(max(z[0] for z in d)+1):
    print(r,sorted((s,v) for (t,s),v in d.items() if t==r))
   print()
