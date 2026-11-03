from itertools import product

def brute(m, x, y=None, xy=1):
 n=m+1+(y is not None); out=[0]*n
 for i in range(m):
  for j in range(i+1,m):out[i]|=1<<j
 for i in range(m):
  if x>>i&1:out[m]|=1<<i
  else:out[i]|=1<<m
 if y is not None:
  for i in range(m):
   if y>>i&1:out[m+1]|=1<<i
   else:out[i]|=1<<(m+1)
  if xy:out[m]|=1<<(m+1)
  else:out[m+1]|=1<<m
 d=[[0]*n for _ in range(1<<n)]
 for i in range(n):d[1<<i][i]=1
 for S in range(1,1<<n):
  for i in range(n):
   if d[S][i]:
    for j in range(n):
     if not S>>j&1 and out[i]>>j&1:d[S|1<<j][j]+=d[S][i]
 return sum(d[-1])
def formula(m,x):
 h=2
 for i in range(m):
  if not x>>i&1:
   if i+1<m and x>>(i+1)&1:h+=1
   for j in range(i):
    if x>>j&1:h+=1<<(i-j-1)
 return h
cases=0
for m in range(2,9):
 for x in range(1<<m):
  if x&1 and not x>>(m-1)&1:
   if brute(m,x)!=formula(m,x):
    raise AssertionError((m,x))
   cases+=1
print('one-apex subset-DP checks:', cases)
m=12;x=2663;y=1740
if brute(m,x,y,1)!=80269:
 raise AssertionError("incorrect witness path count")
out=[0]*14
for i in range(m):
 for j in range(i+1,m):out[i]|=1<<j
for i in range(m):
 for apex,mask in ((12,x),(13,y)):
  if mask>>i&1:out[apex]|=1<<i
  else:out[i]|=1<<apex
out[12]|=1<<13
for i in range(14):
 seen=1<<i
 while True:
  new=seen
  for j in range(14):
   if seen>>j&1:new|=out[j]
  if new==seen:break
  seen=new
 if seen!=(1<<14)-1:
  raise AssertionError("witness is not strongly connected")
print('14-vertex witness: H=80269, strongly connected')
