import json

def need(c,m):
 if not c:raise RuntimeError(m)
def attach(s,w,k):
 n=k+len(s)+1;r=n-1;a=[0]*n
 for u in range(k):
  for v in range(u+1,n-1):a[u]|=1<<v
  a[r]|=1<<u
 for i,row in enumerate(s):
  u=k+i;a[u]|=row<<k
  if w>>i&1:a[u]|=1<<r
  else:a[r]|=1<<u
 return a

def data(a):
 n=len(a);N=1<<n;d=[[0]*n for _ in range(N)];h=[0]*N;h[0]=1
 for u in range(n):d[1<<u][u]=1
 for m in range(1,N):
  for u in range(n):
   if not m>>u&1:continue
   q=m^(1<<u)
   if q:d[m][u]=sum(d[q][v] for v in range(n) if q>>v&1 and a[v]>>u&1)
  h[m]=sum(d[m])
 return h,d

def rooted(a,r):
 h,d=data(a);full=(1<<len(a))-1
 return h[-1],sum(d[m][r]*h[full^m] for m in range(1,full+1))
def strong(a):
 n=len(a)
 for u in range(n):
  seen=1<<u
  while True:
   nxt=seen
   for v in range(n):
    if seen>>v&1:nxt|=a[v]
   if seen==nxt:break
   seen=nxt
  if seen!=(1<<n)-1:return False
 return True

def run():
 hostile=[]
 for s in [(66,68,121,99,107,67,0),(66,68,121,67,107,75,0)]:
  a=attach(s,84,1);h,f=rooted(a,8);need(h==423 and f-h==112,'initial hostile')
  b=list(s);need(b[4]>>5&1,'available named arc');b[4]^=1<<5;b[5]^=1<<4
  aa=attach(b,84,1);hh,ff=rooted(aa,8);need(strong(a) and strong(aa),'hostile strong')
  hostile.append({'before':h,'after':hh,'defect':ff-hh})
 need([v['after'] for v in hostile]==[611,471],'hostile divergence')
 source=attach([sum(1<<j for j in range(i+1,5)) for i in range(5)],21,3)
 need(rooted(source,8)==(157,169),'source physical root')
 t=attach((66,68,121,67,75,91,0),80,1);need(rooted(t,8)==(475,521),'bulk target')
 n=len(t);u=7;r=8;x=n;y=n+1
 out=list(t)+[1<<u,1<<x]
 for v in range(n):
  if v!=u:out[v]|=1<<x
  out[v]|=1<<y
 need(strong(t) and strong(out),'prepared strength')
 need(rooted(out,r)==(1429,1563),'actual births')
 need(all((out[i]&((1<<n)-1))==t[i] for i in range(n)),'old induced graph')
 return {'status':'PASS','hostile':hostile,'physical_chain':[157,475,1429],'physical_root_F':[169,521,1563],'scope':'Independent small literal checks for NOTE.md, not replay of the incoming all-clock native compiler.'}
if __name__=='__main__':print(json.dumps(run(),indent=2,sort_keys=True))
