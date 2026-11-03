import sys

def endpoints(a):
 n=len(a); d=[[0]*n for _ in range(1<<n)]
 for i in range(n): d[1<<i][i]=1
 for mask in range(1,1<<n):
  for i in range(n):
   if d[mask][i]:
    for j in range(n):
     if not (mask>>j&1) and a[i][j]: d[mask|1<<j][j]+=d[mask][i]
 return sum(d[-1]), d[-1]

def add(a,v):
 n=len(a)
 return [row+[int(i==v)] for i,row in enumerate(a)]+[[int(i!=v) for i in range(n)]+[0]]
T=[[0,1,0],[0,0,1],[1,0,0]]
limit=int(sys.argv[1]) if len(sys.argv)>1 else 17
states=[(T,3,())]; hits={3:()}
for n in range(3, limit+4):
 new=[]
 for a,h,w in states:
  for v in range(n):
   b=add(a,v); H=endpoints(b)[0]
   if H <= h: raise RuntimeError('strict growth failed')
   if H<=limit:
    ww=w+(v,); new.append((b,H,ww));hits.setdefault(H,ww)
 print('order',n+1,'states',len(new),'counts',sorted(set(H for _,H,_ in new)),flush=True)
 states=new
 if not states:break
print('hits',sorted(hits.items()));print('missing odd',[i for i in range(3,limit+1,2) if i not in hits])
