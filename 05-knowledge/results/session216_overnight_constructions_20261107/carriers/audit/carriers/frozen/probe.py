import sys,itertools,json
from inherited import *
def add(a,b):
 c=[0]*max(len(a),len(b))
 for i,x in enumerate(a):c[i]+=x
 for i,x in enumerate(b):c[i]+=x
 return c
def mul(a,b):
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return c

def polynomials(rows,v):
 # root=0, v replaced by arbitrary module; skeleton runs indexed by k
 n=len(rows); ext=[i for i in range(n) if i!=v]; N=1<<len(ext); full=N-1
 def paths(allowed,which):
  # count paths across all external subsets, with optional repeated module runs
  dp={}; out=[[0]*(n+3) for _ in range(N)];out[0][0]=int(which=='h')
  for j,u in enumerate(ext):
   if allowed>>j&1 and (which!='s' or u==0):dp[(1<<j,j,0)]=1
  if which!='s':dp[(0,len(ext),1)]=1
  for weight in range(2*n+3):
   for (S,last,k),a in list(dp.items()):
    if S.bit_count()+k!=weight:continue
    u=v if last==len(ext) else ext[last]
    if which!='e' or u==0:out[S][k]+=a
    for j,w in enumerate(ext):
     if allowed>>j&1 and not(S>>j&1) and rows[u]>>w&1:
      key=(S|1<<j,j,k);dp[key]=dp.get(key,0)+a
    if u!=v and rows[u]>>v&1:
     key=(S,len(ext),k+1);dp[key]=dp.get(key,0)+a
  return out
 allp=paths(full,'h'); sp=paths(full,'s');ep=paths(full,'e')
 def convolution(a,b,S):
  out=[0]*(2*n+5);A=S
  while True:
   t=mul(a[A],b[S^A]);out=add(out,t)
   if A==0:break
   A=(A-1)&S
  return out
 rootbit=1<<ext.index(0);D=full^rootbit
 vals=[allp[full],sp[full],ep[full],convolution(ep,allp,full),convolution(sp,allp,full),allp[D],convolution(allp,allp,full),convolution(allp,allp,D)]
 # only k>=1 represent covering the nonempty substituted module
 return [x[1:n+2]+[0]*max(0,n+1-len(x[1:n+2])) for x in vals]

def expand(rows,v,m):
 ex=[i for i in range(len(rows)) if i!=v]; mapping=ex+[v]*m
 return tuple(sum(( (i<j if u==v and w==v else rows[u]>>w&1))<<j for j,w in enumerate(mapping) if i!=j) for i,u in enumerate(mapping))

if __name__=='__main__':
 pairs=[]
 for X,js,r,a,b in PAIRS:pairs.append((X,remark(a,r),remark(b,r)))
 for X,a,r,b,t,w in [(4225,(2,20,41,67,77,27,167,63),6,(6,12,24,49,195,23,47,111),4,'JO'),(4357,(10,20,65,22,37,15,59),3,(2,48,11,67,45,77,23),2,'ZJO')]:
  a,b=remark(a,r),remark(b,t)
  for c in w:a,b=graph_move(a,c),graph_move(b,c)
  pairs.append((X,a,b))
 records=[]
 for X,a,b in pairs:
  banks=[]
  for rows in [a,b]:
   bank=[]
   for v in range(1,len(rows)):
    p=polynomials(rows,v)
    s=p[1]
    if any(s[1:]):continue
    residue=s[0]%3;j=(s[0]-residue)//3
    if residue==0 or not 0<=j<=8:continue
    xp=[h+l-sx for h,l,sx in zip(p[0],p[3],p[1])];xp[0]+=2*j
    key=tuple(tuple(z) for z in [xp,p[2],add(p[4],p[7]),p[5]])
    bank.append((key,v,p))
   banks.append(bank)
  for ka,va,pa in banks[0]:
   for kb,vb,pb in banks[1]:
    if ka==kb:
     rec={'base':X,'a':a,'b':b,'va':va,'vb':vb,'key':ka,'pa':pa,'pb':pb};records.append(rec)
     print(X,va,vb,'nonconstant',any(any(z[1:]) for z in ka),flush=True)
 print(json.dumps(records),file=open('matches.json','w'))
