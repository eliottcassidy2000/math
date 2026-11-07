import itertools,json

def req(b,m):
 if not b:raise RuntimeError(m)
def rows(n,b):
 a=[0]*n;k=0
 for i in range(n):
  for j in range(i+1,n):
   if b>>k&1:a[i]|=1<<j
   else:a[j]|=1<<i
   k+=1
 return a
def tables(a):
 n=len(a);e=[[0]*n for _ in range(1<<n)];s=[[0]*n for _ in range(1<<n)];h=[0]*(1<<n);h[0]=1
 for v in range(n):e[1<<v][v]=s[1<<v][v]=1
 for m in range(1,1<<n):
  for v in range(n):
   if not(m>>v&1):continue
   if m!=(1<<v):
    e[m][v]=sum(e[m^(1<<v)][u] for u in range(n) if a[u]>>v&1)
    s[m][v]=sum(s[m^(1<<v)][u] for u in range(n) if a[v]>>u&1)
  h[m]=sum(e[m])
 return h,s,e
def apex(a,cut):
 n=len(a);b=a.copy()+[cut]
 for i in range(n):
  if not(cut>>i&1):b[i]|=1<<n
 return b
def whole(L,R,cutL,cutR,k):
 l,r=len(L),len(R);n=l+k+r+1;x=n-1;b=[0]*n
 for i in range(n-1):
  for j in range(i+1,n-1):b[i]|=1<<j
 for i in range(l):b[i]=(b[i]&~((1<<l)-1))|L[i]
 for i in range(r):
  z=l+k+i;mask=((1<<r)-1)<<(l+k);b[z]=(b[z]&~mask)|(R[i]<<(l+k))
 for i in range(n-1):
  out=(cutL>>i&1) if i<l else (cutR>>(i-l-k)&1) if i>=l+k else 0
  if out:b[x]|=1<<i
  else:b[i]|=1<<x
 return b
checks=0
for l in range(4):
 for lb in range(1<<(l*(l-1)//2)):
  L=rows(l,lb);hl,sl,el=tables(L);fl=(1<<l)-1
  for r in range(4):
   for rb in range(1<<(r*(r-1)//2)):
    R=rows(r,rb);hr,sr,er=tables(R);fr=(1<<r)-1
    for cl in range(1<<l):
     alpha=sum(hl[fl^q]*sum(sl[q][v] for v in range(l) if cl>>v&1) for q in range(1,fl+1))
     tau=sum(sum(el[p][v] for v in range(l) if cl>>v&1)*sum(sl[fl^p][v] for v in range(l) if cl>>v&1) for p in range(1,fl))
     for cr in range(1<<r):
      beta=hr[-1]+sum(sum(er[p][v] for v in range(r) if not(cr>>v&1))*hr[fr^p] for p in range(1,fr+1))
      d=hl[-1]*tables(apex(R,cr))[0][-1]-hr[-1]*tau
      for k in [1,2]:
       val=tables(whole(L,R,cl,cr,k))[0][-1]
       req(val==alpha*beta*(1<<k)+d,'factor formula')
       checks+=1
print(json.dumps({'status':'PASS','literal_checks':checks}))
