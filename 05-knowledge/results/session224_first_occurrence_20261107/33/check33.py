from itertools import permutations
n=6
pairs=[(i,j) for i in range(n) for j in range(i+1,n)]; inds={p:k for k,p in enumerate(pairs)}
perms=list(permutations(range(n)))
def transform(m,p):
 r=0
 for k,(i,j) in enumerate(pairs):
  a,b=p[i],p[j]; bit=((m>>k)&1)^(a>b)
  if bit:r|=1<<inds[tuple(sorted((a,b)))]
 return r
def rows(m):
 a=[0]*n
 for k,(i,j) in enumerate(pairs):
  if m>>k&1:a[i]|=1<<j
  else:a[j]|=1<<i
 return a
def hp(a):
 d=[[0]*n for _ in range(1<<n)]
 for i in range(n):d[1<<i][i]=1
 for s in range(1,1<<n):
  for i in range(n):
   if d[s][i]:
    t=a[i]&~s
    while t:
     b=t&-t;t-=b;d[s|b][b.bit_length()-1]+=d[s][i]
 return sum(d[-1])
seen=set(); out=[]
for m in range(1<<15):
 if m in seen:continue
 orb={transform(m,p) for p in perms};seen.update(orb)
 if hp(rows(m))==33:
  out.append(dict(mask=m,rows=rows(m),orbit=len(orb),aut=720//len(orb),converse=min(transform(m^32767,p) for p in perms),scores=sorted(x.bit_count() for x in rows(m))))
print(out)
if len(out)!=5: raise ValueError('class count')
if sum(x['orbit'] for x in out)!=2640: raise ValueError('labelled mass')
if sum(x['mask']==x['converse'] for x in out)!=1: raise ValueError('self-converse count')
