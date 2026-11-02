from fractions import Fraction
import random,json

def det(m):
 a=[r[:] for r in m];n=len(a);old=1;sign=1
 for k in range(n-1):
  if not a[k][k]:
   j=next((j for j in range(k+1,n)if a[j][k]),None)
   if j is None:return 0
   a[k],a[j]=a[j],a[k];sign=-sign
  piv=a[k][k]
  for i in range(k+1,n):
   for j in range(k+1,n):
    value=a[i][j]*piv-a[i][k]*a[k][j]
    if value%old:raise RuntimeError('Bareiss divisibility')
    a[i][j]=value//old
   a[i][k]=0
  old=piv
 return sign*a[-1][-1] if n else 1

def strong(a):
 n=len(a);r=[x|(1<<i)for i,x in enumerate(a)]
 for k in range(n):
  for i in range(n):
   if r[i]>>k&1:r[i]|=r[k]
 return all(x==(1<<n)-1 for x in r)

def invariants(a):
 n=len(a);A=[[int(a[i]>>j&1)for j in range(n)]for i in range(n)]
 u=det([[A[i][j]+int(i==j)for j in range(n)]for i in range(n)])
 v=(-1)**n*det(A);D=u+v;s=Fraction(u-v,D)
 K=[[int(i==j)+A[i][j]-A[j][i]for j in range(n)]for i in range(n)]
 if Fraction(det(K),2**(n-1))!=D:raise RuntimeError('determinant identity')
 return {'u':u,'v':v,'disc':D,'s':str(s)}

def main():
 rng=random.Random(137)
 for trial in range(20000):
  n=7;a=[0]*n
  for i in range(n):
   for j in range(i+1,n):
    u,v=(i,j)if rng.randrange(2)else(j,i);a[u]|=1<<v
  if not strong(a):continue
  r=invariants(a)
  if r['v']>0:
   joined=[x|(((1<<n)-1)<<n) for x in a]+[x<<n for x in a]
   q=invariants(joined)
   if not q['disc']<r['disc']**2:raise RuntimeError('expected failure')
   return {'trial':trial,'strong_rows':a,'factor':r,'joined':q,'product_disc':r['disc']**2}
 raise RuntimeError('not found')
if __name__=='__main__':print(json.dumps(main(),indent=2))
