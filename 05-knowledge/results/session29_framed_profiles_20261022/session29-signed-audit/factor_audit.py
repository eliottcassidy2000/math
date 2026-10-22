import itertools as it,math,json
from literal import fixed

def calc(m,q,bits):
 edges=list(it.combinations(range(m),2));ctrl={};b={};T=[0]*(2*m)
 for k,(v,w) in enumerate(edges):
  if not(q>>k&1):v,w=w,v
  ctrl[frozenset((v,w))]=v;b[v,w]=1 if bits>>k&1 else -1
  for s in range(2):
   for t in range(2):
    x,y=2*v+s,2*w+t
    if (1 if s==0 else -1)==-b[v,w]:T[x]|=1<<y
    else:T[y]|=1<<x
 for v in range(m):T[2*v]|=1<<(2*v+1)
 f=0;U=0;H=0
 for p in it.permutations(range(m)):
  z=1;c=0
  for j in range(1,m-1):
   u,v,w=p[j-1:j+2]
   if ctrl[frozenset((u,v))]==v==ctrl[frozenset((v,w))]:z*=1-b[v,u]*b[v,w];c+=1
  if m>1 and ctrl[frozenset(p[-2:])]==p[-1]:z*=1+b[p[-1],p[-2]];c+=1
  f+=z;U+=2**c
  H+=all(ctrl[frozenset((u,v))]==u for u,v in zip(p,p[1:]))
 assert fixed(T)==f
 return f,U,H
checks=0
for m in range(1,5):
 E=m*(m-1)//2
 for q in range(1<<E):
  acc=0
  for bits in range(1<<E):
   f,U,H=calc(m,q,bits);acc+=f;checks+=1
   assert f>=H and f%2==1
   if bits==0:assert f==H
   if m>=4:assert f<U
  assert acc==math.factorial(m)*(1<<E)
Path=__import__('pathlib').Path
Path(__file__).with_name('factor_results.json').write_text(json.dumps({'all_gauge_checks_through_m4':checks,'mean_m_factorial':True,'all_minus_embedding':True,'strict_relaxation_gap_m4':True},indent=2)+'\n')
print('factor checks',checks)
