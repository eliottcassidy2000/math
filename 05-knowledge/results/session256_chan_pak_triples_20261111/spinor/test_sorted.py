import contextlib,io,json,math,runpy
from itertools import product
from pathlib import Path
with contextlib.redirect_stdout(io.StringIO()):
 D=runpy.run_path(str(Path(__file__).with_name('test_spinor.py')))
def poly(p,q,a):
 s=p+q+p*q;t=p+2*q+q*q+2*p*q+p*q*q
 return {'A':(p,s),'B':(t,s),'C':(t,q)}[a]
def inv(N):
 out=[]
 for p in range(1,math.isqrt(N-2)+1):
  if (N-2)%p:continue
  q=(N-2)//p
  for u,v in {(p,q),(q,p)}:
   if (u+1)&u or (v+1)&v:continue
   x=(u+1).bit_length()-1;y=(v+1).bit_length()-1
   if x%2 and math.gcd(x,y)==1:out.append((x,y))
 return sorted(out)
num=direct=0
for k in range(7):
 for tup in product('ABC',repeat=k):
  a=''.join(tup);m,n=D['spin'](a);x,y=m-n,n;p=q=1
  for c in a:p,q=poly(p,q,c)
  assert (p,q)==((1<<x)-1,(1<<y)-1)
  H=2+p*q
  assert H%3!=1
  assert D['count']('O'*x+'I'*y)==H
  if m<=13:assert D['tournament_h']('O'*x+'I'*y)==H;direct+=1
  num+=1
for N in range(3,512,2):
 expected=sorted((x,y) for x in range(1,10,2) for y in range(1,10)
                 if math.gcd(x,y)==1 and 2+((1<<x)-1)*((1<<y)-1)==N)
 assert inv(N)==expected
m,n=D['spin']('CA');morphic=D['count'](D['word']('CA'))
assert (m,n,morphic)==(7,4,73)
print(json.dumps({'addresses_through_depth_6':num,'independent_sorted_graph_counts':direct,
 'odd_inverse_targets_3_through_511':255,'CA':{'m':m,'n':n,'morphic_H':morphic,
 'sorted_H':2+((1<<(m-n))-1)*((1<<n)-1)},
 'sorted_image_through_511':[N for N in range(3,512,2) if inv(N)]},indent=2))
