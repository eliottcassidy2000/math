from itertools import product
from replay_axis import C
Q=(4,2,5,3,7)
E=(12,51,53,52)
def mul(a,b,m):
    return tuple(sum(a[2*i+k]*b[2*k+j] for k in (0,1))%m
                 for i,j in product(range(2),repeat=2))
def clock(Y,m):
    a=(1,0,0,1)
    for q,y in zip(Q,Y):
        e,f,g=C(q)
        a=mul((e,f,y*f,y*g),a,m)
    return a
def st(Y):
 A=mul(clock(Y[1:],63),E,63)
 return ((A[0]+Y[0]*A[1])%63,(A[2]+Y[0]*A[3])%63)
base=[pow(2,r,63) for r in range(6)]
pts={}
for b in product(base,repeat=6):pts[b]=st(b)
lines=[(1,x) for x in range(9)]+[(x,1) for x in range(0,9,3)]
lines2=[(1,x) for x in range(7)]+[(0,1)]
records=[]
for b,z in pts.items():
 delta=[tuple(v-u for u,v in zip(z,pts[b[:i]+(2*b[i]%63,)+b[i+1:]])) for i in range(6)]
 bad3=next((i for i,(a,c) in enumerate(((0,1),(1,0),(1,1),(1,2))) if all((a*x+c*y)%3==0 for x,y in delta)),4)
 bad7=next((i for i,(a,c) in enumerate(lines2) if all((a*x+c*y)%7==0 for x,y in delta)),8)
 nonzero3=any((x%3 or y%3) for x,y in delta)
 nonzero7=any((x%7 or y%7) for x,y in delta)
 if nonzero3 and nonzero7: records.append((z,bad3,bad7))
def cr(a,b):return (28*a+36*b)%63
print('good rec',len(records))
failed=[]
for a,b in lines:
 i3=b%3+1 if a==1 else 0
 for i7,(c,d) in enumerate(lines2):
  R=(cr(a,c),cr(b,d));nz=set()
  for (x,y),A,B in records:
   if A!=i3 and B!=i7: nz.add((R[0]*x+R[1]*y)%63)
   if len(nz)==63:break
  if b%3==0:
   assert not nz
  else:
   assert len(nz)==63,(R,len(nz))
assert all(clock(Y,3)[0]==clock(Y,3)[1]==0 for Y in product((1,2),repeat=5))
print('every residue for all 72 good projective rows, and failure of the other 24 rows')
