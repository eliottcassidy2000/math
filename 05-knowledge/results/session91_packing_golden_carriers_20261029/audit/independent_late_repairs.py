#!/usr/bin/env python3
"""Independent exact checks of late precision and linear-fibre repairs."""
import json, math
N=0
def ck(x,m):
 global N
 N+=1
 if not x: raise RuntimeError(m)
def main():
 for p in (2,3,5,7,11):
  for k in (1,2,3):
   m=p**k
   # Primitive norms: all residues for smaller levels, deterministic corner set at large levels.
   vals=range(m) if m<=125 else range(60)
   minimum=None
   for a in vals:
    for b in vals:
     if a%p==b%p==0:continue
     d=a*a+a*b-b*b;J=m*m//math.gcd(m,d)
     if p in (2,3,7):ck(J==m*m,'inert norm unit')
     if p==5:ck(math.gcd(m,d) in (1,5),'ramified valuation one')
     if p==11:ck(J>=m,'split lower bound')
     minimum=J if minimum is None else min(minimum,J)
   if p==5:ck(minimum==5**(2*k-1),'ramified minimum')
 for n in range(4,15,2):
  W=[s for s in range(1<<n) if all(not(s>>i&1 and s>>((i+1)%n)&1) for i in range(n))]
  legal=set(W);bad={s^(1<<i) for s in W for i in range(n)}-legal
  alts=[sum(1<<i for i in range(par,n,2)) for par in (0,1)]
  for a in alts:
   neighbors={z for z in bad if (a^z).bit_count()==1}
   ck(len(neighbors)==n//2,'alternating diagonal')
   ck(all(sum((s^z).bit_count()==1 for s in W)==1 for z in neighbors),'alternating row isolated')
  # Reader zero-fibre sum plus these two rows recovers all three coefficients.
  for x,y,z in ((1,2,3),(-5,7,0),(0,0,1)):
   out=(x+y+z,(n//2)*y,(n//2)*z)
   yy=out[1]//(n//2);zz=out[2]//(n//2)
   ck((out[0]-yy-zz,yy,zz)==(x,y,z),'linear source recovery')
 print(json.dumps({'status':'PASS','checks':N,'scope':'prime-power norm controls and even-cycle isolated alternating rows n4..14'},indent=2))
if __name__=='__main__':main()
