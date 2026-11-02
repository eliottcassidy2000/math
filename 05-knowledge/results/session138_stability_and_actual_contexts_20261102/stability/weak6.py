from bernstein import *
from reserved_cycles import triangles
n=6;m=15;h,u,v,V=polys(n);t=triangles(n);a=48*(h-u-t)-37*V
for o in representatives(n):
 b=bernstein(a,m,o)
 print(json.dumps({'n':n,'type':'weak_triangle','k':'37/48','orientation':o,'min':int(b.min()),'negatives':int((b<0).sum()),'zero':int((b==0).sum()),'scale':48*2**(n+m),'argmin':int(b.argmin())}),flush=True)
