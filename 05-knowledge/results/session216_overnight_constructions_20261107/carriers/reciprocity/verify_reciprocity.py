from reciprocity_probe import *
import json,hashlib
seen={};checks=0;collisions=0;x=sp.Symbol('x')
for n in range(1,6):
 for bits in range(1<<(n*(n-1)//2)):
  g=[0]*n;z=0
  for i in range(n):
   for j in range(i+1,n):
    if bits>>z&1:g[i]|=1<<j
    else:g[j]|=1<<i
    z+=1
  F,G=polynomials(g)
  need(all(c.q==1 for c in F.all_coeffs()+G.all_coeffs()),'integer coefficients')
  need(sp.expand(F.as_expr().subs(x,1-x)-(-1)**(n-1)*G.as_expr())==0,'rooted reciprocity')
  p=profile(g);need((int(F.eval(1)),int(F.eval(2)),int(F.eval(3)),int(G.eval(1)),int(G.eval(2)))==(p[1],p[4],p[8],p[2],p[3]),'profile identification')
  key=(n%2,p[2],p[3],p[4],p[8])
  if key in seen:need((p[1]-seen[key])%4==0,'same-parity obstruction');collisions+=1
  else:seen[key]=p[1]
  checks+=1
print(json.dumps({'status':'PASS','labeled_rooted_polynomials':checks,'same_parity_collisions':collisions,'max_order':5},indent=2))
