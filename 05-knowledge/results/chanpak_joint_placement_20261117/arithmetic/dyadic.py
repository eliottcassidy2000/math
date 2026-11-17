"""Independent exact polynomial control for the zero-shift collision subsequence."""
from collisions import family,F,H,reverse
import json
checks=0
for k in range(2,9):
 G=2**k;u0,w0,_,_=family(G);u1,w1,_,_=family(G,1)
 U=u1-u0;W=w1-w0;b=max(u0.bit_length(),w0.bit_length());D=2**b
 alpha=256*U*W*D*D
 beta=256*D*(U*w0+W*u0)+142*D*(U+W)+5*(U*(D-reverse(u0,b))+W*(D-reverse(w0,b)))
 gamma=256*u0*w0+142*(u0+w0)+81+5*(F(U)-U+F(u0)-1+F(W)-W+F(w0)-1)
 for j in range(50):
  u,w,_,_=family(G,D*2**j)
  if H(u,w)!=alpha*4**j+beta*2**j+gamma:raise ArithmeticError((k,j))
  checks+=1
print(json.dumps({'exact_polynomial_controls':checks}))
