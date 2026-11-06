"""Optional independent symbolic controls; requires SymPy, not a theorem dependency."""
import sympy as s
S,L,i,z,b=s.symbols('S L i z b')
e=[s.Integer(1)]
for r in range(1,6):
 power=[None]+[s.summation((S-i)**a,(i,1,L)) for a in range(1,r+1)]
 er=s.expand(sum((-1)**(a-1)*e[r-a]*power[a] for a in range(1,r+1))/r)
 e.append(er)
 assert s.Poly(er,S,L).total_degree()<=2*r
 for ell in range(r,8):
  for size in range(-2,10):
   poly=s.prod(z-(size-j) for j in range(1,ell+1))
   assert s.Poly(poly,z).coeff_monomial(z**(ell-r))==(-1)**r*er.subs({S:size,L:ell})
q=s.expand(e[2].subs({S:3*b,L:b})-2*b*(b-1)*e[1].subs({S:3*b-1,L:b-1})+5*b*(b-1)*(b-2)+2*b*(b-1)*(b-2)*(b-3))
expected=b*(b-1)*(3*b*b+17*b+2)/24
assert s.expand(q-expected)==0
assert q.subs(b,2)/s.factorial(2)==2
assert q.subs(b,3)/s.factorial(3)==s.Rational(10,3)
print('PASS symbolic coefficient/degree controls through r5 and real-inventory Q2(b,0)=b(b-1)(3b²+17b+2)/24')
