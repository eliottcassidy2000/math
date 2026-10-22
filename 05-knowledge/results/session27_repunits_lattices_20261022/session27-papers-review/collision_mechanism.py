import sympy as s
z,h=s.symbols('z h'); mass={2:-2,3:4,4:3,5:-8,6:-5,7:12,8:1,9:-8,10:3};P=sum(v*z**k for k,v in mass.items());print('rank-mass polynomial',s.factor(P));print('power sums', [sum(v*k**m for k,v in mass.items()) for m in range(5)])
cA=-2*h*(14*(h-1)**3-37*(h-1)**2+39*(h-1)-15);delta=-10*h*(h-2)*(h-4)
# D5 mass15 on R2, D6 mass6 on R1; D per-rank fourth=A per-rank fourth-delta
print('A proxy c4 difference',sum(v*cA.subs(h,k) for k,v in mass.items()));print('D corrections',-15*delta.subs(h,8)+6*delta.subs(h,10))
print('h fourth moment',sum(v*k**4 for k,v in mass.items()),'leading cA',s.Poly(cA,h).LC())
k=s.symbols('k');a5=s.Rational(4,5)*(k+1)*(42*k**4-175*k**3+301*k*k-245*k+85);d5=s.Rational(16,5)*(k-1)*(336*k**4-2800*k**3+8876*k*k-12705*k+6920)
a5h=s.expand(a5.subs(k,h-1));delta5=s.factor(a5h-d5.subs(k,h/2+1));print('c5 normalized A-D defect',delta5);print('fifth power mass',sum(v*k**5 for k,v in mass.items()));print('A c5 contribution',sum(v*a5h.subs(h,k) for k,v in mass.items()));print('D c5 correction',-15*delta5.subs(h,8)+6*delta5.subs(h,10))
