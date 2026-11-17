from sympy import symbols, resultant, expand, Poly
m,z,w=symbols('m z w')
D=(2*m-1+z)**2-4*(2-z)
P=z**5+z**4-4*z**3-3*z**2+3*z+1
N=resultant(P,D,z)
Q=32*m**5-96*m**4-96*m**3+172*m**2+164*m+23
R=704*m**4+5280*m**3+10384*m**2+7920*m+2112
checks=[expand(Q**2-N-R)==0,
        all(c>0 for c in Poly((N-(Q-1)**2).subs(m,w+19),w).all_coeffs()),
        all(c>0 for c in Poly(Q.subs(m,w+19),w).all_coeffs()),
        3642956**2<N.subs(m,11)<3642957**2]
if not all(checks): raise AssertionError
print('11 norm sieve verified')
