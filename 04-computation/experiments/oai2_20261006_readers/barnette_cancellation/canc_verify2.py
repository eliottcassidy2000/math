from sympy import symbols, expand, cancel, fraction, Poly
p, s, u, w, L_, M_, Zs = symbols('p s u w L M Zs')
x0 = s**2 + u**3
FL = 4*s**2*L_ + (1 + 2*s*x0)*M_
JL = -(1 - 2*s*x0)*L_ + x0**2*M_
xL = x0 + p**2*FL
HL = expand(xL**2*FL - (1 + 2*s*xL)*JL - p**2*JL**2 - p*u)
Q = 2*x0*FL**2 - 2*s*FL*JL - JL**2
C = expand(Q.subs(L_, 0))
num, den = fraction(cancel((Q - C)/L_)); print("Q1 denominator:", den)
Q1 = expand(num/den)
Ls = p*u - p**2*C
num, den = fraction(cancel(HL.subs(L_, Ls)/p**3)); print("H(L*)/p^3 denominator:", den)
num, den = fraction(cancel(-HL.subs(L_, Ls + p**3*Zs)/p**3)); print("W(Z) denominator:", den)
W = expand(num/den)
print("deg_Z W =", Poly(W, Zs).degree(), "; number of terms of W:", len(W.as_ordered_terms()))
# explicit coordinates of T = P[w]/(H + p^3 w) are p,s,u,M,e with e0 = -w - pF^3 - Q1 h; check L = L* + p^3 e and w = W(e) mod (H+p^3w)
h = u - p*Q - p**3*FL**3 - p**2*w
e0 = expand(-w - p*FL**3 - Q1*h)
print("terms in e0:", len(e0.as_ordered_terms()))
