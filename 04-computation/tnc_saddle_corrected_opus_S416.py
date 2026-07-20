import sympy as sp
u=sp.symbols('u')
print("CORRECTED ARGUMENT.  A GENUINE saddle is a root of S(u) := u R'(u) - N R(u) that is")
print("NOT a root of R (else log R is singular there and both sides vanish trivially).")
print("If a genuine saddle u* exists, the saddle value R(u*)/u*^N is NONZERO, so the")
print("Cauchy integral has an m-th-power term that cannot vanish for all m.")
print()
print("SO A TNC VIOLATOR NEEDS: every root of S is a root of R.  Degrees:")
print("   deg R = M+N; leading coeff of S = (M+N)r_d - N r_d = M r_d, so deg S = M+N too.")
print("   'every root of S is a root of R' with equal degrees and leading ratio M")
print("   ==> S = M*R  ==>  u R' - N R = M R  ==>  u R' = (M+N) R  ==>  R = c u^{M+N}.")
print("   But R(0) != 0 (min exponent is exactly -N), so M+N = 0.  CONTRADICTION.")
print("   => A GENUINE SADDLE ALWAYS EXISTS.  The algebraic obstruction is EMPTY.")
print()
print("verify the ODE step and that genuine saddles exist in the test cases:")
for lab,N,R in [("N=1,R=1+u",1,1+u),("N=1,R=(1+u)^2",1,(1+u)**2),
                ("N=2,R=(1+u)^4",2,(1+u)**4),("N=2,R=1+u+u^2+u^3+u^4",2,1+u+u**2+u**3+u**4),
                ("N=2,R=(1-u)^2(1+u)^2",2,(1-u)**2*(1+u)**2),("N=3,R=(1+u)^6",3,(1+u)**6)]:
    Rx=sp.expand(R); S=sp.expand(u*sp.diff(Rx,u)-N*Rx)
    g=sp.gcd(sp.Poly(Rx,u),sp.Poly(S,u)).as_expr()
    quo=sp.simplify(sp.cancel(S/g)); Rquo=sp.simplify(sp.cancel(Rx/g))
    # genuine saddles = roots of S/gcd that are not roots of R
    genuine=[r for r in sp.roots(sp.Poly(quo,u)) if sp.simplify(Rx.subs(u,r))!=0]
    print(f"   {lab:24s} S/gcd = {sp.factor(quo)}   genuine saddles: {genuine}")
print()
print("check the ODE claim: u R' = (M+N) R has only monomial solutions")
c,k=sp.symbols('c k'); f=sp.Function('f')
sol=sp.dsolve(sp.Eq(u*sp.Derivative(f(u),u),k*f(u)),f(u))
print(f"   dsolve(u f' = k f) -> {sol}   (pure power u^k, so R(0)=0 unless k=0)")
