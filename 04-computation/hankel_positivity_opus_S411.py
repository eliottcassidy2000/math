import sympy as sp, numpy as np
from math import factorial
z,zb=sp.symbols('z zb')
def E1(e):
    e=sp.expand(e)
    if e==0: return sp.Integer(0)
    p=sp.Poly(e,z,zb); t=0
    for (a,b),c in zip(p.monoms(),p.coeffs()):
        if a==b: t+=c*factorial(a)
    return sp.expand(t)
print("THE KEY STEP:  if P has all charges >= 0 and P is in the nullcone, then its")
print("charge-0 part P0 must VANISH -- so P is strictly charge-definite.")
print()
print("Proof.  With only nonnegative charges present, the charge-0 part of P^2 can only")
print("come from P_q * P_{-q} with both present, i.e. q = 0.  So  E[P^2] = E[P0^2].")
print("Write P0 = sum_a c_a (z zbar)^a.  Then")
print("      E[P0^2] = sum_{a,b} c_a c_b (a+b)!  =  c^T H c,   H_{ab} = (a+b)!")
print("H is the HANKEL MOMENT MATRIX of the exponential distribution (moments a!),")
print("hence POSITIVE DEFINITE.  So E[P^2] = 0 forces c = 0, i.e. P0 = 0.  QED")
print()
print("(i) H_{ab} = (a+b)! is positive definite -- eigenvalues for sizes 1..7:")
for N in range(1,8):
    H=np.array([[float(factorial(a+b)) for b in range(N)] for a in range(N)])
    ev=np.linalg.eigvalsh(H)
    print(f"    N={N}: min eigenvalue = {ev.min():.6g}   all > 0: {bool(ev.min()>0)}")
print()
print("(ii) check E[P^2] = E[P0^2] on nonneg-charge examples, and that it is > 0:")
for P0c,extra in [([1,-1],z), ([2,-1,0],z**2), ([1,0,-1],z+z**3), ([3,-2,0,0],zb*0+z)]:
    P0=sum(c*(z*zb)**a for a,c in enumerate(P0c))
    P=sp.expand(P0+extra)
    lhs=E1(sp.expand(P**2)); rhs=E1(sp.expand(P0**2))
    quad=sum(P0c[a]*P0c[b]*factorial(a+b) for a in range(len(P0c)) for b in range(len(P0c)))
    print(f"    P0={P0}  E[P^2]={lhs}  E[P0^2]={rhs}  c^T H c={quad}  "
          f"{'MATCH' if lhs==rhs==quad else 'mismatch'}   E[P0]={E1(P0)}")
print()
print("=> CONCLUSION: GMC(2) is PROVED for every P whose charges are all of ONE SIGN")
print("   (weakly).  Such P in the nullcone are strictly charge-definite, and the charge")
print("   lemma then gives E[Q P^m] = 0 for all m > (charge deficit of Q).")
print("   REMAINING: P with charges of BOTH signs -- verified exhaustively, not proved.")
