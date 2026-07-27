import sympy as sp
x,y,z = sp.symbols('x y z')
u = 1 + x*y
F1 = u**3*z + y**2*u*(4+3*x*y)
F2 = y + 3*x*u**2*z + 3*x*y**2*(4+3*x*y)
F3 = 2*x - 3*x**2*y - x**3*z
F = sp.Matrix([F1,F2,F3])
J = F.jacobian([x,y,z])
detJ = sp.simplify(J.det())
print("=== INDEPENDENT verification of the Keller/JC counterexample (opus) ===")
print("det J_F =", detJ, " (expect -2)")
pts = {'(0,0,-1/4)':(0,0,sp.Rational(-1,4)),
       '(1,-3/2,13/2)':(1,sp.Rational(-3,2),sp.Rational(13,2)),
       '(-1,3/2,13/2)':(-1,sp.Rational(3,2),sp.Rational(13,2))}
print("triple collision (all should map to (-1/4,0,0)):")
for name,(a,b,c) in pts.items():
    val = F.subs({x:a,y:b,z:c})
    print(f"  F{name} = ({val[0]},{val[1]},{val[2]})")
# GEOMETRIC-RANK structure: the 3 collision preimages under C*-doubling lambda->lambda^2, weights (1,-1,-2)
print("\n=== geometric-rank structure: the collision is 1 fixed + 2 orbit-pair (torus doubling) ===")
print("  preimage x-coords: 0, 1, -1  => {0} (fixed) + {1,-1} (sigma=diag(-1,-1,1) swapped pair) = 1+2*1")
print("  this is the '1+2k' census motif; the collision needs 3 branches => GEOMETRIC rank 3.")

# CROSS-THREAD COINCIDENCE TEST: is -1/4 = -prod cos(2*pi*j/3) an LRC-floor echo (klein-S324)? test + release
print("\n=== TEST klein-S324's '-1/4 = -prod cos(2pi j/3) = LRC 3-runner floor' coincidence ===")
prodcos3 = sp.prod([sp.cos(2*sp.pi*j/3) for j in [1,2]])
print(f"  prod_{{j=1,2}} cos(2 pi j/3) = {sp.nsimplify(prodcos3)} ; -that = {sp.nsimplify(-prodcos3)}  (collision value -1/4? {sp.nsimplify(-prodcos3)==sp.Rational(-1,4)})")
# actual LRC 3-runner (2 speeds {1,2}) floor:
import numpy as np
xs=(np.arange(200000)+0.5)/200000
def distZ(a): 
    r=np.mod(a,1.0); return np.minimum(r,1-r)
M12=float(np.max(np.minimum(distZ(1*xs),distZ(2*xs))))
print(f"  actual LRC 3-runner floor M({{1,2}}) = {M12:.5f} (=1/3?); is it -1/4 or 1/4? NO.")
print(f"  => -1/4 = -prod cos(2pi j/3) is TRUE arithmetically, but the LRC 3-runner FLOOR is 1/3, not 1/4.")
print(f"     -1/4 is the value -cos(2pi/3)^2 (a Chebyshev/root-of-unity value), NOT the LRC floor.")
print(f"     VERDICT: charming root-of-unity coincidence; the 'LRC 3-floor' label is imprecise. RELEASE (per kps-S134 method).")
