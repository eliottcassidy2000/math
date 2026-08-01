"""Two things: (a) the log substitution turning the factorial functional into an
integral over [0,1]; (b) degree-1 transcendence of int_0^1 e^{Q}."""
from mpmath import mp, mpf, quad, exp, log, factorial as mfac
mp.dps=30
print("(a) THE LOG BRIDGE:  L(g) = int_0^inf g(x) e^{-x} dx = int_0^1 g(-log u) du")
print("    (x = -log u, dx = -du/u, e^{-x} = u)")
for j in range(6):
    lhs = mfac(j)
    rhs = quad(lambda u: (-log(u))**j, [0,1])
    print(f"    j={j}: L(s^j) = {j}! = {int(lhs)}   int_0^1 (-log u)^j du = {mp.nstr(rhs,12)}")
print()
print("    So the FACTORIAL weight and the [0,1] weight are the SAME functional")
print("    in different coordinates -- but the exponent becomes polynomial in")
print("    log u, NOT in u.  int_0^1 e^{Q(t)}dt with Q polynomial is a DIFFERENT")
print("    object; the bridge does not by itself produce it.")
print()
print("(b) DEGREE-1 TRANSCENDENCE, from Lindemann-Weierstrass")
print("    int_0^1 e^{at+b} dt = (e^{a+b} - e^b)/a.")
print("    Suppose it equals an algebraic gamma.  Then")
print("        1*e^{a+b} + (-1)*e^b + (-a*gamma)*e^0 = 0,")
print("    a Qbar-linear relation among e^alpha for the DISTINCT algebraic")
print("    exponents a+b, b, 0.  LW says those exponentials are linearly")
print("    independent over Qbar, forcing the coefficient 1 to vanish --")
print("    contradiction.  Edge cases: b=0 uses exponents a,0; a+b=0 gives")
print("    a*gamma = 1 - e^b with e^b transcendental.  So the value is")
print("    TRANSCENDENTAL for every algebraic a != 0, b.  QED")
print()
for a,b in ((1,0),(2,1),(mpf(1)/3,mpf(2)/5)):
    v=(exp(a+b)-exp(b))/a
    print(f"    a={a}, b={b}: value = {mp.nstr(v,20)}")
