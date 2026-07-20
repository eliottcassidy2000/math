"""opus-2026-07-20-S417 (HYP-8145): the owner's claim-2 {w=0} fiber formula --
verified in full symbolic generality -- and the SHEET DECOMPOSITION THEOREM:
F^{-1}({w=0}) = {x=0} (tame sheet, maps by a plane AUTOMORPHISM) u {T=0}
(exotic sheet, 2:1 branched over the conic v^2 = 16u = the w=0 slice of
kp-THM-1310's Jelonek quartic {L=0}), with the deck involution derived.
Credits: kp THM-1310/1335 (master cubic, Jelonek), S414 anatomy (C = xT)."""
import sympy as sp
x,y,z,u_,v,s_ = sp.symbols('x y z u v s')
U = 1+x*y
A = U**3*z + y**2*U*(4+3*x*y)
B = y + 3*x*U**2*z + 3*x*y**2*(4+3*x*y)
C = 2*x - 3*x**2*y - x**3*z
# claim 1 recap is canon (S414). CLAIM 2, general symbolic with u = (v^2-s^2)/16:
uu = (v**2 - s_**2)/16
p0 = {x:0, y:v, z:uu - 4*v**2}
im0 = [sp.simplify(f.subs(p0)) for f in (A,B,C)]
print("tame point F(0, v, u-4v^2) =", im0, "== (u,v,0):",
      [sp.simplify(im0[0]-uu)==0, sp.simplify(im0[1]-v)==0, im0[2]==0])
res = []
for xp in (2/s_, -2/s_):
    pe = {x:xp, y:v/4 - 3/(2*xp), z:sp.simplify(13/(2*xp**2) - 3*v/(4*xp))}
    ime = [sp.simplify(f.subs(pe)) for f in (A,B,C)]
    res.append([sp.simplify(ime[0]-uu)==0, sp.simplify(ime[1]-v)==0, sp.simplify(ime[2])==0])
print("exotic points map to (u,v,0):", res)
# SHEET DECOMPOSITION: C = x*T with T = 2-3xy-x^2 z
T = 2 - 3*x*y - x**2*z
print("C == x*T:", sp.expand(C - x*T) == 0)
# tame sheet {x=0}: F(0,y,z) = (z+4y^2, y, 0) -- a plane shear automorphism
print("tame sheet map:", [sp.expand(f.subs({x:0})) for f in (A,B,C)],
      " inverse: (u,v) -> (v, u-4v^2)")
# exotic sheet {T=0}: z = (2-3xy)/x^2; derive B and A there
zT = (2 - 3*x*y)/x**2
Bt = sp.simplify(B.subs(z, zT)); At = sp.simplify(A.subs(z, zT))
print("on {T=0}: B =", sp.expand(Bt), " => y = B/4 - 3/(2x):",
      sp.simplify(Bt - (4*y + 6/x)/1) , "check y-law:",
      sp.simplify(sp.solve(sp.Eq(Bt, v), y)[0] - (v/4 - 3/(2*x))) == 0)
Au = sp.simplify(At.subs(y, v/4 - 3/(2*x)))
print("on {T=0} with the y-law: A =", sp.expand(Au),
      " => x^2 = 4/(v^2 - 16A):", sp.simplify(sp.expand(Au) - (v**2 - 4/x**2)/16) == 0)
# deck involution: x -> -x with y -> v/4 + 3/(2x) = y + 3/x
print("deck involution iota(x,y) = (-x, y + 3/x) preserves {T=0} and (A,B):",
      sp.simplify(Bt.subs({x:-x, y:y+3/x}) - Bt) == 0,
      sp.simplify(At.subs({x:-x, y:y+3/x}) - At) == 0)
# fiber census ON the conic: v^2 = 16u => s = 0, exotic escapes; tame remains
pc = {x:0, y:4, z:sp.Rational(1)-64}  # u = 1, v = 4: v^2 = 16u
print("on-conic example (u,v)=(1,4): tame preimage maps to",
      [sp.simplify(f.subs(pc)) for f in (A,B,C)], "(fiber drops 3 -> 1: exotic x = +-2/s escapes as s->0)")
