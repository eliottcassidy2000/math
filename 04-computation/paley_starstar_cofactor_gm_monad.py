#!/usr/bin/env python3
"""
paley_starstar_cofactor_gm_monad.py
monad-explorer-2026-06-07 (deep-research, 11th session)

UNIFICATION of the two THM-438 handoffs on a SINGLE cofactor polynomial.

Let Q_m(s) = sum_{e=m}^{2m-1} R_s(m,e) s^e  be the s-expansion of column m
(Q_m(s) = T_m(x), s=x/(1-x);  Q_m = s^m (1+s)^{m-1} P_m(s/(1+s)) ).

CLAIM (verified m<=4 exact; m=5 modulo c2,c3,c4):
   Q_m(s) = s^m * (1+s) * g_m(s),   deg g_m = m-2,
where g_m(s) = Q_m(s)/(s^m(1+s)) is the COFACTOR.  Then:
   * handoff #1  (deg P_m=m-2)   <=>  (1+s) | Q_m   <=>  g_m is a polynomial
                                 <=>  Q_m(-1)=0.
   * handoff #2  (lead P_m=2^m-1) <=>  g_m(-1) = (-1)^m (2^m - 1)
                                 [ = -((-2)^m) ... actually (-1)^m(2^m-1) = (-2)^m-(-1)^m ].
   * KNOWN:  g_m(0)   = A088368(m)        (the WILD/factorial end, ~ e*m!)
             lead g_m = P_m(1) = R_s(m,2m-1) (the top pole residue).
So the SAME degree-(m-2) polynomial g_m carries BOTH ends:
   g_m(0)  = A088368(m)   [factorial, wild],
   g_m(-1) = (-1)^m(2^m-1) [Mersenne, tame].
The "tame<->wild bridge" (ADD-8 reflection) is literally g_m's two evaluations.

Consequence / sharpened handoff: (1+s)|Q_m is the ALGEBRAIC SHADOW of a
fixed-point-free involution on whatever R_s(m,e) enumerates, flipping e-parity.
This re-localizes the involution search onto the s-expansion (binomial-transform)
objects, with the involution = the factor (1+s).
"""
import sympy as sp

s = sp.symbols('s')
c2,c3,c4 = sp.symbols('c2 c3 c4')
A088368 = {1:1,2:3,3:13,4:69,5:421,6:2867,7:22417}
Q = {
 1: s,
 2: 3*s**2 + 3*s**3,
 3: 13*s**3 + 33*s**4 + 20*s**5,
 4: 69*s**4 + 304*s**5 + 416*s**6 + 181*s**7,
 5: 421*s**5 + 2740*s**6 + (5694+c2)*s**7 + (4852+2*c2+c3)*s**8 + (1477+c2+c3+c4)*s**9,
}

print("="*74)
print("Q_m(s) = s^m (1+s) g_m(s) ;  g_m = cofactor (deg m-2)")
print("="*74)
for m in range(1,5):
    Qm = Q[m]
    # divide by s^m(1+s)
    quo, rem = sp.div(sp.Poly(Qm, s), sp.Poly(s**m*(1+s), s))
    print(f"\nm={m}:  Q_m = {sp.factor(Qm)}")
    print(f"     remainder of Q_m/(s^m(1+s)) = {rem.as_expr()}  (0 <=> handoff#1 holds)")
    g = quo.as_expr()
    print(f"     g_m(s) = {g}   (deg {sp.Poly(g,s).degree() if g!=0 else 0})")
    print(f"       g_m(0)   = {g.subs(s,0)}   (A088368(m)={A088368[m]})")
    print(f"       g_m(-1)  = {g.subs(s,-1)}   (want (-1)^m(2^m-1)={(-1)**m*(2**m-1)})")
    lead = sp.Poly(g,s).all_coeffs()[0] if g!=0 else g
    print(f"       lead g_m = {lead}   (P_m(1)=R_s(m,2m-1))")

print("\n" + "="*74)
print("m=5 (symbolic c2,c3,c4): the handoffs PIN the cofactor")
print("="*74)
Qm = Q[5]
print("Q_5(-1) =", sp.expand(Qm.subs(s,-1)), " => handoff#1 <=> c4 = 0")
# under c4=0:
Q5 = Qm.subs(c4,0)
quo, rem = sp.div(sp.Poly(Q5, s), sp.Poly(s**5*(1+s), s))
print("with c4=0: remainder =", sp.expand(rem.as_expr()), " (0 confirms (1+s)|Q_5)")
g5 = sp.expand(quo.as_expr())
print("g_5(s) =", g5)
print("g_5(0) =", g5.subs(s,0), "(=A088368(5)=421)")
print("g_5(-1) =", sp.expand(g5.subs(s,-1)), " ; handoff#2 wants (-1)^5(2^5-1) = -31")
sol = sp.solve(sp.Eq(g5.subs(s,-1), -31), c3)
print("   handoff#2 => c3 =", sol, " (predicts c3=31 if independent of c2? check:)")
print("   NOTE: g_5(-1) =", sp.expand(g5.subs(s,-1)), "depends on c2,c3 -> ONE eqn, not enough alone.")
print("   lead g_5 = P_5(1) =", sp.expand(sp.Poly(g5,s).all_coeffs()[0]))
