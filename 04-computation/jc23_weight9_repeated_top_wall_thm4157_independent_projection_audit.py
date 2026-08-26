#!/usr/bin/env python3
"""Clean-room exact audit of the exact-M=9 repeated top-edge candidate.

This is deliberately independent of the producing script.  The affine
critical scheme is audited by a primitive linear coordinate Z=X+cT and a
direct resultant of (G_X,G_T), retaining all four universal points rather
than dividing by T or changing to (s,p).
"""

from hashlib import sha256
import json
import sympy as S


X,T,Z,s,p,a,z,Q,tau = S.symbols("X T Z s p a z Q tau")


def G_source(delta, theta, phi, eta):
    P = T + X**2*T**2
    Y = X*T*P
    k = S.Rational(2848,45)-S.Rational(7,6)*delta
    return S.expand(
        -X**2*T/2 - 3*P + S.Rational(8,3)*P**2
        - S.Rational(1376,135)*P**3 + k*Y**2
        + phi*P**2*Y + delta*P**4 + theta*P*Y**2
        + eta*P**3*Y - eta*Y**3)


def poly_digest(poly):
    pp=S.Poly(poly).monic()
    payload="degree=%d;"%pp.degree()+",".join(
        "%s/%s"%(S.Rational(c).p,S.Rational(c).q) for c in pp.all_coeffs())
    return sha256(payload.encode()).hexdigest()


def primitive_projection(delta,theta,phi,eta,c=1):
    g=G_source(delta,theta,phi,eta)
    gx=S.diff(g,X).subs(X,Z-c*T)
    gt=S.diff(g,T).subs(X,Z-c*T)
    rr=S.Poly(S.resultant(gx,gt,T),Z)
    sq=S.gcd(rr,rr.diff())
    pgx=S.Poly(gx,T)
    pgt=S.Poly(gt,T)
    lc_gcd=S.gcd(S.Poly(pgx.LC(),Z),S.Poly(pgt.LC(),Z))
    return g,rr,sq,(pgx.degree(),pgt.degree(),lc_gcd)


CONTROLS = {
    "A": (S.Rational(1),S.Rational(19,11),S.Rational(11,7),S.Rational(23,13),23,
          "41af700e36f2b726e31fb021149c38efa323308163763aa29dddac9b4ceeadaa"),
    "B": (S.Rational(1),S.Rational(-1),S.Rational(11,7),S.Rational(23,13),21,
          "c49fcaa4a6147c33138a06672ba5c051cac88e134bcdf4865a5a504dcf1760a9"),
    "C": (S.Rational(1),S.Rational(-1),S.Rational(0),S.Rational(23,13),19,
          "15f85f45e684b17f239a68a5a64dff205110eb3471edb99abfd67c3a30a866d5"),
    "D": (S.Rational(2048,45),-S.Rational(2048,45),S.Rational(0),S.Rational(23,13),16,
          "5972821cd259985e8c5d979179d010b0073fe91b6837db6e18a4190d24b13e39"),
}


if __name__ == "__main__":
    ledger=[]
    u0=S.Poly(Z**2+6,Z)
    u1=S.Poly((Z+S.Rational(1,6))**2-6,Z)
    universal=u0*u1
    for name,(d,th,ph,e,L,expected_digest) in CONTROLS.items():
        g,r,gg,sidecar=primitive_projection(d,th,ph,e)
        got_digest=poly_digest(r)
        assert r.degree()==L and gg.degree()==0 and got_digest==expected_digest
        assert sidecar[2].degree()==0
        q,rem=S.div(r,universal)
        assert rem.is_zero and q.degree()==L-4
        hess=S.det(S.hessian(g,(X,T)))
        for tv,x2,gv,hv in ((S.Rational(0),-S.Rational(6),S.Rational(0),S.Rational(6)),
                            (-S.Rational(1,6),S.Rational(6),S.Rational(1,2),-S.Rational(6))):
            modulus=S.Poly(X**2-x2,X)
            for expr in (S.diff(g,X),S.diff(g,T),g-gv,hess-hv):
                assert S.rem(S.Poly(S.expand(expr.subs(T,tv)),X),modulus).is_zero
        ledger.append((name,L,q.degree(),sidecar[0],sidecar[1],got_digest))
    payload="\n".join("|".join(map(str,row)) for row in ledger)+"\n"
    print("INDEPENDENT PRIMITIVE-PROJECTION CRITICAL AUDIT")
    print("projection=Z=X+T;eliminate=T from (G_X,G_T)")
    print("infinity_sidecar=gcd of the two T-leading coefficient polynomials in Z")
    for name,L,residual,deg1,deg2,dig in ledger:
        print(f"row={name};degree={L};squarefree=True;universal=4;residual={residual};"
              f"T_degrees=({deg1},{deg2});infinity_gcd=1;monic_sha256={dig}")
    print("universal=T=0,X^2=-6,G=0,Hessian=6;T=-1/6,X^2=6,G=1/2,Hessian=-6")
    print("semantic_sha256="+sha256(payload.encode()).hexdigest())
    print("verdict=PASS")
