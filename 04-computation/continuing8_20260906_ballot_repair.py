"""Exact ballot/Fuss ratio repair and GLOBAL bronze classification.
All identities are elementary factorial cancellations. No producer imports.
"""
from pathlib import Path
from fractions import Fraction as F
from math import comb
import hashlib
import json
import sys
import sympy as S

sys.stdout.reconfigure(newline="\n")
GATES = 0


def need(ok, why):
    global GATES
    GATES += 1
    if not ok:
        raise ArithmeticError(why)


k, a, b, delta, DD = S.symbols("k a b delta DD")
R = (2*k+a)*(2*k+a-1)*(k+b+1)*(k+a-b+1)/((k+b)*(k+a-b)*(2*k+a+2)*(2*k+a+1))
Q = k*k+(a+1-(a-2*b)**2)*k+(a*(a+2)-(2*a+1)*(a-2*b)**2)/4
den = (k+b)*(k+a-b)*(2*k+a+2)*(2*k+a+1)
need(S.cancel(1-R-2*Q/den) == 0, "all-shift exact monic ballot numerator")
qdelta = k*k+(a+1-delta**2)*k+(a*(a+2)-(2*a+1)*delta**2)/4
need(S.expand(16*(qdelta.subs(k,0)+1).subs(a,(DD-1)/2)
              -(DD*DD+2*DD+13-4*DD*delta**2)) == 0, "global metallic Diophantine reduction")
solutions = []
for divisor in (-13, -1, 1, 13):
    d2 = F(divisor+2+F(13,divisor),4)
    need(d2 == (-3 if divisor < 0 else 4), "complete signed divisor alternatives")
    if d2 >= 0:
        aa = (divisor-1)//2
        for dd in (-2, 2):
            bb = (aa-dd)//2
            polynomial = S.expand(Q.subs({a: aa, b: bb}))
            need(S.degree(S.gcd(polynomial, den.subs({a:aa,b:bb})),k) == 0,
                 "no cancellation in metallic alternatives")
            solutions.append((aa, bb, str(polynomial)))
need(sorted((aa,bb) for aa,bb,pq in solutions if aa == 0) == [(0,-1),(0,1)],
     "only positive-metallic alternatives are bronze")

# All valid interior rows in the old finite rectangle, using literal binomials.
for aa in range(-2,7):
    for bb in range(-3,5):
        for kk in range(max(1,1-bb,1-aa+bb),12):
            vals = [comb(2*j+aa,j+bb) for j in (kk-1,kk,kk+1)]
            need(all(vals), "index-valid positive binomial triple")
            need(R.subs({a:aa,b:bb,k:kk}) == F(vals[1]**2,vals[0]*vals[2]),
                 "literal all-shift rational ratio")
            need(den.subs({a:aa,b:bb,k:kk}) > 0, "positive uncancelled denominator")

displayed = [
    (0,0,1/(k*(2*k+1))),
    (0,-1,(k*k-3*k-1)/((k-1)*(k+1)**2*(2*k+1))),
    (1,-1,(k*k-7*k-6)/((k-1)*(k+1)*(k+2)*(2*k+3))),
]
for aa,bb,formula in displayed:
    need(S.cancel(1-R.subs({a:aa,b:bb})-formula) == 0, "old displayed ballot row is correct")
catalan_ratio = 2*(2*k-1)/(k+1)
need(S.cancel(1-catalan_ratio/catalan_ratio.subs(k,k+1)-3/((k+1)*(2*k+1))) == 0,
     "old Catalan row is correct")
need((k*k-7*k-6).subs(k,7) < 0 < (k*k-7*k-6).subs(k,8),
     "another column changes concavity after the displayed short bank")

orders = []
for pp in range(2,9):
    adjacent = S.cancel(S.prod(pp*(k-1)+j for j in range(1,pp+1))
                        /(k*S.prod((pp-1)*(k-1)+j for j in range(1,pp))))
    fuss = S.cancel(adjacent*((pp-1)*(k-1)+1)/((pp-1)*k+1))
    for kk in range(1,10):
        need(adjacent.subs(k,kk) == F(comb(pp*kk,kk),comb(pp*(kk-1),kk-1)),
             "integer-order binomial adjacent ratio")
        f0 = F(comb(pp*(kk-1),kk-1),(pp-1)*(kk-1)+1)
        f1 = F(comb(pp*kk,kk),(pp-1)*kk+1)
        need(fuss.subs(k,kk) == f1/f0, "integer-order Fuss-Catalan adjacent ratio")
    change = S.factor(1-fuss/fuss.subs(k,k+1))
    orders.append({"p":pp,"binomial_adjacent":str(adjacent),"Fuss_adjacent":str(fuss),"Fuss_1_minus_R":str(change)})
    if pp == 4:
        need(S.degree(S.fraction(change)[0],k) == 4, "Fuss numerator need not remain quadratic")

# The coefficient-polynomial reciprocal stratum can maximally alternate.
roots = [1,1,3,3,9,9]
e = [1]
for root in roots:
    e.append(0)
    for j in range(len(e)-1,0,-1):
        e[j] += root*e[j-1]
newton = [F(e[j]**2*comb(6,j-1)*comb(6,j+1),e[j-1]*e[j+1]*comb(6,j)**2) for j in range(1,6)]
expected = [F(65,57),F(4693,4005),F(71289,61009),F(4693,4005),F(65,57)]
need(newton == expected, "canonical reciprocal polynomial normalized Newton ratios")
word = [1 if q>p else -1 if q<p else 0 for p,q in zip(newton,newton[1:])]
need(word == [1,-1,1,-1], "reciprocal antipalindrome is ALSO maximally alternating")
need(word == [-v for v in word[::-1]], "correct antipalindrome normalization")
need(sorted(F(9,r) for r in roots) == roots, "scaled reciprocal root multiset")

result = {"global_Q":str(Q),"denominator":str(den),"metallic_alternatives":solutions,
          "Fuss_orders":orders,"reciprocal_hostile":{"roots":roots,"e":e,
          "normalized_ratios":list(map(str,newton)),"circuit_word":word}}
here = Path(__file__).resolve().parent
dest = here.parent/"05-knowledge/results" if here.name == "04-computation" else here
path = dest/(Path(__file__).stem+"_certificate.json")
path.write_bytes((json.dumps(result,indent=2,sort_keys=True)+"\n").encode("utf-8"))
print("BALLOT all integer shifts; index-valid triple; monic degree-two numerator before cancellation")
print("GLOBAL_BRONZE",solutions,"positive metallic n>=1 only(a,b)=(0,-1),(0,1)")
print("FUSS integer orders all rational; p4 has quartic Newton numerator")
print("RECIPROCAL_HOSTILE ratios",list(map(str,newton)),"circuit",word)
print("CERTIFICATE_SHA256",hashlib.sha256(path.read_bytes()).hexdigest())
print("PASS",GATES,"always-active exact gates; actual LF")
