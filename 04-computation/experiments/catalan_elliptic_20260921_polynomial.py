"""Exact controls for the conditional reading of the user's polynomial formulas.

The integer-point census is bounded, not a complete integral-point theorem.
All decisions remain active under optimized Python. Standard library only.
"""
from fractions import Fraction as Q
from hashlib import sha256
from math import isqrt
from pathlib import Path
import json


def require(test, message):
    if not test:
        raise RuntimeError(message)


def polynomial(x):
    return 2*x**3+4*x*x+1


def rhs(x, coefficients):
    a,b,c=coefficients
    return x**3+a*x*x+b*x+c


def invariants(coefficients):
    a,b,c=coefficients
    cubic_disc=a*a*b*b-4*b**3-4*a**3*c-27*c*c+18*a*b*c
    delta=16*cubic_disc
    c4=16*(a*a-3*b)
    return {"cubic_discriminant":cubic_disc,"elliptic_discriminant":delta,
            "c4":c4,"j":str(Q(c4**3,delta))}


def add(P,R,coefficients=(4,0,4)):
    if P is None:return R
    if R is None:return P
    x,y=P;u,v=R;a,b,_=coefficients
    if x==u and y==-v:return None
    slope=(v-y)/(u-x) if x!=u else (3*x*x+2*a*x+b)/(2*y)
    z=slope*slope-a-x-u
    result=(z,-y+slope*(x-z))
    require(result[1]**2==rhs(result[0],coefficients),"Group-law output off curve")
    return result


def count_points(prime,coefficients):
    points=[(x,y) for x in range(prime) for y in range(prime)
            if (y*y-rhs(x,coefficients))%prime==0]
    residue_counts=[]
    for x in range(prime):
        value=rhs(x,coefficients)%prime
        residue_counts.append(1 if value==0 else 2 if pow(value,(prime-1)//2,prime)==1 else 0)
    require(len(points)==sum(residue_counts),"Independent finite-field point counts disagree")
    return {"count_including_infinity":len(points)+1,"affine_points":points}


def main():
    formula_checks=0
    for numerator in range(-40,41):
        for denominator in range(1,13):
            x=Q(numerator,denominator);b=x+1;c=x*b/2+b*(b+1)/2
            s=2*x**3+3*x*x-x-1+c
            t=-(x*b*b-x)-(b*x*x-b)
            require(c==b*b,"c simplification")
            require(s==x*(2*x*x+4*x+1),"s simplification")
            require(t==-(2*x+1)*(x*x+x-1),"t factorization")
            require(s+t==c and polynomial(x)-s==1-x,"sum or discrepancy")
            require(4*polynomial(x)==rhs(2*x,(4,0,4)),"Weierstrass rescaling")
            formula_checks+=1

    E0=(4,0,4);E4=(109,224,0)
    inv0=invariants(E0);inv4=invariants(E4)
    require(inv0["j"]=="-65536/91" and inv4["j"]=="1408317602329/2153060","j controls")
    require(all(rhs(x,E0)%5!=0 for x in range(5)),"Irreducible cubic mod5")
    require(inv0["cubic_discriminant"]==-1456,"Cubic discriminant control")
    counts={str(p):{"polynomial_curve":count_points(p,E0),"fruit_curve":count_points(p,E4)}
            for p in (3,11,17,19,23)}
    require(counts["11"]["polynomial_curve"]["count_including_infinity"]==16 and
            counts["11"]["fruit_curve"]["count_including_infinity"]==12,"Isogeny hostile count")
    for p in (3,11,17,19,23):
        require(inv0["elliptic_discriminant"]%p!=0 and inv4["elliptic_discriminant"]%p!=0,
                "Point comparison used bad reduction")

    P=(Q(0),Q(2));R=None;multiples=[]
    for n in range(1,13):
        R=add(R,P)
        require(R is not None,"Unexpected finite-order collision in controls")
        multiples.append([str(z) for z in R])
    require(multiples[3]==["20","98"] and multiples[4]==["-24/25","326/125"],
            "Explicit multiple check")
    require(add(add(P,P),add(P,P))==(Q(20),Q(98)),"Independent double-double path")
    for a in range(1,6):
        for b in range(1,6):
            pa=tuple(map(Q,multiples[a-1]));pb=tuple(map(Q,multiples[b-1]))
            require(add(pa,pb)==tuple(map(Q,multiples[a+b-1])),"Addition consistency")

    squares8={y*y%8 for y in range(8)}
    require(all(polynomial(x)%8 not in squares8 for x in (1,3,5,7)),"Odd x obstruction")
    points=[]
    for x in range(-2,200001):
        v=polynomial(x)
        if v>=0 and isqrt(v)**2==v:
            y=isqrt(v);points.append([x,y])
            require(x%2==0 and y%2==1,"Parity of integral point")
            z=x//2;k=(y-1)//2
            require(k*(k+1)==4*z*z*(z+1),"Adjacent-factor model")
    require(points==[[-2,1],[0,1],[10,49]],"Bounded integer-point census changed")
    for q in range(2,101):
        z=q*q-1;A=4*z*q
        require(polynomial(2*z)==A*A+1 and A*A<polynomial(2*z)<(A+1)**2,
                "Square z+1 obstruction")
    output={"status":"PASS","source_sha256":sha256(Path(__file__).read_bytes()).hexdigest(),
            "interpretation":"All pasted juxtaposition/asterisks in b,c,s,t interpreted as multiplication; conditional on intended notation.",
            "universe":{"rational_formula_numerators":[-40,40],"rational_formula_denominators":[1,12],
                        "integer_x_census":[-2,200000],"nonnegative_y_only":True,
                        "finite_field_primes":[3,11,17,19,23],"rational_multiples":[1,12]},
            "formula_controls":formula_checks,"polynomial_curve":inv0,"fruit_curve":inv4,
            "finite_field_counts":counts,"multiples_of_0_2":multiples,
            "integer_points_in_box":points,"complete_integral_point_classification":False,
            "collatz_convergence_proved":False,
            "dependencies":"General elliptic-curve group/reduction theorems are cited in the companion proof note."}
    Path(__file__).with_suffix('.json').write_text(json.dumps(output,indent=2,sort_keys=True)+'\n',
                                                  encoding='utf-8',newline='\n')
    print('PASS: exact polynomial identities, curve invariants, finite-field separation, group operations and bounded integer points')


if __name__=='__main__':main()
