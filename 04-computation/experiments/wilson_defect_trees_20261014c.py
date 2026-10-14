"""Exact factorial remainders and labelled p-adic product-defect trees.

Correcting a residue section changes the source. It is not a Wilson-prime test
for the unchanged section (1,...,p-1). All controls remain active under -O.
"""
from dataclasses import dataclass, replace
from itertools import product
from math import factorial, gcd, isqrt, lcm, prod
import json


def need(condition, message):
    if not condition:
        raise ValueError(message)


def integer(n, minimum=0):
    need(type(n) is int and n >= minimum, "exact integer in domain")
    return n


def is_prime(n):
    integer(n, 2)
    return all(n % d for d in range(2, isqrt(n)+1))


def odd_prime(p):
    integer(p, 3)
    need(p % 2 and is_prime(p), "odd prime required")
    return p


def factorial_mod(n, modulus):
    integer(n)
    integer(modulus, 1)
    value = 1 % modulus
    for a in range(2, n+1):
        value = value*a % modulus
    return value


def factorial_remainder(n):
    """Complete (n-1)! mod n classification for integer n>=2."""
    integer(n, 2)
    if is_prime(n):
        return n-1
    return 2 if n == 4 else 0


@dataclass(frozen=True)
class WilsonReceipt:
    prime: int
    precision: int
    quotient_residue: int


def wilson_receipt(p, precision=1):
    odd_prime(p)
    integer(precision, 1)
    modulus = p**precision
    numerator = factorial_mod(p-1, p*modulus)+1
    need(numerator % p == 0, "Wilson field congruence")
    return WilsonReceipt(p, precision, (numerator//p) % modulus)


def audit_wilson(receipt, expected_prime=None):
    need(type(receipt) is WilsonReceipt, "typed Wilson receipt")
    odd_prime(receipt.prime)
    integer(receipt.precision, 1)
    integer(receipt.quotient_residue)
    need(receipt.quotient_residue < receipt.prime**receipt.precision, "canonical quotient residue")
    need(receipt == wilson_receipt(receipt.prime, receipt.precision), "unchanged factorial source")
    if expected_prime is not None:
        odd_prime(expected_prime)
        need(receipt.prime == expected_prime, "receipt belongs to another supplied prime")
    return receipt


def fermat_quotient_mod(p, a):
    odd_prime(p)
    integer(a, 1)
    need(a % p, "unit source")
    return ((pow(a, p-1, p*p)-1)//p) % p


def product_mod(values, modulus):
    answer = 1
    for value in values:
        answer = answer*value % modulus
    return answer


@dataclass(frozen=True)
class Section:
    prime: int
    depth: int
    residues: tuple


def audit_section(section):
    need(type(section) is Section, "typed labelled section")
    p = odd_prime(section.prime)
    k = integer(section.depth, 1)
    need(type(section.residues) is tuple and len(section.residues) == p-1, "all nonzero residue labels retained")
    modulus = p**k
    for label, value in enumerate(section.residues, 1):
        integer(value, 1)
        need(value < modulus and value % p == label, "canonical source-labelled residue")
    need(product_mod(section.residues, modulus) == modulus-1, "product-minus-one leaf guard")
    return section


def root_section(p):
    odd_prime(p)
    return audit_section(Section(p, 1, tuple(range(1,p))))


def section_defect(section):
    audit_section(section)
    p, k = section.prime, section.depth
    modulus = p**k
    numerator = product_mod(section.residues, p*modulus)+1
    need(numerator % modulus == 0, "retained previous-depth guard")
    return (numerator//modulus) % p


def lift(section, free_digits):
    """Choose labels1..p-2 freely; solve the label p-1 pivot uniquely."""
    audit_section(section)
    p, k = section.prime, section.depth
    need(type(free_digits) is tuple and len(free_digits) == p-2, "complete free-digit address")
    for digit in free_digits:
        integer(digit)
        need(digit < p, "base-p digit")
    used = sum(digit*pow(label,-1,p) for label,digit in enumerate(free_digits,1)) % p
    last = ((section_defect(section)-used)*(p-1)) % p
    digits = free_digits+(last,)
    child = Section(p, k+1, tuple(value+p**k*digit for value,digit in zip(section.residues,digits)))
    return audit_section(child)


def unlift(section):
    audit_section(section)
    need(section.depth >= 2, "root has no parent")
    modulus = section.prime**(section.depth-1)
    parent = audit_section(Section(section.prime, section.depth-1, tuple(a % modulus for a in section.residues)))
    digits = tuple(a//modulus for a in section.residues[:-1])
    need(lift(parent,digits) == section, "unique pivot reconstruction")
    return parent, digits


def encode(p, address):
    odd_prime(p)
    need(type(address) is tuple, "literal finite address")
    section = root_section(p)
    for digits in address:
        section = lift(section, digits)
    return section


def decode(section):
    audit_section(section)
    address = []
    while section.depth > 1:
        section, digits = unlift(section)
        address.append(digits)
    return tuple(reversed(address))


def ordinary_section(p, depth):
    """Returns None on a failed unchanged-section guard; no representative repair."""
    odd_prime(p)
    integer(depth, 1)
    section = Section(p, depth, tuple(range(1,p)))
    if factorial_mod(p-1,p**depth) != p**depth-1:
        return None
    return audit_section(section)


def teichmueller_section(p, depth):
    odd_prime(p)
    integer(depth, 1)
    values = list(range(1,p))
    modulus = p
    for _ in range(1,depth):
        nxt = []
        for value in values:
            defect = (pow(value,p-1,modulus*p)-1)//modulus
            derivative = ((p-1)*pow(value,p-2,p)) % p
            digit = -defect*pow(derivative,-1,p) % p
            nxt.append(value+modulus*digit)
        values, modulus = nxt, modulus*p
    section = audit_section(Section(p,depth,tuple(values)))
    need(all(pow(a,p-1,modulus)==1 for a in values), "root-of-unity branch")
    return section


def inverse_fixed_points(n):
    integer(n,2)
    return tuple(a for a in range(1,n) if gcd(a,n)==1 and a*a % n == 1)


def unit_product(n):
    integer(n,2)
    return product_mod((a for a in range(1,n) if gcd(a,n)==1),n)


def trim(poly):
    out = list(poly)
    while len(out)>1 and out[-1]==0:
        out.pop()
    return tuple(out)


def poly_mul(a,b,p):
    out = [0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            out[i+j] = (out[i+j]+x*y) % p
    return trim(out)


def poly_mod(a,b,p):
    out = list(trim(tuple(x % p for x in a)))
    while len(out)>=len(b) and any(out):
        c = out[-1]*pow(b[-1],-1,p) % p
        shift = len(out)-len(b)
        for i,x in enumerate(b):
            out[shift+i] = (out[shift+i]-c*x) % p
        out = list(trim(out))
    return tuple(out)


def polynomial_wilson_control(f):
    """Finite F3 cubic experiment only; these helpers are not a prime oracle."""
    need(type(f) is tuple and len(f)==4 and f[-1]==1, "literal monic cubic")
    need(all(type(x) is int and 0<=x<3 for x in f), "F3 coefficients")
    need(all(sum(c*x**i for i,c in enumerate(f))%3 for x in range(3)), "cubic irreducibility by no root")
    modulus = poly_mul(f,f,3)
    value = (1,)
    for coefficients in product(range(3),repeat=3):
        if any(coefficients):
            value = poly_mod(poly_mul(value,trim(coefficients),3),modulus,3)
    second = trim(tuple((i*(i-1)*f[i])%3 for i in range(2,len(f))))
    return value, second


CHECKS=0


def check(condition,message):
    global CHECKS
    CHECKS+=1
    if not condition:
        raise RuntimeError(message)


def rejects(fn):
    try:
        fn()
    except ValueError:
        check(True,"malformed source rejected")
        return
    check(False,"malformed source accepted")


def main():
    remainders=[]
    for n in range(2,2001):
        observed=factorial_mod(n-1,n)
        check(observed==factorial_remainder(n),"full remainder classification")
        check((observed==2)==(n in (3,4)),"exact remainder-two iff")
        if observed==2:remainders.append(n)
    check(factorial(3)%4==2 and not is_prime(4),"composite exceptional source")
    checks_wilson=[]
    for p in range(3,2001,2):
        if not is_prime(p):continue
        r=wilson_receipt(p,2)
        check(audit_wilson(r)==r,"source-bound quotient receipt")
        check(r.quotient_residue==((factorial(p-1)+1)//p)%(p*p),"independent full quotient")
        if r.quotient_residue % p==0:
            checks_wilson.append((p,r.quotient_residue//p))
        if p<=199:
            qs=tuple(fermat_quotient_mod(p,a) for a in range(1,p))
            check(sum(qs)%p==r.quotient_residue%p,"Lerch first-defect identity")
            section=teichmueller_section(p,2)
            check(tuple((x-a)//p for a,x in enumerate(section.residues,1))==tuple(a*q%p for a,q in enumerate(qs,1)),"Teich first-digit formula")
    check([p for p,_ in checks_wilson]==[5,13,563],"bounded Wilson census only")
    check(all(digit!=0 for _,digit in checks_wilson),"named examples fail next zero digit")
    # Exhaustive ternary address/section bijection, every node through depth8.
    levels=[root_section(3)]
    level_counts=[1]
    for depth in range(2,9):
        new=[]
        for parent in levels:
            children=tuple(lift(parent,(d,)) for d in range(3))
            check(len(set(children))==3,"three distinct children")
            for d,child in enumerate(children):
                check(unlift(child)==(parent,(d,)),"unique inverse edge")
                check(encode(3,decode(child))==child,"full address inverse")
                new.append(child)
        levels=new
        check(len(set(levels))==3**(depth-1),"complete ternary level count")
        level_counts.append(len(levels))
    # Exhaustive independent tuple enumeration at selected small levels.
    complete=[]
    for p,depth in ((3,2),(3,3),(3,4),(3,5),(5,2)):
        direct=[]
        for digits in product(range(p**(depth-1)),repeat=p-1):
            values=tuple(a+p*d for a,d in enumerate(digits,1))
            if prod(values) % p**depth == p**depth-1:
                section=audit_section(Section(p,depth,values))
                check(encode(p,decode(section))==section,"independent section coverage")
                direct.append(section)
        check(len(direct)==p**((p-2)*(depth-1)),"p-ary free coordinate dimension")
        complete.append((p,depth,len(direct)))
    for p in (3,5,7,11,13,17,19):
        for depth in range(1,9):
            section=teichmueller_section(p,depth)
            check(encode(p,decode(section))==section,"canonical branch lies in lift tree")
            if depth>=2:
                check(unlift(section)[0]==teichmueller_section(p,depth-1),"coherent canonical branch")
    first3=tuple(lift(root_section(3),(d,)).residues for d in range(3))
    check(first3==((1,8),(4,2),(7,5)),"literal ternary children")
    check(ordinary_section(3,2) is None and first3[0]==teichmueller_section(3,2).residues,"source-changing repair hostile")
    for p in (5,13,563):
        check(ordinary_section(p,2) is not None and ordinary_section(p,3) is None,"unchanged path stops at exact depth2")
    check(teichmueller_section(5,2).residues==(1,7,18,24),"aggregate zero does not kill each correction")
    check(audit_wilson(wilson_receipt(13),13).quotient_residue==wilson_receipt(5).quotient_residue==0,"same verdict does not identify the source")
    # Unit involution products versus full-ring factorials and Korselt.
    for n in range(2,501):
        fixed=inverse_fixed_points(n)
        observed=unit_product(n)
        check(observed==product_mod(fixed,n),"inverse-pair cancellation")
        expected=fixed[1] if len(fixed)==2 else 1%n
        check(observed==expected,"elementary-two-group fixedpoint product")
    check(unit_product(561)==1 and factorial_remainder(561)==0,"unit versus full-ring source")
    check(all(pow(a,561,561)==a for a in range(561)),"separate all-ring Korselt control")
    check(unit_product(15)==1 and pow(2,15,15)!=2,"unit product does not imply Korselt")
    check(unit_product(9)==8 and factorial_remainder(9)==0,"unique involution does not imply prime")
    for p in (5,13,563):
        n=p-2
        check(all(pow(a,n,n)==a%n for a in range(n)),"three explicit shifted absolute-Fermat examples only")
    successor_factors=(3,11,17,401,3041)
    check(all(is_prime(p) for p in successor_factors) and prod(successor_factors)==684106401,"marked parent successor factorization")
    check(lcm(*(p-1 for p in successor_factors))==15200 and 15201==27*563,"source-specific 563 arithmetic link")
    # Same formal Wilson definition, genuinely different coefficient domain.
    positive=(2,2,0,1);negative=(1,2,1,1)
    pos=polynomial_wilson_control(positive)
    neg=polynomial_wilson_control(negative)
    check(pos==((2,),(0,)),"F3 cubic Wilson positive")
    check(neg[0]!=(2,) and neg[1]!=(0,),"F3 cubic non-Wilson negative")
    check(polynomial_wilson_control((1,2,0,1))[0]==(2,),"second derivative-zero prime")
    field_census=[]
    for lower in product(range(3),repeat=3):
        f=lower+(1,)
        if all(sum(c*x**i for i,c in enumerate(f))%3 for x in range(3)):
            value,second=polynomial_wilson_control(f)
            check((value==(2,))==(second==(0,)),"complete F3 irreducible cubic comparison")
            field_census.append((f,value==(2,)))
    check(len(field_census)==8 and sum(flag for _,flag in field_census)==2,"finite polynomial Wilson census")
    base=root_section(3);child=lift(base,(0,));wr=wilson_receipt(5)
    hostiles=(lambda:factorial_remainder(True),lambda:factorial_remainder(1),
        lambda:odd_prime(9),lambda:wilson_receipt(2),lambda:wilson_receipt(5,1.0),
        lambda:audit_wilson(replace(wr,prime=13),5),lambda:audit_wilson(replace(wr,quotient_residue=False)),
        lambda:audit_section(replace(base,prime=3.0)),lambda:audit_section(replace(base,depth=True)),
        lambda:audit_section(replace(base,residues=(1,2.0))),lambda:audit_section(replace(base,residues=(2,1))),
        lambda:audit_section(Section(3,2,(1,2))),lambda:lift(base,(True,)),lambda:lift(base,(3,)),
        lambda:lift(base,()),lambda:unlift(base),lambda:encode(3,[(0,)]),
        lambda:audit_section(replace(child,residues=(1,5))),lambda:fermat_quotient_mod(5,5),
        lambda:polynomial_wilson_control((0,0,0,1)))
    for fn in hostiles:rejects(fn)
    print(json.dumps({"status":"PROVED scoped classifications and reversible section tree; global integer Wilson distribution OPEN",
        "checks":CHECKS,"remainder_universe":[2,2000],"remainder_two_sources":remainders,
        "wilson_prime_bound":2000,"wilson_prime_and_next_digit":checks_wilson,
        "ternary_level_counts_depth1_through8":level_counts,"ternary_depth2_sections":first3,
        "independent_full_section_counts":complete,"teichmueller5_mod25":teichmueller_section(5,2).residues,
        "unit_product_universe":[2,500],"f3_positive":positive,"f3_negative":negative,
        "f3_irreducible_cubic_census":field_census,"shifted_wilson_sources":[3,11,561],
        "successor_marker":{"source":684106401,"lambda":15200,"lambda_plus_one":"27*563"},
        "malformed_controls":len(hostiles),
        "scope":"Corrected sections vary representatives; they do not prove a fixed prime Wilson, a fixed Collatz source native, or ROOT coverage."},indent=2,sort_keys=True))


if __name__=="__main__":
    main()
