"""Primitive prime-power Hecke support, Cartan frames and exact s7 traces.

Finite checks support the proofs in the companion note. No routine infers
truncated-series supercongruences from the trace. Stdout is the JSON receipt.
"""
from dataclasses import dataclass, asdict, replace
from fractions import Fraction
from math import gcd, isqrt
import json


def integer(value, minimum=None):
    if type(value) is not int or (minimum is not None and value < minimum):
        raise ValueError("exact integer outside domain")
    return value


def prime(p):
    integer(p, 3)
    if p % 2 == 0 or any(p % d == 0 for d in range(3, isqrt(p) + 1, 2)):
        raise ValueError("odd prime required")
    return p


def parameters(p, r):
    prime(p)
    integer(r, 1)
    return p ** r


def val(n, p):
    integer(n)
    prime(p)
    if not n:
        return None
    answer = 0
    while n % p == 0:
        n //= p
        answer += 1
    return answer


def rational_val(x, p):
    x = Fraction(x)
    return None if not x else val(x.numerator, p) - val(x.denominator, p)


@dataclass(frozen=True)
class Budget:
    p: int
    r: int
    branch_counts: tuple
    degree: int
    product_valuation: int
    root_valuation: Fraction
    trace_power: int


def budget(p, r):
    parameters(p, r)
    counts = (1,) + tuple(p ** j - p ** (j - 1) for j in range(1, r)) + (p ** r,)
    degree = p ** (r - 1) * (p + 1)
    product = 2 * sum(j * count for j, count in enumerate(counts))
    slope = Fraction(product, degree)
    return Budget(p, r, counts, degree, product, slope, -(-slope.numerator // slope.denominator))


def audit_budget(packet):
    if type(packet) is not Budget:
        raise ValueError("Budget required")
    integer(packet.degree, 1)
    integer(packet.product_valuation, 1)
    integer(packet.trace_power, 1)
    if type(packet.branch_counts) is not tuple or any(type(x) is not int for x in packet.branch_counts):
        raise ValueError("exact branch counts required")
    if type(packet.root_valuation) is not Fraction or packet != budget(packet.p, packet.r):
        raise ValueError("budget is not bound to its prime and depth")
    return packet


def primitive_matrices(p, r):
    modulus = parameters(p, r)
    return tuple((p ** (r-j), b, p ** j)
                 for j in range(r+1) for b in range(p ** j)
                 if gcd(gcd(p ** (r-j), b), p ** j) == 1)


def all_matrices(p, r):
    parameters(p, r)
    return tuple((p ** (r-j), b, p ** j)
                 for j in range(r+1) for b in range(p ** j))


def pair(z):
    if type(z) is not tuple or len(z) != 2:
        raise ValueError("an exact coefficient pair is required")
    for x in z:
        integer(x)
    return z


def quadratic(p, r, d):
    modulus = parameters(p, r)
    integer(d)
    if pow(d % p, (p-1)//2, p) != p-1:
        raise ValueError("unramified nonsplit quadratic algebra required")
    return modulus


def multiply(z, w, d, modulus):
    a, b = pair(z)
    c, e = pair(w)
    return ((a*c + d*b*e) % modulus, (a*e+b*c) % modulus)


def projective(z, p, r):
    modulus = parameters(p, r)
    a, b = pair(z)
    a, b = a % modulus, b % modulus
    if a % p:
        return (1, b*pow(a, -1, modulus) % modulus)
    if b % p:
        return (a*pow(b, -1, modulus) % modulus, 1)
    raise ValueError("primitive vector required")


def projective_points(p, r):
    modulus = parameters(p, r)
    return tuple((1, b) for b in range(modulus)) + tuple((a, 1) for a in range(0, modulus, p))


def circle(z, p, r, d):
    """[z] -> z / conjugate(z), on the unramified quadratic algebra."""
    modulus = quadratic(p, r, d)
    a, b = projective(z, p, r)
    inverse_norm = pow((a*a-d*b*b) % modulus, -1, modulus)
    square = multiply((a,b), (a,b), d, modulus)
    return tuple(x*inverse_norm % modulus for x in square)


def circle_inverse(w, p, r, d):
    """Two-chart Hilbert 90, retaining the missing w=-1 chart."""
    modulus = quadratic(p, r, d)
    a, b = pair(w)
    if not (0 <= a < modulus and 0 <= b < modulus) or (a*a-d*b*b) % modulus != 1:
        raise ValueError("normalized norm-one point required")
    z = ((1+a) % modulus, b)
    if all(x % p == 0 for x in z):
        z = multiply((0,1), ((1-a) % modulus, -b % modulus), d, modulus)
    return projective(z, p, r)


def series(values, size):
    return list(values[:size]) + [0] * max(0, size-len(values))


def mul(a, b, size):
    out = [0]*size
    for i, x in enumerate(a[:size]):
        if x:
            for j, y in enumerate(b[:size-i]):
                out[i+j] += x*y
    return out


def inv(a, size):
    if not a or a[0] != 1:
        raise ValueError("normalized unit series required")
    out = [1]+[0]*(size-1)
    for n in range(1, size):
        out[n] = -sum(a[j]*out[n-j] for j in range(1, min(n+1, len(a))))
    return out


def power(a, exponent, size):
    integer(exponent, 0)
    out = [1]+[0]*(size-1)
    while exponent:
        if exponent & 1:
            out = mul(out, a, size)
        exponent //= 2
        if exponent:
            a = mul(a, a, size)
    return out


def dilate(a, exponent, size):
    return [a[n//exponent] if n % exponent == 0 and n//exponent < len(a) else 0 for n in range(size)]


def logarithmic_derivative(a):
    size = len(a)
    return mul([n*x for n,x in enumerate(a)], inv(a,size), size)


def cyclic_norm(a, degree, size):
    """Product F(zeta_m^b q^(1/m)); exact recurrence, no complex roots."""
    integer(degree, 1)
    integer(size, 1)
    if len(a) <= degree*(size-1):
        raise ValueError("insufficient original-series precision")
    if any(type(x) is not int for x in a):
        raise ValueError("integer series required")
    derivative = logarithmic_derivative(a)
    out = [1]+[0]*(size-1)
    for n in range(1,size):
        numerator = sum(derivative[degree*j]*out[n-j] for j in range(1,n+1))
        out[n], remainder = divmod(numerator,n)
        if remainder:
            raise ArithmeticError("cyclic norm lost integrality")
    return out


def primitive_norm(a, p, r, size):
    """The exact primitive product divided by F(q)^psi(p^r)."""
    parameters(p,r)
    norms = {j: cyclic_norm(a,p**j,size) for j in range(r+1)}
    product = mul(dilate(a,p**r,size), norms[r], size)
    for j in range(1,r):
        product = mul(product,dilate(norms[j],p**(r-j),size),size)
        product = mul(product,inv(dilate(norms[j-1],p**(r-j-1),size),size),size)
    return mul(product,inv(power(series(a,size),budget(p,r).degree,size),size),size)


def theta7(size):
    out = [0]*size
    bound = isqrt(2*size)+2
    for x in range(-bound,bound+1):
        for y in range(-bound,bound+1):
            n = x*x+x*y+2*y*y
            if n < size:
                out[n] += 1
    return out


def parameter7(theta, size):
    eta = [1]+[0]*(size-1)
    for step in (1,7):
        for n in range(step,size,step):
            for _ in range(3):
                for k in range(size-1,n-1,-1):
                    eta[k] -= eta[k-n]
    ratio = mul(eta,inv(power(theta,3,size),size),size)
    return [0]+ratio[:size-1]


def polynomial_in_t(a, t, degree):
    size = len(a)
    remainder = a[:]
    term = [1]+[0]*(size-1)
    coefficients = []
    for k in range(degree+1):
        coefficients.append(remainder[k])
        remainder = [x-coefficients[-1]*y for x,y in zip(remainder,term)]
        term = mul(term,t,size)
    if any(remainder):
        raise ArithmeticError("supplied degree bound failed at retained precision")
    return tuple(coefficients)


def evaluate(coefficients, x):
    answer = Fraction(0)
    for c in reversed(coefficients):
        answer = answer*x+c
    return answer


def s7_trace(p, r):
    """Exact q-series reconstruction; arithmetic conclusion uses note hypotheses."""
    parameters(p,r)
    if p == 7:
        raise ValueError("bad level prime")
    data = budget(p,r)
    top_degree = (2*data.degree)//3
    size = top_degree+5
    needed = p**r*(size-1)+1
    theta = theta7(needed)
    form = mul(theta,theta,needed)
    reciprocal = inv(form,needed)
    top = primitive_norm(form,p,r,size)
    trace_parts = inv(dilate(form,p**r,size),size)
    for j in range(1,r+1):
        stride = p**j
        selected = [reciprocal[stride*n] for n in range(size)]
        block = dilate(selected,p**(r-j),size)
        for n in range(size):
            trace_parts[n] += p**(3*j)*block[n]
        if j < r:
            selected = [reciprocal[p**(j-1)*n] for n in range(size)]
            block = dilate(selected,p**(r-j-1),size)
            for n in range(size):
                trace_parts[n] -= p**(3*j-1)*block[n]
    trace = mul(series(form,size),trace_parts,size)
    numerator = mul(top,trace,size)
    t = parameter7(theta,size)
    top_poly = polynomial_in_t(top,t,top_degree)
    trace_poly = polynomial_in_t(numerator,t,(2*(data.degree-1))//3)
    return top_poly,trace_poly


def verify():
    checks = 0
    def check(condition):
        nonlocal checks
        checks += 1
        if not condition:
            raise AssertionError(checks)
    def rejects(call):
        try:
            call()
        except (ValueError,ArithmeticError):
            check(True)
        else:
            check(False)

    parameters_checked = []
    for p in (3,5,7,11):
        d = next(d for d in range(2,p) if pow(d,(p-1)//2,p)==p-1)
        for r in range(1,4):
            b = audit_budget(budget(p,r))
            check(sum(b.branch_counts)==b.degree)
            check(b.root_valuation == 2*r-Fraction(2*(p**r-1),(p-1)*b.degree))
            check(2*r-1 < b.root_valuation < 2*r)
            check(b.trace_power == 2*r)
            matrices = primitive_matrices(p,r)
            check(len(matrices)==b.degree)
            check(sum(2*val(d0,p) for _,_,d0 in matrices)==b.product_valuation)
            full = all_matrices(p,r)
            content_hist = {}
            for a,b0,d0 in full:
                content = val(gcd(gcd(a,b0),d0),p)
                content_hist[content] = content_hist.get(content,0)+1
            for j,count in content_hist.items():
                k = r-2*j
                check(count==(1 if k==0 else budget(p,k).degree))
            points = projective_points(p,r)
            images = tuple(circle(z,p,r,d) for z in points)
            check(len(set(images)) == len(points) == b.degree)
            for z,w in zip(points,images):
                check(circle_inverse(w,p,r,d)==z)
                check((w[0]*w[0]-d*w[1]*w[1])%(p**r)==1)
            # Actual nonsplit Cartan action from any fixed source is bijective.
            for source in points[:3]:
                orbit = {projective(multiply(source,z,d,p**r),p,r) for z in points}
                check(orbit==set(points))
                for z in points:
                    check(circle(multiply(source,z,d,p**r),p,r,d)==
                          multiply(circle(source,p,r,d),circle(z,p,r,d),d,p**r))
            check(circle_inverse((p**r-1,0),p,r,d)==(0,1))
            if r>1:
                fibers = {}
                for z in points:
                    reduced = projective(z,p,r-1)
                    fibers[reduced]=fibers.get(reduced,0)+1
                check(set(fibers.values())=={p})
            parameters_checked.append([p,r,b.degree,str(b.root_valuation)])

    # Independent resultant route for integral cyclic norms.
    import sympy as sp
    x,q=sp.symbols('x q')
    for seed in ((1,-1),(1,2,-1,3),(1,0,2,1)):
        f=sum(c*x**i for i,c in enumerate(seed))
        for m in (1,3,5,9):
            resultant=sp.Poly(sp.resultant(x**m-q,f,x),q)
            # Normalize the resultant's convention by the known constant 1.
            expected=sp.Poly(resultant.as_expr()/resultant.nth(0),q)
            size=7
            actual=cyclic_norm(series(seed,m*(size-1)+1),m,size)
            check(actual==[int(expected.nth(i)) for i in range(size)])
            if m>1:
                p=3 if m in (3,9) else 5
                check(all((actual[i]-(seed[i] if i<len(seed) else 0))%p==0 for i in range(size)))
        for p,r in ((3,1),(3,2),(5,2)):
            size=7
            actual=primitive_norm(series(seed,p**r*(size-1)+1),p,r,size)
            check(actual[0]==1 and all(c%p==0 for c in actual[1:]))

    # Composed one-step paths retain p+1 scalar backtracking copies.
    for p in (3,5,7):
        size=8
        h=[(n*n+3*n+1)*(-1)**n for n in range(p*p*(size-1)+1)]
        def one_step(values, retained):
            return [((values[n//p] if n%p==0 else 0)+p**3*values[p*n])
                    for n in range(retained)]
        composed=one_step(one_step(h,p*(size-1)+1),size)
        cyclic=[(h[n//(p*p)] if n%(p*p)==0 else 0)
                +p**3*(h[n] if n%p==0 else 0)-p**2*h[n]
                +p**6*h[p*p*n] for n in range(size)]
        check(composed==[cyclic[n]+(p+1)*p*p*h[n] for n in range(size)])
        check((p+1)**2==budget(p,2).degree+(p+1))

    traces=[]
    for p,r in ((3,1),(3,2),(3,3),(5,1),(5,2),(11,1)):
        top,num=s7_trace(p,r)
        check(top[0]==1)
        check(all(c%p==0 for c in top[1:]))
        cm_rows=[]
        for t0,disc in ((Fraction(-1),-7),(Fraction(1,2),-12),(Fraction(1,8),-19)):
            if t0.denominator%p==0 or disc%p==0 or pow(disc%p,(p-1)//2,p)!=p-1:
                continue
            den=evaluate(top,t0)
            value=evaluate(num,t0)/den
            check(rational_val(den,p)==0)
            check(value==0 or rational_val(value,p)>=2*r)
            cm_rows.append({"t":str(t0),"discriminant":disc,"trace":str(value),"valuation":rational_val(value,p)})
        traces.append({"p":p,"r":r,"top":top,"trace_numerator":num,"CM":cm_rows})
        if p==3 and r==2:
            primitive_trace=evaluate(num,Fraction(-1))/evaluate(top,Fraction(-1))
            full_trace=primitive_trace+p*p
            check(rational_val(primitive_trace,p)>=4)
            check(rational_val(full_trace,p)==2)
            imprimitive={"p":p,"r":r,"primitive_trace":str(primitive_trace),"full_trace":str(full_trace),"full_valuation":2}

    # A generic integral parameter does not inherit inert-CM valuation.
    top,num=s7_trace(3,2)
    check(rational_val(evaluate(num,Fraction(0))/evaluate(top,Fraction(0)),3)==0)
    # The same combinatorial valuation bill at p=2,r=2 is 3, not 4.
    check(Fraction(2*(1+2*4),1+1+4)==3)
    check(Fraction(2*(1+2*4),1+1+4)<4)

    malformed=(lambda:budget(True,1),lambda:budget(9,1),lambda:budget(3,False),
               lambda:audit_budget(replace(budget(3,2),r=2.0)),
               lambda:audit_budget(replace(budget(3,2),trace_power=True)),
               lambda:projective((3,3),3,2),lambda:projective((True,0),3,1),
               lambda:circle((1,0),3,1,1),lambda:circle((1,0),3,1,3),
               lambda:circle_inverse((0,0),3,1,2),lambda:circle_inverse((-1,0),3,1,2),
               lambda:cyclic_norm([1,2],3,3),lambda:s7_trace(7,1))
    for call in malformed:
        rejects(call)
    return {"checks":checks,"scope":"primitive support and q-series finite checks; no truncated-series inference",
            "Cartan_parameters":parameters_checked,"s7_traces":traces,"imprimitive_hostile":imprimitive,
            "type_and_guard_hostiles":len(malformed)}


if __name__=='__main__':
    print(json.dumps(verify(),indent=2,sort_keys=True))
