"""Exact consecutive power rows, successor recurrences, and bounded root receipts.

Run python -B [-O] this_file.py. No numerical root finder or global higher-power
nonexistence claim is used. Output is deterministic JSON.
"""
from dataclasses import dataclass, replace
from fractions import Fraction
from math import comb, gcd, isqrt
import json


def integer(value, minimum=0):
    if type(value) is not int or value < minimum:
        raise ValueError("exact integer outside domain")
    return value


def power_sums(n, degree):
    """S_j(n)=sum(i**j, i=1..n), via exact telescoping; O(degree**2) operations."""
    integer(n)
    integer(degree)
    values = [n]
    for j in range(1, degree + 1):
        top = (n + 1)**(j + 1) - 1 - sum(comb(j + 1, r)*values[r] for r in range(j))
        value, rem = divmod(top, j + 1)
        if rem:
            raise ArithmeticError("power-sum divisibility")
        values.append(value)
    return tuple(values)


def defect(centre, n, power):
    integer(centre)
    integer(n, 1)
    integer(power, 1)
    moments = power_sums(n, power)
    return centre**power - 2*sum(comb(power, j)*centre**(power-j)*moments[j]
                                  for j in range(1, power + 1, 2))


def direct_defect(centre, n, power):
    return sum((centre-i)**power for i in range(n+1)) - sum((centre+i)**power for i in range(1,n+1))


@dataclass(frozen=True)
class RootDecision:
    n: int
    power: int
    lower: int
    upper: int
    lower_defect: int
    upper_defect: int


def locate(n, power):
    """Exact unique positive centre: equal endpoints or a strict adjacent bracket.

    For power>=3 the initial integer interval has width power, independently n.
    A zero endpoint is an exact integer (and hence rational) balanced centre.
    """
    integer(n, 1)
    integer(power, 1)
    u = n*(n+1)
    if power <= 2:
        c = power*u
        return RootDecision(n, power, c, c, 0, 0)
    moments = power_sums(n, power)
    def evaluate(c):
        return c**power - 2*sum(comb(power,j)*c**(power-j)*moments[j]
                                for j in range(1,power+1,2))
    lo, hi = power*u, power*(u+1)
    dlo, dhi = evaluate(lo), evaluate(hi)
    if not dlo < 0 < dhi:
        raise ArithmeticError("uniform bracket theorem")
    while hi-lo > 1:
        mid = (lo+hi)//2
        value = evaluate(mid)
        if value == 0:
            return RootDecision(n,power,mid,mid,0,0)
        if value < 0:
            lo,dlo = mid,value
        else:
            hi,dhi = mid,value
    return RootDecision(n,power,lo,hi,dlo,dhi)


def audit_decision(receipt):
    if type(receipt) is not RootDecision:
        raise ValueError("typed root decision required")
    for value in (receipt.n,receipt.power,receipt.lower,receipt.upper):
        integer(value,1)
    for value in (receipt.lower_defect,receipt.upper_defect):
        if type(value) is not int:
            raise ValueError("exact signed defects required")
    if receipt != locate(receipt.n,receipt.power):
        raise ValueError("source-bound root decision does not replay")
    return receipt


@dataclass(frozen=True)
class Row:
    n: int
    power: int
    start: int
    centre: int
    total: int


def row(n, power):
    decision = locate(n,power)
    if decision.lower != decision.upper:
        return None
    c = decision.lower
    start = c-n
    if start < 1:
        raise ArithmeticError("positive-window theorem")
    total = power_sums(c,power)[power] - power_sums(start-1,power)[power]
    return Row(n,power,start,c,total)


def audit_row(receipt):
    if type(receipt) is not Row:
        raise ValueError("typed balanced row required")
    for value in (receipt.n,receipt.power,receipt.start,receipt.centre,receipt.total):
        integer(value,1)
    expected = row(receipt.n,receipt.power)
    if expected is None or expected != receipt:
        raise ValueError("balanced row does not replay")
    return receipt


def successor(receipt):
    r = audit_row(receipt)
    if r.power not in (1,2):
        raise ValueError("successor rule proved only for powers one and two")
    out = row(r.n+1,r.power)
    if out.start != r.start + 2*r.power*r.n + 2*r.power-1:
        raise ArithmeticError("successor address")
    return out


def split_start(left, right, power):
    """Complete positive-integer start test for arbitrary adjacent counts, p=1,2."""
    integer(left,1)
    integer(right,1)
    integer(power,1)
    if power not in (1,2):
        raise ValueError("arbitrary-count compiler is scoped to powers one and two")
    if left <= right:
        return None
    d,n = left-right,right
    linear = d*d-d-2*n*n
    if power == 1:
        numerator,denominator = -linear,2*d
    else:
        constant = 2*power_sums(left-1,2)[2]-power_sums(left+right-1,2)[2]
        discriminant = linear*linear-4*d*constant
        if discriminant < 0:
            return None
        root = isqrt(discriminant)
        if root*root != discriminant:
            return None
        numerator,denominator = -linear+root,2*d
    if numerator < denominator or numerator % denominator:
        return None
    return numerator//denominator


def square_reflection(n):
    """Rational Householder map taking the square row to its zero-padded right side."""
    integer(n,1)
    v = (-(n+1),)*n + (2*n*(n+1),)
    norm = sum(x*x for x in v)
    return tuple(tuple(Fraction(int(i==j))-Fraction(2*v[i]*v[j],norm)
                       for j in range(n+1)) for i in range(n+1))


def matvec(matrix, vector):
    return tuple(sum((a*b for a,b in zip(line,vector)),Fraction(0)) for line in matrix)


def cube_pairs(value):
    """Inherited two-cube divisor criterion, restricted to positive unordered pairs."""
    integer(value,1)
    divisors = set()
    for d in range(1,isqrt(value)+1):
        if value % d == 0:
            divisors.update((d,value//d))
    pairs=[]
    for d in sorted(divisors):
        top = d**3-value
        if top % (3*d):
            continue
        xy=top//(3*d)
        discriminant=d*d-4*xy
        if discriminant < 0:
            continue
        e=isqrt(discriminant)
        if e*e==discriminant and (d-e)%2==0 and d>e:
            pairs.append(((d-e)//2,(d+e)//2))
    return tuple(sorted(pairs))


def differences(values):
    return tuple(b-a for a,b in zip(values,values[1:]))


CHECKS=0


def check(condition,message):
    global CHECKS
    CHECKS+=1
    if not condition:
        raise RuntimeError(message)


def rejects(call):
    try:
        call()
    except ValueError:
        check(True,"malformed receipt rejected")
        return
    check(False,"malformed receipt accepted")


def main():
    # Independent direct summation of moments and central defects.
    for n in range(0,33):
        moments=power_sums(n,16)
        for p,value in enumerate(moments):
            check(value==sum(i**p for i in range(1,n+1)),"telescoped moments")
    for n in range(1,17):
        for p in range(1,11):
            for c in (n,n+1,p*n*(n+1),p*n*(n+1)+1):
                check(defect(c,n,p)==direct_defect(c,n,p),"odd-moment identity")
    # Complete source-owned root decisions in a declared finite rectangle.
    rational_rows=[]
    first_failed_unit=[]
    for p in range(1,33):
        bad=None
        for n in range(1,129):
            r=locate(n,p)
            check(audit_decision(r)==r,"root decision replay")
            u=n*(n+1)
            if r.lower==r.upper:
                check(defect(r.lower,n,p)==0 and r.lower>n,"literal balanced row")
                check(r.lower%2==0,"integer-root parity")
                rational_rows.append((n,p,r.lower))
            else:
                check(r.upper==r.lower+1 and r.lower_defect<0<r.upper_defect,"strict unique-root bracket")
                check(p*u<=r.lower<r.upper<=p*(u+1),"uniform exponent-width bound")
            if p>=3 and defect(p*u+1,n,p)<=0 and bad is None:
                bad=n
        if bad is not None:
            first_failed_unit.append((p,bad))
    check(all(p<=2 for _,p,_ in rational_rows),"declared finite higher-power census")
    check(first_failed_unit[0]==(15,3),"first failure in declared rectangle")
    bad15=defect(181,3,15)
    bad16=defect(33,1,16)
    check(bad15==-61355600116348100287252240367 and bad16==-2148849686515840130559,"historical bracket hostiles")
    check(locate(3,15).lower==181 and locate(1,16).lower==33,"roots lie beyond former unit bracket")
    # Arbitrarily large input height is addressed without looping through the row.
    large=locate(10**100,17)
    check(large.lower_defect<0<large.upper_defect and large.upper-large.lower==1,"large-height exact arithmetic")
    # Cubic and quartic no-rational-row proofs have explicit sign polynomials.
    for n in range(1,129):
        u=n*(n+1)
        check(defect(3*u,n,3)==-u*u//2,"cubic lower")
        check(2*defect(3*u+1,n,3)==17*u*u+12*u+2,"cubic upper")
        check(defect(4*u,n,4)==-8*u**3,"quartic lower")
        check(defect(4*u+1,n,4)==(4*u+1)*(14*u*u+8*u+1),"quartic upper")
    # Actual successor rows and the complete quadratic finite-difference bill.
    samples=[]
    for p in (1,2):
        r=row(1,p)
        for n in range(1,65):
            check(r==row(n,p) and audit_row(r)==r,"successor source")
            check(r.start==p*n*n+(p-1)*n,"start polynomial")
            check(sum(i**p for i in range(r.start,r.centre+1))==r.total,"left total")
            check(sum(i**p for i in range(r.centre+1,r.centre+n+1))==r.total,"right total")
            if n<=5:samples.append(r.__dict__)
            r=successor(r)
        starts=tuple(row(n,p).start for n in range(1,33))
        check(set(differences(differences(starts)))=={2*p},"constant second difference")
        check(set(differences(differences(differences(starts))))=={0},"third difference vanishes")
    sums1=tuple(n*(n+1)*(2*n+1)//2 for n in range(13))
    sums2=tuple(n*(n+1)*(2*n+1)*(12*n*n+12*n+1)//6 for n in range(13))
    check(set(differences(differences(differences(sums1))))=={6},"linear row sum cubic")
    fifth=sums2
    for _ in range(5):fifth=differences(fifth)
    check(set(fifth)=={480},"square row sum quintic")
    check(set(differences(fifth))=={0},"sixth sum difference vanishes")
    for n in range(1,65):
        check(row(n,1).total==3*power_sums(n,2)[2],"square-pyramid shadow")
        check(row(n,2).total==n*(n+1)*(2*n+1)*(12*n*n+12*n+1)//6,"square common total")
    # Arbitrary split counts: compare every tested start with the iff compiler.
    split_examples=[]
    for left in range(1,11):
        for right in range(1,11):
            for p in (1,2):
                found=split_start(left,right,p)
                if found is not None and left-right!=1:
                    split_examples.append((left,right,p,found))
                for start in range(1,101):
                    equal=sum(i**p for i in range(start,start+left))==sum(i**p for i in range(start+left,start+left+right))
                    check(equal==(start==found),"arbitrary-count exact membership")
    for left,right,start in ((17,8,18),(35,10,4)):
        check(split_start(left,right,2)==start,"unequal-width square split")
        check(sum(i*i for i in range(start,start+left))==sum(i*i for i in range(start+left,start+left+right)),"literal unequal-width square identity")
        split_examples.append((left,right,2,start))
    # The inherited raw-block coincidence is unique among all four side types.
    coincidences=[]
    for n in range(1,65):
        q=row(n,2)
        qb=((q.start,q.centre),(q.centre+1,q.centre+n))
        for m in range(1,65):
            f=row(m,1)
            fb=((f.start,f.centre),(f.centre+1,f.centre+m))
            for i in range(2):
                for j in range(2):
                    if qb[i]==fb[j]:coincidences.append((n,i,m,j,qb[i]))
    check(coincidences==[(3,0,4,1,(21,24))],"unique raw shared block")
    # Rational quadratic isometry; it does not preserve the integer positive cone.
    for n in range(1,9):
        r=row(n,2)
        left=tuple(range(r.start,r.centre+1))
        right=tuple(range(r.centre+1,r.centre+n+1))+(0,)
        H=square_reflection(n)
        check(matvec(H,left)==right and matvec(H,right)==left,"marked norm reflection")
        for j in range(n+1):
            basis=tuple(int(i==j) for i in range(n+1))
            check(matvec(H,matvec(H,basis))==basis,"reflection involution")
    check(matvec(square_reflection(1),(1,0))==(Fraction(3,5),Fraction(4,5)),"integer-lattice hostile")
    check(matvec(square_reflection(1),(0,1))==(Fraction(4,5),Fraction(-3,5)),"positive-cone hostile")
    check(gcd(3,gcd(4,5))==1 and 3*3+4*4==5*5,"primitive consecutive root")
    check(5*5+12*12==13*13 and 12-5!=1,"not every primitive triangle is a tower root")
    # Cube norm and Gaussian norm keep different arithmetic predicates.
    check(cube_pairs(1729)==((1,12),(9,10)),"taxicab divisor reconstruction")
    check(all(a**3+b**3==1729 for a,b in cube_pairs(1729)),"literal cube equalities")
    squares7={i*i%7 for i in range(7)}
    check([(a,b) for a in squares7 for b in squares7 if (a+b)%7==0]==[(0,0)],"minus-one inert at seven")
    check(1729%7==0 and 1729%49!=0,"1729 cannot be two integer squares")
    check(not any(a*a+b*b==1729 for a in range(42) for b in range(42)),"independent square census")
    # Public receipt types and exact source binding are checked before equality.
    r=row(2,2);d=locate(3,15)
    hostiles=(lambda:locate(True,2),lambda:locate(2,2.0),lambda:locate(0,2),
              lambda:power_sums(2,-1),lambda:defect(3,1,True),
              lambda:audit_row(replace(r,n=2.0)),lambda:audit_row(replace(r,start=r.start+1)),
              lambda:audit_row(replace(r,power=3)),lambda:audit_decision(replace(d,lower=float(d.lower))),
              lambda:audit_decision(replace(d,upper_defect=True)),lambda:audit_decision(replace(d,n=4)),
              lambda:split_start(True,1,1),lambda:split_start(3,2,3),lambda:square_reflection(0),
              lambda:cube_pairs(Fraction(1729)))
    for hostile in hostiles:rejects(hostile)
    report={"status":"PROVED scoped exact classification/recursions; higher-power finite census only",
            "checks":CHECKS,"root_rectangle":{"n":[1,128],"power":[1,32],"integer_rows":len(rational_rows),
                "higher_power_integer_rows":0,"first_unit_bracket_failures":first_failed_unit},
            "uniform_bracket":"p*n*(n+1) < c < p*(n*(n+1)+1), p>=3",
            "historical_hostiles":{"p15_n3_c181":str(bad15),"p16_n1_c33":str(bad16)},
            "sample_rows":samples,"arbitrary_split_examples":split_examples,
            "arbitrary_split_universe":{"left":[1,10],"right":[1,10],"start":[1,100],"power":[1,2]},
            "raw_block_coincidences":coincidences,"taxicab_pairs":cube_pairs(1729),
            "malformed_controls":len(hostiles),"scope":"No all-p>=5 nonexistence, Collatz, JC, or plane-colouring conclusion."}
    print(json.dumps(report,indent=2,sort_keys=True))


if __name__=="__main__":
    main()
