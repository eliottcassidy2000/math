"""Exact source-weighted reset controls. Standard library; no removable asserts."""
from fractions import Fraction as F
from collections import defaultdict
import json

SCALE = 1 << 100


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def valuation(n, p):
    k = 0
    while n % p == 0:
        n //= p
        k += 1
    return k


def weight(j):
    return F(1, j*(j+1))


def reset(j):
    require(j >= 2, 'active integer source')
    k = (j & -j).bit_length()-1
    ell = j >> k
    y = 3**k * ell
    a = 1 + valuation(3*y-1, 2)
    target = ((3*y-1) // (1 << (a-1)) + 1) // 2
    return k, y, a, target


def point_mass(y, z):
    if y < 3 or y % 2 == 0:
        return F(0)
    return sum(z**k * weight((1 << k)*(y//3**k)) for k in range(valuation(y,3)+1))


def literal_controls():
    for j in range(2, 20002):
        k,y,a,target = reset(j)
        n = 2*j-1
        for i in range(k):
            require(valuation(3*n+1,2) == 1, 'exact ones prefix')
            n = (3*n+1)//2
            require((n+1)//2 == (3**(i+1))*(j//(1 << (i+1))), 'integer clock')
        require((n+1)//2 == y, 'marked reset source')
        require(valuation(3*n+1,2) == a and a >= 2, 'first reset valuation')
        require(((3*n+1)//(1 << a)+1)//2 == target, 'literal reset target')
        require(target % 3 != 2, 'forbidden endpoint phase')
        require(target <= F(3*y+1,4), 'reset size bound')
    checks = 0
    for z in (F(1), F(3,4)):
        literal = defaultdict(F)
        for j in range(2,502):
            k,y,_,_ = reset(j)
            literal[y] += z**k * weight(j)
        for y in range(1,502):
            require(literal[y] == point_mass(y,z), 'exact inverse-fibre law')
            checks += 1
    # The explicit root family, independently tested by the actual odd map.
    roots = []
    for k in range(5):
        for t in range(4):
            exponent = (2*t+1)*3**k
            j = (1 << k)*((1 << exponent)+1)//3**(k+1)
            if j == 1:
                continue
            kk,_,a,target = reset(j)
            require(kk == k and a == exponent+1 and target == 1, 'root-family arrival')
            if k<=2 and t<=1:
                roots.append({'k':k,'t':t,'j':j,'odd_source':2*j-1,'last_valuation':a})
    # Exhaustive finite iff, not only forward examples.
    for j in range(2,20002):
        k,_,a,target = reset(j)
        legal = (a-1) % (2*3**k) == 3**k
        expected = (1 << k)*((1 << (a-1))+1)//3**(k+1) if legal else -1
        require((target == 1) == (legal and j == expected), 'root-family converse')
    return {'sources':20000,'point_fibres':checks,'root_examples':roots,
            'parent_descent_hostile':{'odd_source':27,'odd_endpoint':2*reset(14)[3]-1}}


def tail_controls(limit=1000000):
    thresholds = (10,100,1000,10000)
    zs = (F(1),F(3,4))
    sums = {z:{J:[0,0,0,0,0] for J in thresholds} for z in zs}
    powers = {z:[z**k for k in range(limit.bit_length()+2)] for z in zs}
    for j in range(2,limit+1):
        k,y,a,target = reset(j)
        den = j*(j+1)
        for z in zs:
            v = powers[z][k]
            prefix = SCALE*v.numerator//(v.denominator*den)
            v *= z
            full = SCALE*v.numerator//(v.denominator*den)
            for J in thresholds:
                if y > J:
                    sums[z][J][0] += prefix
                if target > J:
                    sums[z][J][1] += full
                    sums[z][J][2+target%3] += full
    rows = []
    for z in zs:
        for J in thresholds:
            bounds=[]
            for index,amount in enumerate(sums[z][J]):
                low=F(amount,SCALE)
                high=low+F(limit,SCALE)+F(1,limit+1)
                if index==4:
                    require(amount==0,'exact forbidden phase')
                    low=high=F(0)
                bounds.append([str(low),str(high)])
            require(F(sums[z][J][0],SCALE) <= F(2,J)/(1-3*z/4), 'independent prefix envelope')
            require(F(sums[z][J][1],SCALE) <= 6*z/((1-3*z/4)*(4*J-1)), 'independent full-block envelope')
            require(sums[z][J][4] == 0, 'no discarded forbidden-phase mass')
            rows.append({'z':str(z),'cutoff_J':J,'exact_tail_intervals':bounds,
                         'J_times_midpoints_diagnostic':[float(J*(F(a)+F(b))/2) for a,b in bounds],
                         'predicted_leading_coefficients':[str(1/(2*(1-3*z/4))),
                             str(z/(4-3*z)),str(z/(5*(4-3*z))),str(4*z/(5*(4-3*z))),'0']})
    return {'source_cutoff':limit,'fixed_point_scale_bits':100,
            'rounding_error_per_sum':str(F(limit,SCALE)),
            'omitted_original_source_mass':str(F(1,limit+1)),
            'interval_order':['prefix','full','full phase0','full phase1','full phase2'],
            'rows':rows}


def root_mass(K=4,T=12):
    mass=0
    count=0
    first_omitted=(1 << (K+1))*((1 << (3**(K+1)))+1)//3**(K+2)
    # Both increasing t and increasing k at t=0 at least double j.
    tail=F(9,4)*weight(first_omitted)
    for k in range(K+1):
        for t in range(T+1):
            exponent=(2*t+1)*3**k
            j=(1 << k)*((1 << exponent)+1)//3**(k+1)
            if j==1:
                continue
            mass += SCALE//(j*(j+1))
            count += 1
        exponent=(2*T+3)*3**k
        j=(1 << k)*((1 << exponent)+1)//3**(k+1)
        # Consecutive t gives at least doubling of j and weight ratio <=1/3.
        tail += F(3,2)*weight(j)
    lo=F(mass,SCALE)
    hi=lo+F(count,SCALE)+tail
    # Compact outward rational bound, avoiding a huge-denominator printout.
    q=10**12
    a=lo.numerator*q//lo.denominator
    b=-((-hi.numerator*q)//hi.denominator)
    require(F(a,q)<=lo<=hi<=F(b,q), 'outward root-mass enclosure')
    return {'k_max':K,'t_max':T,'included_sources':count,
            'first_block_terminal_mass_interval':[str(F(a,q)),str(F(b,q))],
            'already_terminal_plus_first_block_interval':[str(F(1,2)+F(a,q)),str(F(1,2)+F(b,q))]}


if __name__=='__main__':
    print(json.dumps({'literal':literal_controls(),'tails':tail_controls(),
                      'root_mass':root_mass()},indent=2,sort_keys=True))
