#!/usr/bin/env python3
"""Independent exact checks of Session31 structural repairs; no assertions.

The proofs are in AUDIT.md. This finite checker uses literal dynamics, a binary
valuation-word constructor, piecewise polynomial integration, and rational
exponential certificates. It never imports either audited checker.
"""
from fractions import Fraction as Q
from math import comb, factorial
from pathlib import Path
import hashlib
import json

COUNT = {}
def need(condition, label):
    COUNT[label] = COUNT.get(label, 0) + 1
    if not condition:
        raise RuntimeError(label)

def valuation(n):
    if n <= 0:
        raise ValueError('valuation requires n > 0')
    k = 0
    while not n % 2:
        n //= 2
        k += 1
    return k

def odd_step(x):
    if x <= 0 or x % 2 != 1:
        raise ValueError('positive odd start required')
    raw = 3*x+1
    k = valuation(raw)
    return raw//2**k, k

def run(x):
    start = x
    word = []
    trace = [x]
    while True:
        z = 3*x+1
        trace.append(z)
        k = 0
        while z % 2 == 0:
            z //= 2
            k += 1
            trace.append(z)
        word.append(k)
        if k > 1:
            return {'start':start, 'end':z, 'word':word, 'trace':trace}
        x = z

def realize_word(word):
    # 2^E*x_L = 3^L*x_0 + B. Prescribe x_L odd modulo 2.
    E = 0
    B = 0
    for k in word:
        B = 3*B + 2**E
        E += k
    modulus = 2**(E+1)
    residue = ((2**E-B)*pow(3**len(word), -1, modulus)) % modulus
    return residue if residue else modulus

def check_collatz():
    examples = []
    for r in range(1, 17):
        for s in range(1, 17):
            for t in range(1, 17):
                target = [1]*(r-1)+[s+1]+[1]*(t-1)+[2]
                x = realize_word(target)
                z = x
                actual = []
                for _ in target:
                    z,k = odd_step(z)
                    actual.append(k)
                need(actual == target, 'binary_word_realization')
                # Separately inspect the manuscript audit's mod-3 construction.
                mod = 3**r
                b = ((2**s-1)*pow(2**(s+t), -1, mod)) % mod
                if not b % 2:
                    b += mod
                a = (2**(s+t)*b-2**s+1)//mod
                need(a > 0 and a % 2 == 1 and b > 0 and b % 2 == 1,
                     'mod3_construction_positive_odd')
                z = a*2**r-1
                R = run(z)
                next_R = run(R['end'])
                need(R['word'] == [1]*(r-1)+[s+1] and len(next_R['word']) == t,
                     'mod3_construction_literal_runs')
                need(R['end'] == 2**t*b-1, 'mod3_construction_endpoint')
    for r,s,t in [(20,20,20),(1,80,1),(80,1,80),(73,61,47),(32,64,96)]:
        target = [1]*(r-1)+[s+1]+[1]*(t-1)+[2]
        x = realize_word(target)
        z = x
        for k in target:
            z, actual = odd_step(z)
            need(k == actual, 'asymmetric_large_word')
        examples.append({'r':r,'s':s,'t':t,'start':str(x)})
    for m in range(1, 81):
        for q in (1,3,17):
            x = q*2**(3*m+1)-5
            initial = x
            for j in range(m):
                R = run(x)
                expected = q*9**(j+1)*2**(3*(m-j-1)+1)-5
                need(R['word'] == [1,2] and R['end'] == expected,
                     'fixed_growth_cycle')
                need(min(R['trace']) == x and R['end'] > x >= initial,
                     'no_original_baseline_crossing')
                x = R['end']
    for r in range(3,64,2):
        for k in range(2,65,2):
            y = 2**r-1
            x = (2**k*y-1)//3
            z, v = odd_step(x)
            need(z == y and v == k and len(run(z)['word']) == r,
                 'inverse_regeneration')
    for n in range(3,164,4):
        N = 2**n-2
        R = run(N//2)
        need(R['word'] == [1]*(n-2)+[4] and R['end'] == (3**(n-1)-1)//8,
             'fixed_four_halving_cascade')
        if n >= 11:
            need(R['end'] > N, 'noncontracting_complete_cycle')
    # Independent memoized census, rather than replaying and hashing each orbit.
    cache = {1:(1,False,False,0)}
    for start in range(1,10001):
        n = start
        stack = []
        while n not in cache:
            stack.append(n)
            n = n//2 if n%2 == 0 else 3*n+1
        peak,hit1822,hit9232,steps = cache[n]
        for n in reversed(stack):
            peak=max(peak,n); hit1822 |= n == 1822; hit9232 |= n == 9232; steps+=1
            cache[n]=(peak,hit1822,hit9232,steps)
    counts=[sum(cache[n][1] for n in range(1,10001)),
            sum(cache[n][2] for n in range(1,10001)),
            sum(cache[n][0] == 9232 for n in range(1,10001))]
    need(counts == [3997,4011,1579], 'independent_memoized_census')
    return {'triple_universe':'r,s,t = 1..16', 'triple_count':16**3,
            'large_word_examples':examples,'growth_prefixes':'m = 1..80; q in {1,3,17}',
            'census':counts}

def kernel_knots(r,j):
    # Pairing from the RIGHT gives alternating unit ramps and plateaus.
    return Q((r-j+1)//2)

def kernel_moment(r,p):
    total = Q(0)
    for j in range(r):
        left,right=kernel_knots(r,j),kernel_knots(r,j+1)
        slope = right-left
        intercept = left-j*slope
        total += intercept*Q((j+1)**(p+1)-j**(p+1),p+1)
        total += slope*Q((j+1)**(p+2)-j**(p+2),p+2)
    return total

def check_kernel():
    c = Q(2,3)
    for r in range(1,65):
        for j in range(r):
            a,b=kernel_knots(r,j),kernel_knots(r,j+1)
            need(a>=0 and b>=0,'kernel_nonnegative_knots')
            for theta in (Q(1,3),Q(1,2),Q(2,3)):
                u=j+theta
                direct=sum((-1)**(r+d)*max(Q(0),d-u) for d in range(1,r+1))
                need(direct == a+(b-a)*theta,'kernel_independent_linear_interpolation')
        mass = 2*kernel_moment(r,0)
        fourth = kernel_moment(r,2)
        need(mass == Q(r*(r+1),2),'kernel_mass')
        need(fourth == Q(r*(r+1)*(r*(r+1)-1),24),'kernel_fourth_remainder')
        for p in range(13):
            direct=sum((-1)**(r+d)*(c+d)**p for d in range(-r,r+1))-c**p
            integral=Q(0)
            if p>=2:
                for j in range(0,p-1,2):
                    integral += 2*p*(p-1)*comb(p-2,j)*c**(p-2-j)*kernel_moment(r,j)
            need(direct == integral,'polynomial_peano_identity')
        lam=r*(r+1)
        expected=(1,lam,lam*lam-lam,lam**3-3*lam**2+3*lam,
                  lam**4-6*lam**3+17*lam**2-17*lam)
        for m,value in enumerate(expected):
            direct=sum((-1)**(r+d)*d**(2*m) for d in range(-r,r+1))
            need(value == direct,'explicit_euler_moments')
    return {'radii':'1..64','polynomial_degrees':'0..12','center':'2/3',
            'method':'Exact polynomial integration of linearly interpolated paired kernels'}

def exp_interval(x, terms=120):
    if x < 0 or x >= terms+2:
        raise ValueError('unsupported exponential certificate')
    term=Q(1)
    total=term
    for j in range(1,terms+1):
        term *= x/j
        total += term
    next_term=term*x/(terms+1)
    upper=total + next_term/(1-x/Q(terms+2))
    return total,upper

def ceil_decimal(q, places=24):
    scale=10**places
    z=(q.numerator*scale+q.denominator-1)//q.denominator
    return f'{z//scale}.{z%scale:0{places}d}'

def check_rational_certificates():
    # Certify log(10) with EXPONENTIAL Taylor bounds, not an atanh log routine.
    log10_lower=Q(2302585092994045684017991,10**24)
    log10_upper=log10_lower+Q(1,10**24)
    need(exp_interval(log10_lower)[1] < 10 < exp_interval(log10_upper)[0],
         'log10_interval_by_rational_exp')
    expected={10:'0.002916802930051020125517',100:'0.000004988801113793516135',
              1000:'0.000000007248352833489898',10000:'0.000000000009544628508874'}
    tails={}
    for j in range(1,5):
        N=10**j; L=j*log10_upper
        I4=(L/3+Q(1,9))/N**3
        I6=(L**2/5+2*L/25+Q(2,125))/N**5
        U=3*(1+Q(1,N+1))*(I4+I6)
        need(U < Q(expected[N]),'published_tail_decimal_outward')
        need(ceil_decimal(U) == expected[N],'independent_tail_decimal_match')
        tails[str(N)]=ceil_decimal(U)
    # Entirely small rational alternative to the audit's Galois magnitude check.
    exp_partial=sum(Q(12,5)**j/factorial(j) for j in range(9))
    need(exp_partial > 11,'log11_less_than_12_over_5')
    need(Q(5,4)**11 > 11,'root11_less_than_5_over_4')
    upper=Q(5,4)*((2*Q(12,5)-3)/11**3+(Q(12,5)-1)**2/11**4)
    need(upper == Q(136,73205) < Q(1,500),'curvature11_small_rational_bound')
    need(6**13 < 13*5**13,'root13_greater_than_6_over_5')
    need(Q(3,4)>Q(16,25),'cos_pi_over_13_greater_than_4_over_5')
    need(Q(1,500)-Q(9,5)*Q(6,5) == -Q(1079,500), 'conjugate_real_part_bound')
    need((11*10*12*4*12)%13 != 0,'compositum_degree_prime_to_13')
    # e > 5/2 and log(3/2) > 2/5 prove exp(1/e)<3/2.
    need(sum(Q(1,factorial(j)) for j in range(4))>Q(5,2),'e_greater_than_5_over_2')
    need(exp_interval(Q(2,5),32)[1] < Q(3,2),'exp_2_over_5_less_than_3_over_2')
    return {'tail_upper_decimals':tails, 'alternate_curvature11_upper':'136/73205',
            'conjugate_modulus_strict_lower':'1079/500',
            'degree_divisor':11*10*12*4*12,
            'warning':'Field-existence and all-order statements are proofs in AUDIT.md, not finite computations.'}

def main():
    root=Path(__file__).resolve().parent.parent
    targets=['session31-collatz/AUDIT.md','session31-collatz/verify_audit.py',
             'session31-square-blocks/AUDIT.md','session31-square-blocks/verify_square_blocks.py']
    def resolve_source(relative):
        first, rest = relative.split('/', 1)
        folder = first.removeprefix('session31-')
        portable = root / (rest if folder == 'root' else folder + '/' + rest)
        return portable if portable.is_file() else root / relative
    sources={relative:resolve_source(relative) for relative in targets}
    report={'status':'PASS: independent finite exact checks; see proofs for all-order scope',
            'collatz':check_collatz(),'kernel':check_kernel(),'rational':check_rational_certificates(),
            'check_counts':COUNT,
            'audited_sha256':{relative:hashlib.sha256(p.read_bytes()).hexdigest() for relative,p in sources.items()}}
    target=Path(__file__).with_name('independent_receipt.json')
    target.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps(report,indent=2,sort_keys=True))
if __name__=='__main__':
    main()
