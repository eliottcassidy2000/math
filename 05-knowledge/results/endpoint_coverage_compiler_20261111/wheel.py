"""Exact small-prime wheel for 99*2^k+23; no primality inference from survival."""
from fractions import Fraction
from itertools import combinations
from math import gcd,lcm,isqrt
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent


def prime(n):
    return n >= 2 and all(n%d for d in range(2,isqrt(n)+1))


def order2(p):
    value, t = 2%p,1
    while value != 1:
        value = value*2%p
        t += 1
    return t


def merge(a,m,b,n):
    g = gcd(m,n)
    if (b-a)%g:
        return None
    t = ((b-a)//g)*pow(m//g,-1,n//g) % (n//g)
    return (a+m*t)%lcm(m,n),lcm(m,n)


def main():
    primes = [p for p in range(2,100) if prime(p)]
    active, inactive = [], []
    for p in primes:
        if p == 2 or 99%p == 0 or p == 23:
            inactive.append(p)
            continue
        period = order2(p)
        residues = [a for a in range(period) if (99*pow(2,a,p)+23)%p == 0]
        if residues:
            if len(residues) != 1:
                raise RuntimeError('One discrete-log class expected.')
            active.append(dict(p=p,period=period,residue=residues[0]))
        else:
            inactive.append(p)
    density = Fraction(1)
    for size in range(1,len(active)+1):
        for subset in combinations(active,size):
            state = (0,1)
            for item in subset:
                state = merge(*state,item['residue'],item['period'])
                if state is None:
                    break
            if state is not None:
                density += Fraction((-1)**size,state[1])
    period = lcm(*(a['period'] for a in active))
    if (period*density).denominator != 1:
        raise RuntimeError('Wheel density must give an integral survivor count.')
    if not prime(419) or 99*2**2+23 != 419:
        raise RuntimeError('Large-prime hostile failed.')
    modulus = lcm(period,order2(419))
    for k in range(1,5001):
        literal = all((99*pow(2,k,p)+23)%p for p in primes)
        compact = all(k%a['period'] != a['residue'] for a in active)
        if literal != compact:
            raise RuntimeError('Prime wheel mismatch.')
    for k in range(5,101):
        if bin(99*2**k+23)[2:] != '1100011'+'0'*(k-5)+'10111':
            raise RuntimeError('Binary regular language mismatch.')
    if any((2-a['residue'])%a['period'] == 0 for a in active):
        raise RuntimeError('419 must survive all small primes.')
    result = dict(primes_below_100=primes,active=active,inactive=inactive,
                  wheel_period=period,survivor_residues=int(period*density),
                  density_numerator=density.numerator,density_denominator=density.denominator,
                  order2_mod419=order2(419),composite_survivor_exponents=f'2+{modulus}*t, t>=1',
                  direct_wheel_controls=5000,binary_language_controls=96)
    (HERE / 'wheel.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
