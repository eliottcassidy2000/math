"""Prime-divisor certificates for complementary regular norms.

Exact finite checks of an analytically proved, unbounded integer theorem.
Run normally or with -O. All validity gates remain active.
"""
from fractions import Fraction
from math import factorial
import json
from pathlib import Path
import sys

sys.stdout.reconfigure(newline='\n')
gates = 0

def require(ok, label):
    global gates
    gates += 1
    if not ok:
        raise ArithmeticError(label)

def primes(n):
    out = []
    for k in range(2, n+1):
        if all(k % p for p in out if p*p <= k):
            out.append(k)
    return out

def valuation(n, p):
    if not n:
        raise ArithmeticError('zero valuation not allowed')
    result = 0
    while n % p == 0:
        result += 1
        n //= p
    return result

def vq(q, p):
    return valuation(q.numerator, p)-valuation(q.denominator, p)

def row(H, z):
    return [Fraction(factorial(H)*factorial(z+2*H),
                     factorial(j)*factorial(z+2*j)*factorial(3*H-3*j))
            for j in range(H+1)]

def mod(q, p):
    require(q.denominator % p != 0, 'p-integral rational reduction')
    return q.numerator*pow(q.denominator, -1, p) % p

def multiply_mod(a, b, phi, p):
    H = len(phi)-1
    product = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            product[i+j] = (product[i+j]+x*y) % p
    while len(product)>H:
        c = product.pop()
        for j in range(H):
            product[len(product)-H+j] = (product[len(product)-H+j]-c*phi[j]) % p
    return product+[0]*(H-len(product))

def determinant(a, p):
    a = [list(r) for r in a]
    value = 1
    for j in range(len(a)):
        pivot = next((i for i in range(j, len(a)) if a[i][j] % p), None)
        if pivot is None:
            return 0
        if pivot != j:
            a[j], a[pivot] = a[pivot], a[j]
            value = -value
        value = value*a[j][j] % p
        inv = pow(a[j][j], -1, p)
        for i in range(j+1, len(a)):
            c = a[i][j]*inv % p
            for k in range(j, len(a)):
                a[i][k] = (a[i][k]-c*a[j][k]) % p
    return value % p

def certificate(H, z, p):
    S = z+2*H
    e = valuation(2*S-1, p)
    require(p>4*H and e>6*H//p, 'analytic prime-divisor hypotheses')
    phi, psi = row(H,z), row(2*H,2*z)
    require(phi[-1] == psi[-1] == 1, 'monic factorial rows')
    for c in phi:
        require(vq(c,p) == 0, 'every Phi coefficient is a unit')
    for j,c in enumerate(psi[:-1]):
        require(vq(c,p) == e-(6*H-3*j)//p, 'exact doubled coefficient valuation')
        require(vq(c,p)>0, 'every nonleading Psi coefficient vanishes')
    pp, qq = [mod(c,p) for c in phi], [mod(c,p) for c in psi]
    require(qq == [0]*(2*H)+[1], 'actual reduced doubled row is a monomial')
    # Build the actual response multiplication matrix; do not infer its norm
    # merely from a coefficient statistic.
    tmod = [0,1] if H>1 else [-pp[0] % p]
    response = [1]+[0]*(H-1)
    for _ in range(2*H):
        response = multiply_mod(response,tmod,pp,p)
    columns = []
    for j in range(H):
        columns.append(response)
        response = multiply_mod(response,tmod,pp,p)
    actual = determinant([[-columns[j][i] % p for j in range(H)]
                          for i in range(H)],p)
    predicted = (-1)**H*pow(pp[0],2*H,p) % p
    require(actual == predicted != 0, 'actual characteristic norm residue')
    return {'H':H,'z':z,'p':p,'e':e,'norm_mod_p':actual}

bank = []
prime_bank = primes(301)
for H in range(1,25):
    for z in range(1,102,2):
        for p in prime_bank:
            M = 2*z+4*H-1
            if M%p == 0 and p>4*H and valuation(M,p)>6*H//p:
                bank.append(certificate(H,z,p))
require(any(r['H']>=7 for r in bank), 'new complementary heights occur')
require(any(r['p']<=6*r['H'] and r['e']>=2 for r in bank), 'square-prime extension exercised')
require(any(r['H']==1 and r['z']==11 and r['p']==5 and r['e']==2 for r in bank),
        'small square-prime positive control')

# Periodic classes retain the prime residue coordinate; z is not bounded by p.
periodic = []
for H,p in ((7,43),(12,73),(20,127),(31,191)):
    z0 = ((p+1)//2-2*H) % p
    if z0 == 0:
        z0 += p
    for k in (0,1,2,5):
        periodic.append(certificate(H,z0+k*p,p))

# This hostile falsifies the monomial reduction when its exponent threshold
# is deleted. It does not falsify noncancellation over Q.
phi,psi = row(1,1),row(2,2)
require(phi == [Fraction(1),Fraction(1)], 'minimal hostile Phi')
require(psi == [Fraction(1),Fraction(10),Fraction(1)], 'minimal hostile Psi')
require([mod(c,5) for c in psi] == [1,0,1], 'e=1, p<=6H is not a monomial')
require(sum(c*(-1)**j for j,c in enumerate(psi)) == -8,
        'failed modular ansatz still has a nonzero rational response')

destination = Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).with_name(
    'continuing6_20260906_norm_units_certificate.json')
payload = {'status':'FINITE-EXACT checks of proved analytic prime-divisor theorem',
           'universe':'H=1..24; odd z=1..101; all prime divisors through301 satisfying p>4H and v_p(2z+4H-1)>floor(6H/p)',
           'rows':bank,'periodic_controls':periodic,
           'rule':'N_H(z) mod p = (-1)^H Phi_(H,z)(0)^(2H) != 0'}
destination.write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8',newline='\n')
print('finite_prime_divisor_cases',len(bank))
print('new_height_cases',sum(r['H']>=7 for r in bank))
print('square_prime_extension_cases',sum(r['p']<=6*r['H'] for r in bank))
print('periodic_controls',len(periodic))
print('largest_complementary_height',max(r['H'] for r in bank+periodic))
print('minimal_hostile_Psi_mod5 t^2+1; rational_response -8')
print('always_active_gates',gates)
print('PASS')
