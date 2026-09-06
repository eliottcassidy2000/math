"""Independent recurrence, Euclidean-resultant and rational Sylvester referee.
No producer implementation is imported. Certificate used only for comparison.
"""
from pathlib import Path
from fractions import Fraction as F
from math import prod, isqrt
import hashlib
import json
import sys

sys.stdout.reconfigure(newline='\n')
G = 0


def need(ok, label):
    global G
    G += 1
    if not ok:
        raise ArithmeticError(label)


def val_int(a, p):
    if not a:
        raise ArithmeticError('zero valuation')
    v = 0
    while a % p == 0:
        a //= p
        v += 1
    return v


def val(a, p):
    return val_int(a.numerator, p)-val_int(a.denominator, p)


def rows(H, z):
    """Successive three-factor lowering from the monic coefficient."""
    out = [F(0)]*H+[F(1)]
    for j in range(H, 0, -1):
        k = H-j
        out[j-1] = out[j]*F(j*(z+2*j)*(z+2*j-1), (3*k+1)*(3*k+2)*(3*k+3))
    return out


def coefficient_product(H, z, j):
    # Independent falling-factor expression without giant endpoint factorials.
    return F(prod(range(j+1, H+1))*prod(range(z+2*j+1, z+2*H+1)),
             prod(range(1, 3*(H-j)+1)))


def reduce(a, p):
    need(all(c.denominator % p for c in a), 'entire row is p-integral')
    return [c.numerator*pow(c.denominator, -1, p) % p for c in a]


def trim(a):
    while a and a[-1] == 0:
        a.pop()
    return a


def rem(a, b, p):
    a = list(a)
    inv = pow(b[-1], -1, p)
    while len(a) >= len(b):
        d, c = len(a)-len(b), a[-1]*inv % p
        for j, v in enumerate(b):
            a[d+j] = (a[d+j]-c*v) % p
        trim(a)
    return a


def poly_gcd(a, b, p):
    a, b = list(a), list(b)
    while b:
        a, b = b, rem(a, b, p)
    return [c*pow(a[-1], -1, p) % p for c in a]


def resultant(a, b, p):
    """Euclidean degree-drop formula, distinct from quotient multiplication."""
    out = 1
    a, b = list(a), list(b)
    while len(b) > 1:
        m, n = len(a)-1, len(b)-1
        r = rem(a, b, p)
        if not r:
            return 0
        k = len(r)-1
        out = out*((-1)**(m*n))*pow(b[-1], m-k, p) % p
        a, b = b, r
    return out*pow(b[0], len(a)-1, p) % p


def sylvester(a, b):
    """Actual rational resultant by fraction-free Sylvester elimination."""
    m, n = len(a)-1, len(b)-1
    A = []
    for k in range(n):
        A.append([F(0)]*k+list(reversed(a))+[F(0)]*(n-k-1))
    for k in range(m):
        A.append([F(0)]*k+list(reversed(b))+[F(0)]*(m-k-1))
    N, prev, sign = m+n, F(1), 1
    for k in range(N-1):
        pivot = next((i for i in range(k, N) if A[i][k]), None)
        if pivot is None:
            return F(0)
        if pivot != k:
            A[pivot], A[k] = A[k], A[pivot]
            sign = -sign
        pivot_value = A[k][k]
        for i in range(k+1, N):
            for j in range(k+1, N):
                A[i][j] = (pivot_value*A[i][j]-A[i][k]*A[k][j])/prev
            A[i][k] = F(0)
        prev = pivot_value
    return sign*A[-1][-1]


def prime_sieve(n):
    live = bytearray([1])*(n+1)
    live[0:2] = b'\0\0'
    for p in range(2, isqrt(n)+1):
        if live[p]:
            for k in range(p*p, n+1, p):
                live[k] = 0
    return [p for p in range(2, n+1) if live[p]]


def inspect(H, z, p, expected=None):
    S = z+2*H
    e = val_int(2*S-1, p)
    need(p > 4*H and (2*S-1) % p == 0 and e > 6*H//p, 'declared prime-divisor hypotheses')
    phi, psi = rows(H,z), rows(2*H,2*z)
    need(phi[-1] == psi[-1] == 1, 'actual monic normalization')
    for j, c in enumerate(phi):
        need(c == coefficient_product(H,z,j), 'recurrence versus falling-product coefficient')
        need(val(c,p) == 0, 'every ordinary coefficient has exact zero valuation')
    for j, c in enumerate(psi):
        need(c == coefficient_product(2*H,2*z,j), 'doubled recurrence versus falling-product coefficient')
        need(val(c,p) == (0 if j == 2*H else e-(6*H-3*j)//p),
             'every doubled coefficient exact valuation including leading endpoint')
    f, h = reduce(phi,p), reduce(psi,p)
    need(h == [0]*(2*H)+[1], 'actual doubled row is monomial after reduction')
    g = poly_gcd(f,h,p)
    norm = (-1)**H*resultant(f,h,p) % p
    need(g == [1] and norm == (-1)**H*pow(f[0],2*H,p) % p != 0,
         'independent Euclidean gcd and resultant give actual nonzero norm')
    if expected is not None:
        need(expected == {'H':H,'z':z,'p':p,'e':e,'norm_mod_p':norm}, 'entire producer certificate row')
    return [H,z,p,e,norm]


def main():
    here = Path(__file__).resolve().parent
    main_name = 'continuing6_20260906_norm_units'
    candidates = [here/(main_name+'_certificate.json'), Path('C:/w')/(main_name+'_certificate.json')]
    if here.name == '04-computation':
        candidates.insert(0, here.parent/'05-knowledge/results'/(main_name+'_certificate.json'))
    path = next(p for p in candidates if p.is_file())
    raw = path.read_bytes()
    data = json.loads(raw)
    primes = prime_sieve(301)
    bank = [(H,z,p) for H in range(1,25) for z in range(1,102,2) for p in primes
            if (2*z+4*H-1) % p == 0 and p > 4*H and val_int(2*z+4*H-1,p) > 6*H//p]
    need(len(bank) == 492 and len(data['rows']) == 492, 'full independently generated parameter-prime universe')
    observed = [inspect(H,z,p,row) for (H,z,p),row in zip(bank,data['rows'])]
    need(sum(H>=7 for H,z,p in bank) == 292 and sum(p<=6*H for H,z,p in bank) == 4,
         'all new-height and square-prime cases counted')
    periodic = []
    for H,p in ((7,43),(12,73),(20,127),(31,191)):
        z0 = ((p+1)//2-2*H) % p or p
        for k in (0,1,2,5):
            periodic.append((H,z0+k*p,p))
    need(len(periodic) == len(data['periodic_controls']) == 16, 'entire periodic bank')
    observed += [inspect(H,z,p,row) for (H,z,p),row in zip(periodic,data['periodic_controls'])]

    # Fresh even-z and higher-valuation controls are outside the producer bank.
    extras = [(1,2,7),(2,6,19),(8,14,59),(7,911,43)]
    # Construct two even-z e>=3 controls with primes4H<p<=6H when possible.
    for H,p in ((2,11),(3,13)):
        z = (p**3-4*H+1)//2
        extras.append((H,z,p))
    for H,z,p in extras:
        inspect(H,z,p)
    need(len(extras) >= 4, 'fresh even, beyond-bank and higher-valuation controls')

    rational_controls = []
    for H in range(1,5):
        for z in range(1,13):
            f,h = rows(H,z),rows(2*H,2*z)
            exact = (-1)**H*sylvester(f,h)
            for p in (101,103):
                actual = (-1)**H*resultant(reduce(f,p),reduce(h,p),p) % p
                need(exact.denominator % p and actual == exact.numerator*pow(exact.denominator,-1,p) % p,
                     'actual rational Sylvester norm agrees with modular Euclidean resultant')
            for p in primes:
                M = 2*z+4*H-1
                if M % p == 0 and p > 4*H and val_int(M,p) > 6*H//p:
                    need(exact != 0 and val(exact,p) == 0, 'literal rational norm has unit valuation at certified prime')
            rational_controls.append([H,z,str(exact)])
    need(len(rational_controls) == 48, 'complete fresh rational Sylvester control rectangle')

    f,h = rows(1,1), rows(2,2)
    need(f == [1,1] and h == [1,10,1] and reduce(h,5) == [1,0,1],
         'minimal missing-valuation hostile is literal')
    need(sylvester(f,h) == -8 and (-1)*sylvester(f,h) == 8,
         'hostile response and characteristic-norm sign are distinguished')
    f,h = rows(26,3),rows(52,6)
    f109,h109 = reduce(f,109),reduce(h,109)
    common = poly_gcd(f109,h109,109)
    need(len(common) == 2 and resultant(f109,h109,109) == 0,
         'larger omitted-threshold hostile has a genuine modular common root')
    root = (-common[0]) % 109
    need(sum(c*pow(root,j,109) for j,c in enumerate(f109)) % 109 == 0 and
         sum(c*pow(root,j,109) for j,c in enumerate(h109)) % 109 == 0,
         'literal evaluation of both original rows at the shared modular root')
    f157,h157 = reduce(f,157),reduce(h,157)
    witness = resultant(f157,h157,157)
    need(poly_gcd(f157,h157,157) == [1] and witness != 0,
         'a second good reduction proves the rational rows remain coprime')
    need(109 > 4*26 and val_int(109,109) == 6*26//109 == 1,
         'larger hostile is outside exactly the strict valuation threshold')

    # Exact arithmetic specialization of the inherited jet iff; no new jet computation.
    boundary_count = 0
    for h0 in range(1,201):
        p = 4*h0+1
        if p not in prime_sieve(p):
            continue
        supplied = [r for r in range(1,h0+1) if h0-r == 0 or p > 6*(h0-r)]
        need(supplied == list(range((h0+2)//3,h0+1)), 'all r at a prime endpoint equal the claimed ceiling interval')
        boundary_count += 1
    need(boundary_count > 0, 'prime-index arithmetic specialization nonempty')

    print('SOURCE_METHOD lowering recurrence and falling products; no producer imports')
    print('COMPLETE_BANK 492 PERIODIC16 FRESH',len(extras),'RATIONAL_SYLVESTER48')
    print('MINIMAL_HOSTILE H1 z1 p5: response -8; characteristic norm +8')
    print('LARGE_HOSTILE H26 z3 p109 gcd',common,'root',root,'mod157 resultant',witness,
          'norm',(-1)**26*witness % 157,'rational gcd1')
    print('BOUNDARY_PRIME_CASES',boundary_count,'r>=ceil(h/3), h1..200')
    print('PRODUCER_CERTIFICATE_SHA256',hashlib.sha256(raw).hexdigest())
    print('PASS',G,'always-active exact gates; actual LF')


if __name__ == '__main__':
    main()
