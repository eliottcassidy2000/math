"""Independent exact referee for the nearby-endpoint resultant identity.
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


def endpoint(H,S):
    a=[F(0)]*H+[F(1)]
    for k in range(H):
        a[H-k-1]=a[H-k]*F((H-k)*(S-2*k)*(S-2*k-1),(3*k+1)*(3*k+2)*(3*k+3))
    return a


def residual(H,m):
    a=[F(0)]*m+[F(1)]
    for k in range(m):
        a[m-k-1]=a[m-k]*F((2*H-k)*(2*m+1-2*k)*(2*m-2*k),(3*k+1)*(3*k+2)*(3*k+3))
    return a


def evaluate(a,t):
    out=0
    for c in reversed(a):out=out*t+c
    return out


def test(H,m,z,p,record=None):
    need(0<=m<2*H and p>6*H and (2*z+4*H-2*m-1)%p==0,'complete integer endpoint hypotheses')
    f,h=endpoint(H,F(z+2*H)),endpoint(2*H,F(2*z+4*H))
    need(f==rows(H,z) and h==rows(2*H,2*z),'independent endpoint and z-index coefficient recurrences agree')
    a,b=reduce(f,p),reduce(h,p)
    R=reduce(residual(H,m),p)
    need(all(a) and all(R) and b==[0]*(2*H-m)+R,'full rational row reduction including all zero coefficients')
    half=endpoint(H,F(2*m+1,2))
    need(a==reduce(half,p),'auxiliary half-endpoint polynomial is the same reduction')
    actual=(-1)**H*resultant(a,b,p)%p
    small=resultant(R,a,p)
    predicted=(-1)**H*pow(a[0],2*H-m,p)*small%p
    need(actual==predicted,'actual full Euclidean resultant equals degree-m resultant formula')
    need(poly_gcd(a,b,p)==poly_gcd(a,R,p),'entire common-root obstruction preserved after zero-block removal')
    if record is not None:
        need(record=={'H':H,'m':m,'z':z,'p':p,'norm':actual,'residual_norm':small},'every field of the frozen producer row')
    if m==1:
        need(R==[2*H,1] and small==evaluate(a,-2*H)%p,'rank-one residual and evaluation normalization')
        need(reduce(endpoint(H,F(z+p+2*H)),p)==a and reduce(endpoint(2*H,F(2*(z+p)+4*H)),p)==b,
             'actual positive-z periodic lift, both rows')
    return [H,m,z,p,actual,small]


def main():
    here=Path(__file__).resolve().parent
    stem='continuing7_20260906_norm_rank_drop'
    choices=[here/(stem+'_certificate.json'),Path('C:/w/continuing7_20260906_root')/(stem+'_certificate.json')]
    if here.name=='04-computation':choices.insert(0,here.parent/'05-knowledge/results'/(stem+'_certificate.json'))
    path=next(p for p in choices if p.is_file())
    raw=path.read_bytes();data=json.loads(raw)
    primes=prime_sieve(40000)
    targets=[]
    for H in range(1,19):
        ps=[p for p in primes if 6*H<p<=6*H+44][:2]
        need(len(ps)==2,'complete prime choice at each declared height')
        for m in range(min(4,2*H-1)+1):
            for p in ps:targets.append((H,m,(p+2*m+1)//2-2*H,p))
    need(len(targets)==len(data['records'])==172,'complete producer finite universe')
    for args,row in zip(targets,data['records']):test(*args,record=row)

    full_rank=[]
    for H in range(1,11):
        p=next(p for p in primes if p>6*H)
        for m in range(2*H):
            full_rank.append(test(H,m,(p+2*m+1)//2-2*H,p))
    need(len(full_rank)==110 and any(m>=H for H,m,z,p,a,b in full_rank),
         'complete residual-rank range includes identities without dimension decrease')

    # A separate literal Sylvester determinant checks both orientation signs.
    rational_count=0
    for H in range(1,5):
        p=next(p for p in primes if p>6*H)
        for m in range(2*H):
            z=(p+2*m+1)//2-2*H
            phi,psi=endpoint(H,F(z+2*H)),endpoint(2*H,F(2*z+4*H))
            Q=residual(H,m)
            actual=(-1)**H*sylvester(phi,psi)
            small=F(1) if m==0 else sylvester(Q,phi)
            predicted=(-1)**H*phi[0]**(2*H-m)*small
            need(actual.denominator%p and predicted.denominator%p,
                 'rational Sylvester values have legal reductions')
            av=actual.numerator*pow(actual.denominator,-1,p)%p
            pv=predicted.numerator*pow(predicted.denominator,-1,p)%p
            need(av==pv,'literal full and residual Sylvester determinants agree after reduction')
            rational_count+=1
    need(rational_count==20,'all low-height residual ranks checked by rational matrices')

    need(len(data['linear_evaluations'])==80,'complete declared exceptional evaluation bank')
    for H,row in enumerate(data['linear_evaluations'],1):
        phi=endpoint(H,F(3,2));A=evaluate(phi,-2*H)
        terms=[phi[H-k]/(2*H)**k for k in range(H+1)]
        need(terms[0]==1 and terms[1]==F(1,16) and all(a>0 for a in terms),'positive alternating coefficients from independent recurrence')
        for k in range(1,H):
            ratio=F(H-k,2*H)*F((F(2*k)-F(3,2))*(F(2*k)-F(1,2)),(3*k+1)*(3*k+2)*(3*k+3))
            need(terms[k+1]/terms[k]==ratio and 0<ratio<F(1,2),'exact consecutive ratio and its positive upper bound')
        normalized=A/((-2*H)**H)
        need(normalized==sum((-1)**k*a for k,a in enumerate(terms)), 'literal exceptional value equals alternating sum')
        need(F(15,16)<=normalized<1 and (normalized==F(15,16))==(H==1),
             'sharp lower endpoint and nonzero auxiliary value')
        need(row=={'H':H,'numerator':A.numerator,'denominator':A.denominator}, 'every exact exceptional rational value')

    A1=evaluate(endpoint(1,F(3,2)),-2)
    A2=evaluate(endpoint(2,F(3,2)),-4)
    need(A1==F(-15,8) and A2==F(9601,640) and 9601 in primes and 28807 in primes,
         'minimal-height hostile and both exact prime certificates')
    phi,psi=endpoint(2,F(4802)),endpoint(4,F(9604))
    f,b=reduce(phi,9601),reduce(psi,9601)
    bad_gcd=poly_gcd(f,b,9601)
    need(b==[0,0,0,4,1] and bad_gcd==[4,1] and resultant(f,b,9601)==0,
         'the dropped residual-unit guard produces a genuine common linear factor')
    f2,b2=reduce(phi,28807),reduce(psi,28807)
    good=resultant(f2,b2,28807)
    need(good==20135 and poly_gcd(f2,b2,28807)==[1], 'same rational rows separated by another legal prime')
    need(data['hostile']=={'H':2,'z':4798,'p':9601,'other_prime':28807,'other_norm':good}, 'complete hostile record')
    new=test(7,1,9,43)
    need(new[-2]==16 and 43>42 and 4*11-1==43 and 7+4==11,
         'actual new carried-boundary dictionary and unit')
    need(all(p<=28 for p in primes if 45%p==0), 'old endpoint45 has no eligible prime')
    need(data['new_carried']=={'h':11,'r':4,'H':7,'m':1,'z':9,'p':43,'norm':16,'residual_norm':new[-1]},
         'entire new carried-boundary consumer record')

    print('METHOD coefficient recurrences, full Euclidean resultants/gcds, and20 rational Sylvester controls')
    print('PRODUCER172 FULL_RANK_RANGE110 EXCEPTIONAL_VALUES80')
    print('HOSTILE H2 z4798: gcd_mod9601=t+4, norm_mod28807=20135; rational gcd1')
    print('NEW_BOUNDARY h11 r4 H7 z9 norm16mod43; exact order3 by inherited jet iff')
    print('SCOPE dimension decreases only for m<H; modular unit iff, rational nonvanishing sufficient')
    print('PRODUCER_CERTIFICATE_SHA256',hashlib.sha256(raw).hexdigest())
    print('PASS',G,'always-active exact gates; actual LF')


if __name__=='__main__':main()
