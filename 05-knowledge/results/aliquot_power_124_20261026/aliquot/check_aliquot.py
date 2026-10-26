#!/usr/bin/env python3
"""Complete finite inverse fibres and controls for the unique-6 theorem.

Two independent sum implementations, explicit exceptions, standard library.
The all-perfect-number theorem is proved in README, not inferred from a scan.
"""
from pathlib import Path
from math import isqrt
import argparse
import json

HERE=Path(__file__).resolve().parent
checks=0

def require(value,message):
    global checks
    checks+=1
    if not value:
        raise RuntimeError(message)

def prime(n):
    if n<2:return False
    return all(n%d for d in range(2,isqrt(n)+1))

def s_factor(n):
    """Multiplicative sigma from prime powers, then subtract n."""
    rest=n; sigma=1; p=2
    while p*p<=rest:
        term=1;power=1
        while rest%p==0:
            rest//=p;power*=p;term+=power
        sigma*=term
        p+=1
    if rest>1:sigma*=1+rest
    return sigma-n

def s_divisors(n):
    """Direct paired divisors; no factorization implementation imported."""
    if n==1:return 0
    total=1
    for d in range(2,isqrt(n)+1):
        if n%d==0:
            total+=d
            if d*d!=n:total+=n//d
    return total

LIMIT=60
BOUND=(LIMIT-1)**2
frozen=json.loads((HERE/'INHERITED_FIBRES.json').read_text())
fibres={m:[] for m in range(2,LIMIT+1)}
for n in range(1,BOUND+1):
    m=s_factor(n)
    require(m==s_divisors(n),'independent divisor implementations disagree')
    require((m==0)==(n==1),'zero fibre')
    require((m==1)==prime(n),'one fibre equals primes')
    if m>1:
        require(n<=(m-1)**2,'sharp finite inverse bound')
        square_prime=isqrt(n)**2==n and prime(isqrt(n))
        require(square_prime or 4*n<=(m-1)**2,'stronger non-prime-square bound')
        require((n==(m-1)**2)==square_prime,'equality boundary is exactly prime squares')
    if m in fibres:fibres[m].append(n)

for m in range(2,LIMIT+1):
    # Independent target-wise path: direct divisors in the smaller universe,
    # plus the exceptional prime square. No s_factor values are reused.
    candidates={n for n in range(2,(m-1)**2//4+1) if s_divisors(n)==m}
    if prime(m-1):candidates.add((m-1)**2)
    require(sorted(candidates)==fibres[m],'target-wise finite search differs')
    require(fibres[m]==frozen['fibres'][str(m)],'inherited exact prefix differs')

small={2:[],3:[4],4:[9],5:[],6:[6,25],7:[8],8:[10,49],9:[15],10:[14],11:[21],12:[121],13:[27,35]}
for m,expected in small.items():require(fibres[m]==expected,'user list or small correction')
untouchables=[m for m in fibres if not fibres[m]]
require(untouchables==[2,5,52],'complete untouchables through60')
require([m for m in untouchables if prime(m)]==[2,5],'finite prime untouchables through60')

# A source cutoff equal to the target cutoff cannot certify the fibres.
require(s_divisors(59*59)==60 and 59*59>60,'source60 cutoff misses3481 in fibre60')
# Prime-square and semiprime guards are different; repeated primes do not work.
require(s_divisors(3*3)==4 and 3+3+1==7,'equal-prime Goldbach hostile')
require(s_divisors(4*4)==15 and 4+1==5,'composite-square hostile')
require(s_divisors(28)==28 and s_divisors(27*27)!=28,'perfect28 has no naive square predecessor')

goldbach_controls=[]
for p in [3,5,7,11,13,17,19]:
    for q in [3,5,7,11,13,17,19]:
        if p>=q:continue
        require(s_divisors(p*q)==p+q+1,'distinct-prime semiprime identity')
        goldbach_controls.append([p,q,p*q,p+q+1])

square_to_perfect=[]
for p in range(2,1001):
    if not prime(p):continue
    require(s_factor(p*p)==p+1,'prime-square identity')
    if s_divisors(p+1)==p+1:square_to_perfect.append([p,p*p,p+1])
require(square_to_perfect==[[5,25,6]],'finite control for unique6 theorem')

perfect_controls=[]
for r in [2,3,5,7,13,17,19,31]:
    mersenne=2**r-1;N=2**(r-1)*mersenne
    require(prime(r) and prime(mersenne),'declared Mersenne control')
    require(s_factor(N)==N,'Euclid-Euler positive control')
    require((N-1)%3==0 if r>2 else N==6,'all later even-perfect predecessors fail mod3')
    if r>2:require(N-1>3,'divisibility3 is a genuine composite obstruction')
    perfect_controls.append({'r':r,'N':N,'candidate_square_base':N-1,'candidate_base_mod3':(N-1)%3})

binary_controls=[]
for r in range(1,21):
    n=2**r;target=n-1
    require(s_divisors(n)==target,'binary divisor chain')
    binary_controls.append([r,n,target])
require([s_divisors(n) for n in [8,7,1]]==[7,1,0],'8 to7 to1 to0')
require([s_divisors(n) for n in [25,6]]==[6,6],'25 to6 fixed point')

out={'status':'PASS','checks':checks,'universe':{'target_max':LIMIT,'complete_composite_source_bound':BOUND,'prime_square_perfect_control_prime_max':1000,'binary_exponents':[1,20]},'fibres':{str(m):v for m,v in fibres.items()},'exceptional_fibres':{'0':[1],'1':'all primes; infinite, proved in the note'},'untouchables_through60':untouchables,'perfect_prime_square_controls':square_to_perfect,'euclid_euler_controls':perfect_controls,'binary_controls':binary_controls,'distinct_goldbach_controls':goldbach_controls,'scope':'The inverse-bound proof certifies all listed fibres. The general unique6 theorem uses Euclid-Euler and mod3, not the finite perfect controls. Odd-untouchable uniqueness and the complete list of prime untouchables remain OPEN.'}
parser=argparse.ArgumentParser();parser.add_argument('--output');args=parser.parse_args()
encoded=json.dumps(out,sort_keys=True,indent=2)+'\n'
if args.output:Path(args.output).write_text(encoded)
else:print(encoded,end='')
