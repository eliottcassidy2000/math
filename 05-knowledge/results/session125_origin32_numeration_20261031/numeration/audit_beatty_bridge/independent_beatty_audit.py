"""Independent exact-Q(sqrt(5)) audit; no floating point or disabled asserts."""
import hashlib
import importlib.util
import itertools
import json
import random
import re
import sys
from fractions import Fraction
from pathlib import Path

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'beatty_bridge.py'
spec=importlib.util.spec_from_file_location('beatty_under_test',SOURCE)
b=importlib.util.module_from_spec(spec)
spec.loader.exec_module(b)
checks=0

def require(ok,why):
    global checks
    checks+=1
    if not ok:
        raise RuntimeError(why)

# Sign of a+b*sqrt(5), by rational arithmetic and monotonicity of squaring.
def sign(a,c):
    if c == 0:
        return (a>0)-(a<0)
    if a == 0:
        return (c>0)-(c<0)
    if (a>0) == (c>0):
        return 1 if a>0 else -1
    magnitude=(a*a>5*c*c)-(a*a<5*c*c)
    return magnitude if a>0 else -magnitude

# Exact quadratic field pairs, used independently of the integer-shift code.
def add(x,y):
    return x[0]+y[0],x[1]+y[1]

def mul(x,y):
    return x[0]*y[0]+5*x[1]*y[1], x[0]*y[1]+x[1]*y[0]

BETA=(Fraction(3,2),Fraction(1,2))
BETA_INV=(Fraction(3,2),Fraction(-1,2))
TAU=(Fraction(-1,2),Fraction(1,2))

def subtract(x,y):
    return x[0]-y[0],x[1]-y[1]

def legal(w):
    return re.search('21*2',w) is None

def fractional(w):
    # w is MSB-first fractional digit string.
    x=(Fraction(0),Fraction(0))
    for d in reversed(w):
        x=mul(add((int(d),0),x),BETA_INV)
    return x

def independent_word(n):
    # All Fibonacci terms constructed first, not the even-subsequence recurrence.
    fs=[0,1]
    while fs[-1] <= max(n,3):
        fs.append(fs[-1]+fs[-2])
    weights=fs[2::2]
    while weights and weights[-1]>n:
        weights.pop()
    if not weights:
        return '0'
    out=[]
    for a in reversed(weights):
        d,n=divmod(n,a)
        out.append(str(d))
    return ''.join(out)

def main():
    require(mul(BETA,BETA_INV)==(1,0),'conjugate reciprocal')
    limit=250000
    # Monotone exact oracle: each positive-n next ceiling is previous+2 or +3;
    # compare the candidate to beta*n in Q(sqrt(5)), without isqrt or floats.
    shifts=[0]
    for n in range(1,limit+2):
        candidate=shifts[-1]+2
        if sign(2*candidate-3*n,-n)<0:
            candidate+=1
        require(sign(2*candidate-3*n,-n)>0,'oracle strict upper ceiling')
        require(sign(2*(candidate-1)-3*n,-n)<0,'oracle strict lower ceiling')
        shifts.append(candidate)
    cursor=0
    for n in range(limit+1):
        while shifts[cursor+1] <= n:
            cursor+=1
        require(b.shift(n)==shifts[n],'isqrt ceiling vs exact field comparison')
        require(b.parent(n)==cursor,'isqrt parent vs interval oracle')
        require(b.parent(n)==3*n-shifts[n],'parent/shift exact linear relation')
        require(b.divide(n)==(cursor,n-shifts[cursor]),'quotient/remainder oracle')
        word=independent_word(n)
        require(b.encode(n)==word,'division vs independent weight greedy')
        locked=word.rstrip('1').endswith('2')
        require(b.free_state(n)==(not locked),'state from last non-1 digit')
        require(shifts[n+1]-shifts[n]==(2 if locked else 3),'exact branching gap')
        # rho=(2P-3n-n*sqrt(5))/2; certify all threshold sides directly.
        require(sign(2*shifts[n]-3*n,-n)>=0,'rho nonnegative')
        require(sign(2*shifts[n]-3*n-2,-n)<0,'rho strictly below one')
        threshold=sign(2*shifts[n]-3*n+1,-n-1)
        require(threshold!=0 and (threshold<0)==(not locked),'strict tau threshold')

    fs=[0,1]
    for _ in range(24):
        fs.append(fs[-1]+fs[-2])
    weights=fs[2::2]
    words=0
    for length in range(11):
        for ds in itertools.product('012',repeat=length):
            w=''.join(ds)
            if not legal(w):
                continue
            words+=1
            n=sum(int(d)*weights[i] for i,d in enumerate(reversed(w)))
            p=sum(int(d)*weights[i+1] for i,d in enumerate(reversed(w)))
            rho=fractional(w[::-1])
            require(rho==subtract((p,0),mul(BETA,(n,0))),'fractional expansion exact identity')
            require(sign(*rho)>=0 and sign(*subtract(rho,(1,0)))<0,'strict finite fractional bounds')
            if legal('2'+w):
                t=fractional(w)
                require(sign(*subtract(t,TAU))<0,'locked fractional upper bound')
            require(b.decode(w)==n,'field/weights/decoder agreement')
            require(b.shift(n)==p,'legal weighted append zero')

    rng=random.Random(125_2618_09)
    large=0
    for bits,count in [(512,100),(4096,16),(16384,2)]:
        for _ in range(count):
            n=rng.getrandbits(bits) or 1
            s,p=b.shift(n),b.parent(n)
            require(sign(2*s-3*n,-n)>0 and sign(2*s-3*n-2,-n)<0,'large exact ceiling')
            require(sign(2*n-3*p,-p)>=0 and sign(2*n-3*(p+1),-p-1)<0,'large exact parent')
            require(b.shift(s)==3*s-n,'large double shift recurrence')
            require(b.decode(b.encode(n))==n,'large exact encode/decode')
            large+=1

    intervals=0
    for _ in range(10000):
        start=rng.randrange(limit//2)
        length=rng.randrange(1,limit//2)
        count=b.free_count(start,length)
        require(count==shifts[start+length]-shifts[start]-2*length,'telescoping count')
        # |count-tau*length|<1, with tau=(-1+sqrt5)/2.
        twice_rational=2*count+length
        require(sign(twice_rational-2,-length)<0 and sign(twice_rational+2,-length)>0,'strict balanced discrepancy')
        intervals+=1

    orbit_cases=100000
    child_cases=0
    for n in range(1,orbit_cases+1):
        s=shifts[n]
        t=3*s-n
        C=n*n-3*n*s+s*s
        require(b.orbit_norm(n)==C and C>0,'norm positive exact quadratic form')
        require(b.orbit_norm(s)==C,'norm shift invariance')
        require(n*s+C==(s-n)**2 and n*t+C==s*s and s*t+C==(t-s)**2,'all D(C) pair equations')
        word=independent_word(n)
        primitive=word.rstrip('0')
        depth=len(word)-len(primitive)
        # Recover primitive value by the already independent integer interval oracle.
        core=n
        for _ in range(depth):
            core=3*core-shifts[core]
        require(b.orbit_core(n)==(core,depth),'core and depth from trailing zeros')
        require(independent_word(core)==primitive,'primitive canonical word')
        require(b.iterate_shift(core,depth)==n,'core reconstructs input')
        require(b.orbit_core(s)==(core,depth+1),'shift increments depth exactly')
        require(b.orbit_norm(core)==C,'primitive core retains norm')
        require((C==1)==(primitive=='1') and (C==4)==(primitive=='2'),'norm one/four classification')
        D=2*s-3*n
        require(D*D-5*n*n==4*C and D%2==n%2,'Pell norm and parity identity')
        for digit in range(shifts[n+1]-shifts[n]):
            child=s+digit
            require(b.orbit_norm(child)==C+digit*D+digit*digit,'legal child norm increment')
            require(digit==0 or digit*D+digit*digit>=4,'later nonzero digit increment')
            child_cases+=1
    require(b.orbit_norm(5)==b.orbit_norm(7)==11,'norm collision')
    require(b.orbit_core(5)==(5,0) and b.orbit_core(7)==(7,0),'norm is not complete orbit invariant')
    require(b.orbit_core(0)==(0,0) and b.orbit_norm(0)==0,'zero orbit convention')
    require(b.orbit_norm(5)!=b.orbit_norm(b.parent(5)),'general quotient does not preserve norm')

    # Zero is a mandatory special case for both reciprocal floor and Wythoff.
    require(b.shift(0)==0 and b.parent(0)==0 and b.divide(0)==(0,0),'zero boundary')
    require(b.free_state(0) and b.gap(0)==3,'zero Free first gap')
    require(b.decode('')==0 and b.decode('000')==0,'empty/leading-zero boundary')
    require(b.encode(0)=='0','display zero')
    require(b.free_count(0,0)==0,'empty interval')
    require(b.shift(b.shift(1))==8 and sign(14-7,-3)>0 and sign(12-7,-3)<0,'iterated ceiling not power ceiling')

    report={
      'status':'PASS','optimized_python':not __debug__,'explicit_checks':checks,
      'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
      'exact_field_integer_cases':limit+1,'legal_fractional_words_through_length':10,
      'legal_fractional_words_checked':words,'large_integer_cases':large,'largest_bit_width':16384,
      'random_balanced_intervals':intervals,'orbit_norm_core_cases':orbit_cases,'legal_child_norm_cases':child_cases,
      'oracle':'Exact Q(sqrt(5)) rational sign comparisons and independent Fibonacci weights; no isqrt in oracle, no floats.',
      'scope':'Canonical arithmetic only, no graph existence inference; no source or remote edits.'
    }
    text=json.dumps(report,indent=2)+'\n'
    print(text,end='')
    if len(sys.argv)>1:
        Path(sys.argv[1]).write_text(text)

if __name__=='__main__':
    main()
