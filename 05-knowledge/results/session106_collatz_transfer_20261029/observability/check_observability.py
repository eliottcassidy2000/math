"""Exact source-faithful survival-rank controls; no asymptotic convergence claim."""
from fractions import Fraction
import json
from pathlib import Path
from hashlib import sha256

def require(ok, message):
    if not ok: raise ValueError(message)

def step(n,b=3):
    require(n>0 and n%2==1,'positive odd input')
    x=b*n+1
    while x%2==0: x//=2
    return x

def predecessor(y):
    require(y>0 and y%2==1 and y%3,'admissible reverse endpoint')
    choices=(2,4) if y%3==1 else (3,5)
    for a in choices:
        num=(1<<a)*y-1
        require(num%3==0,'reverse integrality')
        x=num//3
        if x%3 and x>y:
            require(x>y and x%2==1 and step(x)==y,'exact proper inverse')
            return x,a
    raise ValueError('two exponents unexpectedly both divisible by nine')

def rank(a):
    a=[[Fraction(x) for x in row] for row in a]; r=0
    for c in range(len(a[0])):
        p=next((i for i in range(r,len(a)) if a[i][c]),None)
        if p is None: continue
        a[r],a[p]=a[p],a[r]
        v=a[r][c]; a[r]=[x/v for x in a[r]]
        for i in range(len(a)):
            if i!=r:
                v=a[i][c]; a[i]=[x-v*y for x,y in zip(a[i],a[r])]
        r+=1
    return r

def main():
    xs=[1]; exponents=[]
    for _ in range(80):
        x,a=predecessor(xs[-1]); xs.append(x); exponents.append(a)
    for j,x in enumerate(xs):
        orbit=x
        for t in range(j):
            require(orbit>1,'no premature target')
            orbit=step(orbit)
        require(orbit==1,'exact hitting clock')
        require(3**j*x<=32**j,'exponential witness size bound')
    rows=[]
    for n in (1,2,3,5,10,20,40):
        m=[[int(j>i) for j in range(1,n+1)] for i in range(n)]
        require(rank(m)==n,'survival response rank')
        rows.append({'horizon':n,'rank':n,'largest_witness':str(xs[n])})
    # Independent evaluation of the triangular matrix using literal trajectories.
    for j in range(1,41):
        v=xs[j]
        for i in range(40):
            require(int(v!=1)==int(j>i),'literal survival matrix entry')
            if v!=1: v=step(v)
    # Density of long-lived sources does not imply a fixed source survives forever.
    mersenne=[]
    for k in range(1,41):
        n=2**(k+1)-1; v=n
        for i in range(k):
            require(v==3**i*2**(k+1-i)-1 and v>1,'long-lived source')
            v=step(v)
        require(v==2*3**k-1 and v>1,'survives k steps')
        mersenne.append(n)
    return {'status':'PASS','scope':'Exact inverse chain through80; independent survival matrices through40; proofs all-order in note',
            'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
            'chain_first_16':xs[:16], 'inverse_exponents_first_15':exponents[:15],
            'rank_rows':rows,'max_chain_digits':len(str(xs[-1])),'long_survival_controls':len(mersenne)}

if __name__=='__main__': print(json.dumps(main(),indent=2,sort_keys=True))
