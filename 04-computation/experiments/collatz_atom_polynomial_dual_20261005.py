"""Signed polynomial lower certificates for one discrete probability atom.

All arithmetic is exact. Moment intervals are premises, not invented oracle
answers. Actual Collatz controls use an explicitly replayed finite ROOT head.
"""
from fractions import Fraction as F
from math import factorial
import json


def natural(n, least=0):
    if type(n) is not int or n < least:
        raise ValueError('exact natural required')
    return n


def rational(x):
    if type(x) not in (int,F):
        raise ValueError('exact rational required')
    return F(x)


def node(m):
    natural(m)
    return F(1,1 << m)


def multiply(a,b):
    c = [F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            c[i+j] += x*y
    return tuple(c)


def evaluate(poly,x):
    value=F(0)
    for c in reversed(poly):
        value=value*x+c
    return value


def selector(m,d):
    natural(m);natural(d)
    center=node(m)
    poly=(F(0),)*d+(center**(-d),)
    for j in tuple(range(m))+(m+1,):
        x=node(j)
        poly=multiply(poly,(-x/(center-x),1/(center-x)))
    return poly


def tail_constant(m):
    natural(m)
    c=F(1)
    for j in range(1,m+1):
        c /= 1-node(j)
    return c


def conditioning(m,d):
    natural(m);natural(d)
    c=F(2**(m*d)*(2**(m+1)+1))
    for j in range(m):
        c *= F(2**j+1)/(1-node(m-j))
    return c


def moment_interval_bounds(m,d,intervals):
    """Conditional bounds from truthful intervals for moments0..degree.

    This validates exact types and elementary probability bounds, not existence
    of a law realizing an arbitrary supplied packet. The caller must certify
    that all intervals refer to the SAME probability on {2^-j:j>=0}.
    """
    poly=selector(m,d)
    if type(intervals) is not tuple or len(intervals) != len(poly):
        raise ValueError('one exact interval per required moment')
    packet=[]
    for interval in intervals:
        if type(interval) is not tuple or len(interval)!=2:
            raise ValueError('interval must be an exact pair')
        lo,hi=map(rational,interval)
        if not 0 <= lo <= hi <= 1:
            raise ValueError('probability moment interval outside [0,1]')
        packet.append((lo,hi))
    if packet[0] != (F(1),F(1)):
        raise ValueError('normalized probability moment0 must be exactly1')
    lower=upper=F(0)
    for c,(lo,hi) in zip(poly,packet):
        lower += c*(lo if c>=0 else hi)
        upper += c*(hi if c>=0 else lo)
    error=tail_constant(m)/4**d
    return {'lower':lower,'upper':upper+error,'selector_expectation_upper':upper,
            'tail_error':error,'coefficient_norm':sum(map(abs,poly),F(0))}


def geometric_moment(ell,missing=None):
    natural(ell)
    value=F(2**ell,2**(ell+1)-1)
    if missing is None:
        return value
    natural(missing)
    mass=node(missing+1)
    return (value-mass*node(missing)**ell)/(1-mass)


def moment_readout(poly,oracle):
    return sum((c*oracle(j) for j,c in enumerate(poly)),F(0))


def determinant(matrix):
    a=[list(map(F,row)) for row in matrix]
    out=F(1)
    for i in range(len(a)):
        pivot=next((j for j in range(i,len(a)) if a[j][i]),None)
        if pivot is None:return F(0)
        if pivot!=i:
            a[i],a[pivot]=a[pivot],a[i];out=-out
        value=a[i][i];out*=value
        for j in range(i+1,len(a)):
            ratio=a[j][i]/value
            for k in range(i,len(a)):
                a[j][k]-=ratio*a[i][k]
    return out


def root_word(n,cap=512):
    natural(n,1);natural(cap)
    if n%2!=1:raise ValueError('positive odd source')
    out=[]
    for _ in range(cap+1):
        if n==1:return tuple(out)
        if len(out)==cap:break
        raw=3*n+1;a=(raw&-raw).bit_length()-1
        out.append(a);n=raw >> a
    raise ValueError('declared finite control cap exhausted')


def replay(n,word):
    natural(n,1)
    if n%2!=1 or type(word) is not tuple:raise ValueError('exact source/word required')
    for a in word:
        natural(a,1)
        if n==1:raise ValueError('root padding')
        raw=3*n+1
        if (raw&-raw).bit_length()-1 != a:raise ValueError('false valuation')
        n=raw >> a
    if n!=1:raise ValueError('not a completed ROOT word')


def injection_atom(m):
    natural(m)
    n=6*m+3
    word=root_word(n)
    replay(n,word)
    t=len(word)-1;k=sum((a-1)//2 for a in word)
    return F(2*factorial(k)*factorial(t+1),factorial(t+k+2)),word


def head_moment_intervals(head,degree):
    """Certified finite head plus probability normalization, no guessed tail."""
    natural(degree)
    if type(head) is not tuple or not head:raise ValueError('nonempty exact head')
    head=tuple(map(rational,head))
    if any(x<0 for x in head) or sum(head)>1:raise ValueError('invalid subprobability head')
    remainder=1-sum(head)
    out=[(F(1),F(1))]
    for ell in range(1,degree+1):
        lower=sum((mass*node(j)**ell for j,mass in enumerate(head)),F(0))
        out.append((lower,lower+remainder*node(len(head))**ell))
    return tuple(out)


def main():
    checks=0
    def check(ok,label):
        nonlocal checks
        checks+=1
        if not ok:raise ValueError(label)

    examples=[]
    for m in range(13):
        check(tail_constant(m)<=F(32,9),'uniform product bound')
        previous=None
        for d in range(13):
            poly=selector(m,d);error=tail_constant(m)/4**d
            check(len(poly)==m+d+2,'declared polynomial degree')
            check(sum(map(abs,poly),F(0))==conditioning(m,d)==abs(evaluate(poly,F(-1))),
                  'exact coefficient conditioning')
            for j in range(65):
                value=evaluate(poly,node(j))
                if j==m:check(value==1,'target cardinal value')
                elif j<m or j==m+1:check(value==0,'finite exact zeros')
                else:check(-error<=value<=0,'whole-tail sign and finite controls')
                if previous is not None:
                    check(value>=evaluate(previous,node(j)),'pointwise increasing lower selectors')
            actual=node(m+1)
            lower=moment_readout(poly,geometric_moment)
            check(lower<=actual<=lower+error,'geometric exact dual interval')
            erased=moment_readout(poly,lambda ell:geometric_moment(ell,m))
            check(erased<=0<=erased+error,'missing-atom hostile')
            eps=F(1,2**(m+d+20))
            packet=tuple((F(1),F(1)) if ell==0 else
                         (max(F(0),geometric_moment(ell)-eps),min(F(1),geometric_moment(ell)+eps))
                         for ell in range(len(poly)))
            bounds=moment_interval_bounds(m,d,packet)
            check(bounds['lower']<=actual<=bounds['upper'],'certified moment uncertainty')
            check(lower-bounds['lower']<=eps*conditioning(m,d),'conditioning error budget')
            if m in (0,1,4) and d in (1,4,8):
                examples.append({'m':m,'d':d,'actual_atom':str(actual),
                                 'exact_lower':str(lower),'tail_error':str(error),
                                 'coefficient_norm':str(conditioning(m,d))})
            previous=poly

    hankel=[]
    for size in range(1,6):
        value=determinant([[geometric_moment(i+j,1) for j in range(size)] for i in range(size)])
        check(value>0,'strict Hankel positivity despite missing atom1')
        hankel.append(str(value))

    # Actual examples are representation checks of an already certified finite
    # head. They are NOT claims of a new ROOT source.
    certified=[injection_atom(m) for m in range(16)]
    head=tuple(x[0] for x in certified)
    actual_controls=[]
    for target in range(5):
        first=None
        for d in range(33):
            bounds=moment_interval_bounds(target,d,head_moment_intervals(head,target+d+1))
            check(bounds['lower']<=head[target]<=bounds['upper'],'actual lambda conditional interval')
            if bounds['lower']>0:
                first=d;break
        check(first is not None,'finite actual head re-encodes a positive certificate')
        actual_controls.append({'source':6*target+3,'index':target,'first_tested_d':first,
                                'actual_atom':str(head[target]),'lower_bound':str(bounds['lower']),
                                'word':list(certified[target][1]),
                                'scope':'already certified head; no additional ROOT source'})

    # A retained head omitting target16 admits a completion at node18 instead.
    # Every dual lower bound from those moment intervals must remain nonpositive.
    for d in range(13):
        bounds=moment_interval_bounds(16,d,head_moment_intervals(head,17+d))
        check(bounds['lower']<=0,'same finite head cannot invent next atom')

    bad=[lambda:selector(True,1),lambda:selector(1,1.0),lambda:selector(-1,0),
         lambda:tail_constant(False),lambda:geometric_moment(1,True),
         lambda:root_word(True),lambda:replay(1.0,()),
         lambda:moment_interval_bounds(0,0,((F(1),F(1)),(0.2,0.3))),
         lambda:moment_interval_bounds(0,0,((F(0),F(1)),(F(0),F(1)))),
         lambda:moment_interval_bounds(0,0,((1,1),(F(1),F(0)))),
         lambda:head_moment_intervals((F(2),),1)]
    for job in bad:
        try:job()
        except ValueError:check(True,'invalid exact premise rejected')
        else:raise ValueError('invalid premise accepted')

    print(json.dumps({'status':'PROVED signed dual and error bounds; universal Collatz OPEN',
        'polynomial_universe':{'target_indices':'0..12','damping_degrees':'0..12','tested_nodes':'0..64'},
        'geometric_examples':examples,'missing_atom1_Hankel_determinants':hankel,
        'actual_lambda_head':{'indices':'0..15','sources':'3,9,...,93','literal_cap':512},
        'actual_positive_controls':actual_controls,'no_new_source_controls':13,
        'invalid_inputs':len(bad),'checks':checks},indent=2,sort_keys=True))
    print('PASS: a positive signed-moment lower bound certifies the selected atom; independent moment premises remain essential.')


if __name__=='__main__':main()
