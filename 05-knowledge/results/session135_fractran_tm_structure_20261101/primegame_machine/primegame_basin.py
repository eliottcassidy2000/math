#!/usr/bin/env python3
"""Total basin decider for the exact ordered Conway 1987 PRIMEGAME.

No search cutoff is used. A priori bound is proved in THEORY.md.
Spectator factors are retained in integer input/output but never factored.
"""
from math import comb, prod

PRIMES = (2,3,5,7,11,13,17,19,23,29)
FRACTIONS = ((17,91),(78,85),(19,51),(23,38),(29,33),(77,29),
             (95,23),(77,19),(1,17),(11,13),(13,11),(15,2),(1,7),(55,1))
LETTERS = 'ABCDEFGHIJKLMN'
CONTROLS = (0,11,13,17,19,23,29)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def valuations(n):
    if not isinstance(n,int) or isinstance(n,bool) or n < 1:
        raise ValueError('positive integer required')
    out=[]
    for p in PRIMES:
        e=0
        while n%p == 0:
            n//=p; e+=1
        out.append(e)
    return tuple(out),n


def integer(x, spectator=1):
    return spectator*prod(p**e for p,e in zip(PRIMES,x))


def vectors():
    return tuple((valuations(u)[0],valuations(v)[0]) for u,v in FRACTIONS)


VECTORS = vectors()
GUARDS = tuple(tuple(i for i,e in enumerate(b) if e) for a,b in VECTORS)
DELTAS = tuple(tuple(u-v for u,v in zip(a,b)) for a,b in VECTORS)


def step(x):
    for i,guard in enumerate(GUARDS):
        if all(x[j]>0 for j in guard):
            return tuple(e+delta for e,delta in zip(x,DELTAS[i])), LETTERS[i]
    raise RuntimeError('unreachable: denominator 1')


def integer_step(n):
    for i,(u,v) in enumerate(FRACTIONS):
        if n%v==0:
            return n//v*u,LETTERS[i]
    raise RuntimeError('unreachable: denominator 1')


def measures(x):
    a,b,c,d,e,f,g,h,i,j=x
    return a+c+i,b+d+g+h+i+j,e+f+g+h+i+j  # T,Q,k


def periodic(x):
    """Return least period 2/6 if already periodic, else 0.

    S129 proves 2 and 6 are the only possible least periods.
    """
    y=x; p2=False
    for t in range(1,7):
        y,_=step(y)
        if t==2: p2=(y==x)
    if y!=x: return 0
    return 2 if p2 else 6


def onehot_escape(x):
    """Iff this k<=1 input ultimately joins the spectator-scaled seed-2 orbit."""
    a,b,c,d,*ctrl=x
    require(sum(ctrl)<=1,'onehot predicate outside its proved domain')
    q=next((p for p,e in zip(PRIMES[4:],ctrl) if e),0)
    if q==0: return a+b>0
    if q in (11,13): return a+c>0 and b+d>0
    if q in (17,19,29): return a+c>0
    if q==23: return True
    raise RuntimeError('invalid control')


def state_bound(x):
    T,Q,k=measures(x)
    return (T+1)*comb(Q+2,2)*comb(k+6,6)


def classify_vector(x):
    """Return a certificate for basin 'prime', 'cycle2', or 'cycle6'.

    'prime' means eventual exact coalescence with the standard seed-2 orbit,
    not that the arbitrary input itself is a pure-power output. The witness
    is either a onehot state with the proved predicate or an actual cycle.
    """
    if len(x)!=10 or any(not isinstance(e,int) or isinstance(e,bool) or e<0 for e in x):
        raise ValueError('ten nonnegative integer exponents required')
    x=tuple(x); initial=x; bound=state_bound(x); steps=0
    while True:
        T,Q,k=measures(x)
        if k<=1:
            return dict(basin='prime' if onehot_escape(x) else 'cycle2',
                        witness=list(x), witness_kind='onehot',
                        transitions=steps, bound=bound)
        p=periodic(x)
        if p:
            return dict(basin='cycle'+str(p),witness=list(x),
                        witness_kind='periodic',transitions=steps,bound=bound)
        x,_=step(x); steps+=1
        # This is a theorem-backed consistency check, NOT a heuristic cutoff.
        require(steps<=bound,'proved finite-state bound violated')
        T0,Q0,k0=measures(initial)
        T1,Q1,k1=measures(x)
        require(T1==T0 and Q1<=Q0 and k1<=k0,'positive-control budget violated')


def classify_integer(n):
    x,s=valuations(n)
    result=classify_vector(x)
    result['spectator']=s
    result['input']=n
    return result


if __name__=='__main__':
    import argparse,json
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('integer',type=int)
    args=parser.parse_args()
    print(json.dumps(classify_integer(args.integer),sort_keys=True,indent=2))
