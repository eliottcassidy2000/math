from itertools import combinations, product
from functools import cache

@cache
def h(rows):
    @cache
    def d(s,v):
        r=s^(1<<v)
        if not r: return 1
        return sum(d(r,w) for w in range(len(rows)) if rows[v]&r&(1<<w))
    return sum(d((1<<len(rows))-1,v) for v in range(len(rows))) or 1

def sub(rows, s):
    v=[i for i in range(len(rows)) if s&(1<<i)]
    return tuple(sum(1<<j for j,w in enumerate(v) if rows[u]&(1<<w)) for u in v)

@cache
def Q(rows):
    f=(1<<len(rows))-1
    return sum(h(sub(rows,s))*h(sub(rows,f^s)) for s in range(f+1))

def join(a,b):
    c=((1<<len(b))-1)<<len(a)
    return tuple(x|c for x in a)+tuple(x<<len(a) for x in b)

def glued(A,D):
    b=join(join((0,),A),join((0,),D))
    v=len(b)
    U=1|(1<<(len(A)+1))
    return tuple(x | ((1<<v) if not U&(1<<i) else 0) for i,x in enumerate(b))+(U,)

def sources(n):
    e=list(combinations(range(n),2))
    for bits in product((0,1), repeat=len(e)):
        out=[0]*n
        for (u,v),k in zip(e,bits):
            if k: u,v=v,u
            out[u]|=1<<v
        yield tuple(out)

count=0
bank=[a for n in range(4) for a in sources(n)]
for A,D in product(bank,repeat=2):
    if A and D:
        want=(2*Q(A)+h(A))*Q(D)-(Q(A)-h(A))*h(D)
        if bool(want%3) != bool(len(A)%2 and not len(D)%2 and h(A)%3 and h(D)%3):
            raise AssertionError('mod 3')
    elif A:
        want=Q(A)+h(A)
    elif D:
        want=3*Q(D)-h(D)
    else:
        want=1
    if h(glued(A,D))!=want:
        raise AssertionError((A,D,want,h(glued(A,D))))
    count+=1
print(count,'literal ordered-block two-contact checks pass.')
