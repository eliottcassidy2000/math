from orbital import rows7, induced_table
from itertools import product
from random import Random
from functools import cache

G=tuple((a,b) for a in range(7) for b in (1,2,4))
small=tuple(s for s in range(128) if s.bit_count()<=2)

@cache
def move(s,a,b):
    return sum(1<<((a+b*u)%7) for u in range(7) if s>>u&1)

@cache
def orbit(x,y):
    return min((move(x,a,b),move(y,a,b)) for a,b in G)

def subs(s):
    t=s
    while True:
        yield t
        if t==0: break
        t=(t-1)&s

def lookup(H):
    reps={orbit(a,b) for a,b in product(small,repeat=2)}
    data={p:H(*p) for p in sorted(reps)}
    c={}
    for a,b in sorted(reps):
        c[a,b]=sum((-1)**((a^x).bit_count()+(b^y).bit_count())*data[orbit(x,y)]
                   for x,y in product(subs(a),subs(b)))
    return c

@cache
def full_orbits():
    return tuple(sorted({orbit(x,y) for x,y in product(range(128),repeat=2)}))

def predict(c,x,y):
    return sum(c[orbit(a,b)] for a,b in product(small,repeat=2) if a&x==a and b&y==b)

def inverse(c,N):
    return tuple((x,y) for x,y in full_orbits() if predict(c,x,y)==N)

def graph(x,y,extra=False):
    U=rows7((1,2,4))
    rows=[a|((not(x>>j&1))<<7)|((not(y>>j&1))<<8) for j,a in enumerate(U)]
    rows += [x|256, y]
    if extra:
        # z dominates U and w, and loses to v.
        rows[7]|=512
        rows.append(127|256)
    return rows

if __name__=='__main__':
    benchmark=tuple(tuple(len({orbit(x,y) for x,y in product(small,repeat=2)
                                if x.bit_count()==i and y.bit_count()==j})
                          for j in range(3)) for i in range(3))
    assert benchmark==((1,1,1),(1,3,7),(1,7,21)) and len(full_orbits())==824
    rng=Random(827)
    for extra in (False,True):
        H=lambda x,y:induced_table(graph(x,y,extra))[-1]
        c=lookup(H)
        assert len(c)==43 and all(n%2==0 for p,n in c.items() if p!=(0,0))
        for x,y in [(36,63),(96,115)]+[(rng.randrange(128),rng.randrange(128)) for _ in range(24)]:
            assert predict(c,x,y)==H(x,y)
        for a in range(128):
            for x,y in ((a,0),(a,127),(0,a),(127,a)):
                assert predict(c,x,y)%3==0
        if not extra:
            assert orbit(36,63) in inverse(c,2471)
            assert orbit(96,115) in inverse(c,2621)
    print('43 invariant probes recover 824 mask-pair orbits, and exact target inverses replay')
