from fractions import Fraction as Q
from itertools import combinations_with_replacement
from random import Random

def coords(path, base):
    lo, h = Q(0), Q(1)
    for side, length in path:
        w = h / base**length
        if side:
            lo += h-w
        h = w
    return lo+h/2,h

def cosh(path, other, base):
    x, y = coords(path, base)
    u, v = coords(other, base)
    return ((x-u)**2+y*y+v*v)/(2*y*v)

def distance(path, other):
    j = 0
    while j < min(len(path),len(other)) and path[j] == other[j]:
        j += 1
    return sum(edge[1] for edge in path[j:]+other[j:])

def threshold(base, alpha, k):
    s = alpha * base**k
    return (s+1/s)/2

def run(base, gap, types):
    rng = Random(127)
    nodes = [()]
    for j in range(6):
        new = []
        for p in nodes[2**j-1:2**(j+1)-1]:
            a,b = rng.choice(types)
            if rng.randrange(2):
                a,b=b,a
            new.extend((p+((0,a),),p+((1,b),)))
        nodes.extend(new)
    gamma = gap**2
    assert base*gamma > 1
    assert Q(5,4)/gamma < base
    for p,q in combinations_with_replacement(nodes,2):
        C = cosh(p,q,base)
        D = distance(p,q)
        assert (D==0 or C>=threshold(base,gamma,D))
        assert C<threshold(base,gamma,D+1)
        # independently verify distance-only recovery without ancestor data
        lo, hi = 0,1
        while C>=threshold(base,gamma,hi):
            lo,hi=hi,2*hi
        while hi-lo>1:
            k=(lo+hi)//2
            if C>=threshold(base,gamma,k):
                lo=k
            else:
                hi=k
        assert lo==D
        if D:
            comparable = p==q[:len(p)] or q==p[:len(q)]
            straight=threshold(base,Q(1),D)
            assert C!=straight and (C>straight)==comparable
    print(base,len(nodes),len(nodes)*(len(nodes)+1)//2)

run(3,Q(53,81),[(2,2),(2,4),(1,4),(3,4)])
run(5,Q(3,5),[(1,1)])
