"""Finite controls for valuation buckets, trail contraction and translation-fixed overlays."""
from itertools import permutations
from math import gcd
from pathlib import Path
import sys


def oddpart(n):
    return n // (n & -n)


def jump(n):
    return oddpart(3*n-1)


def label(n):
    q = 3*oddpart(n)-1
    return (q & -q).bit_length()-1


def terminal_class(word):
    # Invert the odd steps symbolically: x_0=(2^A*y+C)/3^m.
    a = 0
    c = 0
    three = 1
    for b in word:
        assert b > 0
        c = 3*c+(1<<a)
        a += b
        three *= 3
    residue = (-c*pow(1<<a,-1,three)) % three
    return a, c, three, residue


def exponent(word, anchor):
    """Return a final halving exponent linking the prescribed word to the anchor."""
    _, _, modulus, residue = terminal_class(word)
    top = 3*modulus
    target = ((3*residue-1)*pow(anchor,-1,top)) % top
    power = 1
    period = 2
    mod = 3
    if target%3 == 1:
        power = 2
    while mod < top:
        mod *= 3
        power = next(k for k in (power,power+period,power+2*period)
                     if pow(2,k,mod)==target%mod)
        period *= 3
    assert pow(2,power,top)*anchor%top == (3*residue-1)%top
    return power


def explicit_seed(word, anchor, e):
    value, rem = divmod((1<<e)*anchor+1,3)
    assert not rem
    for b in reversed(word):
        value, rem = divmod((value<<b)+1,3)
        assert not rem
    return value


def euler(rows, omit=None):
    """An Euler trail in balanced rows or with prescribed endpoints after path deletion."""
    n=len(rows)
    rest=rows[:]
    if omit:
        for u,v in zip(omit,omit[1:]):
            rest[u] ^= 1<<v
        start=omit[-1]
    else:
        start=0
    stack=[start]
    rev=[]
    while stack:
        u=stack[-1]
        if not rest[u]:
            rev.append(stack.pop())
        else:
            bit=rest[u]&-rest[u]
            rest[u]^=bit
            stack.append(bit.bit_length()-1)
    return rev[::-1]


def cylinder(word, sign=-1):
    a,c,three,_=terminal_class(word)
    modulus=1<<(a+1)
    return ((1<<a)-sign*c)*pow(three,-1,modulus)%modulus


def euclid_bounds():
    earlier=Path(__file__).resolve().parents[1]/"chanpak_odd_square_spines_20261122"
    sys.path.insert(0,str(earlier))
    from near_witnesses import (A,D,EDGE,WINNING_SLICES,flipped,literal_receiver,nontags,
                                ordinal,star_rows)
    receiver, recv_tags=literal_receiver(nontags(4))
    tags=12|(35<<len(A))|(recv_tags<<(len(A)+len(D)))
    for j,(*_, record) in enumerate(WINNING_SLICES):
        x,y,n,offset=record
        for orientation,middle in enumerate((flipped(D,*EDGE),D)):
            rows=star_rows(ordinal(ordinal(A,middle),receiver),tags,x,y)
            degrees=[r.bit_count() for r in rows]
            pieces=sum(max(0,2*v+1-n) for v in degrees)
            visits=sorted((max(v,n-1-v) for v in degrees),reverse=True)
            weight=sum(i*k for i,k in enumerate(visits))
            lower=((916917,916914),(916001,915998))[j][orientation]
            assert 1<<(weight+lower-1)>=3**(lower-1) and 1<<(weight+lower)<3**lower
            print("E61 midpoint offset",offset,"orientation",orientation,"pieces",pieces,
                  "minimum valuation excess",weight,"odd period lower bound",lower)


def check_interchanges():
    labels=(1,7,3,2)
    before={(1,7),(1,3),(7,3),(3,2),(1,2),(2,7)}
    after=before ^ {(1,7),(7,1)}
    def words(edges):
        return {p for p in permutations(labels)
                if all(edge in edges for edge in zip(p,p[1:]))}
    low,high=words(before),words(after)
    assert low-high=={(1,7,3,2)}
    assert high-low=={(7,1,3,2),(3,2,7,1),(2,7,1,3)}
    assert len(low)==3 and len(high)==5
    for x,p in ((541,(1,7,3,2)),(555,(7,1,3,2))):
        assert x==cylinder(p)
        for b in p:
            assert label(x)==b
            x=jump(x)
        assert x==5
    for i in range(4):
        for u in range(1,8):
            for v in range(u+1,60):
                w=[4]*i+[u,v,5,1]
                w2=[4]*i+[v,u,5,1]
                _,c,m,r=terminal_class(w)
                _,c2,m2,r2=terminal_class(w2)
                same=(v-u)%(2*3**(i+1))==0
                assert m==m2 and (r==r2)==same
                if same:
                    assert abs(c-c2)//m == (1<<(4*i))*abs((1<<v)-(1<<u))//(3**(i+2))
    print("3-to-5 flip seed residues modulo 16384:",sorted(map(cylinder,low)),"to",sorted(map(cylinder,high)))


def check():
    families=Path(__file__).resolve().parents[3]/"session212/cycle-levels/classes7.rows"
    for number,line in enumerate(families.read_text().splitlines(),1):
        rows=list(map(int,line.split()))
        if any(row.bit_count()!=3 for row in rows):
            continue
        paths=[p for p in permutations(range(7))
               if all(rows[u]>>v&1 for u,v in zip(p,p[1:]))]
        arcs={(u,v) for u in range(7) for v in range(7) if rows[u]>>v&1}
        for p in paths:
            completion=euler(rows,p)
            tour=list(p)+completion[1:]
            assert len(tour)==22 and tour[-1]==tour[0]
            assert set(zip(tour,tour[1:]))==arcs
        print("regular order-seven row",number,"paths",len(paths))
    for differences in ((1,2,3),(1,2,4)):
        allowed=set(differences)
        paths=sum(all((v-u)%7 in allowed for u,v in zip(p,p[1:]))
                  for p in permutations(range(7)))
        assert paths == (175 if differences==(1,2,3) else 189)
        fixed=0
        # Anchor the first difference for cyclic rotation of the base word.
        for tail in permutations(differences[1:]):
            pattern=(differences[0],)+tail
            trail=[0]
            for b in pattern*7:
                trail.append((trail[-1]+b)%7)
            if len(set(zip(trail,trail[1:])))==21:
                fixed+=1
        assert fixed == (2 if gcd(sum(differences),7)==1 else 0)
        print("connection set",differences,"full-translation-fixed tours",fixed)
    for n in range(2,7):
        words=list(permutations(range(1,n+1)))
        residues={cylinder(p) for p in words}
        assert len(residues)==len(words)
        modulus=1<<(1+n*(n+1)//2)
        for p in words:
            x=cylinder(p)
            assert x%2 and (-x)%modulus==cylinder(p,+1)
            for b in p:
                assert label(x)==b
                x=jump(x)
    word=[3,4,5,3]
    for anchor in (1,5,7,17):
        e=exponent(word,anchor)
        start=explicit_seed(word,anchor,e)
        x=start
        for b in word:
            assert label(x)==b
            step=(3*x-1)//2
            assert jump(x)==oddpart(step) and label(step)==label(jump(x))
            x=jump(x)
        assert label(x)==e and jump(x)==anchor
        print("buckets",word,"anchor",anchor,"final exponent",e,"seed",start)
    assert [label(n) for n in (7,21)]==[2,1]
    check_interchanges()
    euclid_bounds()
    assert 2**(108+63)<3**108 and 2**(107+63)>=3**107


if __name__=="__main__":
    check()
