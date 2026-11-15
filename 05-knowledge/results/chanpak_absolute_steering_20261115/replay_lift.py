"""Literal subset-DP and exact rational normal-form tests."""
from collections import deque
from fractions import Fraction
from functools import lru_cache
from itertools import product


def check(test):
    if not test:
        raise ArithmeticError("nested recontact audit failed")


def move(x, ch):
    g, s, d = x
    if ch == "E": return g+s-d, g+s, g-d
    if ch == "A": return g+2*s, g+2*s, g
    if ch == "B": return g+s-d, 2*g+2*s, 0
    if ch == "C": return 5*g+9*s-d, 6*g+12*s, 2*g
    if ch == "D": return 5*g+6*s-2*d, 9*g+12*s, g
    raise ValueError(ch)


def apply(x, w):
    for ch in w:
        x = move(x, ch)
    return x


def ordinal(*parts):
    out = []
    n = sum(map(len, parts))
    for R in parts:
        i = len(out)
        high = ((1 << n)-1) ^ ((1 << (i+len(R)))-1)
        out.extend((r << i) | high for r in R)
    return tuple(out)


def apex(T, marks):
    n = len(T)
    return tuple(r | ((not (marks >> i & 1)) << n)
                 for i, r in enumerate(T)) + (marks,)


def literal(T, w):
    u = len(T)-1
    for ch in w:
        R, marks = {"E": ((), 0), "A": ((0,), 0), "B": ((0,), 1),
                    "C": ((2,4,1), 1), "D": ((2,4,1), 3)}[ch]
        T = apex(ordinal(T, R), (1 << u) | (marks << len(T)))
        u = len(T)-1
    return T


def read(T, marks):
    n = len(T)
    full = (1 << n)-1

    @lru_cache(None)
    def start(m, i):
        rest = m ^ (1 << i)
        if not rest:
            return 1
        return sum(start(rest, j) for j in range(n)
                   if rest >> j & T[i] >> j & 1)

    @lru_cache(None)
    def H(m):
        return sum(start(m, i) for i in range(n) if m >> i & 1) if m else 1

    S = d = Q = 0
    for m in range(full+1):
        g = H(full^m)
        Q += H(m)*g
        t = sum(start(m, i) for i in range(n) if m >> i & marks >> i & 1)*g
        S += t
        d += t if m.bit_count()%2 else -t
    return H(full), S, d, Q


def strong(T):
    for i in range(len(T)):
        seen = 1 << i
        while True:
            new = seen
            for j, row in enumerate(T):
                if seen >> j & 1:
                    new |= row
            if new == seen:
                break
            seen = new
        if seen != (1 << len(T))-1:
            return False
    return True


def tournaments(n):
    edges = [(i,j) for i in range(n) for j in range(i)]
    for bits in product((0,1), repeat=len(edges)):
        rows = [0]*n
        for (i,j), b in zip(edges,bits):
            if b: rows[i] |= 1 << j
            else: rows[j] |= 1 << i
        yield tuple(rows)


def read_charge():
    checks = 0
    for k in range(4):
        for D in tournaments(k):
            for marks in range(1<<k):
                C = apex(D, marks)
                if not strong(C):
                    continue
                c, *_ = read(C, 1<<k)
                a, SD, dD, _ = read(D, marks)
                for L in [(), (0,), (2,4,1)]:
                    for R in [(), (0,), (2,4,1), ordinal((0,), (2,4,1), (0,))]:
                        l, *_ = read(L, 0)
                        r, _, _, qr = read(R, 0)
                        U = ordinal(L,D,R)
                        Z = (marks << len(L)) | (((1<<len(R))-1) << (len(L)+k))
                        h,S,d,_ = read(U,Z)
                        N, *_ = read(apex(U,Z),0)
                        check((h,N,S+h,d) == (l*a*r,l*c*r,l*qr*(SD+a), h if R else l*dD))
                        checks += 1
    print("fixed-baseline literal checks", checks)


def reach(alphabet, X):
    seed = (1,1,1)
    witnessed = {seed: ""}
    queue = deque([seed])
    counts = {}
    while queue:
        x = queue.popleft()
        w = witnessed[x]
        counts.setdefault(x[0],w)
        for ch in alphabet:
            y = move(x,ch)
            check(y[0]>=x[0])
            if y[0]<=X and y not in witnessed:
                witnessed[y] = w+ch
                queue.append(y)
    return counts


def main():
    words = ["".join(w) for w in product("EAB",repeat=4)
             if all(bad not in "".join(w) for bad in ("BE","AA","BB"))]
    check(len(words)==24)
    sums = []
    for s,d in product((1,4),(0,1)):
        x=(1,s,d)
        gains = [Fraction(sum(apply(x,w)[:2]),1+s) for w in words]
        z=sum(1/g for g in gains)
        check(z < Fraction(7,8))
        sums.append(round(float(z),9))
        check(apply(x,"BE")==apply(x,"EA"))
        g=1
        check(apply(x,"AA")==tuple((g+2*s)*k for k in (3,3,1)))
        check(apply(x,"BB")==tuple((3*g+3*s-d)*k for k in (1,2,0)))
    check(256*7**50 < 8**50)
    print("normal-form reciprocal block bounds", sums)
    seed=(0,)
    for n in range(4):
        for w in product("EABCD",repeat=n):
            T=literal(seed,w)
            check(read(T,1<<(len(T)-1))[:3]==apply((1,1,1),w))
    for w,g in [("C",13),("EC",23),("EAC",69)]:
        T=literal(seed,w)
        check(read(T,1<<(len(T)-1))[0]==g)
        print(w,g,"order",len(T),"strong",strong(T))
    old=reach("EAB",1000)
    new=reach("EABCD",1000)
    check((len(old),len(new))==(141,267))
    check(all(g not in old and g in new for g in (13,23,69)))
    check(not {11,19,33}&new.keys())
    print("count palettes through 1000", len(old),len(new))
    read_charge()
    xz = (3,1)   # exposed singleton followed by A, state (3,3,1)
    macros = ("EA","BA","DA")
    m = {"EA":(5,2,2),"BA":(9,2,4),"DA":(51,4,21)}
    count = 0
    for i in range(6):
        for w in product(macros, repeat=i):
            g,_,y = apply((3,3,1),"".join(w))
            target = (g,(g-y)//2)
            phases = []
            while target != xz:
                X,Z = target
                val = (Z & -Z).bit_length()-1
                j = {1:"EA",2:"BA",0:"DA"}[val]
                a,b,c = m[j]
                check(Z%c==0)
                new_x = Z//c
                check(new_x%2==1 and (X-a*new_x)%b==0)
                target = (new_x,(X-a*new_x)//b)
                check(target[0]<X and target[1]>0)
                phases.append(j)
            check(tuple(reversed(phases))==w)
            count += 1
    print("exact two-coordinate macro inverse controls", count)
    for T,expected in [((6,12,24,17,3),{41,43}),((6,12,16,5,11),{33,43})]:
        check(strong(T) and read(T,16)[:3]==(15,21,3))
        check({read(apex(T,16|(1<<j)),0)[0] for j in range(4)}==expected)
    print("two-contact split of equal exposed states", [41,43], [33,43])



if __name__=="__main__":
    main()
