"""Standalone exact selector controls; pure Python, standard library only."""
from itertools import product
from fractions import Fraction
from math import isqrt

def p(address):
    total = 0
    while address:
        high, bit = divmod(address, 2)
        if bit == 0:
            total += high + (high & 1)
        address = high
    return total


def flipped(rows, i, j):
    out = list(rows)
    out[i] ^= 1 << j
    out[j] ^= 1 << i
    return tuple(out)


def star_rows(core, tags, x, y):
    """Ordinal left flank, marked core, right flank, and last vertex as apex."""
    left = bin(x)[2:] if x else ''
    right = bin(y)[2:][::-1].translate(str.maketrans('01', '10')) if y else ''
    tape = left + ''.join(str((tags >> i) & 1) for i in range(len(core))) + right
    n = len(tape)
    mask = (1 << n) - 1
    rows = [mask ^ ((1 << (i + 1)) - 1) for i in range(n)]
    width = len(core)
    core_mask = ((1 << width) - 1) << len(left)
    for i, row in enumerate(core):
        index = len(left) + i
        rows[index] = (rows[index] & ~core_mask) | (row << len(left))
    apex = 0
    for i, bit in enumerate(tape):
        if bit == '1':
            apex |= 1 << i
        else:
            rows[i] |= 1 << n
    return tuple(rows + [apex])


def subpaths(rows, starts=None):
    n = len(rows)
    size = 1 << n
    end = [[0] * n for _ in range(size)]
    for i in range(n):
        if starts is None or starts >> i & 1:
            end[1 << i][i] = 1
    for subset in range(1, size):
        for i, count in enumerate(end[subset]):
            if count:
                choices = rows[i] & ~subset
                while choices:
                    bit = choices & -choices
                    end[subset | bit][bit.bit_length() - 1] += count
                    choices ^= bit
    return [sum(row) for row in end]


def profile(core, tags):
    n = len(core)
    full = (1 << n) - 1
    paths = subpaths(core)
    paths[0] = 1
    h = paths[-1]
    q = sum(paths[s] * paths[full ^ s] for s in range(full + 1))

    def eo(graph, marked):
        rooted = subpaths(graph, marked)
        return tuple(sum(rooted[s] * paths[full ^ s] for s in range(1, full + 1)
                         if s.bit_count() % 2 == parity) for parity in (0, 1))

    reverse = tuple(row ^ full ^ (1 << i) for i, row in enumerate(core))
    k = subpaths(star_rows(core, tags, 0, 0))[-1]
    return h, q, eo(core, tags), eo(reverse, full ^ tags), k


def count_from_profile(data, x, y):
    h, q, (e, o), (end_e, end_o), k = data
    return (k + h * (p(x) + p(y)) + q * x * y + (end_e + end_o) * x +
            (e + o) * y - (end_e - end_o) * (x & 1) - (e - o) * (y & 1))

def strong(rows):
    full = (1 << len(rows)) - 1
    for root in range(len(rows)):
        seen = 1 << root
        while True:
            new = seen
            for i, row in enumerate(rows):
                if seen >> i & 1:
                    new |= row
            if new == seen:
                break
            seen = new
        if seen != full:
            return False
    return True


def signed_divisors(m):
    for d in range(1,isqrt(abs(m))+1):
        if m%d:continue
        for f in {d,abs(m)//d}:
            yield f
            yield -f


def adjusted(D,ex,ey):
    h,q,(e,o),(ep,op),k=D
    return h,q,e+o,ep+op,k-ex*(ep-op)-ey*(e-o)


def selector(L,U,N,ex,ey):
    """All two-positive-address solutions within a specified parity."""
    h,q,s,r,k=adjusted(L,ex,ey)
    H,Q,S,R,K=adjusted(U,ex,ey)
    a,b,c,t=H*q-h*Q,H*s-h*S,H*r-h*R,(H-h)*N-2*h-H*k+h*K
    d=a*t+b*c
    cap=(N-k)//q
    if cap<1:return set()
    if a and d:
        def candidates():
            for f in signed_divisors(d):
                x,j=divmod(f-b,a)
                y,l=divmod(d//f-c,a)
                if not j and not l:
                    yield x,y
    else:
        def candidates():
            for x in range(1,cap+1):
                den=a*x+b
                rhs=t-c*x
                if den:
                    y,r=divmod(rhs,den)
                    if not r:yield x,y
                elif not rhs:
                    for y in range(1,cap//x+1):yield x,y
    found={(x,y) for x,y in candidates() if x>0 and y>0 and x%2==ex and y%2==ey
             and count_from_profile(L,x,y)==N}
    assert all(count_from_profile(U,x,y)==N+2 for x,y in found)
    return found


def append(v,b,k):return (v<<k)+b


def reverse_width(v,k):return int(format(v,f'0{k}b')[::-1],2)


def gamma(v,k):
    q=1<<k
    return Fraction(p(v)+((v&1) and v<q//2))-Fraction((q-1-reverse_width(v,k))*v,q-1)


def gamma_by_runs(v,k):
    s=list(map(int,format(v,f'0{k}b')))
    cycles=sum(s[i]==1 and s[(i+1)%k]==0 for i in range(k))
    decay=sum(Fraction((1-s[(i-d)%k])*2**(k-d),2*(2**k-1))
              for i in range(k) if s[i] for d in range(1,k+1))
    return cycles-decay


def check():
    lower=(10,4,9,18,7)
    upper=flipped(lower,0,1)
    a,b=profile(lower,2),profile(upper,2)
    assert a==(13,104,(5,9),(27,38),23)
    assert b==(11,100,(9,11),(32,40),29)
    h,q,s,r,k=adjusted(a,0,0)
    H,Q,S,R,K=adjusted(b,0,0)
    N=1221
    A,B,C,T=H*q-h*Q,H*s-h*S,H*r-h*R,(H-h)*N-2*h-H*k+h*K
    assert (A,B,C,T)==(-156,-106,-221,-2344)
    assert (4,2) in selector(a,b,N,0,0)
    actual=tuple(star_rows(z,2,4,2) for z in (lower,upper))
    assert tuple(subpaths(z)[-1] for z in actual)==(1221,1223)
    assert all(strong(z) for z in actual)
    assert sum((u^v).bit_count() for u,v in zip(*actual))==2
    for L,U in [(a,b),(b,a)]:
        for x,y in product(range(1,9),repeat=2):
            h=count_from_profile(L,x,y)
            assert ((x,y) in selector(L,U,h,x%2,y%2)) == (count_from_profile(U,x,y)==h+2)
    num=0
    for k in range(1,10):
        q=1<<k
        for v in range(q):
            g=gamma(v,k)
            assert g==gamma_by_runs(v,k)
            assert (g>0)==(0<v<q-1)
            assert g==gamma((q-1)^reverse_width(v,k),k)
            w=[13]
            for j in range(5):w.append(append(w[-1],v,k))
            # A fixed suffix 101 changes exponential and constant, not the linear drift.
            pc=[p(append(r,5,3)) for r in w]
            # Remove first increment; parity is stationary after the first fixed block.
            slopes=[pc[j+1]-pc[j] for j in range(1,5)]
            assert all(Fraction(q*slopes[j]-slopes[j+1],q-1)==g for j in range(3))
            num+=1
    print('nonpersistent pair 1221/1223; exact selector crosscheck 128 queries; periodic blocks',num)


if __name__=='__main__':check()
