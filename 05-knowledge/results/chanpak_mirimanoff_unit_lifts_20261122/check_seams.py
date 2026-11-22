"""Construct the marked tournaments, reproduce local certificates, and steer gaps."""

import argparse
from math import gcd, isqrt



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


CASES = {
    8: ((2,12,25,129,43,79,31,119),51,(3,5),0,20),
    9: ((6,20,40,19,37,139,63,95),106,(5,6),0,74),
    10: ((14,84,56,146,97,75,141,55),79,(3,6),0,119),
    12: ((4,17,10,19,37,207,31,95),4,(5,6),0,44964),
    # Append the tagged singleton to the deletion before attaching the suffix.
    14: ((68,81,74,115,101,71,0),65,(3,5),0,1191170),
    16: ((4,17,34,7,13,27,63,127),86,(2,3),294,5699),
}
def strong(rows):
    target=(1<<len(rows))-1
    for u in range(len(rows)):
        mask=1<<u
        while True:
            new=mask
            for i,r in enumerate(rows):
                if mask>>i&1:new|=r
            if new==mask:break
            mask=new
        if mask!=target:return False
    return True

def verify():
    for exponent,(D,tags,(u,v),x,y) in CASES.items():
        expected=(1<<(exponent-1))*((1<<exponent)-1)
        reverse=flipped(D,u,v)
        pair=tuple(star_rows(d,tags,x,y) for d in (D,reverse))
        delta=x.bit_length()
        assert pair[1]==flipped(pair[0],u+delta,v+delta)
        assert all(strong(T) for T in pair)
        data=tuple(profile(d,tags) for d in (D,reverse))
        actual=tuple(count_from_profile(d,x,y) for d in data)
        assert actual==(expected-1,expected+1),(exponent,actual)
        if exponent<=11:assert tuple(subpaths(T)[-1] for T in pair)==actual
        print(exponent,len(pair[0]),actual,data)
def receiver(n):
    return tuple(sum(1<<j for j in range(n)
          if j==i-1 or j>=i+2) for i in range(n))

def receiver_replay(N=11):
    h,c,z,k=1,1,-1,1
    Q=[1,2,4]
    D=[0,1,1]
    for n in range(3,N+1):
        Q.append(2*Q[-1]+Q[-2]+2*Q[-3])
        D.append(2*D[-1]+D[-2]+2*D[-3])
    for n in range(1,N+1):
        a,tag=receiver(n),1<<(n-1)
        hh,q,(e,o),(te,to),kk=profile(a,tag)
        assert (h,c,z,k)==(hh,e+o,e-o,kk)
        assert subpaths(a,tag)[-1]==1
        assert q==Q[n] and te+to==q-h-D[n] and te-to==1-h
        assert n<3 or strong(a)
        h,c,z,k=k,h+c,-h-z,k+c-z
    print("receiver recurrence PASS")

if __name__=='__main__':
    verify()
    receiver_replay()

