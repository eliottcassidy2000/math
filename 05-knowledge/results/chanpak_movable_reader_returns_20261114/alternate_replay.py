#!/usr/bin/env python3
"""Literal and target-first tests for a reciprocal source/reader return."""
from functools import cache

A = (42, 124, 121, 112, 97, 64, 1)
B = (40, 109, 121, 112, 99, 64, 1)
C = (112, 97, 99, 39, 110, 64, 8)
D = (248, 97, 227, 230, 238, 128, 160, 2)

def revcomp(w):
    return ''.join('1' if x == '0' else '0' for x in reversed(w))

def h_word(x):
    h = 1
    while x:
        q, e = divmod(x, 2)
        if e == 0:
            h += q + q % 2
        x = q
    return h

def profile(rows, t, full_profile=False):
    n = len(rows)
    full = (1 << n) - 1
    for i, row in enumerate(rows):
        assert row & (1 << i) == 0 and row <= full
        for j in range(i):
            assert ((row >> j) & 1) != ((rows[j] >> i) & 1)
    @cache
    def start(i, mask):
        mask ^= 1 << i
        if not mask:
            return 1
        available = mask & rows[i]
        result = 0
        while available:
            bit = available & -available
            available ^= bit
            result += start(bit.bit_length() - 1, mask)
        return result
    @cache
    def count(mask):
        if not mask:
            return 1
        available = mask
        result = 0
        while available:
            bit = available & -available
            available ^= bit
            result += start(bit.bit_length() - 1, mask)
        return result
    H = count(full)
    ret = H
    if full_profile:
        back = tuple(full ^ (1 << i) ^ r for i, r in enumerate(rows))
        ending = profile_from_root(back, t)
        P = sum(start(t,m)*count(full ^ m) for m in range(1,full+1) if m>>t&1)
        ret = (H, P, ending)
    start.cache_clear()
    count.cache_clear()
    return ret

def profile_from_root(rows, t):
    @cache
    def start(i, mask):
        if not mask:
            return 1
        return sum(start(j,mask^(1<<j)) for j in range(len(rows)) if rows[i]&mask&(1<<j))
    return start(t, ((1 << len(rows))-1) ^ (1<<t))

def strong(rows):
    for i in range(len(rows)):
        reach = 1 << i
        for _ in rows:
            for j, r in enumerate(rows):
                if reach >> j & 1:
                    reach |= r
        if reach != (1 << len(rows)) - 1:
            return False
    return True

def word_rows(w):
    m=len(w)
    full=(1<<m)-1
    return tuple((full ^ ((1<<(i+1))-1)) | ((x=='0')<<m)
                 for i,x in enumerate(w)) + (sum((x=='1')<<i for i,x in enumerate(w)),)

def glue(source, w):
    n=len(source)
    rhs=word_rows(w)
    m=len(rhs)
    full=(1<<m)-1
    apex=n+m-1
    L=[r | (full<<n) for r in source]
    L[-1] ^= 1 << apex
    R=[(full^(1<<i)^r)<<n for i,r in enumerate(rhs)]
    R[-1] |= 1<<(n-1)
    return tuple(L+R)

def decode_word(M, reciprocal=True):
    """Exact inverse with final 00; reciprocal mode requires initial 11."""
    if M < 1 or not M%2:
        return None
    for k in range(4 if reciprocal else 3, (M-1).bit_length()+2):
        a=k//2 if k>=4 else 1
        b=k-a
        # At width 3 the only possible one-sided word is 100.
        if k==3:
            if M==5:
                return "100"
            continue
        a,b=b,a
        ymin, ymax = 3<<(b-2), 1<<b
        xmin, xmax = (3<<(a-2) if reciprocal else 1<<(a-1)), 1<<a
        for x in range(xmin,min(xmax,(M-1)//ymin+1)):
            hx=h_word(x)
            def score(y):
                return x*y + hx + h_word(y) - 1 + (x&y&1)
            lo,hi=ymin,ymax
            while lo<hi:
                mid=(lo+hi)//2
                if score(mid)<M:
                    lo=mid+1
                else:
                    hi=mid
            if lo<ymax and score(lo)==M:
                return format(x,'b')+revcomp(format(lo,'b'))
    return None

def decode(N):
    """Return two source/reader pairs or reject this exact framed construction."""
    if N%2==0 or (N+9)%31:
        return None
    W=decode_word((N+9)//31-1)
    if W is None:
        return None
    return (A,W[:-1]), (B,revcomp(W)[:-1])

if __name__=='__main__':
    from math import isqrt
    def prime(p):
        return p>1 and all(p%d for d in range(2,isqrt(p)+1))
    for R in (A,B):
        assert strong(R)
        assert profile(R,6,True)==(31,31,9)
    assert sorted(v.bit_count() for v in A)!=sorted(v.bit_count() for v in B)
    words={h_word(int('11'+format(y,f'0{(k-4)}b')+'00',2))
           for k in range(5,11) for y in range(1<<(k-4))}
    words.add(h_word(12))
    for m in range(1,257,2):
        assert bool(decode_word(m))==(m in words)
    for N,W,p in ((2719,'1100100',2719),(92743,'110000010000',13249)):
        assert prime(p)
        assert 31*(1+h_word(int(W,2)))-9==N
        pairs=((A,W[:-1]),(B,revcomp(W)[:-1]))
        for L,w in pairs:
            T=glue(L,w)
            assert strong(T)
            assert profile(T,len(T)-1)==N
        Tv=[glue(L,w) for L,w in pairs]
        d=[sorted(t.bit_count() for t in rows) for rows in Tv]
        assert d[0]!=d[1] and d[0]!=sorted(len(Tv[1])-1-v for v in d[1])
        assert decode(N) is not None
        print('literal pair:',N,W,revcomp(W),len(Tv[0]),flush=True)
    assert decode_word(19) is None and decode_word(19, False)=='10100'
    for W,N0,p in (('11001000',11083,1583),('10100',1171,167)):
        assert prime(N0) and prime(p) and 7*p==N0-2
        for source, E, N in ((C,9,N0),(D,11,N0-2)):
            assert strong(source) and profile(source,len(source)-1,True)==(59,59,E)
            M=(N+E)//59-1
            assert h_word(int(W,2))==M and decode_word(M,False) is not None
            T=glue(source,W[:-1])
            assert strong(T) and profile(T,len(T)-1)==N
    print('profiles, degree separation, framed inverse and eight literal outputs passed')
