from binary_cycle import create,cycles,literal
from cycle_digit_gap import delete
from math import isqrt
from functools import cache

def contact_rows(word):
    n=len(word)
    T=[sum(1<<j for j in range(i+1,n))|((1<<n) if not int(b) else 0)
       for i,b in enumerate(word)]
    T.append(sum(1<<i for i,b in enumerate(word) if b=='1'))
    return tuple(T)

if __name__=='__main__':
    for K in range(1,45):
        D=bin(K)[2:]
        T,u=create(K,False)
        n=len(D)
        assert cycles(T,u)==K
        L=delete(T,u)
        S=delete(L,n) # index n is a
        # Same literal counts as the two words in their standard
        # reordered b-before-a transitive spine.
        x=literal(S)[-1];y=literal(L)[-1]
        assert x==literal(contact_rows(D+'0'))[-1]
        assert y==literal(contact_rows(D+'00'))[-1]
        assert y-x==2*K and K+1<=x and y<=4*K+1
        assert len(S)==n+2 and len(L)==n+3
        if K&(K-1):assert y<=4*K-1
    print('direct induced path-gap and binary cycle relation checked for 44 targets')
