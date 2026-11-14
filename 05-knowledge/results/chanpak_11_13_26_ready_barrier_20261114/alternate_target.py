"""Literal two-hub witness at the target unavailable to terminal seven."""
from bilateral import path_count
from math import isqrt

def two_hub(word):
    m=len(word)
    rows=[sum(1<<j for j in range(i+1,m)) for i in range(m)]+[1<<(m+1),0]
    for i,c in enumerate(word):
        for h in (0,1):
            if int(c)>>(1-h)&1:rows[m+h]|=1<<i
            else:rows[i]|=1<<(m+h)
    return tuple(rows)

def strong(rows):
    for i in range(len(rows)):
        R=Q=1<<i
        while Q:
            b=Q&-Q;Q^=b
            z=rows[b.bit_length()-1]&~R
            R|=z;Q|=z
        if R!=(1<<len(rows))-1:return False
    return True

if __name__=='__main__':
    w="1002303022110"
    a=two_hub(w)
    assert len(a)==15 and strong(a)
    assert path_count(a)==180131==7*25733
    assert all(25733%d for d in range(2,isqrt(25733)+1))
    print('literal alternate for 180131:',w,a)
