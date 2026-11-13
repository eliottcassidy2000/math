from collections import Counter
from functools import cache
from near_twins import path_count

def slot_transfer(A,root,t):
    n=len(A); S=(1<<n)-1^(1<<root)
    @cache
    def paths(T,last,k):
        if not T:return int(k==0 or (k==1 and (last<0 or A[last]>>root&1)))
        out=0
        if k and (last<0 or A[last]>>root&1):
            out=paths(T,root,k-1)
        U=T if last<0 else T&A[last]
        while U:
            b=U&-U;U-=b
            out+=paths(T^b,b.bit_length()-1,k)
        return out
    @cache
    def slots(S,t):
        if not t:return ((0,1),) if not S else ()
        V=Counter()
        T=S
        while True:
            for j in range(T.bit_count()+2):
                h=paths(T,-1,j)
                if h:
                    for k,z in slots(S^T,t-1):V[j+k]+=h*z
            if not T:break
            T=(T-1)&S
        return tuple(V.items())
    return dict(sorted((j,z) for j,z in slots(S,t) if j))


def expand(A,L):
    mask=(1<<(L+4))-16
    R=[(a&15)|(mask if a>>4&1 else 0) for a in A[:4]]
    for i in range(4,L+4):
        R.append(A[4]|((1<<(L+4))-(1<<(i+1))))
    return R

A=(22,24,18,5,8); B=(20,25,18,5,8)
got=[slot_transfer(G,4,t) for G in (A,B) for t in (1,2)]
expected=[{1:11,2:6},{1:100,2:98,3:30},{1:13,2:6},{1:104,2:100,3:28}]
if got!=expected: raise ValueError(got)
for l in range(1,7):
    h=(path_count(expand(A,l)),path_count(expand(B,l)))
    if h!=(6*2**l-1,6*2**l+1):
        raise ValueError(h)
print('opposed gate slot rows:',got)
