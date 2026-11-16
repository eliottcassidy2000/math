from math import isqrt
def coeff(a,w,S):
    out=[0]*5
    def rec(L,last,r):
        if not L:
            if r:out[r]+=1
            if last!=w and a[last]>>w&1:out[r+1]+=1
            return
        if last!=w and (last<0 or a[last]>>w&1):
            rec(L,w,r+1)
        choices=L if last<0 else L&a[last]
        while choices:
            b=choices&-choices;choices-=b
            rec(L^b,b.bit_length()-1,r)
    rec(S^(1<<w),-1,0)
    return tuple(out[1:])
A=(2,8,11,97,15,87,151,63)
assert [coeff(A,4,S) for S in (255,255^4,255^128)]==[(79,36,0,0),(43,18,0,0),(37,18,0,0)]
for t in (2,10):
    T=2**t
    h,p,q=(36*T+7,18*T+7,18*T+1)
    assert h==p+q-1 and all(x%d for x in (h,p,q) for d in range(2,isqrt(x)+1))
for e in [-1,1]:
    for f in [-1,1]:
        p=(2+(e-1)*(-f-1))%3;q=(2+(-e-1)*(f-1))%3
        assert (p!=0 and q!=0)==(e==f)
        if e==f:assert (p+q-1)%3==0
print('Asymmetric J coefficients, two prime hosts, and exact parity lock of the other template verified.')
