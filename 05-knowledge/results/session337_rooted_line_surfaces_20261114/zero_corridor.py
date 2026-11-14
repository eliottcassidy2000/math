from math import isqrt
from terrace import pathcount
from itertools import product

def rows():
    earlier=[(1,1,1),(1,2,2),(3,3,4)]
    for row in earlier:yield row
    while True:
        row=tuple(map(sum,zip(*earlier)));yield row
        earlier=earlier[1:]+[row]
def exact(N):
    if N<=0:return
    for q,(e,f,g) in enumerate(rows()):
        if e>N:break
        for s in range(2,N.bit_length()):
            d,r=divmod(N-e-g*(1<<s),f)
            if r or d<=0:continue
            D=d*d-(4<<s)
            if D<0:continue
            t=isqrt(D)
            if t*t!=D or (d-t)%2:continue
            x,y=(d-t)//2,(d+t)//2
            if all(a>=2 and a&(a-1)==0 for a in (x,y)):
                yield x.bit_length()-1,q,y.bit_length()-1
if __name__=='__main__':
    for q,(e,f,g) in enumerate(rows()):
        if q==13:break
        for a,b in product(range(1,5),repeat=2):
            assert pathcount((a,)+(0,)*q+(b,))==e+f*(2**a+2**b)+g*2**(a+b)
    assert list(exact(573))==[(2,3,4)]
    assert not list(exact(385))
    assert list(exact(53))==[(1,1,3),(1,2,2)]
    print('573:',list(exact(573)),'385:',list(exact(385)),'53:',list(exact(53)))
