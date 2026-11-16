from itertools import product
from math import isqrt

def count(w):
    h,x,p=1,0,0
    for b in w:
        if b=='0': h+=x+p
        p=int(b);x=2*x+p
    return h

# Rows: 3k+1 and 3k+2 both at 2k+1; 3k+3 at 2k+2.
def row(n):
    q,r=divmod(n-1,3)
    return 2*q+1+(r==2)

def expected_low(h):
    out=set()
    for r in range(1,h.bit_length()+1):
        for m,g in ((1,0),(3,2),(5,2),(9,4)):
            z=m<<r
            if h==2*z+1-g:
                out.add((3*z+3-g,3+g))
    return out

def main():
    for j in range(10000):
        for i in range(3):
            n=6*j+2*i+1
            dest=(3*n+1)//2
            assert n==row(dest) and dest%3==2
            assert row(n)==4*j+i+1
            assert row(dest)-row(n)==2*j+i
        m=3*(j+1)
        assert row(2*m)==2*row(m)
        assert row(2*(3*j+1))==4*j+1
        assert row(2*(3*j+2))==4*j+3
    for j in range(20):
        s=3+(1<<j)
        assert row(2*s)-row(2*s-3)==2
        assert row(2*s)-row(s)+1==3+2*((1<<(j+1))//3-(1<<j)//3)
    ks=set(); seen=set()
    for n in range(2,16):
        for w in product('01',repeat=n-2):
            s='1'+''.join(w)+'0'
            h=count(s)
            all_child=set(); same_gap=False
            for j in range(n+1):
                pair=set()
                for b in '01':
                    out=count(s[:j]+b+s[j:])
                    if out>h:
                        d=2*out-3*h
                        assert d>=3 or (h,out)==(3,5)
                        if h>=23 and d<=7:
                            assert (out,d) in expected_low(h)
                    pair.add(out-2*h)
                    all_child.add(out-2*h)
                same_gap|={-1,1}<=pair
            if {-1,1}<=all_child:
                assert same_gap
                seen.add(h)
                k=isqrt(h)
                if k*k!=h: ks.add(h)
    assert ks=={5,19,23,79}
    for n in range(1,101):
        for m,g in ((1,0),(3,2),(5,2),(9,4)):
            z=m<<n
            v=format(z,'b')[::-1].translate(str.maketrans('01','10'))
            p='1'+v
            c='10'+v
            assert count(p)==2*z+1-g and count(c)==3*z+3-g
    for p,j,children in [('100',1,{9,11}),('1000',2,{17,19}),('11010',2,{37,39}),
                         ('11000',3,{45,47}),('1010000',4,{157,159})]:
        assert {count(p[:j]+b+p[j:]) for b in '01'}==children
    for n in range(1,101):
        p='1'+'0'*n+'1000'+'1'*n+'0'
        assert count(p)==(2**(n+2)+3)**2
    print("Tower formulas: checked 30,000 odd starts and 30,000 doubling starts.")
    print("All strong words through length 15: growth, small offsets and shared-gap forks verified.")
    print("Observed fork counts:", sorted(seen))
if __name__=='__main__': main()
