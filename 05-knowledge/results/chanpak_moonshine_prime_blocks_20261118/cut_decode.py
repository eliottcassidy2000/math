"""An exact outside-in decoder for single-contact collisions of rooted words."""
from functools import lru_cache
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE.parent/x) for x in ('reflected_words','ans_inverse','ratio_entry','exact_carry','joint_growth','half_reflection','figurate_overlap','tape_carry')]
from reflected_return import debit
from ans_replay import run,bits_lsb
from fixtures import state,forward
from ratio_entries import BASE
from reflected_return import SHORT
from exact_product_carry import even_pair,make

def low_bits(D,d,p,T,top):
    """Return one p-bit x with specified top bit and D*x+d*rev_p(x)=T."""
    if T<0:return None
    if not p:return 0 if T==0 else None
    cap=D+d
    @lru_cache(None)
    def scan(i,u,v):
        if not (0<=u<=cap and 0<=v<=cap):return None
        j=p-1-i
        if i>j:return 0 if u==v else None
        if i==j:
            return next((x<<i for x in ((top,) if i==0 else (0,1))
                         if u+cap*x-(T>>i&1)==2*v),None)
        for x in (0,1):
            for y in ((top,) if i==0 else (0,1)):
                r=u+D*x+d*y-(T>>i&1)
                if r&1:continue
                z=scan(i+1,r//2,2*v-D*y-d*x+(T>>j&1))
                if z is not None:return z|(x<<i)|(y<<j)
        return None
    return scan(0,0,T>>p)

def cuts(h,e,t,s,d,p):
    """Raw target is coincidence; supports differ only by the extra p in the lower-t reader."""
    D=s-t
    assert D>0 and p>=0
    target=(t<<p)-(e if p==0 else 0)
    coef=(D<<(p+1))+d
    for y in range(target//coef+1):
        rem=target-coef*y-d*(y&1)
        for b in ((0,1) if p else (0,)):
            x=low_bits(D,d,p,rem-d*b,b)
            if x is None:continue
            W1=x+(y<<(p+1));W0=W1+(1<<p)
            H0=h+t*W0-e*(W0&1)-d*debit(W0)
            H1=h+s*W1-e*(W1&1)-d*debit(W1)
            assert H0==H1
            yield W0,W1,H0

def revk(y,k):
    return int(f'{y:0{k}b}'[::-1],2) if k else 0
def anti_height(k,y):
    L=1<<k;assert 0<=y<L
    if not k:return 5
    return 8*L**2-6*L+3+6*revk(y,k)+y*y-(6*L-2)*y-2*debit(y)-(y&1)
def address(k,y):
    L=1<<k
    return L**2+L*(L-1-y)+revk(y,k)
def select(N):
    p=N.bit_length()-1
    if p%2==0:p-=1
    if p<1:return None
    k=(p-1)//2
    lo,hi=0,1<<k
    while lo<hi:
        i=(lo+hi)//2
        if anti_height(k,i)>N:lo=i+1
        else:hi=i
    if lo==(1<<k) or anti_height(k,lo)!=N:return None
    W=address(k,lo)
    return k,lo,W+(1<<p),W
def actual():
    A=(6,4,0)
    assert state(A,1)==(1,0,2,1) and state(A,0)==(1,0,4,1)
    for k in range(10):
        L=1<<k
        h=[anti_height(k,i) for i in range(L)]
        assert all(a-b>=2*L for a,b in zip(h,h[1:]))
        for i,N in enumerate(h):
            W=address(k,i);p=2*k+1
            assert 1+4*W-debit(W)==N
            assert ((W<<1)+revk(W,p)+1)==1<<(p+1)
            assert select(N)==(k,i,W+(1<<p),W)
            if k<=3:
                for r,Z in ((1,W+(1<<p)),(0,W)):
                    assert state(forward(A,r,bits_lsb(Z,Z.bit_length())),r)[0]==N
        print('strict monotone exact cut, width',2*k+1,'size',len(h),'extremes',h[-1],h[0])
    assert select(53) is None
    T=((126,124,88,112,96,68,0),3);U=((0,1,35,7,15,27,63),3)
    assert state(*T)==(5,0,18,3) and state(*U)==(5,0,24,3)
    for (A,r),w in ((T,3),(U,2)):
        assert state(forward(A,r,bits_lsb(w,2)),r)[0]==53
    print('other roots give single-contact witness at prime 53')

def main():
    actual()
    for p in range(1,9):
        for D,d in ((2,1),(186,3),(120,1),(1,3)):
            for w in range(1<<p):
                R=int(f'{w:0{p}b}'[::-1],2);T=D*w+d*R
                z=low_bits(D,d,p,T,w>>(p-1))
                assert z is not None and z>>(p-1)==w>>(p-1)
                assert D*z+d*int(f'{z:0{p}b}'[::-1],2)==T
    for label,par in (('863',(863,3,1023,1209,3)),('1063',(1063,1,1847,1967,1)),('29',(29,0,66,78,3))):
        for p in range(41):
            result=list(cuts(*par,p))
            if result:
                print(label,p,result[:2])
                break
        else:print(label,'no decoded coincidence through p=40')
if __name__=='__main__':main()
