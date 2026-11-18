"""Independent reversed continuations and their exact scalar / marked inverses."""
from pathlib import Path
from itertools import combinations
from math import isqrt
import sys
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE.parent/x) for x in ('ans_inverse','exact_carry','half_reflection','joint_growth','figurate_overlap','tape_carry')]
from reflection_search import find
from ans_replay import run,bits_lsb
from fixtures import state,forward
from exact_product_carry import even_pair,make

def reverse(w):
    return int(format(w,'b')[::-1],2)
def debit(w):
    I=[i for i in range(w.bit_length()) if w>>i&1]
    return sum((1<<(j-i-1))+int(j==i+1) for i,j in combinations(I,2))

class Return:
    def __init__(self,a,b,w,d):
        self.a,self.b,self.w,self.d=a,b,w,d
        self.k=w.bit_length()
        self.R=reverse(w)
        assert w&1 and a*w==b*self.R and w%b==0
        self.v=w//b;self.K=debit(w);self.C=a*w
        assert self.K==debit(self.R)
    def macro(self,H,u,e,j):
        return H+self.C*(u<<j)-self.d*self.K-(e if j==0 else 0), (u<<(self.k+j))-self.d*self.v,self.d
    def bitwords(self,j):
        return (bits_lsb(z<<j,self.k+j) for z in (self.w,self.R))
    def inv(self,H,u,e):
        assert u&1 and e==self.d
        x=u+self.d*self.v
        split=(x&-x).bit_length()-1
        j=split-self.k
        if j<0:return None
        prev=x>>split
        h=H-self.C*(prev<<j)+self.d*self.K+(self.d if j==0 else 0)
        assert self.macro(h,prev,self.d,j)==(H,u,e)
        return h,prev,self.d,j

SHORT=Return(11,13,13,3)
BASE=(863,93,3)

def shell(T):
    K=143*93*(1<<(T-4))
    return (99*K+99)//100,107*K//100

def select(N):
    """Exact restricted scalar decoder. Returns schedules, not all witnesses."""
    if N==BASE[0]:
        yield ()
        return
    N-=BASE[0]
    T=4
    while shell(T)[1]<N:T+=1
    if N<shell(T)[0]:return
    def scan(left,H,u,e,J):
        if not left:
            if H==N+BASE[0]:yield J
            return
        for i in range(4,left+1):
            if 0<left-i<4:continue
            h,v,f=SHORT.macro(H,u,e,i-4)
            yield from scan(left-i,h,v,f,J+(i-4,))
    yield from scan(T,*BASE,())

def test():
    assert find(11,13)[0]==13 and find(7,19)[0] is None
    A=((0,41,67,5,15,29,59),2)
    B=((0,65,83,7,11,31,57),6)
    assert state(*A)==(29,0,66,3) and state(*B)==(29,0,78,3)
    K=SHORT.macro(29,6,0,0);assert K==BASE
    S=[]
    for (T,r),bits,a in zip((A,B),SHORT.bitwords(0),(11,13)):
        U=forward(T,r,bits)
        assert state(U,r)==(863,3,93*a,3)
        S.append((U,r))
    print('two literal order-11 states at 863',*(state(*U) for U in S))
    assert all(863%i for i in range(2,isqrt(863)+1))
    for J in ((0,), (1,), (0,0),(1,0,2),(0,1,3,8),(1,12,6,0)):
        H,u,e=BASE
        for j in J:
            Q=SHORT.macro(H,u,e,j)
            assert SHORT.inv(*Q)==(H,u,e,j)
            H,u,e=Q
        assert J in set(select(H))
        print('restricted scalar and marked inverse',J,H)
    # Exact product-carry seed, independent bit reflection.
    w,n=find(1847,1967);assert (w,n)==(51214312118786901,1111)
    T=Return(1847,1967,w,1)
    assert reverse(w)==48089900601626541
    A,B=(make((0,),*p) for p in even_pair((1,1,1,0,0,0,0),1))
    r=len(A)-1
    assert state(A,r)==(1063,1,1847,1) and state(B,r)==(1063,1,1967,1)
    H,u,e=1063,1,1
    for j in (0,2,1,0,7):
        h,v,f=T.macro(H,u,e,j)
        assert T.inv(h,v,f)==(H,u,e,j)
        for b,bits in zip((T.a,T.b),T.bitwords(j)):
            assert run((H,e,b*u,1),bits)==(h,f,b*v,1)
        H,u,e=h,v,f
    h,_,_=T.macro(1063,1,1,0)
    for Q,bits in zip((A,B),T.bitwords(0)):
        Z=forward(Q,r,bits)
        assert all((Z[i]>>j&1)+(Z[j]>>i&1)==1 for i,j in combinations(range(len(Z)),2))
    print('56-letter reflection returns 1063 to',h,'in both order-68 tournaments')

if __name__=='__main__':test()
