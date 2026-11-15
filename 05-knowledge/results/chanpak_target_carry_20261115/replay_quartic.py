# Load the published predecessor packet (or its local working name).
from pathlib import Path
import sys
base=Path(__file__).resolve().parent.parent
for directory in ("chanpak_absolute_steering_20261115","chanpak_absolute_steering"):
    if (base/directory).is_dir():
        sys.path.insert(0,str(base/directory))
from math import gcd
from random import Random
from replay_pair import pair_compiler, prescribe_pair, rev
from replay_even_pair import compile_even, prescribe
from replay_steering import B,Q,T,S,p,output

def need(x):
    if not x:
        raise AssertionError('quartic shell assertion failed')

def setup(M,q,u):
    m=M//(M&-M)
    a=(M&-M).bit_length()-1
    g=gcd(45,M)
    k=max(a+1,(g-1).bit_length())
    A=1<<k
    K=2*q*u+q-1
    comp=compile_even(M) if a else pair_compiler(M)
    compiler=prescribe if a else prescribe_pair
    D=1<<comp[-1]
    need(q//2%D==0)
    def core(v,d):
        X=q*u+v
        Y=K-X
        x,y=A*X+d,A*Y+d
        C=(A-1-rev(d,k))*K+2*p(d)+(d<A//2)+u*(q-1)+2*p(u)+(u&1)
        return B+Q*x*y+T*x+S*y+45*C
    lo=core(q-1,0)+Q*A*A*q*q//16
    hi=core(q//2,0)-Q*A*A*q*q//16
    def place(N):
        need(lo<=N<=hi and N%2)
        d=(N-B-S*A*K)*pow(T+S,-1,g)%g if g>1 else 0
        H=B+Q*d*(A*K+d)+T*d+S*(A*K+d)+45*((A-1-rev(d,k))*K+2*p(d)+(d<A//2))
        need((N-H)%(2*g)==0)
        z=(N-H)//(2*g)*pow(45//g,-1,M//g)%(M//g) if M//g>1 else 0
        def program(j):
            X0=(q*u)//D+j; Y0=(K+1)//D-1-X0
            tape=compiler(comp,X0,Y0,(p(X0)+p(Y0))//2,0,z)
            X=D*X0+int(tape,2)
            v=X-q*u
            F=output(A*X+d,A*(K-X)+d)
            need(F==core(v,d)+45*(p(v)+p(q-1-v)))
            need((F-N)%(2*M)==0)
            return F
        i,j=(q//2)//D,q//D-1
        U,L=program(i),program(j)
        need(L<=N<=U)
        while j-i>1:
            zmid=(i+j)//2
            Hmid=program(zmid)
            if Hmid>=N: i,U=zmid,Hmid
            else: j,L=zmid,Hmid
        need(L<=N<=U and U-L<=(2*D*(Q*A*A+A*abs(T-S))+45)*q)
        return L,U
    return lo,hi,place

if __name__=='__main__':
    rng=Random(940301)
    n=0
    for M in (1,2,3,4,5,6,9,10,12,15,16,18,32,64,128):
        comp=compile_even(M) if M%2==0 else pair_compiler(M)
        for r in range(max(11,2*comp[-1]+2),max(11,2*comp[-1]+2)+2):
            q=1<<r
            for u in ((q+199)//200,q//100):
                lo,hi,place=setup(M,q,u)
                for _ in range(3):
                    N=(rng.randrange(lo+1,hi-1)|1)
                    place(N)
                    n+=1
    print('PASS: target-chosen common phase, literal tape formulas, quartic modulus brackets:',n)
