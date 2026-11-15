from replay_recontact_energy import st, ordinal
from random import Random
def tr(r):
    return tuple((1<<r)-(1<<(i+1)) for i in range(r))
def main():
    R=(16,1,3,7,14)
    G=Random(615)
    for _ in range(120):
        L=G.randrange(1,6)
        b=G.randrange(1<<L)
        bm=int(f'{b:0{L}b}'[::-1],2)
        m=G.randrange(32)
        h,_,S,_=st(ordinal(R,tr(L)),m|(bm<<5))
        if S*pow(h,-1,1<<L)%(1<<L)!=b:raise AssertionError((m,L,b,S))
    print('120 literal terminal-word recoveries modulo the exact tail scale')
    for j,c in [(3,2),(4,4),(5,8)]:
        h=15; S=2*(17*(1<<j)-h); u=(S-c)//h;a=2*u
        L=a.bit_length()
        def pick(t):
            binary=int(f'{t:0{L}b}'[::-1],2)
            b=(1<<L)-1-binary
            return (S<<L)+h*b
        T,U=pick(a),pick(a+1)
        if abs(T-U)!=h*(1<<(L-1)) or max(T,U)>=2*min(T,U):raise AssertionError((j,T,U))
        print('Fold endpoints, height',j,':',T,U)
if __name__=='__main__':main()
