"""Independent literal certification of a four-count window with free u."""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]
                       / 'chanpak_phase_windows'))
from replay import attach, paths, profile, g

BANK = [
    ((2,8,3,21,7),4,2,3,2),
    ((2,8,3,21,7,31,63),77,1,1,1),
    ((2,8,3,21,7),22,1,1,5),
    ((2,8,3,21,7,31),19,2,3,1),
]
OFFSETS = (752,754,756,758)

def reverse(x,k):
    return int(f'{x:0{k}b}'[::-1], 2)

def params(spec):
    rows,U,k,w,b=spec
    h,f,P,Q,S,T,K=profile(attach(rows,U,0,0,0,0))
    M=1<<k
    D=(f*b+Q)*M+h*(M-1-reverse(w,k))
    J=h*(1-(w>>(k-1)))
    C=K+P*b+Q*w+f*w*b+h*(g(w)+g(b)-3)+S*(b&1)+T*(w&1)
    return D,h,J,C

def graph(u,spec):
    rows,U,k,w,b=spec
    x,y=(u<<k)+w,b
    return attach(rows,U,x.bit_length(),y.bit_length(),x,y)

def actual(u,spec):
    R=graph(u,spec)
    return sum(paths(R)[-1])

def prove_target(N):
    for c,spec in zip(OFFSETS,BANK):
        if (N-c)%2==0: continue
        D, h = 908, 9
        lo=max(1,(N-c-h+D+h-1)//(D+h))
        hi=(N-c-h)//D
        if hi<lo: continue
        M=2
        while D*M<=h*(hi//M+1):M*=2
        for r in range(min(M,hi+1)):
            a=max(0,(lo-r+M-1)//M)
            b=(hi-r)//M
            while a<=b:
                t=(a+b)//2
                u=M*t+r
                value=D*u+h*g(u)+c
                if value==N:
                    yield u,c,graph(u,spec)
                    break
                if value>N: b=t-1
                else: a=t+1

if __name__=='__main__':
    for offset,spec in zip(OFFSETS,BANK):
        assert params(spec)==(908,9,0,offset)
        for u in (1,2,3,4,6,9,15,16,22,33):
            N=908*u+9*g(u)+offset
            assert actual(u,spec)==N
            assert (u,offset,graph(u,spec)) in set(prove_target(N))
    assert not list(prove_target(2027)), 'hostile must miss the whole bank'
    if len(sys.argv)>1:
        N=int(sys.argv[1])
        for u,c,R in prove_target(N):
            print('WITNESS',N,u,c,R)
    print('PASS symbolic coefficient signatures, 40 literal counts, inverse.')
