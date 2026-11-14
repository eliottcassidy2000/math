#!/usr/bin/env python3
from local_11 import J, S0, mat, mm, mv, g
THREE=(
    (319,(2,4,9,3,15,31),46,('OOOZ',)),
    (321,(0,1,3,7,15,31,63),89,('O','OZ')),
    (323,(0,1,3,7,15,31,63),89,('OOOZ',)),
)
LINEAR=(
    (752,(126,124,56,112,96,64,4),19,2),
    (754,(254,252,248,112,224,192,128,8),107,1),
    (756,(62,28,56,48,32,2),45,5),
    (758,(254,252,248,112,224,192,128,8),211,1),
)
def attach(r,m,x,y):
    left=bin(x)[2:]; right=''.join(str(1-int(b)) for b in bin(y)[2:][::-1])
    s=len(left); N=s+len(r)+len(right);F=(1<<N)-1
    C=left+''.join(str(m>>j&1) for j in range(len(r)))+right
    R=[F^((1<<(i+1))-1) for i in range(N)]
    for i,a in enumerate(r):R[s+i]=(a<<s)|(F^((1<<(s+len(r)))-1))
    return tuple(a|((1-int(C[i]))<<N) for i,a in enumerate(R))+(sum(int(b)<<i for i,b in enumerate(C)),)

def dp(r,seed=None):
    n=len(r);F=(1<<n)-1
    e=[[0]*n for _ in range(F+1)]
    for i in range(n) if seed is None else (seed,):e[1<<i][i]=1
    for s,V in enumerate(e):
        for i,a in enumerate(V):
            if a:
                q=r[i]&~s
                while q:
                    b=q&-q;q^=b;e[s|b][b.bit_length()-1]+=a
    return [1]+list(map(sum,e[1:]))

def rooted(r,u):
    h=dp(r);b=dp(r,u);F=len(h)-1
    f=d=0
    for s in range(1,F+1):
        if s>>u&1:
            v=b[s]*h[F^s];f+=v;d+=v if s.bit_count()&1 else -v
    return h[F],f,d

def word(r,u,ws):
    for w in ws:
        for ch in 'J'+w:
            N=len(r)
            b=0 if ch=='O' else 1<<u
            r=tuple(a|((1-(b>>i&1))<<N) for i,a in enumerate(r))+(b,)
            if ch=='J':u=N
    return r,u

def main():
    assert mm(mat('OZ'),mat('O'))[:3]==(9,9,-3)
    assert mat('OOOZ')[:3]==(9,9,-1)
    for d,a,U,ws in THREE:
        for x,y in ((2,2),(2,4),(4,2)):
            r=attach(a,U,x,y)
            s=rooted(r,len(r)-1)
            o=s
            for w in ws:o=mv(mat(w),o)
            assert o[0]==3456*x*y+459*x+2430*y+d
            if x==y==2:
                R,u=word(r,len(r)-1,ws)
                assert dp(R)[-1]==o[0]
    for u in (1,2,3,5,9):
        for i,(delta,r,U,b) in enumerate(LINEAR):
            T=attach(r,U,u,b)
            h=908*u+9*g(u)+delta
            z=1344*u+1099+6*(i==3)
            s=rooted(T,len(T)-1)
            assert s==(h,z,9)
            N=2252*u+9*g(u)+delta+1090+6*(i==3)
            R,_=word(T,len(T)-1,('',))
            assert dp(R)[-1]==mv(J,s)[0]==N
    print('PASS: literal changed-reader windows and all three transported free-height offsets')
if __name__=='__main__':main()
