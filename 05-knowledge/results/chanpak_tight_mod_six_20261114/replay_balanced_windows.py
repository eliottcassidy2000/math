from math import isqrt

def g(x):
    h = 1
    while x:
        x, b = divmod(x, 2)
        if not b:
            h += x + x % 2
    return h

def bilateral(x,y):
    return x*y+g(x)+g(y)-1+(x%2)*(y%2)

def geometry(x,y):
    z = bin(x)[2:] + ''.join(str(1-int(b)) for b in bin(y)[2:][::-1])
    n=len(z)
    return tuple(((1<<n)-(1<<(i+1)))|((int(z[i])^1)<<n) for i in range(n))+(sum(int(b)<<i for i,b in enumerate(z)),)

def count(row):
    n=len(row)
    m=(1<<n)-1
    P=[[0]*n for _ in range(m+1)]
    for i in range(n):
        P[1<<i][i]=1
    for s in range(1,m+1):
        for i in range(n):
            if not P[s][i]:
                continue
            t=(m^s)&row[i]
            while t:
                b=t&-t
                P[s|b][b.bit_length()-1]+=P[s][i]
                t^=b
    return sum(P[m])

for x in range(6):
    for u in range(6):
        C={ (x*v+g(8*u+v)-g(u)+(x%2)*(v%2))%3 for v in range(8) }
        assert (len(C)<3)==(u==[5,2,1,4,3,0][x])
        assert len({bilateral(x,16*u+v)%6 for v in range(16)})==3

for k in range(2,12):
    M=2**k
    for u in [0,1,2,15,8991]:
        B=M*M*u*(u+1)+M*M//4-u+2*g(u)+(u%2)-3
        for q in range(M//2):
            a,b=M//2-1-q,M//2+q
            H=bilateral(M*u+a,M*u+b)
            expect=B+2*g(q)-q*q-2*q-2*(q>>(k-2))+(q%2)
            assert H==expect
            if k>=4 and q<=2:
                assert H==B+2-2*q
                if k<=5 and u<=2:
                    assert count(geometry(M*u+a,M*u+b))==H
print("PASS: 36 nibble states, exact complementary chart and literal small-order DP")

def E(q): return 2*g(q)-q*q-2*q+q%2
for v in range(4000):
    assert [(E(4*v+i)-E(4*v))%3 for i in range(4)] == [
        0, 1, (2-2*(v%2))%3, (1-2*(v%2))%3]
    for j in range(1, min(12, v+2)):
        assert E(2**j-1)==4-4**j and E(2**j)==2-4**j
for k in range(4,9):
    M=2**k
    for u in range(8):
        HH=[bilateral(M*u+M//2-1-q,M*u+M//2+q) for q in range(M//2)]
        assert all(x>y for x,y in zip(HH,HH[1:]))
        assert all(isqrt(H)//M==u for H in HH)
print("PASS: aligned mod-6 and power-index mod-4 windows; exact target inverse")

def select(N):
    if N <= 1 or N % 2 == 0:
        return None
    for k in range(4, (N+1).bit_length()):
        M=2**k
        u=isqrt(N)//M
        B=M*M*u*(u+1)+M*M//4-u+2*g(u)+(u%2)-3
        lo,hi=0,M//2-1-(u==0)
        while lo<=hi:
            q=(lo+hi)//2
            H=B+E(q)-2*(q>>(k-2))
            if H>N: lo=q+1
            elif H<N: hi=q-1
            else:
                x=M*u+M//2-1-q
                y=M*u+M//2+q
                assert x>0 and y>0 and bilateral(x,y)==N
                return u,k,q,geometry(x,y)
    return None

for N in (61,63,65,573,575,577):
    witness=select(N)
    assert witness and count(witness[-1])==N

if __name__=='__main__':
    import json,sys
    if len(sys.argv)>1:
        witness=select(int(sys.argv[1]))
        print(json.dumps(witness))
