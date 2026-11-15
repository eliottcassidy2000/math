from heisenberg import block, mm, inv, pp, filt
from math import gcd
def comm(a,b,m):
    return mm(mm(mm(a,b,m),inv(a,m),m),inv(b,m),m)
def setup(m,L):
    half=pow(2,-1,m)
    ts=(2,3*half%m,half)
    X,Y=block('AA'+'B'*(L-2),ts,m),block('A'+'B'*(L-1),ts,m)
    V=mm(X,inv(Y,m),m)
    C=[V]
    for _ in range(3): C.append(mm(mm(X,C[-1],m),inv(X,m),m))
    lam=[pow(t,e*L,m) for t in ts[:2] for e in (-1,1)]
    E=[filt(C,lam,j,m) for j in range(2 if m%5==0 else 4)]
    Z=comm(E[0],E[1],m)
    if m%5:return X,Z,comm(E[2],E[3],m)
    return X,Z,pp(X,4,m)
def read(u,m):
    return (1+u[-1]+u[0][-1])%m,(2*u[0][1]+pow(2,-1,m)*u[0][1]*(u[1][-1]+u[2][-1]))%m
def main():
    for m in (7,13,19,37,49,91):
        L=1+sum(gcd(i,m)==1 for i in range(m))
        x,z,t=setup(m,L)
        if not (gcd(z[0][-1],m)==gcd(t[1][-1],m)==1 and not t[0][-1] and not t[2][-1]):raise AssertionError(m)
        print('long-tape centers',m,L)
    for m in (5,25,125):
        x,z,d=setup(m,5)
        st={read(x,m):(x,0,0)}
        r=1
        while r<m:
            new={}
            for u,i,j in st.values():
                for a in range(5):
                    for b in range(5):
                        v=mm(u,mm(pp(z,a*r,m),pp(d,b*r,m),m),m)
                        h,q=read(v,m)
                        h0,q0=read(u,m)
                        if h%r!=h0%r or q%r!=q0%r:raise AssertionError(m)
                        new[(h%(r*5),q%(r*5))]=(v,i+a*r,j+b*r)
            r*=5
            if len(new)!=r*r:raise AssertionError((m,r,len(new)))
            st=new
        print('5-primary profile plane',m,len(st))
if __name__=='__main__':main()
