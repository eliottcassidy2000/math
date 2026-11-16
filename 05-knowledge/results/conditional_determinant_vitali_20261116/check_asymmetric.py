from itertools import product

PATCH='0001120'
def rev(v,n):
    return int(f'{v:0{n}b}'[::-1],2)
def rc(v,n):
    return (1<<n)-1-rev(v,n)
def addr(w,row,back=False):
    w=w[::-1] if back else w
    return int(''.join(str(((int(z)>>row)&1)^back) for z in w),2)
def gap(w):
    a,b=zip(*(divmod(int(z),2) for z in w))
    def f(j,i):
        return (1-a[i])*b[j]-(1-b[i])*a[j]
    return (b[0]-a[0]+b[-1]-a[-1]+sum(f(i+1,i) for i in range(len(w)-1))
            +sum((1<<(i-j-1))*f(j,i) for i in range(len(w)) for j in range(i)))
def choices(a):
    options=[]
    for x,y in product(range(4),repeat=2):
        u,v=a+5*x,a+5*y
        if u//16!=v//16:continue
        u0=-10*rev(u%16,4)+8*rev(4*x,4)
        v0=-10*rev(v%16,4)+8*rev(4*y,4)
        options.append((abs(u0-v0),x,y,u//16))
    return max(options)
assert min(choices(i)[0] for i in range(17))==54
A,B,C,D=(addr(PATCH,r,c) for r,c in [(1,False),(0,False),(1,True),(0,True)])
assert B-A==10 and D-C==8 and gap(PATCH)==8
def plan(d,L):
    assert d%2==0
    s=(-abs(d)//2)%4
    T0=(abs(d)+8+10*s)//8
    opt=[]
    c=0
    for j in range(L):
        a=((T0>>(4*j))&15)+c
        g,x,y,c=choices(a)
        opt.append((x,y))
    return T0,s,opt
def main():
    for d in list(range(-120,121,2))+[0,2,10024,2**38,2**80-2]:
        for L in range(1,5):
            T0,S0,opts=plan(d,L)
            ds=set();lengths=set()
            for es in product(*opts):
                t=3*16**L+sum(e*16**i for i,e in enumerate(es))
                T,S=5*t+T0,4*t+S0
                l=T.bit_length();r=4*(L+1)
                u=''.join('30'[b=='0'] for b in f'{T:0{l}b}')
                v=''.join('30'[b=='0'] for b in f'{rc(S,r):0{r}b}')
                w=u+PATCH+v
                la,lb,ra,rb=(addr(w,j,k) for j,k in [(1,False),(0,False),(1,True),(0,True)])
                d0=lb*ra-la*rb
                assert gap(w)==-abs(d)
                assert d0==(2**(l+r)*(128*(10*S-8*T)+10*C-8*A)
                             +10*rc(T,l)*2**r-8*rc(S,r)*2**l)
                assert w[0]=='3' and w[-2:]=='00'
                ds.add(d0);lengths.add(len(w))
            assert len(ds)==2**L and len(lengths)==1
    print('All 17 carry inputs, distinct determinants, signed gap and strong boundary controls passed.')
if __name__=='__main__':main()
