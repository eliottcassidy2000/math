from itertools import product

def check(test, record):
    if not test:
        raise AssertionError(record)

def tt(n):
    return tuple((1<<n)-(1<<(i+1)) for i in range(n))

def hostsig(word):
    b, a, t = 0, 1, 0
    for bit in word:
        s = int(bit)
        a += (1-s)*(b+t)
        b = 2*b+s
        t = s
    return b, a

def rev(b, k):
    if not k:
        return 0
    return int(f'{b:0{k}b}'[::-1], 2)

def binary(b, k=None):
    if k == 0:
        return ''
    return bin(b)[2:] if k is None else f'{b:0{k}b}'

def arm(m, h, q, v):
    x, k = (1<<m)-1, len(v)
    b, a = hostsig(v)
    return x*q*(1<<k)+h*((1<<k)-rev(b,k)), x*q*(b+1)+h*(a+1-x)

def glue(word, W, m):
    n,k=len(word), len(W)
    # R: transitive spine then outlet. Then W, then the transitive U.
    v=n
    allw=((1<<k)-1)<<(n+1)
    allu=((1<<m)-1)<<(n+1+k)
    R=tuple(r|allw|allu|((1<<v) if x=='0' else 0) for r,x in zip(tt(n),word))
    R+=(allw+sum(1<<i for i,x in enumerate(word) if x=='1'),)
    return R+tuple((r<<(n+1))|allu for r in W)+tuple((r<<(n+1+k))|(1<<v) for r in tt(m))

def dp(rows, profile=False):
    n=len(rows); zz=1<<n
    H=[[0]*n for _ in range(zz)]
    T=[0]*zz if profile else None
    if profile: T[0]=1
    for mask in range(1,zz):
        f=mask
        while f:
            p=f&-f; f-=p; i=p.bit_length()-1
            old=mask-p
            if old==0: H[mask][i]=1;continue
            q=rows[i]&old
            while q:
                b=q&-q; q-=b
                H[mask][i]+=H[old][b.bit_length()-1]
        if profile:T[mask]=sum(H[mask])
    if profile:return T[-1],sum(T[k]*T[zz-1-k] for k in range(zz))
    return sum(H[-1])

def strong(rows):
    mask=1
    for i in range(len(rows)):
        for j in range(len(rows)):
            if mask>>j&1:mask|=rows[j]
    h=1
    for i in range(len(rows)):
        for j in range(len(rows)):
            if rows[j]&h:h|=1<<j
    return mask==h==(1<<len(rows))-1

triples = [
    (1281,1279,1,1,'111111110',1,(82,108,9,1,78,29,44),(81,1200)),
    (6241,6093,2,4,'1111100',2,(54,188,128,5,12,28,191,57),(73,2056)),
    (16393,15369,1,10,'1110',4,(98,28,121,1,9,90,26),(73,1088)),
]
for D,E,m,s,v,m2,W,hq in triples:
    h,q=hq
    check(dp(W,True)==hq and (strong(W) or D==6241),W)
    check(arm(m,1,1<<s,v)==(D,E),(D,E))
    check(arm(m2,h,q,'')==(D,E+2),(D,E))
    g=1
    # The larger transitive arm in the third pair is independently
    # evaluated by the fast host recurrence at 24 gap lengths below.
    if D!=16393:
        check(dp(glue('1'+'0'*g+v,tt(s),m))==(D<<g)+E,(D,g))
    check(dp(glue('1'+'0'*g,W,m2))==(D<<g)+E+2,(D,g))
    for g in range(1,25):
        for j,hh,qq,u,low in [(m,1,1<<s,v,E),(m2,h,q,'',E+2)]:
            b,A=hostsig('1'+'0'*g+u)
            x=(1<<j)-1
            check(x*qq*(b+1)+hh*(A+1-x)==(D<<g)+low,(D,g,j))
    print('literal pair:',D,E,E+2)

records={}
for m,s,k in product(range(1,10),range(1,10),range(8)):
    for b in range(1<<k):
        D,E=arm(m,1,1<<s,binary(b,k))
        L=D.bit_length()
        m0=next(j for j in range(1,L) if (1<<L)-(1<<(L-j))<D<=(1<<L)-(1<<(L-j-1)))
        k0=L-m0-1
        b0=rev((1<<L)-(1<<(k0))-D,k0)
        check(m0==m and arm(m0,1,2,binary(b0,k0))==(D,E),(m,s,k,b,D,E))
        check(arm(m0,1,2,binary(b0,k0+1))==(2*D,E),(D,E))
        check(D not in records or records[D]==E,(D,E))
        records[D]=E
print('transitive high coefficients checked:',len(records))
for a in range(1,10001):
    p=binary(a)
    A=hostsig(p)[1]
    b1,A1=hostsig(p+'111111110')
    minus=2*(b1+1)+A1
    plus=1200*(a+1)+81*A
    check((minus,plus)==(1279+1280*a+A-1,1281+1200*a+81*(A-1)),a)
    check((abs(plus-minus)==2)==(a>1 and a&(a-1)==0),a)
print('arbitrary-prefix exact gap checked through:',10000)
