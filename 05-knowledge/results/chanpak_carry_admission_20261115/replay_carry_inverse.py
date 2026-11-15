from itertools import product

def p(b):
    v=0
    while b:
        z=b//2
        if b%2==0:v+=z+(z%2)
        b=z
    return v

def cv(w):
    return int(w,2) if w else 0

def sig(U,W,v):
    u,Q=U
    h,q=W
    x=Q-u; d=x*q; t=u*h
    k=len(v); r=cv(v[::-1]); b=cv(v)
    return (d+t)*2**k-t*r, d-x*h+2*t+d*b+t*p(b), t

def inv(U,W,lam,c,k):
    u,Q=U; h,q=W; t=u*h; d=(Q-u)*q
    r,rem=divmod((d+t)*2**k-lam,t)
    if rem or r<0 or r>=2**k or (k and not r%2):return None
    v=f'{r:0{k}b}'[::-1] if k else ''
    if sig(U,W,v)[:2]!=(lam,c):return None
    return v

def tt(n):return tuple((1<<n)-(1<<(i+1)) for i in range(n))

def glue(word,W,U):
    k=len(word); n=k+1; m=len(W)
    r=list(tt(k))+[0]
    for i,b in enumerate(word):
        if b=='1': r[k]|=1<<i
        else:r[i]|=1<<k
    wm=((1<<m)-1)<<n; um=((1<<len(U))-1)<<(n+m)
    return tuple(a|wm|(um if i!=k else 0) for i,a in enumerate(r))+tuple((a<<n)|um for a in W)+tuple((a<<(n+m))|(1<<k) for a in U)

def dp(r):
    n=len(r); full=(1<<n)-1
    end=[[0]*n for i in range(full+1)]
    for s in range(1,full+1):
        z=s
        while z:
            bit=z&-z; i=bit.bit_length()-1; t=s^bit
            a=r[i]&t
            if not t: end[s][i]=1
            while a:
                j=a&-a; end[s][i]+=end[t][j.bit_length()-1];a^=j
            z^=bit
    return sum(end[-1])

def main():
    profs=((1,2),(1,8),(3,12),(5,32),(19,256),(345,5640))
    count=0
    for U,W in product(profs,repeat=2):
        for k in range(10):
            for b in (range(2**(k-1),2**k) if k else (0,)):
                v=f'{b:0{k}b}' if k else ''
                lam,c,_=sig(U,W,v)
                if not 1<=c<=lam or not (lam%2 and c%2) or inv(U,W,lam,c,k)!=v:
                    raise AssertionError((U,W,v,lam,c))
                for g in (1,3,8):
                    q=1<<g
                    N=q*lam+c
                    if not (N+q)//(q+1)<=lam<=(N-1)//q:
                        raise AssertionError((N,lam))
                    count+=1
    print('Canonical carried profile/target controls:',count)
    for u in (3,5):
        for l in range(1,4):
            a=u<<l
            for m,W,v,x,want in [(3,tt(1),'',2*a,30*a+7),(1,(2,4,1),'0',a,30*a+9)]:
                word=bin(x)[2:]+v
                got=dp(glue(word,W,tt(m)))
                if got!=want:raise AssertionError((a,word,got,want))
    print('12 independent literal shifted-ray counts')
    check=dp(glue('100',tt(1),tt(1)))
    lam,c,_=sig((1,2),(1,2),'')
    if check!=(4*lam+c): raise AssertionError(check)
    # Prefix 10, suffix 110 separately: c exceeds 2 and the carry is nonzero.
    lam,c,_=sig((1,2),(1,2),'110')
    actual=dp(glue('10'+'110',tt(1),tt(1)))
    if not c>=2 or actual!=2*lam+c or inv((1,2),(1,2),lam,c,3)!='110':
        raise AssertionError((lam,c,actual))
    print('Literal nonzero carry:',(lam,c,c//2,actual))
    lam,c,_=sig((1,2),(1,2),'111111110')
    z=dp(glue('10'+'111111110',tt(1),tt(1)))
    if (lam,c,z)!=(1281,1279,3841):raise AssertionError((lam,c,z))
    print('Literal carry 639:',(lam,c,c//2,z))

if __name__=='__main__':main()
