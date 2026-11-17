from math import gcd

def ck(x,y):
    if x!=y:raise AssertionError((x,y))
def F(u):
    z=1
    while u:
        x,b=divmod(u,2)
        if not b:z+=x+(x&1)
        u=x
    return z
def M(s):
    u=int(s or '0',2);D=1<<len(s);v=int(s[::-1] or '0',2)
    return F(u),u+1,D-v,D
def det(m):
    a,b,c,d=m;return a*d-b*c
def grow(A,root,word):
    n=len(A);s=n+len(word);out=list(A)+[0]*(len(word)+1)
    out[root]|=1<<s;out[s]=((1<<n)-1)^(1<<root)
    for j,x in enumerate(word):
        i=n+j;out[i]=((1<<n)-1)+sum(1<<v for v in range(i+1,s))
        if x=='1':out[s]|=1<<i
        else:out[i]|=1<<s
    return tuple(out),s
def tab(a):
    n=len(a);d=[[0]*n for _ in range(1<<n)];h=[1]+[0]*(len(d)-1)
    for i in range(n):d[1<<i][i]=1
    for s in range(1,len(d)):
        h[s]=sum(d[s])
        for i,c in enumerate(d[s]):
            if not c:continue
            z=a[i]&~s
            while z:
                b=z&-z;z-=b;d[s|b][b.bit_length()-1]+=c
    return h,d
def state(a,r):
    h,t=tab(a);f=len(h)-1
    # full starts at root, independently by the converse
    at=tuple(sum(((x>>i)&1)<<j for j,x in enumerate(a)) for i in range(len(a)))
    h2,u=tab(at)
    return h[-1],sum(t[s][r]*h[f^s] for s in range(f+1)),u[-1][r],t[-1][r]
def iso(a,b,m):
    ck(tuple(sum(((a[i]>>j)&1)<<m[j] for j in range(len(a))) for i in sorted(range(len(a)),key=lambda x:m[x])),b)
def gauge(a,r,A,B):
    X,x=grow(a,r,'0'+A);Y,y=grow(X,x,B)
    U,u=grow(a,r,A);V,v=grow(U,u,B+'1')
    n=len(a);f=len(Y)-1
    iso(Y,V,tuple(i if i<n or i==f else f-1 if i==n else i-1 for i in range(f+1)))
    ck(y,v)
def act(m,h,r):
    a,b,c,d=m;return a*h+b*r,c*h+d*r
def rc(s):return ''.join('1' if b=='0' else '0' for b in reversed(s))
def read_singleton(I,r):
    n=len(I);q=n+1
    T=tuple(sum(((x>>i)&1)<<j for j,x in enumerate(I)) for i in range(n))
    return tuple(x|(1<<n)|((i!=r)<<q) for i,x in enumerate(T))+(1<<q,1<<r)

def norm_controls():
    from math import isqrt
    for w,p in [('10',1669),('100',6709),('1010',25309)]:
        a,b='0'+w+'1','0'+rc(w)+'1'
        A,B,C,D=M(a);ck(M(b),(A,C,B,D))
        h,R=act(M(a),1,2);f,g=act(M(b),h,R)
        ck(f+2*g,h*h+R*R);ck(f+2*g,p)
        ck(all(p%d for d in range(2,isqrt(p)+1)),True)
        if w=='10':
            X,r=grow((2,0),1,a);ck(state(X,r),(h,R,0,1))
            Y,r=grow(X,r,b);ck(state(Y,r),(f,g,0,1))
            ck(tab(read_singleton(Y,r))[0][-1],p)
    for k in range(1,40):
        a=2**k;w='1'*k+'0'
        h,R=act(M('0'+w+'1'),1,2)
        ck((h,R),(9*a-3,18*a+2));ck(gcd(h,R),1)
    print('39 primitive norm formulas; three primes and one literal final-host DP verified.')

def main():
    ctr=0
    for k in range(1,14):
        vals=[]
        for u in range(1<<k):
            s=format(u,f'0{k}b')
            w=s.lstrip('0').rstrip('1')
            g=len(s)-len(s.lstrip('0'))
            if not w:
                ck(det(M(s)),0);continue
            t=len(s)-g-len(w);A,B,C,D=M(w)
            ck(M(s),(A,B<<t,C<<g,D<<(g+t)))
            a=w.count('1');b=w.count('0')
            p=det(M(w))
            ck(p%2,1)
            ck(p>=2**a+2**b-1,True)
            ck(p==2**a+2**b-1,w=='1'*a+'0'*b)
            if s==w:vals.append(p)
            ctr+=1
        if k>=2:ck(min(vals),2**(k//2)+2**((k+1)//2)-1)
    for a,b in [('010','01010100'),('1','0'),('0010',''),('10','11')]:
        gauge((2,0),1,a,b)
    A,r=grow((2,0),1,'0110');ck(state(A,r),(19,42,0,1))
    B,s=grow(A,r,'010');ck(state(B,s),(183,450,0,1))
    ck((gcd(19,42),gcd(183,450)),(1,3))
    norm_controls()
    print(ctr,'active determinant/padding checks; literal rooted gauge and primitive loss verified.')
if __name__=='__main__':main()
