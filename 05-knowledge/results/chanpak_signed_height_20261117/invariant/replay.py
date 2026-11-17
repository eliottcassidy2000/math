from math import isqrt,gcd
from random import Random
from itertools import product

def ck(a,b):
    if a!=b:raise AssertionError((a,b))
def F(u):
    c=1
    while u:
        x,b=divmod(u,2)
        if not b:c+=x+(x&1)
        u=x
    return c

def migration(k,u):
    return F(u),u+1,(1<<k)-int(format(u,f'0{k}b')[::-1],2),1<<k

def apply(M,h,R):
    a,b,c,d=M
    return a*h+b*R,c*h+d*R

def edge(j,i,limit=None):
    if not (j>=2 and 2<=i<=j+1):return None
    q=1<<i;Q=1<<j;a=q+1;b=Q+1
    t=1+Q//2-q//4;k=2*i*t
    if limit is not None and k>limit:return None
    D=1<<k;T=(D-1)//a
    m=((b*D)+(D//2+1+T))//a
    return m,T,k

def val(n):
    if n<1:return -1
    return (n&-n).bit_length()-1

def inverse(last_i,h):
    path=[last_i]
    while h!=1:
        if h<1:return None
        v=val(h-1)
        i=(v+1)//2
        if i<2 or path[-1]!=i:return None
        q=1<<i;a=q+1
        C=[]
        if v&1:
            js= [(i-1,False)] if i>=3 else []
        else:
            Y=(a*a*h-a)//q
            e=val(Y-1)
            js=[((e+1)//2,False)] if (e+1)//2>=i else []
            w=e+i+1
            if w%(2*i)==0:
                j0=val(w//(2*i)-1+q//4)+1
                if j0>=i and w//(2*i)-1+q//4==(1<<(j0-1)):js.append((j0,True))
        for j,seed in js:
            data=edge(j,i,h.bit_length()+1)
            if data is None:continue
            m,T,_=data
            u,r=divmod(h+T,m)
            if r or u<1 or u>=h or (seed and u!=1):continue
            if u>1 and (val(u-1)+1)//2!=j:continue
            if not (v&1) and not seed and u==1:continue
            C.append((j,u))
        if not C:return None
        ck(len(C),1)
        j,h=C[0];path.append(j)
    return tuple(reversed(path))

def inverse5(h):
    path=[]
    while h!=1:
        x=291*h-51
        if h<1:return None
        v=val(x)
        t=(v-4)//2 if not v%2 else (v-3)//2
        if t<0:return None
        for _ in range(t):
            h,r=divmod(h+51,292)
            if r or h<1:return None
            path.append('A')
        if not v%2:return ''.join(reversed(path)) if h==1 else None
        h,r=divmod(h+1023,5976)
        if r or h<1:return None
        path.append('B')
    return ''.join(reversed(path))

def grow(A,root,k,u):
    n=len(A);s=n+k;out=list(A)+[0]*(k+1)
    out[root]|=1<<s
    out[s]=((1<<n)-1)^(1<<root)
    for j in range(k):
        i=n+j;out[i]=((1<<n)-1)+sum(1<<v for v in range(i+1,s))
        if u>>(k-1-j)&1:out[s]|=1<<i
        else:out[i]|=1<<s
    return tuple(out),s
def table(A):
    n=len(A);d=[[0]*n for _ in range(1<<n)]
    for j in range(n):d[1<<j][j]=1
    for m in range(1,len(d)):
        for j,p in enumerate(d[m]):
            if not p:continue
            z=A[j]&~m
            while z:
                b=z&-z;z-=b;d[m|b][b.bit_length()-1]+=p
    return d
def prof(a,r):
    d=table(a);f=len(d)-1
    h=[1]+[sum(j) for j in d[1:]]
    return h[f],sum(h[s]*h[f^s] for s in range(f+1)),sum(d[s][r]*h[f^s] for s in range(f+1)),d[f][r]
def strong(a):
    f=(1<<len(a))-1
    for i in range(len(a)):
        s=1<<i
        for _ in a:
            for j,x in enumerate(a):
                if s>>j&1:s|=x
        if s!=f:return False
    return True

def main():
    A=(174,124,120,240,65,80,1,118);B=(174,124,120,240,65,16,33,118)
    ck(tuple(x^y for x,y in zip(A,B)),(0,0,0,0,0,64,32,0))
    ck((prof(A,0)[:2],prof(B,0)[:2]),((187,3736),(197,3728)))
    T,r=grow((6,4,0),2,6,50)
    ck((prof(T,r)[0],prof(T,r)[2:]),(241,(301,1)));ck(strong(T),True)
    ck(187*241+3736*301,1169603)
    rng=Random(98370);n=0
    for j in range(2,9):
        for i in range(2,j+2):
            q=1<<i;Q=1<<j
            m,T,k=edge(j,i);u=T-1
            a,b,c,d=migration(k,u)
            ck((a+(Q+1)*b,b),(m,T));ck(apply((a,b,c,d),1,Q),(m-T,(q+1)*(m-T)-1))
            ck(val(m),i)
            ck(inverse(i,m-T),(j,i));n+=1
    for h in (1,5,17):
        ck(apply(migration(8,50),h,5*h-1),(292*h-51,5*(292*h-51)-1))
        b,g=apply(migration(8,84),h,5*h-1)
        ck(apply(migration(4,2),b,g),(5976*h-1023,5*(5976*h-1023)-1))
    count5=0
    for d in range(8):
        for w in product('AB',repeat=d):
            h=1
            for x in w:h=292*h-51 if x=='A' else 5976*h-1023
            ck(inverse5(h),''.join(w));count5+=1
    distinct=set()
    for z in range(2,7):
        a,T,_=edge(z,z);b,t,_=edge(z,z+1);c,s,_=edge(z+1,z)
        core=(a-1)-T;ck(val(core),2*z)
        V=(c*b-1)*T-(a-1)*(c*t+s);ck(val(abs(V)),3*z)
        for d in range(7):
            for w in product('AB',repeat=d):
                h=1;path=[z]
                for x in w:
                    if x=='A':h=a*h-T;path.append(z)
                    else:h=c*(b*h-t)-s;path.extend((z+1,z))
                ck(inverse(z,h),tuple(path))
        distinct.add(a)
    ck(len(distinct),5)
    for _ in range(300):
        j=rng.randrange(2,7);h=1;jumps=[j]
        for z in range(rng.randrange(1,12)):
            i=rng.randrange(2,min(j+2,8));m,T,_=edge(j,i)
            h=m*h-T;j=i;jumps.append(j)
        ck(inverse(j,h),tuple(jumps))
        ck(inverse(j,h+2),None)
    print(count5,'independent slope-five binary words pass.')
    print(n,'explicit edge matrices, 635 binary-loop words, 300 arbitrary token paths and adjacent ghosts pass.')
if __name__=='__main__':main()
