from pathlib import Path
import sys
from math import isqrt
root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/'bank'))
from steer import run,params
sys.path.insert(0,str(root/'arithmetic'))
from selector import F,count,scale,ceilsqrt

def ck(a,b=True):
    if a!=b:raise ArithmeticError((a,b))

def adjust(M,N,X,Y,xr,sr):
    s,n,d,C,V,W,G,K=params(M)
    K-=2 # M is taken divisible by4; remove the two-bit initial toggle.
    D=1<<K
    q,z,f=X+Y+1,X+1,F(X)+F(Y)-Y
    ck(f,F(X)+F(Y)-q+z)
    ck(q%4,0); ck(z%2,0)
    Q=q<<K
    dy=1<<s
    ck(Q%dy,0)
    for x2 in (1,3):
        S=Q-1
        F2=(N-81-142*S-256*x2*(S-x2))*pow(5,-1,dy)%dy
        f2=(F2+x2+1)%dy
        if f2%4==f%4:break
    else:raise ArithmeticError("parity selector")
    tail=format(x2,'0'+str(s)+'b')
    _,zs,fs=run((0,0,0),tail)
    need2=(f2-fs)%dy
    word=G
    q,z,f=run((q,z,f),word)
    zr=(xr+1)%n
    fr=(sr-Q+zr)%n
    z0=(zr-zs)*pow(2,-len(tail),n)%n
    def even(v):return v+n*(v%2)
    _,_,df=run((q%n,even(z0),0),tail)
    f0=(fr-df)%n
    dest=q*pow(2,-1,n)%n
    step=C*(2*(z-dest)%n)
    q,z,f=run((q,z,f),step);word+=step
    ck(z%n,dest)
    onward=C*(2*(dest-z0)%n)
    _,_,ff=run((q%n,even(dest),0),onward)
    step=V*((f+ff-f0)%n)+onward
    q,z,f=run((q,z,f),step);word+=step
    ck((z%n,f%n),(z0,f0))
    how=((f-need2)//4)*pow(n,-1,dy//4)%(dy//4)
    ck((f-need2)%4,0)
    word+=W*(n*how)+G
    pad=K-s-len(word)
    ck(pad>=0);ck(pad%d,0)
    word+='1'*pad+tail
    x=(X<<K)+int(word,2);y=(Y<<K)+D-1-int(word,2)
    ck(x+y,Q-1)
    ck(x%n,xr%n)
    ck((F(x)+F(y))%n,sr%n)
    ck(count(x,y)%M,N%M)
    return x,y

def decode(N,M):
    s=2
    while M%(1<<(s+1))==0:s+=1
    n=M//(M&-M);MM=(1<<s)*n
    K=params(MM)[-1]-2
    D=1<<K
    A=scale(N)
    L=1
    while (2*L)**3<=A:L*=2
    ck(L>=4*D and L>=2**s)
    xr,sr=(-2*(N-81))%n,(57*(N-81))%n
    base=((xr+1)*pow(L,-1,n)%n)*L-1
    total=base+L*n*((A-base+L*n-1)//(L*n))
    b=(total-1)//2
    a=b-ceilsqrt(2*total*(total-A+1))
    ck(0<a<b)
    ck(count(a,total-a)<N<=count(b,total-b))
    while b-a>1:
        j=(a+b)//2
        if count(j,total-j)<N:a=j
        else:b=j
    X=(a//D)|1; Y=(total+1)//D-1-X
    ck(X>0 and Y>0)
    x,y=adjust(MM,N,X,Y,xr,sr)
    H=count(x,y)
    ck(H%M,N%M)
    ck(abs(H-N)**3<=(100000*D*n)**3*N)
    return H-N

if __name__=='__main__':
    trials=0
    for M in (1,5,7,8,12,24,32,40,80):
        o=M//(M&-M)
        s=max(2,(M&-M).bit_length()-1)
        K=params((1<<s)*o)[-1]-2
        for j in (11,19,27):
            N=(1<<(6*K+120))+j
            decode(N,M);trials+=1
    for M in (4,8,12,24,36,60,72,100,160):
        s,n,_,_,_,_,_,K=params(M);K-=2
        q=(3+4+1)<<K
        for xr,sr in ((1,2),(4,9),(3,3)):
            y=q-1-xr
            t=(81+142*(q-1)+256*xr*y+5*sr)%n
            t+=n*((1-t)%2)
            adjust(M,t,3,4,xr,sr)
    print(f'{trials} one-third absolute+phase controls, 27 larger-clock route controls pass.')
