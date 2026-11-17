from replay import F, calc, profile, check, D, M as masks
from math import gcd,isqrt
from fractions import Fraction
from itertools import product

def run(s,bits):
    q,z,f=s
    for b in bits:
        if b=='0':
            # Last x parity is z-1 mod 2; y parity is complementary.
            f+=2*z-q+((z-1)&1)-2
            z=2*z-1
        else:
            f+=1-((z-1)&1)
            z=2*z
        q*=2
    return q,z,f

def params(m):
    s=0;n=m
    while n%2==0:n//=2;s+=1
    d=1
    while pow(2,d,n)!=1%n:d+=1
    C='01'+'1'*((-2)%d)
    U='1'+C+'1'*((-1)%d)
    E='001'+'1'*((-3)%d)
    V=E+U*((-3)%n)
    r=next(i for i in range(max(s,2),max(s,2)+d) if (i+3)%d==0)
    W='1'*r+'001'
    G='1'*(d*((max(s,2)+d-1)//d))
    # Prefix, guards, up to 2n translations and n f adjustments;
    # at most 2^s neutral dyadic macrosteps. Round to legal length.
    L=2+2*len(G)+2*n*len(C)+n*len(V)+(1<<max(s-2,0))*n*len(W)+s
    L+=((2+s-L)%d)
    return s,n,d,C,V,W,G,L

def suffix(m,X,Y,xr,sr,x2=1):
    s,n,d,C,V,W,G,L=params(m)
    suf=format(x2,'0'+str(s)+'b') if s else ''
    q0=X+Y+1;z0=X+1;f0=F(X)+F(Y)-q0+z0
    check((q0%4,z0%2),(2,0))
    qr=q0*pow(2,L,n)%n
    zr=(xr+1)%n
    fr=(sr-qr+zr)%n
    modulation=1<<s
    fr2=(sr+x2+1)%modulation
    drift=run((0,0,0),suf)[2]
    need2=(fr2-drift)%modulation
    prefix=next(t for t in ('11','01') if (run((q0,z0,f0),t)[2]-need2)%min(4,modulation)==0)
    word=prefix+G
    q,z,f=run((q0,z0,f0),word)
    # The final suffix takes z via z -> 2^s*z + zs; work backward.
    qs,zs,fs=run((0,0,0),suf)
    z_goal=(zr-zs)*pow(2,len(suf)*-1,n)%n
    dq,dz,df=run((q%n,z_goal+n*(z_goal&1),0),suf)
    check(dq%n,qr)
    check(dz%n,zr)
    f_goal=(fr-df)%n
    dest=q*pow(2,-1,n)%n
    shift=2*(z-dest)%n
    part=C*shift
    q,z,f=run((q,z,f),part)
    word+=part
    check(z%n,dest)
    shift2=2*(dest-z_goal)%n
    _,_,c=run((q%n,z%n+n*((z%n)&1),0),C*shift2)
    part=V*((f+c-f_goal)%n)+C*shift2
    q,z,f=run((q,z,f),part)
    word+=part
    check((z%n,f%n),(z_goal,f_goal))
    reps=next(i for i in range(1<<max(0,s-2)) if (f-4*n*i-need2)%modulation==0)
    word+=W*(n*reps)+G
    count=L-len(suf)-len(word)
    check((count>=0,count%d),(True,0))
    word+='1'*count+suf
    x=(X<<L)+int(word,2)
    y=(Y<<L)+((1<<L)-1-int(word,2))
    check(x%n,xr%n);check(x%modulation,x2%modulation)
    check((x+y+1),(X+Y+1)<<L)
    q,z,f=run((q0,z0,f0),word)
    check((q,z),(x+y+1,x+1))
    check(f,F(x)+F(y)-q+z)
    check((F(x)+F(y))%m,sr%m)
    return x,y

def decoder(m,N):
    s,n,d,_,_,_,_,L=params(m)
    h,P,l,r,lm,rm,k=profile(D,masks[0])
    S2=1<<s
    V=1<<max(s+1,2,(N.bit_length()-1)//10)
    A=P*V*V
    B=r*V+h*(V//2-1)
    C=(l+P)*V+h*(V-1)
    j=0
    if n%3==0:j=next(j for j in range(3) if (k+2*j-h+r-rm-N)%3==0)
    K=k+2*j-h+r-rm
    # Replay uses brute force for the local residue tuple instead of CRT.
    for a,b in product(range(n),repeat=2):
        rhs=(N-A*a*b-B*a-C*b-K)%n
        g=gcd(h,n)
        if rhs%g==0:break
    else:raise AssertionError((m,N))
    z=(rhs//g*pow(h//g,-1,n//g))%(n//g)
    z2=((N-K)*pow(h,-1,S2)-1)%S2
    z+=n//g*((z2-z)*pow(n//g,-1,S2)%S2)
    DD=1<<L
    targetS0= ((a+b+1)*pow(DD,-1,n)-1)%n
    targetS0+=n*((1-targetS0)*pow(n,-1,4)%4)
    dd=4*n
    E=Fraction(K)-Fraction(B*C,A)
    # Approximate the threshold from below then correct by at most two steps.
    Sguess=isqrt((4*(N-E)//A))
    S0=targetS0+dd*((Sguess//DD-targetS0)//dd)
    quad=lambda S: A*(S+Fraction(B+C,A))**2/4+E
    while quad(DD*(S0+1)-1)<N:S0+=dd
    S=DD*(S0+1)-1
    value=lambda X: A*X*(S-X)+B*X+C*(S-X)+K
    lo,hi=0,(S*A+B-C)//(2*A)+1
    while hi-lo>1:
        x=(lo+hi)//2
        if value(x)>=N:hi=x
        else:lo=x
    X=(hi//DD)|1
    Y=S0-X
    check(X>0 and Y>0,True)
    x,y=suffix(m,X,Y,a,z)
    check((x+y)%n,(a+b)%n)
    u,w=V*x+1,V*y
    p=(h,P,l,r,lm,rm,k+2*j)
    H=calc(p,u,w)
    check(H% m,N%m)
    check(H,A*x*y+B*x+C*y+K+h*(F(x)+F(y)))
    # A deliberately generous coefficient, evaluated with exact integers.
    check(abs(H-N)**5 <= (1000000*DD*DD*m*m)**5*N*N,True)
    return L,j

def main():
    # All allowed dyadic sums and all odd states for a modest composite wheel.
    checked=0
    for m in (1,3,5,7,8,9,10,12,13,20,24,36,40,52,72):
        s,n,_,_,_,_,_,L=params(m)
        for xr,sr in product(range(n),range(0,m,2 if s else 1)):
            X=3+2*(xr%3);Y=(2-X)%4+3
            # Find a positive even Y with X+Y=1 mod4.
            Y=(1-X)%4+4
            suffix(m,X,Y,xr,sr)
            checked+=1
        for a in (1,3,23):
            N=(10**(2*(L+100)))+a
            decoder(m,N)
    print(f'{checked} complementary routes and 45 simultaneous absolute/phase checks passed.')
if __name__=='__main__':main()
