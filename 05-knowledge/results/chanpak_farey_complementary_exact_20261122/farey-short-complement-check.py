from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'chanpak_decimal_kernel_20261122'))
from math import gcd
from random import Random
from local_word_inverse import steer
from hyperbola_controls import chart

def f(x):
    out=0
    while x:
        q=x>>1
        if not x&1:out+=q+(q&1)
        x=q
    return out

def dig(x,w):return format(x,f"0{w}b")
def rev(x,w):return int(dig(x,w)[::-1],2)
def eg(a,b):
    if not b:return abs(a),1 if a>=0 else -1,0
    d,x,y=eg(b,a%b);return d,y,x-a//b*y
def bank(a,M):
    u=(5-M)%8;v=u+M;b=(v+1).bit_length();R=[]
    for U,V,w in ((u,v,b),(u+1,v+1,b),(0,M,M.bit_length())):
        U<<=a+1;V<<=a+1;w+=a+1;B=1<<w
        def K(W):return f(W)+f(B-1-W)+(W>>(w-1))
        r=rev(V,w)-rev(U,w)
        assert r&1 and (V-U)%(2*(M<<a))==0
        R.append((U,V,w,r,K(V)-K(U),(f(V)-f(U))//2))
    if M%3==0:
        assert (R[0][4]-R[1][4])%3
        if not sum(t[4] for t in R)%3:
            i=next(i for i,t in enumerate(R) if t[4]%3)
            U,V,w,r,c,d=R[i];R[i]=(V,U,w,-r,-c,-d)
    return R

def steer_both(D,x0,y0,desx,dsum):
    a=(D&-D).bit_length()-1;M=D>>a;period=2*D;width=(period-1).bit_length()
    x,y=x0,y0
    options=[];minus=0
    def adv(word):
        nonlocal x,y
        n=len(word);v=int(word,2)
        x=(x<<n)+v;y=(y<<n)+(1<<n)-1-v
    if M>1:
        R=bank(a,M)
        g=M;coef=[]
        for _,_,_,r,_,_ in R:
            gg,u,v=eg(g,-4*r);coef=[c*u for c in coef]+[v];g=gg
        assert g in (1,3)
        for j in range((M-1).bit_length()):
            s=x+y;const=0
            for _,_,w,r,c,_ in R:
                s=((s+1)<<width)-1
                const+=s*r+c
                s=((s+1)<<w)-1
            step=1<<j
            if (step-const)%g:step=-step;minus-=step
            assert (step-const)%g==0
            words=[[],[]]
            for z,(U,V,w,r,_,dx) in zip(coef,R):
                d0=(z*((step-const)//g))%M
                d2=(dx*pow(r,-1,1<<a))%(1<<a)
                d=d2+(1<<a)*((d0-d2)*pow(1<<a,-1,M)%M)
                load=dig((2*d-(x<<width))%period,width)
                adv(load)
                assert (x//2-d)%D==0 and not x%2 and y%2
                words[0].extend((load,dig(U,w)));words[1].extend((load,dig(V,w)))
                adv(dig(U,w))
            options.append(("".join(words[0]),"".join(words[1]),step<0))
    tail=steer(D,(x%period,(f(x)//2)%D),(desx,0))
    adv(tail)
    rank=(dsum-f(x)-f(y)+minus)%M
    use=[]
    for j,(u,v,negative) in enumerate(options):
        use.append(v if bool(rank>>j&1)^negative else u)
    word="".join(use)+tail
    x,y=x0,y0;adv(word)
    assert x%period==desx%period and (f(x)//2)%(1<<a)==0
    assert (f(x)+f(y)-dsum)%M==0
    assert x+y==((x0+y0+1)<<len(word))-1
    return x,y,len(word)

if __name__=="__main__":
    r=Random(61);checks=0
    for n in (4,8):
        h,q,c,t,e,rbit,k=chart(n)
        for D in list(range(1,100))+[1215,73872,6*5**12,12*3**10*5**5,64*(2**41-1)]:
            a=(D&-D).bit_length()-1;M=D>>a
            s=(1<<(a+4))-1;x0=r.randrange(s+1);y0=s-x0
            _,_,length=steer_both(D,x0,y0,1,0)
            S=((s+1)<<length)-1
            g=gcd(h,M)
            T=k+rbit+t*S+g*(1+r.randrange(2*D))
            if not T&1:T+=g
            C=(f(y0)-f(x0))//2 + x0//2
            def phase(x):
                return k+rbit+c*S+2*h*C+h+(t-c+q*S-h)*x-q*x*x
            x2=1
            for i in range(1,a):
                if (phase(x2)-T)%(1<<(i+1)):x2+=1<<i
            # Require odd x, a fixed dyadic root, and x=S mod M.
            dy=1<<max(1,a)
            des=(x2+dy*((S-x2)*pow(dy,-1,M)%M))%(2*D)
            dsum=((T-k-rbit-t*S)//g)*pow(h//g,-1,M//g)%(M//g)
            x,y,ln=steer_both(D,x0,y0,des,dsum)
            assert ln==length and x+y==S
            assert (k+q*x*y+t*x+c*y+rbit*(x&1)+e*(y&1)+h*(f(x)+f(y))-T)%D==0
            checks+=1
    print(checks,"literal complementary chart phase checks")
