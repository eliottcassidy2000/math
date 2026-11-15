from math import gcd, isqrt
from random import Random

from replay_prefix import advance, need, p, period

B,Q,T,S,h = 2855,5760,4034,3914,45


def complement(s):
    return ''.join('1' if c=='0' else '0' for c in s)


def pair_advance(x,y,c,s,m):
    for b in s:
        t = x if b=='0' else y
        c += (t+(t&1))//2
        x,y = 2*x+(b=='1'),2*y+(b=='0')
        x,y,c = x%(2*m),y%(2*m),c%m
    return x,y,c


def pair_compiler(m):
    need(m%2==1)
    o = period(m)
    D = "1"*((-2)%o)+"01"
    R = "1"*((-1)%o) + D + "1" + D*((-2)%m) if m>1 else ""
    need(len(R)%o==0)
    l=(2*m-1).bit_length()
    return m,o,R,l,2*l+(m-1)*len(R)


def prescribe_pair(comp,x,y,c,goal_x,goal_c):
    m,o,R,l,L=comp
    def bits(v): return format(v%(2*m),f'0{l}b')
    w=((x+y+1)<<l)%(2*m)
    tt=(-1+w//2)%m
    if tt%2==0: tt+=m
    init=bits(tt-(x<<l))
    X,Y,C=pair_advance(x,y,c,init,m)
    need(X%2==1 and (X+Y+1)%m==w%m)
    final=bits(goal_x-(X<<l))
    _,_,C_out=pair_advance(X,Y,C,final,m)
    a=2*(goal_c-C_out)%m
    tape=init+R*a+"1"*((m-1-a)*len(R))+final
    need(len(tape)==L)
    xx,yy,cc=pair_advance(x,y,c,tape,m)
    need(xx%m==goal_x%m and cc==goal_c%m)
    return tape


def rev(v,r):
    return int(format(v,f'0{r}b')[::-1],2)


def bracket(N,m):
    comp=pair_compiler(m)
    D=1<<comp[-1]
    a=isqrt((N-B)//Q)
    r=max(1,a.bit_length()//5)
    q=1<<r
    g=gcd(45,m)
    v=((N-B+S*q)*pow(T+S,-1,g))%g if g>1 else 0
    need(v<q)
    pmin=2*a//q+6
    P=((pmin+D*m)//(D*m))*D*m-1
    C=B+Q*v*(v+q*P)+T*v+S*(q*P+v)+45*(P*(q-1-rev(v,r))+2*p(v)+(v<q//2))
    need((N-C)%g==0)
    qgoal=((N-C)//g)*pow(90//g,-1,m//g)%(m//g) if m//g>1 else 0
    def F0(u):
        return C+Q*q*q*u*(P-u)+(T-S)*q*u
    def first_below(target):
        lo,hi=(P+1)//2,P
        need(F0(lo)>=target and F0(hi)<target)
        while hi-lo>1:
            v=(lo+hi)//2
            if F0(v)>=target: lo=v
            else: hi=v
        return hi
    up=first_below(N)//D-1
    down=(first_below(N-45*P)+D-1)//D
    outs=[]
    for u0 in (down,up):
        v0=(P-(D-1))//D-u0
        need(u0>=0 and v0>=0)
        code=prescribe_pair(comp,u0,v0,(p(u0)+p(v0))//2,0,qgoal)
        u=D*u0+int(code,2)
        z=D*v0+int(complement(code),2)
        x=q*u+v
        y=q*z+v
        H=B+Q*x*y+T*x+S*y+45*(p(x)+p(y))
        need((H-N)%m==0 and (H-N)%2==0)
        need(H==F0(u)+45*(p(u)+p(z)))
        outs.append(H)
    need(outs[0]<=N<=outs[1])
    # q~a^(1/5), with constant losses from aligning programmable blocks.
    need(max(abs(H-N) for H in outs) <= 100000*(D*m+1)**2*(a//q+1))
    return outs


def main():
    rng=Random(901281)
    t=0
    for m in (1,3,5,7,9,11,15,21,25,27,33,45):
        comp=pair_compiler(m)
        D=1<<comp[-1]
        for j in range(40):
            x,y=rng.getrandbits(j),rng.getrandbits(j)
            gx,gc=rng.randrange(m),rng.randrange(m)
            tape=prescribe_pair(comp,x,y,(p(x)+p(y))//2,gx,gc)
            X,Y=D*x+int(tape,2),D*y+int(complement(tape),2)
            need(X+Y==D*(x+y+1)-1 and X%m==gx and (p(X)+p(Y))//2%m==gc)
            t+=1
    for m in (1,3,5,7,9,15,21):
        L=pair_compiler(m)[-1]
        for j in range(10):
            bits=max(100,5*(L+2)+j)
            N=(1<<bits)+(rng.getrandbits(bits)|1)
            bracket(N,m)
            t+=1
    print("PASS: literal paired addresses, arbitrary initial sums, and two-sided target brackets;",t,"controls")

if __name__=="__main__":
    main()
