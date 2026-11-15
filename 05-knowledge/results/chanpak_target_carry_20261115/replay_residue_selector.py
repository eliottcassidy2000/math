from fractions import Fraction
from random import Random
from replay_carry_selector import B,Q,T,S,output,reverse,check,cmp,extract

def route(s,L,U,N,M,target,side):
    r=s.bit_length()
    mod=2*M
    if (side==1 and target<0) or (side==-1 and target>=(1<<r)):
        return None
    target=max(0,min((1<<r)-1,target))
    # carry, lower and upper address comparisons, carry comparison,
    # x address, reversed zero suffixes, combined surcharge, last bits
    states={(0,0,0,0,0,0,0,0,1,1):(0,0)}
    for i in range(r):
        bL,bU,bs,bt=((L>>i)&1,(U>>i)&1,(s>>i)&1,(target>>(r-1-i))&1)
        powbit=pow(2,i,mod)
        nxt={}
        for (c,d,e,k,a,Rx,Ry,P,lx,ly),(cval,x) in states.items():
            for bx in (0,1):
                for by in (0,1):
                    z=bx+by+c-bs
                    if z not in (0,2):
                        continue
                    cp=z//2
                    key=(cp,cmp(bx,bL) or d,cmp(bx,bU) or e,k or cmp(cp,bt),
                         (a+bx*powbit)%mod,(2*Rx+1-bx)%mod,(2*Ry+1-by)%mod,
                         (P+bx*(Rx+1-lx)+by*(Ry+1-ly))%mod,bx,by)
                    cand=(2*cval+cp,x+(bx<<i))
                    if key not in nxt or side*(cand[0]-nxt[key][0])>0:
                        nxt[key]=cand
        states=nxt
    eligible=[]
    for (c,d,e,k,x,Rx,Ry,P,lx,ly),ans in states.items():
        y=(s-x)%mod
        if c==0 and d>=0 and e<=0 and k*side<=0 and (B+Q*x*y+T*x+S*y+45*P-N)%mod==0:
            eligible.append(ans)
    return max(eligible,key=lambda a: side*a[0]) if eligible else None

def steer(s,L,U,N,M):
    r=s.bit_length()
    K=1<<(r+1)
    t0=Fraction(L+U-s,2)
    w=Fraction(U-L,2)
    rho=Fraction(reverse(s,r),1<<r)
    center=B+Q*(Fraction(s*s,4)-t0*t0)+Fraction((T+S)*s,2)+(T-S)*t0+45*s*(1-rho/2)
    E=Q*(2*abs(t0)*w+w*w)+abs(T-S)*w+45*(abs(t0)+w)+90*(r+1)
    target=K*(center-N)/Fraction(135*s,2)
    n,d=target.numerator,target.denominator
    a=route(s,L,U,N,M,n//d,1)
    b=route(s,L,U,N,M,-(-n//d),-1)
    live=[z for z in (a,b) if z is not None]
    if not live:
        return None,target,E
    c,x=min(live,key=lambda z:(abs(z[0]-target),z[0]))
    model=center-Fraction(135*s*c,2*K)
    cnt=output(x,s-x)
    check((cnt-N)%(2*M)==0 and abs(cnt-model)<=E)
    R=abs(N-model)
    if R+E<2*M:
        check(cnt==N)
    return (c,x,R),target,E

if __name__=='__main__':
    rng=Random(301118)
    k,empty,good,ruled=0,0,0,0
    for s in range(1,221):
        for i in range(14):
            L=rng.randrange(s+1)
            U=rng.randrange(L,s+1) if i%3==0 else rng.randrange(L,min(L+8,s)+1)
            x0=rng.randrange(L,U+1)
            N=output(x0,s-x0)+(0 if i%2==0 else 2*rng.randint(-100,100))
            M=rng.choice((1,2,3,4,5,9,16,31,127,1009))
            res,target,E=steer(s,L,U,N,M)
            allowed={extract(x,s-x,s.bit_length()) for x in range(L,U+1)
                     if (output(x,s-x)-N)%(2*M)==0}
            if res is None:
                check(not allowed)
                empty+=1
            else:
                c,x,R=res
                check(L<=x<=U and extract(x,s-x,s.bit_length())==c)
                check(abs(c-target)==min(abs(z-target) for z in allowed))
                if R>E:
                    check(all(output(x,s-x)!=N for x in range(L,U+1)))
                    ruled+=1
                if R+E<2*M:
                    good+=1
            k+=1
    print('PASS: terminal-residue carry windows:',k,'no path:',empty,
          'envelope exclusions:',ruled,'certified exact:',good)
