from math import gcd, isqrt
from random import Random
from replay_pair import B,Q,T,S,complement, pair_advance, pair_compiler, rev
from replay_prefix import need, p


def compile_even(M):
    a = (M&-M).bit_length()-1
    m=M>>a
    _,o,R,*_=pair_compiler(m)
    l=(2*M-1).bit_length()
    pad=a+1+(-a-1)%o
    K=a+2+(-a-2)%o
    L0=pad
    resets=[]
    for j in range(a-1):
        d=j+2
        top=L0+m*K
        zeros={top,top+d}|{L0+i*K for i in range((-1-pow(2,d,m))%m)}
        raw=''.join('0' if i in zeros else '1' for i in range(top+d,-1,-1))
        raw += '1'*(-len(raw)%o)
        resets.append(raw*m)
    L=2+2*l+(m-1)*len(R)+pad+sum(map(len,resets))
    return M,a,m,o,R,l,pad,resets,L


def prescribe(comp,x,y,q,G,qgoal):
    M,a,m,o,R,l,pad,resets,L=comp
    P=((x+y+1)<<L)-1
    X,Y,C=pair_advance(x,y,q,'11',M)
    inv=(C-Y//2)%2
    y4=2*((qgoal-inv)%2)
    def choose_x(raw,base):
        # least entry congruent to raw mod the dyadic base, and G mod m
        return raw+base*((G-raw)*pow(base,-1,m)%m) if m>1 else raw
    dest=choose_x((P-y4)%4,4)
    w=((X+Y+1)<<l)%(2*M)
    def bits(z): return format(z%(2*M),f'0{l}b')
    d=1<<(a+1)
    # Normalize x to -1+w/2 mod odd m, -1 mod dyadic side.
    xt=d-1+d*((w//2-d)*pow(d,-1,m)%m) if m>1 else d-1
    init=bits(xt-(X<<l))
    X,Y,C=pair_advance(X,Y,C,init,M)
    need(X==xt and X%2==1)
    final=bits(dest-(xt<<l))
    _,_,drift=pair_advance(xt,Y,0,final,M)
    want=(qgoal-drift)%M
    oddrep=2*(want-C)%m
    code='11'+init+R*oddrep+'1'*((m-1-oddrep)*len(R)+pad)
    X,Y,C=pair_advance(x,y,q,code,M)
    need(X==xt and Y%d==0 and want%m==C%m and want%2==C%2)
    bitcode=((want-C)//2*pow(-m,-1,1<<(a-1)))%(1<<(a-1))
    for j,word in enumerate(resets):
        if bitcode>>j&1:
            code+=word
    code+='1'*(L-l-len(code))+final
    need(len(code)==L)
    X,Y,C=pair_advance(x,y,q,code,M)
    need(X%m==G%m and C==qgoal%M)
    return code


def bracket(N,M):
    comp=compile_even(M)
    a,m=comp[1:3]
    D=1<<comp[-1]
    b=isqrt((N-B)//Q)
    r=max(a+1,b.bit_length()//5)
    radix=1<<r
    g=gcd(45,m)
    v=((N-B+S*radix)*pow(T+S,-1,g))%g if g>1 else 0
    need(v<radix)
    low_sum=2*b//radix+6
    P=((low_sum+D*m)//(D*m))*D*m-1
    C=B+Q*v*(v+radix*P)+T*v+S*(radix*P+v)+45*(P*(radix-1-rev(v,r))+2*p(v)+(v<radix//2))
    need((N-C)%(2*g)==0)
    qgoal=((N-C)//(2*g))*pow(45//g,-1,M//g)%(M//g) if M//g>1 else 0
    def base(u): return C+Q*radix*radix*u*(P-u)+(T-S)*radix*u
    def first_below(H):
        lo,hi=(P+1)//2,P
        need(base(lo)>=H and base(hi)<H)
        while hi-lo>1:
            mid=(hi+lo)//2
            if base(mid)>=H: lo=mid
            else: hi=mid
        return hi
    above=first_below(N)//D-1
    below=(first_below(N-45*P)+D-1)//D
    values=[]
    for u0 in (below,above):
        v0=(P+1)//D-1-u0
        need(u0>=0 and v0>=0)
        s=prescribe(comp,u0,v0,(p(u0)+p(v0))//2,0,qgoal)
        u,z=D*u0+int(s,2),D*v0+int(complement(s),2)
        x,y=radix*u+v,radix*z+v
        H=B+Q*x*y+T*x+S*y+45*(p(x)+p(y))
        need((H-N)%(2*M)==0 and H==base(u)+45*(p(u)+p(z)))
        values.append(H)
    need(values[0]<=N<=values[1])
    need(max(abs(H-N) for H in values)<=100000*(D*m+1)**2*(b//radix+1))
    return values


def main():
    rng=Random(93511)
    t=0
    for M in (2,4,6,8,10,12,16,18,20,24,28,32,36,40,64,72):
        comp=compile_even(M)
        m=comp[2]
        L=comp[-1]
        D=1<<L
        for j in range(35):
            x,y=rng.getrandbits(j),rng.getrandbits(j)
            G,qgoal=rng.randrange(m),rng.randrange(M)
            tape=prescribe(comp,x,y,(p(x)+p(y))//2,G,qgoal)
            X,Y=D*x+int(tape,2),D*y+int(complement(tape),2)
            need(X+Y==D*(x+y+1)-1 and X%m==G and (p(X)+p(Y))//2%M==qgoal)
            t+=1
    for M in (2,4,6,8,10,12,16,18,24,32,40,64):
        L=compile_even(M)[-1]
        for j in range(8):
            bits=max(100,5*(L+2)+j)
            N=(1<<bits)+(rng.getrandbits(bits)|1)
            bracket(N,M)
            t+=1
    print("PASS: dyadic parity-coupling and odd CRT, exact addresses and brackets;",t,"controls")

if __name__=="__main__":
    main()
