"""Direct subset path weights on the 8-core plus exact tape/block replay."""
from typing import NamedTuple

class P(NamedTuple):
    h:int; Q:int; n:int; s:int; d:int; t:int; f:int

# Rows encode outgoing neighbors; bit i is vertex i.
C=(4,1,130,7,143,159,63,75)
MASK=35
Z="000100011001"
W="010111000111"
A=4096

def need(test):
    if not test:
        raise AssertionError("literal slide replay mismatch")

def table(rows):
    v=len(rows); H=[1]+[0]*((1<<v)-1)
    S=[[0]*v for _ in H]
    for u in range(1, len(H)):
        for i in range(v):
            if u>>i&1:
                rest=u^(1<<i)
                S[u][i]=sum(S[rest][j] for j in range(v) if rows[i]>>j&1) if rest else 1
        H[u]=sum(S[u])
    return S,H

def local(rows,m):
    f=(1<<len(rows))-1
    S,H=table(rows)
    R,_=table(tuple(f^v^(1<<i) for i,v in enumerate(rows)))
    out=S[0]; end=R[0]
    for i in range(1,f+1):
        for p in range(len(rows)):
            out[p]+=S[i][p]*H[f^i]
            end[p]+=R[i][p]*H[f^i]
    marked=[i for i in range(len(rows)) if m>>i&1]
    unmarked=[i for i in range(len(rows)) if not m>>i&1]
    adj=tuple(v|(1<<len(rows) if i in unmarked else 0) for i,v in enumerate(rows))+(m,)
    return P(H[-1],sum(H[i]*H[f^i] for i in range(f+1)),table(adj)[1][-1],
             sum(out[i] for i in marked),sum(R[f][i] for i in marked),
             sum(end[i] for i in unmarked),sum(S[f][i] for i in unmarked))

def glue(U,V):
    a,q,n,s,d,t,f=U
    b,r,o,c,e,v,g=V
    return P(a*b,q*r,b*n+(o-b)*a+v*s+g*d,r*s+a*c,a*e,q*v+b*t,b*f)

BIT=(P(1,2,1,0,0,1,1),P(1,2,1,1,1,0,0))
def word(s):
    need(s)
    D=BIT[int(s[0])]
    for bit in s[1:]:
        D=glue(D,BIT[int(bit)])
    return D
def bank(k,j):
    need(0<=j<=k)
    return glue(glue(word("1"+W*j), CORE),word(Z*(k-j)+"0"))

CORE=local(C,MASK)
def main():
    need(CORE==P(45,1440,265,517,31,577,25))
    need(word(Z)==P(1,4096,115,281,1,1655,1))
    need(word(W)==P(1,4096,165,1479,1,453,1))
    need(glue(BIT[1],CORE)==P(45,2880,867,1957,31,1154,0))
    z=glue(CORE,word(Z));w=glue(word(W),CORE)
    need(z[:2]==w[:2] and z.s==w.s and z.t==w.t)
    for k in range(33):
        B=A**k; D=(B-1)//(A-1)
        S=1957*B+12645*D
        sigma=1957*D+12645*((D-k)//(A-1))
        T=5760*B-1726-2880*2440*D
        N=898+S+1655*sigma+5175*k
        for j in range(k+1):
            need(bank(k,j)==P(45,5760*B,N-2*j,2*S,0,T,0))
        need(T%45==29)
    # Tape z/w can vary while core's unsiged profiles stay matched at each slide.
    all_tiles=[]
    for z in range(A):
        a=(A-1)*CORE.s+CORE.h*z
        if a%CORE.Q: continue
        w=a//CORE.Q
        if not 0<=w<A: continue
        zs,ws=f"{z:012b}",f"{w:012b}"
        if CORE.Q*int(zs[::-1],2)!=(A-1)*(CORE.Q-CORE.h-CORE.t)+CORE.h*int(ws[::-1],2):
            continue
        if (ws[-1],zs[0]) != ("1","0"): continue
        left=glue(glue(BIT[1],CORE),word(zs+"0"))
        right=glue(glue(word("1"+ws),CORE),BIT[0])
        need(left[:2]==right[:2] and left[3:]==right[3:])
        all_tiles.append((zs,ws,left.n-right.n))
    need(all_tiles==[
       ("000100011001","010111000111",2),
       ("001100011001","010111010111",-48),
       ("010100011001","010111100111",34),
       ("011100011001","010111110111",-16)])
    def surcharge(t):
        if not t: return 0
        u=t//2
        return surcharge(u)+(u+(u%2) if t%2==0 else 0)
    for k in (0,1,2):
        P0=bank(k,0)
        for j in range(k+1):
            Pj=bank(k,j)
            for x in range(16):
                for y in range(16):
                    out=Pj
                    if x:
                        out=glue(word(bin(x)[2:]),out)
                    if y:
                        out=glue(out,word(''.join('1' if a=='0' else '0' for a in bin(y)[2:][::-1])))
                    expect=P0.n-2*j+Pj.Q*x*y+Pj.t*x+Pj.s*y+45*(surcharge(x)+surcharge(y))
                    need(out.n==expect)
    for k in range(4):
        D=bank(k,0)
        print(f"height={k} order={11+12*k} odd interval=[{D.n-2*k},{D.n}] S={D.s} T={D.t}")
    print("8-core and all four compatible width-12 tiles verified; symbolic profile checked for 561 slide positions through height 32.\n1536 exact two-sided placements tested.")

if __name__ == "__main__":
    main()
