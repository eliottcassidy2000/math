"""Literal and symbolic controls for the S331 rooted and directed returns."""
from pathlib import Path
from math import gcd
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent /
                      "session328_fermat_gluing_and_phase_repair_20261114"))
from fermat_profiles import ends, glue, rooted, strong


def verify(ok, label):
    if not ok:
        raise RuntimeError(label)


def delete(A, v):
    return tuple(r & ((1 << v)-1) | (r >> (v+1) << v)
                 for i, r in enumerate(A) if i != v)


def total(A):
    return sum(ends(A)[-1]) if A else 1


def state(A, r):
    return (total(A), total(delete(A,r)), list(rooted(A))[r][1], ends(A)[-1][r])


def step(S, j, k):
    H,h,L,E=S
    d=1<<j
    if k==0:
        return H+d*L,h,2*d*L-h,h
    if k==j+1:
        return H+d*L,H,H+d*L,H
    return (H+d*L, H+d//2*L if j>1 else H+L-E,
            2**(j-k+1)*H+(2**k)*3**(j-k)*L-(E if k==1 else 0),
            2**(k-1)*L)


def realize(A,r,word):
    for j,k in word:
        n=len(A)
        A=glue(A,r,j)
        if k:r=n+k-1
    return A,r


def compose(H,e,prefix,bit,guard=False):
    S=step((H,e,H,e),2,prefix)
    verify(gcd(S[0],S[2])==1, "safe prefix")
    J=(2,4) if bit==0 else (4,2)
    for j in J:S=step(S,j,j)
    k=4 if guard and (S[0]+2*S[2])%3==0 else 1
    T=step(S,k,k+1)
    return (T[0],T[1]), [(2,prefix),*[(j,j) for j in J],(k,k+1)]


def gate(mods):
    A=[]
    last=None
    for P,s,t in mods:
        n=len(A)
        A.extend((r << n)|((1<<n)-1) for r in P)
        if last is not None:
            A[last]^=1<<(n+s)
            A[n+s]^=1<<last
        last=n+t
    return tuple(A)


def fixed(A,s,t):
    n=len(A)
    dp=[[0]*n for _ in range(1<<n)]
    dp[1<<s][s]=1
    for X,row in enumerate(dp):
        for v,c in enumerate(row):
            if c:
                Y=A[v]&~X
                while Y:
                    b=Y&-Y
                    Y^=b
                    dp[X|b][b.bit_length()-1]+=c
    return dp[-1][t]


def main():
    for A,r in [((0,),0),((2,4,1),1),((2,4,9,19,7),0),((2,4,9,3),2)]:
        S=state(A,r)
        for j in (1,2,4):
            for k in range(j+2):
                B,t=realize(A,r,[(j,k)])
                verify(state(B,t)==step(S,j,k),("port",A,r,j,k,state(B,t),step(S,j,k)))
    print("literal corrected ports, including first-interior subtraction",flush=True)

    A=(2,4,9,19,7)
    for sizes,L in [((7,),1403),((2,2,3),1203)]:
        B,t=realize(A,0,[(r-1,0) for r in sizes])
        verify(state(B,t)[:3]==(713,5,L), "one versus three")
    A23=(2,4,25,19,35,15)
    A31=(2,12,33,21,7,27)
    ordinal=tuple(r|((63)<<6) for r in A23)+tuple(r<<6 for r in A31)
    verify(total(ordinal)==total(A23)*total(A31)==713, "ordinal 23 times 31")
    values={}
    for x in (0,1):
        for y in (0,1):
            breaks=[3]*x+[6]*y+[10]
            bs=[b-a for a,b in zip([0]+breaks,breaks)]
            B,_=realize(A,0,[(z-1,0) for z in bs])
            values[x,y]=total(B)
    verify(2*(values[1,1]-values[1,0]-values[0,1]+values[0,0])==-5*8,
           "mixed difference")
    print("three 713 constructions and mixed-window test",flush=True)

    start,root=realize((0,),0,[(2,3)])
    verify(state(start,root)==(5,1,5,1),"ready seed")
    for bit in (0,1):
        predicted,word=compose(5,1,0,bit)
        B,t=realize(start,root,word)
        verify(strong(B) and (total(B),total(delete(B,t)))==predicted,
               "literal return at common order")
        print("literal primitive return",len(B),predicted,flush=True)
    # A literal ready source where old prefix is bad and the interior prefix repairs.
    start2,t2=realize((0,),0,[(4,5)])
    verify(state(start2,t2)==(17,1,17,1),"forced interior seed")
    B,t=realize(start2,t2,[(2,1)])
    verify(state(B,t)==step(state(start2,t2),2,1) and gcd(total(B),state(B,t)[2])==1,
           "interior safe-prefix control")

    for guard in (False,True):
        layer={(5,1)}
        for d in range(1,7):
            next_layer=set()
            for H,e in layer:
                for q in (0,1):
                    if (8*H-e if q==0 else 10*H-e)%5:
                        for bit in (0,1):
                            pair,_=compose(H,e,q,bit,guard)
                            h,v=pair
                            verify(gcd(h,v)==1 and (not guard or h%3!=0), "phase guard")
                            verify(pair not in next_layer,"injective profile replay")
                            next_layer.add(pair)
            layer=next_layer
            verify(len(layer)==4*3**(d-1), "ternary growth")
        print("guard",guard,"depth six",len(layer),flush=True)

    A=((6,2),(5,0))
    M0=((2201,-212),(709,-68))
    M1=((2321,-272),(685,-80))
    def matvec(M,z,n):
        return tuple(sum(a*b for a,b in zip(row,z))%n for row in M)
    orbit=set()
    v=(5,1)
    for k in range(48):
        orbit.add(v)
        v=matvec(M0,v,7)
    verify(len(orbit)==48 and (0,0) not in orbit and v==(5,1),
           "complete mod-seven orbit")
    S=((1,0),(0,1))
    for k in range(48):
        u,v=[matvec(M0,z,49) for z in zip(*S)]
        S=tuple(zip(u,v))
    verify(S==((43,14),(35,1)), "lift matrix")
    a,b=M1[0]; c,d=M1[1]
    inv=tuple(tuple(pow((a*d-b*c)%7,-1,7)*k%7 for k in row)
              for row in ((d,-b),(-c,a)))
    verify(matvec(A,(5,1),7)==(4,4)
           and matvec(M1,matvec(A,matvec(inv,(5,1),7),7),7)==(0,1),
           "independent tangent pulses")
    for word in ("110000","111011"):
        H,e=5,1
        for bit in map(int,word):
            H,e=compose(H,e,0,bit)[0]
            verify(gcd(H,165)==1,"fixed-order 99 guard")
        verify(e*pow(H,-1,11)%11==2, "projective guarded return")
    x=(5,1)
    for bit in map(int,"0001011011111"):
        x=matvec((M0,M1)[bit],x,49)
    verify(x==(5,1), "coprime loop control")
    for M in (M0,M1):
        S=((1,0),(0,1))
        for i in range(1,337):
            u,v=[matvec(M,z,49) for z in zip(*S)]
            S=tuple(zip(u,v))
            verify((S==((1,0),(0,1)))==(i==336), "matrix order modulo 49")
    H,e=5,1
    for l in range(1,337):
        H,e=(2201*H-212*e)%49,(709*H-68*e)%49
        verify(((H,e)==(5,1))==(l==336),"exact orbit mod 49")
    verify((2321*H-272*e)%49==14,"unbounded exact valuation one")

    C=(34,4,137,19,39,78,159,59)
    verify(strong(C) and fixed(C,0,7)==2,"literal directed interlaced cell")
    mods=[(((2,4,9,3),0,3),((2,4,25,19,35,79,31),6,0)),
          (((2,8,3,21,7),4,0),((2,4,25,19,35,15),5,0))]
    for m,H in zip(mods,[10641,16883]):
        V=gate([((0,),0,0),*m,((0,),0,0)])
        verify(strong(V) and fixed(V,0,12)==6 and total(V)==H,"same marked product")
        print("marked 6 vs unmarked",H,flush=True)


if __name__=="__main__":
    main()
