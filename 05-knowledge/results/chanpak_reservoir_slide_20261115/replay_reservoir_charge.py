"""Literal small DP and prescribed-target unit slide inversion."""
from math import isqrt
from replay_unit_slide import bank,local,need,table

SINGLE=(0,)
C3=(2,4,1)
C5H=(6,12,8,1)  # four vertices, count 5
C37=(30,44,56,48,34,1)
C37b=(22,44,56,49,34,1)

def ordinal(a,b):
    n=len(a)
    return tuple(r|(((1<<len(b))-1)<<n) for r in a)+tuple(r<<n for r in b)

def apex(rows,mask):
    n=len(rows)
    return tuple(r|((1-(mask>>i&1))<<n) for i,r in enumerate(rows))+(mask,)

def mirror(rows,mask):
    F=(1<<len(rows))-1
    return tuple(F^r^(1<<i) for i,r in enumerate(rows)),F^mask

def tape(m):
    return tuple(((1<<m)-1)^((1<<(i+1))-1) for i in range(m))

def add_tail(R,mask,T):
    return ordinal(R,T),mask|(((1<<len(T))-1)<<len(R))

def literal(A,u,B,v,T=(),U=()):
    AA,x=add_tail(A,u,T)
    BB,y=add_tail(B,v,U)
    L=apex(AA,x)
    R,z=mirror(BB,y)
    Y=apex(ordinal(L,R),(1<<len(AA))|(z<<len(L)))
    return len(Y),table(Y)[1][-1]

def prime(n):
    return n>=2 and all(n%p for p in range(2,isqrt(n)+1))

def ceildiv(a,b):
    return -((-a)//b)

def decode(N,k):
    if N%2==0 or N%37==0 or k%3==0:
        return []
    D=bank(k,0); eta=45+D.s; charge=704*eta
    m0=next(m for m in range(36) if charge*pow(2,m,37)%37==N%37)
    lo=max(1,ceildiv(N-37*D.n,charge))
    hi=(N-37*(D.n-2*k))//charge
    if hi<lo:
        return []
    m_min=(lo-1).bit_length()
    out=[]
    for m in range(m0+36*max(0,ceildiv(m_min-m0,36)),hi.bit_length(),36):
        j,remainder=divmod(37*D.n+charge*2**m-N,74)
        need(remainder==0 and 0<=j<=k)
        out.append((m,j))
    return out

def main():
    need((local(C37,0).h,local(C37,0).Q)==(37,352))
    need((local(C37b,0).h,local(C37b,0).Q)==(37,376))
    need((local(C5H,0).h,local(C5H,0).Q)==(5,32))
    need(literal(C3,1,C3,1)==(8,57))
    need(literal(C3,1,C3,1,tape(1),tape(1))==(10,169))
    need(literal(C3,1,C3,1,(),tape(3))==(11,313))
    need(literal(SINGLE,1,C3,1,C37)==(12,4409) and prime(4409))
    need(literal(SINGLE,1,SINGLE,1,C5H)==(8,133) and prime(19) and 133==7*19)
    need(literal(SINGLE,1,C3,1,C37b)==(12,4697))
    need(literal(SINGLE,1,C3,1,C37,C5H)==(16,136093))
    need(pow(2,18,37)==36 and pow(2,12,37)==26)
    for k in range(22):
        D=bank(k,0); eta=45+D.s
        need(eta%37==(0,4,34)[k%3])
        if k%3==0:
            continue
        need({704*eta*pow(2,m,37)%37 for m in range(36)}==set(range(1,37)))
        for m in range(0,60,11):
            for j in (0,k):
                N=37*(D.n-2*j)+704*eta*2**m
                need((m,j) in decode(N,k))
    need(decode(11721300267,1)==[(0,1)])
    need(decode(11721300269,1)==[])
    print("Literal DP: mirrored charges 57,169,313; prime 4409; seven-prime 133")
    print("Same-count tags (37,352)/(37,376); crossed tags (37,352)/(5,32) -> 136093")
    print("Fixed-height 37 decoder: 11721300267 -> [(0,1)], next odd not in that bank")

if __name__=="__main__":
    main()
