"""Physical one-apex recurrence, ternary lift, and extremal-run checks."""
from itertools import product

def g(x):
    h=1
    while x:
        x,b=divmod(x,2)
        if not b:
            h+=x+(x&1)
    return h

def ver(x,s,t):
    z=((3*x+1)*4**t-1)//3
    h0=g(2**s*x)
    A=(3*2**s-2)*(3*x+1)
    C=9*h0-A
    H=(A*4**t+6*t+C)//9
    assert (A*4**t+6*t+C)%9==0
    assert g(2**s*z)==H
    return H

def literal(P,n):
    R=[]
    for i in range(n):
        row=((1<<n)-1)^((2<<i)-1)
        if not (P>>(n-1-i))&1:row|=1<<n
        R.append(row)
    R.append(int(f'{P:0{n}b}'[::-1],2))
    E=[[0]*(n+1) for _ in range(1<<(n+1))]
    for i in range(n+1):E[1<<i][i]=1
    for S in range(1,1<<(n+1)):
        for j in range(n+1):
            if (S>>j)&1:
                E[S][j]+=sum(E[S^(1<<j)][i] for i in range(n+1) if (R[i]>>j)&1)
    return sum(E[-1])

case=0
for n in range(2,17):
    seen={}
    for x in range(1<<(n-1),1<<n,2):
        s=f'{x:b}';t=s.count('10')
        H=g(x)
        assert (H-1-t-(x%3==2))%3==0
        seen[t]=max(seen.get(t,0),H)
        if n<=7:assert literal(x,n)==H
        case+=1
    for t,H in seen.items():
        m=t-1;a=(n-2*m)//2;b=n-2*m-a
        Z=int('1'*a+'01'*m+'0'*b,2)
        assert g(Z)==H
        if m:
            D=(4**m-1)//3
            K=(4**m-1+6*m)//9
            assert H==K+D*(2**a+2**b-2)+4**m*(2**a-1)*(2**b-1)+2

for x,s in product(range(1,30,2),range(1,7)):
    h0=g(x<<s);A=(3*2**s-2)*(3*x+1);alpha=(A+2)//3
    images={(h0+alpha*t+t*(t-1)//2)%3 for t in range(3)}
    assert len(images)==2
    assert {ver(x,s,t)%3 for t in range(3) if (alpha-1+t)%3}==images
    for t,j,z in product(range(9),range(1,5),range(3)):
        assert (ver(x,s,t+z*3**j)-ver(x,s,t)-z*3**j*(alpha-1+t))%3**(j+1)==0

for b,j in product(range(1,7),range(1,5)):
    M=2**b*3**j;rows=set()
    period=2**(b-1)*3**j
    for s in [1,2]:
        rows.update(ver(1,s,t)%M for t in range((b+1)//2,(b+1)//2+period))
    assert rows==set(range(1,M,2))
print(f'PASS: {case} normalized addresses, literal controls, and joint CRT through 2^6·3^4')

for n in range(3,17):
    smallest=None
    for P in range(1<<(n-1),1<<n,2):
        H=g(P);u=f'{P:b}'
        largest=(1<<n)-H
        L=n-(largest).bit_length()
        if L>0:
            assert u.startswith('1'*L) and u.endswith('0'*L)
            Z=u[L:n-L];w=(1<<L)-1
            if Z:
                Pz=int(Z,2);Qz=(1<<len(Z))-1-int(Z[::-1],2)
                assert H==g(Pz)+w*(Pz+Qz)+w*w*2**len(Z)+1-int(Z[0])+int(Z[-1])
            else:
                assert H==w*w+2
        D=H-1-(1<<(n-2))
        assert D>0
        L=n-3-(D.bit_length()-1)
        if L>0:
            assert u.startswith('1'+'0'*L) and u.endswith('1'*L+'0')
            Z=u[1+L:n-1-L]
            Pz=int(Z,2) if Z else 0
            Qz=(1<<len(Z))-1-int(Z[::-1],2) if Z else 0
            assert H==(1<<(n-2))+((2+Pz+Qz)<<L)+g(Pz)
        if smallest is None or H<smallest:smallest=H
    if n>=4:assert smallest==(1<<(n-2))+(1<<((n-2)//2))+(1<<((n-1)//2))+1
print('PASS: both absolute middle decoders and exact lower edge through width 16')


# Exact lower-collar formula, including an empty middle.
def bitrev(z,l):
    return int(f'{z:0{l}b}'[::-1],2) if l else 0
for L,l in product(range(1,6),range(6)):
    for z in range(2**l):
        Z=f'{z:0{l}b}' if l else ''
        P=int('1'+'0'*L+Z+'1'*L+'0',2)
        n=2+2*L+l
        assert g(P)==2**(n-2)+(2+z+(2**l-1-bitrev(z,l)))*2**L+g(z)
# First possible nonexact dyadic matches for the fixed missing-on-ray target.
def Hr(t,A,C,M):
    r=(A*pow(4,t,9*M)+6*t+C)%(9*M)
    assert r%9==0
    return r//9

for b in range(10,15):
    M=2**b
    j=1 if b%2 else 2
    t=(278+j*(M//2))//3
    assert Hr(t,16,11,M)==63 and Hr(t+1,40,5,M)==63
    assert all(Hr(s,16,11,M)!=63 for s in range(t))
    assert all(Hr(s,40,5,M)!=63 for s in range(t+1))
print("PASS: lower collars and exact first hits of the 63 dyadic coset")

for n in range(4,17):
    assert min(g(x) for x in range(1<<(n-1),1<<n,2))==1+2**(n-2)+2**((n-2)//2)+2**((n-1)//2)
print("PASS: sharp lower height through width sixteen")

for L,l in product(range(1,6),range(6)):
    B,M=2**(2*L+l),2**L
    for w in range(2**l):
        W=f'{w:0{l}b}' if l else ''
        Z=''.join('0' if b=='1' else '1' for b in W)
        q=2+w+(2**l-1-bitrev(w,l))
        assert (l%2==0 or q%3==0)
        top=int('1'*L+Z+'0'*L,2)
        bot=int('1'+'0'*L+W+'1'*L+'0',2)
        assert g(top)==B-q*M+g(w)+2
        assert g(bot)==B+q*M+g(w)
print("PASS: exact reflected height charts and ternary carry")

for l in range(1,13):
    groups={}
    for w in range(2**l):
        delta=w-bitrev(w,l)
        eps=tuple(((w>>i)&1)-((w>>(l-1-i))&1) for i in range(l//2))
        assert groups.setdefault(delta,eps)==eps
print("PASS: coarse coefficient uniquely identifies asymmetric mirror bits")
