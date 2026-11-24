from itertools import product
def rows(gaps):
    n=1+len(gaps)+sum(gaps)
    a=[((1<<n)-1)^((1<<(i+1))-1) for i in range(n)]
    i=0
    for g in gaps:
        j=i+1+g;a[i]^=1<<j;a[j]^=1<<i;i=j
    return a
def hp(a):
    n=len(a);D=[[0]*n for _ in range(1<<n)]
    for i in range(n):D[1<<i][i]=1
    for s in range(1,1<<n):
        for i,v in enumerate(D[s]):
            if not v:continue
            w=a[i]&~s
            while w:
                b=w&-w;w-=b;D[s|b][b.bit_length()-1]+=v
    return sum(D[-1])
def state(gaps):
    A,O,E,P=1,0,0,0
    for g in gaps:
        if not g:A,O,E,P=A+2*E+P,A+E,O,P
        else:A,O,E,P=A+2*E+P,0,0,(1<<g)*(A+O+E+P)
    return A+2*E+P
for m in range(7):
 for g in product(range(3),repeat=m):
    if 1+len(g)+sum(g)<=10:assert hp(rows(g))==state(g),(g,hp(rows(g)),state(g))
def K(t,a,b):
 z=0
 for x in product(range(2),repeat=t):
    s=(a,)+x+(b,);v=1;j=0
    while j<len(s):
      if not s[j]:j+=1;continue
      i=j
      while j<len(s) and s[j]:j+=1
      if i and j<len(s):v*=2 if (j-i)%2==0 else 0
    z+=v
 return z
assert [[K(3,a,b) for b in (0,1)] for a in (0,1)]==[[5,6],[6,7]]
print('verified exact strip height and unimodular bare port')
from fractions import Fraction
def mm(A,B,m=None):
    c=[[sum(A[i][k]*B[k][j] for k in (0,1)) for j in (0,1)] for i in (0,1)]
    return c if m is None else [[x%m for x in row] for row in c]
def inv(A,m=None):
    a,b=A[0];c,d=A[1];e=a*d-b*c
    if m is None:return [[Fraction(d,e),Fraction(-b,e)],[Fraction(-c,e),Fraction(a,e)]]
    return [[x*pow(e,-1,m)%m for x in row] for row in [[d,-b],[-c,a]]]
def chain(A,m=None):
    R=[[1,0],[0,1]]
    for V in A:R=mm(R,V,m)
    return R
D=[[1,0],[0,2]]
K1,K3,K4=[[[K(t,a,b) for b in (0,1)] for a in (0,1)] for t in (1,3,4)]
assert chain([inv(K4),K3,inv(K1),D,K1,D])==[[1,1],[0,1]]
for m in [3,5,7,9,11,15,27,77]:
    A,B,C,E=[mm(k,d,m) for k,d in [(K1,D),(K1,mm(D,D)),(K3,D),(K4,D)]]
    d=mm(inv(A,m),B,m)
    k1,k3,k4=[mm(z,inv(d,m),m) for z in [A,C,E]]
    assert chain([inv(k4,m),k3,inv(k1,m),d,k1,d],m)==[[1,1],[0,1]]
print('verified exact dyadic identities and modular transfers')

assert K4==[[9,11],[11,13]]
U=[[1,1],[0,1]]; L=[[1,0],[1,1]]
assert chain([inv(D),U,D,inv(U)])==U
assert chain([D,L,inv(D),inv(L)])==L

# Residue search on only actual moves, including widths large enough for
# reset at the two-primary clock. Each reached state is the literal row.
for a,m in [(1,3),(2,3),(3,3),(3,5),(4,9),(5,15),(6,21),(3,1)]:
    n=2**a*m
    g1=next(g for g in range(a,a+2*m) if pow(2,g,m)==2%m)
    g2=next(g for g in range(a,a+2*m) if pow(2,g,m)==4%m)
    gens=[mm(k,[[1,0],[0,pow(2,g,n)]],n)
          for k in (K1,K3,K4) for g in (g1,g2)]
    A,B=gens[-2]
    seed=(1,pow(2,g1,n))
    seed2=((A[0]+2*B[0])%n,(A[1]+2*B[1])%n)
    todo=[seed,seed2]; seen=set(todo); spectrum=set()
    for x,y in todo:
        spectrum.add((x+y)%n)
        if len(spectrum)==n//2: break
        for (A,B) in gens:
            z=((x*A[0]+y*B[0])%n,(x*A[1]+y*B[1])%n)
            if z not in seen:
                seen.add(z);todo.append(z)
    assert spectrum==set(range(1,n,2))
print("checked simultaneous mixed-modulus scalar spectra using only physical steps")
