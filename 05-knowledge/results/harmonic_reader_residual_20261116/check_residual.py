from fractions import Fraction
from itertools import combinations
for n in range(1,15):
    D=1<<n; peak=Fraction(-1); hits=[]
    for b in range(D):
        word=f'{b:0{n}b}'
        R=int(word[::-1],2)
        ix=[i for i in range(n) if b>>i&1]
        J=sum((1<<(j-i-1))+(j==i+1) for i,j in combinations(ix,2))
        E=b*R-D*J
        assert E>=0 and bool(E)==bool(b)
        if E>peak:peak=E;hits=[word]
        elif E==peak:hits.append(word)
    k=(n+1)//2
    assert Fraction(peak,D)==Fraction(2*k,3)-Fraction(2,9)*(1-Fraction(1,4**k))
    assert len(hits)==1+(n%2==0)
    assert all(all(u!=v for u,v in zip(w,w[1:])) and w.count('1')==k for w in hits)
print("Sharp signed residual and equality masks checked through n=14.")
for n in range(1,11):
    D=1<<n
    total=0
    for b in range(D):
        ix=[i for i in range(n) if b>>i&1]
        J=sum((1<<(j-i-1))+(j==i+1) for i,j in combinations(ix,2))
        total+=b*int(f'{b:0{n}b}'[::-1],2)-D*J
    assert Fraction(total,D*D)==Fraction(n,8)+Fraction(1,4*D)
print("Exact mean signed residual checked through n=10.")
