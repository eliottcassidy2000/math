from math import gcd
def v(n,p):
    assert n
    k=0
    while n%p==0:
        n//=p
        k+=1
    return k
def primes(n):
    return [p for p in range(3,n,2) if all(p%a for a in range(2,int(p**.5)+1))]
def act(x):
    j=v(3*x+1,2)
    return (3*x+1)>>j,j
for r in range(1,501):
    n=2*r+1; m=2**n-1; E=2**(2*r)*m
    x,y=E-1,E+1
    for i in range(2*r):
        x,j=act(x)
        assert j==(2 if i==2*r-1 else 1)
    assert x==(3**(2*r)*m-1)//2
    for i in range(r):
        y,j=act(y)
        assert j==(2+v(3**r-1,2) if i==r-1 else 2)
    assert y==(3**r*m+1)>>v(3**r-1,2)
    assert v(E-1,3)==2+v(r,3)+v(n,3)
    assert v(E+1,7)==(n%6==5)
    assert v(E+1,1093)==(n%364==149)
for p in primes(80):
    if p==3:continue
    q=p*p
    fq=lambda a:(pow(a,p-1,q)-1)//p
    for j in range(1,12):
        for a in range(p):
            b=((3*a+1)*pow(2,-j,p))%p
            l=(3*pow(a,p,q)+1-pow(2,j,q)*pow(b,p,q))%q
            assert l%p==0
            carry=(pow(1+3*a,p,q)-1-pow(3*a,p,q))%q//p
            assert l//p==(j*(3*a+1)*fq(2)-3*a*fq(3)-carry)%p
for p,d,expected in [(1093,182,[(47,74,389),(500,None,650)]),(3511,1755,[(1815,None,1803),(3452,None,661)])]:
    q=p*p
    assert pow(2,p-1,q)==1 and pow(2,p-1,q*p)!=1
    h={pow(4,i,p):i for i in range(d)}
    assert len(h)==d and pow(4,d,p)==1
    assert [(a,h.get(a),((pow(a,p,q)+pow(1-a,p,q)-1)%q)//p)
            for a in range(p) if (2*a*a-a+1)%p==0]==expected
print("500 exact rail/exit checks; all affine defects at 11 guards and primes 5..79; both Wieferich subgroup/carry checks pass")
