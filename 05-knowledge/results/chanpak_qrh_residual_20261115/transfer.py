"""Exact modular address selectors (no conjectural assumptions)."""
def p(x):
    z=0
    while x:
        if x%2==0:
            z+=x//2+(x//2)%2
        x//=2
    return z

def ord2(m):
    t,s=1,2%m
    while s!=1:
        s=s*2%m;t+=1
    return t

def macros(m):
    ell=ord2(m)
    A='1'; Ai='1'*(ell-1)
    C='01'+Ai+Ai
    # Apply A^{-1}, C, A, then C^{-2}; C^m=Id modulo odd m.
    E=Ai+C+A+C*((-2)%m)
    return ell,A,C,E

def steer(x,y,m,a,w=None):
    """An odd x seed; optional even y seed for complementary extension."""
    assert m>=3 and m%2 and x%2==1
    assert y is None or y%2==0
    ell,A,C,E=macros(m)
    X=x%m
    W=(p(x)+(0 if y is None else p(y)))%m
    S=(x+(0 if y is None else y))%m
    at=0 if y is None else S*pow(2,-1,m)%m
    j=2*(X-at)%m
    k=2*(at-a)%m
    # The C trajectory is x -> x - 1/2, and (paired) y -> y + 1/2.
    def cshift(z,n):
        # Sum the initial affine increment along n steps.
        c0=z+1 if y is None else 2*z-S+3*pow(2,-1,m)
        step=pow(2,-1,m) if y is None else 1
        return n*c0-n*(n-1)//2*step
    inc=cshift(X,j)+cshift(at,k)
    v=0 if w is None else ((W+inc-w)*(1 if y is None else 2))%m
    suffix=C*j+E*v+C*k
    # All choices can be padded to the SAME exact suffix length.
    K=(2*(m-1)*len(C)+(m-1)*len(E)+ell-1)//ell*ell
    assert len(suffix)%ell==0 and (K-len(suffix))%ell==0
    return suffix+A*(K-len(suffix))

def append(x,y,s):
    for bit in s:
        b=int(bit)
        x=2*x+b
        if y is not None:y=2*y+(1-b)
    return x,y

def count(x,y,j):
    return 432*x*y+129*x+258*y+15*(x%2)+30*(y%2)+45*(p(x)+p(y))+(135,137,139)[j]

def count_selector(N,m,x=1):
    assert m>=3 and m%2
    j=next(j for j,k in enumerate((135,137,139)) if (N-303-k)%3==0)
    # At y=1 and x odd, H = 561x+45p(x)+303+k.
    # (561,45)/3 = (187,15), and 187*13-15*162=1.
    z=(N-303-(135,137,139)[j])//3
    s=steer(x,None,m,13*z,-162*z)
    X,_=append(x,None,s)
    assert count(X,1,j)%m==N%m
    return X,1,j

def suffix_length(m):
    ell,A,C,E=macros(m)
    return (2*(m-1)*len(C)+(m-1)*len(E)+ell-1)//ell*ell

def joint_selector(N,m,x,y,phase):
    """Adjust an existing complementary-phase point, preserving x+y."""
    from math import gcd
    assert m>=3 and gcd(m,30)==1 and (x+y)%phase==phase-1
    ellK=suffix_length(m)
    assert 1<<(ellK+1)<=phase
    u,v=append(x>>(ellK+1),y>>(ellK+1),"1")
    s=(x+y)%m
    a=x%m;b=(s-a)%m
    j=2*N%3
    const=432*a*b+129*a+258*b+15+(135,137,139)[j]
    w=(N-const)*pow(45,-1,m)%m
    tail=steer(u,v,m,a,w)
    assert len(tail)==ellK
    X,Y=append(u,v,tail)
    assert X+Y==x+y and count(X,Y,j)%(3*m)==N%(3*m)
    return X,Y,j
