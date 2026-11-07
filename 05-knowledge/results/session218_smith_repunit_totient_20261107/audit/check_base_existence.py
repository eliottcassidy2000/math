"""Independent constructive all-base Smith check via smallest-prime-factor sieve."""
L=200000
spf=list(range(L+1))
for p in range(2,int(L**0.5)+1):
    if spf[p]==p:
        for n in range(p*p,L+1,p):
            if spf[n]==n: spf[n]=p

def factors(n):
    out=[]
    while n>1: out.append(spf[n]); n//=spf[n]
    return out

def ds(n,b):
    s=0
    while n: s+=n%b; n//=b
    return s

def smith(n,b): return ds(n,b)==sum(ds(p,b) for p in factors(n))
exceptions={6,8,9,12}
count=0
for n in range(4,L+1):
    if spf[n]==n: continue
    F=sum(factors(n))
    if n in exceptions:
        bases=[d+1 for d in range(1,n-F+1) if (n-F)%d==0]
        assert all(not smith(n,b) for b in bases),(n,bases)
    elif n==4:
        assert all(smith(n,b)==(b>=5) for b in range(2,101))
    else:
        b=n-F+1
        assert b>=2 and smith(n,b),(n,F,b)
        if len(factors(n))==2 and spf[n]==2:
            assert ds(n,b)==4
        else:
            assert b>n/2 and n//b==1 and n%b==F-1 and max(factors(n))<b
        count+=1
print('PASS constructive base n−F(n)+1 for',count,'composites<=200000; exception sieve6,8,9,12; special4.')
