from collections import Counter
from functools import cache
from hashlib import sha1
from math import gcd
from pathlib import Path

def check(ok):
    if not ok:
        raise RuntimeError("control failed")

@cache
def p(x):
    if x<2: return 0
    y,b=divmod(x,2)
    return p(y)+(y+y%2)*(1-b)

def rows():
    local=Path(__file__).parent/'forms_orders3_5.txt'
    if not local.exists():
        local=Path(__file__).parent/'../chanpak_guard_depth_20261115/forms_orders3_5.txt'
    raw=local.read_bytes()
    header=f'blob {len(raw)}\0'.encode()
    check(sha1(header+raw).hexdigest()=='f58b1e2679edd0740a12247acfef423af1170b43')
    return [tuple(map(int,line.split('|')[0].split())) for line in raw.decode().splitlines() if line and not line.startswith('#')]

def base(c,x,y):
    a,Q,s,g,t,f,n=c
    return Q*x*y+t*x+s*y+n+f*(x%2)+g*(y%2)

def value(c,x,y):
    return base(c,x,y)+c[0]*(p(x)+p(y))

def main():
    bank=rows()
    check(len(set(bank))==168)
    for c in bank:
        a,Q,s,g,t,f,n=c
        check(a>0 and Q>=4*a and min(c)>=0 and f<=a and g<=a)
        check(Q%2==0 and (s+g)%2==0 and (t+f)%2==0 and n%2==1)
        for x in range(1,49):
            for y in range(1,x+1):
                check(value(c,x,y+1)>value(c,x,y))
                check(value(c,y+1,x)>value(c,y,x))
    print("168 coefficient and parity rows; both shorter-coordinate cuts")
    cases=0
    for r in range(1,11):
        V=1<<r
        for b in range(V):
            j=int(f'{b:0{r}b}'[::-1],2)
            for u in (0,1,2,3,5,32,997):
                check(p(V*u+b)==p(u)+p(b)+(V-1-j)*u+(b<V//2)*(u%2))
                cases+=1
    print(f"{cases} split controls")
    c=(11,100,0,0,89,11,11)
    check(c in bank)
    for M in (1,2,3,5,11,22,33,121):
        d=gcd(11,M)
        hist=Counter(base(c,x,y)%(2*d) for x in range(2*d) for y in range(2*d))
        for b in range(1,2*M,2):
            check(hist[b%(2*d)]==(4 if d==1 else 4*(21 if b%11==0 else 10)))
    print("single-carrier multiplication law checked at eight moduli")

if __name__=='__main__':
    main()
