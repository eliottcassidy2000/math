from fractions import Fraction as F
from itertools import product
from replay_carry_inverse import p

def dyadic(eta,c,r,X):
    q=1<<r
    j=(1-eta)*q
    if j.denominator!=1:raise ValueError(eta)
    j=int(j)
    l=1
    while q*(l+2+abs(c))>=2**(l-1):
        l+=1
    low=1<<(l-1)
    finite={a for a in range(1,min(low,X+1)) if p(a)-eta*a==c}
    rays=[]
    for v in (j,j-1):
        if not 0<=v<q:continue
        s=int(f'{v:0{r}b}'[::-1],2)
        delta=1-(s>>(r-1))
        b=p(s)-eta*s-c
        if v==j:
            if b.denominator!=1 or b<0:continue
            bounds=2*int(b)+2*delta+3
            zs=[z for z in [0,*range(1,bounds,2)] if z-p(z)-delta*(z%2)==b]
            for u in range(1,bounds,2):
                if u-p(u)-1==b:rays.append((q*u,s))
        else:
            b=-b
            if b.denominator!=1 or b<0:continue
            bounds=2*int(b)+1
            zs=[z for z in range(0,bounds,2) if p(z)==b]
            for u in range(0,bounds,2):
                if p(u)+delta==b:rays.append((q*(u+1),s-q))
        finite.update(a for z in zs if low<=(a:=q*z+s)<=X)
    for a,b in rays:
        v=a*2+b
        while v<=X:
            if v>=low:finite.add(v)
            v=2*(v-b)+b
    return finite,rays

def main():
    X=2048
    n=0
    for r in (1,2,3,4):
        q=1<<r
        for j,c in product(range(q+1),range(-7,8)):
            eta=1-F(j,q)
            got,_=dyadic(eta,c,r,X)
            want={a for a in range(1,X+1) if p(a)-eta*a==c}
            if got!=want:raise AssertionError((eta,c,sorted(got^want)[:10]))
            n+=1
    print('Exact dyadic line decompositions checked:',n)
    for eta,c in [(F(2,5),F(8,5)),(F(1,3),F(-2,3)),(F(1,2),-2)]:
        got=set()
        for L in range(1,X.bit_length()+1):
            u=1<<L
            ng=u*(1-eta)
            rad=2*(L+2+abs(c))
            for v in range(max(0,-(-(ng-rad).numerator//(ng-rad).denominator)),min(u-1,int(ng+rad))+1):
                a=int(f'{v:0{L}b}'[::-1],2)
                if u//2<=a<=X and p(a)-eta*a==c:got.add(a)
        want={a for a in range(1,X+1) if p(a)-eta*a==c}
        if got!=want:raise AssertionError((eta,c,got,want))
    print('Non-dyadic grid pruning and exact numerical reconstruction checked')
if __name__=='__main__':main()
