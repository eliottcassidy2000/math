from fractions import Fraction as F
from itertools import product

def stats(w):
    l=len(w)
    A=2**l
    B=int(''.join(map(str,w)),2)
    C=sum(2**i for i,b in enumerate(w) if b==0)
    theta=[F(int(''.join(map(str,w[i:]+w[:i])),2),A-1) for i in range(l)]
    mu=sum((w[i-1]-theta[i] for i,b in enumerate(w) if b==0),F(0))
    E=sum(w[i]==w[i-1] for i in range(l))
    return A,B,C,mu,E

def transfer(h,N,S,d,w):
    for b in w:
        N,S,d=N+(1-b)*(S+d),2*S+b*h,b*h
    return N,S,d

def main():
    z=0
    for l in range(1,13):
        lo=F((2**(l//2)+2**((l+1)//2)-2),2**l-1)
        stats_l=[]
        for w in product(range(2),repeat=l):
            A,B,C,mu,E=stats(w)
            h,N,S,d=15,163,242,0
            T=F(C)*(S+F(h*B,A-1))/(A-1)
            for k in (0,1,2,5):
                val=transfer(h,N,S,d,w*k)[0]
                expect=N+T*(A**k-1)+h*k*mu+(d-h*w[-1] if k and not w[0] else 0)
                if val!=expect:raise AssertionError((w,k,val,expect))
                z+=1
            if len(set(w))==2:
                if not 0<mu<=F(l,3)-F(2*E,9):raise AssertionError((w,mu,E))
                balanced=(E==l-2 and abs(2*sum(w)-l)<=1)
                if (mu==lo)!=balanced or mu<lo:raise AssertionError(('min',w,mu,lo))
                stats_l.append((mu,E))
            elif mu:raise AssertionError(w)
        if stats_l:
            m=max(mu for mu,_ in stats_l)
            if any((mu==m)!=(E==l%2) for mu,E in stats_l):raise AssertionError(('max',l,m))
            print('length',l,'exact drift extrema',lo,m)
    print('Literal repeated-block transfer arithmetic checks:',z)
if __name__=='__main__':main()
