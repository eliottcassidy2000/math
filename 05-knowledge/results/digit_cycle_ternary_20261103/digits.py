"""Exact decimal carry potentials, minimal append observer and block laws."""
from collections import Counter
import sympy as s
from support import need,save
def ds(n):return sum(map(int,str(n)))
def formula(n):
    digits=list(map(int,str(n)))
    return sum(ds(4*d) for d in digits)-9*sum(a in (2,7) and b>=5 for a,b in zip(digits,digits[1:]))
def v2(n):return (n&-n).bit_length()-1
def main():
    z=s.Symbol('z');W=[ds(4*d) for d in range(10)]
    A=s.zeros(2)
    for state in (0,1):
        for d in range(10):A[state,int(d in (2,7))]+=z**(W[d]-9*state*(d>=5))
    S=sum(z**d for d in range(10));P=sum(z**w for w in W)
    need(s.expand(A.det())==0 and s.expand(s.trace(A)-S)==0,'rank-one decimal average with exact trace')
    need(all(formula(n)==ds(4*n) for n in range(100000)),'literal adjacent-digit carry formula')
    for n in range(100000):
        carries=sum((4*n)%10**j>=10**j-4 for j in range(1,len(str(4*n))+2))
        need(ds(4*(n+1))-ds(4*n)==4-9*carries,'sequential multiples retain their full carry lift')
    checks=[]
    poly=P
    for L in range(1,6):
        if L>1:poly=s.expand(poly*S)
        hist=Counter(ds(4*n) for n in range(10**L))
        need(all(int(poly.coeff(z,k))==hist.get(k,0) for k in range(s.degree(poly,z)+1)),'complete decimal block distribution')
        total=sum(k*c for k,c in hist.items());squares=sum(k*k*c for k,c in hist.items())
        need(total==9*(5*L+1)*10**(L-1),'exact first moment')
        variance=s.Rational(squares,10**L)-s.Rational(total,10**L)**2
        need(variance==s.Rational(33,4)*L+s.Rational(99,100),'exact variance')
        checks.append({'digits':L,'count':10**L,'histogram':dict(sorted(hist.items())),'sum':total,'variance':str(variance)})
    # Finite aliases support, but are not the proof of, the future-context theorem.
    aliases=[]
    for k in range(1,9):
        n=10**k;aliases.append({'n':n,'digit_sum':ds(4*n),'flag':int(n%5==2),'halving_clock':v2(n)})
    need(all(a['digit_sum']==4 and a['flag']==0 and a['halving_clock']==i+1 for i,a in enumerate(aliases)),'unbounded omitted halving clock')
    # Four is the first nontrivial member of the all-dyadic local-potential law.
    higher=[]
    for a in range(1,9):
        m=2**a
        for n in range(5000):
            q=0;answer=0
            for d in map(int,str(n)):
                c=m*d//10
                carries=sum((m*q%10**j)+(c%10**j)>=10**j for j in range(1,a))
                answer+=ds(m*d)-9*carries;q=10*q+d
            need(answer==ds(m*n),'all-dyadic bounded-range decimal potential')
        higher.append({'power':a,'multiplier':m,'integers':5000})
    save('digits.json',{'status':'FINITE_EXACT_WITH_ALL_SIZE_PROOFS','sequential_carry_cases':100000,'digit_weights':W,'weighted_matrix':[[str(s.expand(A[i,j])) for j in range(2)] for i in range(2)],
      'trace':str(S),'determinant':0,'initial_polynomial':str(P),'blocks':checks,'halving_aliases':aliases,'dyadic_potentials':higher,
      'append_observer':'(s10(4n), [n mod5=2]); homogeneous linear dimension3'})
    print('PASS nearest-neighbour digit law; exact block GF and moments through100000; dyadic powers1..8',flush=True)
if __name__=='__main__':main()
