"""Exact overlap-kernel controls and a lawful cancellation-depth counterexample."""
from fractions import Fraction as F
from collections import Counter
import json


def need(ok, message):
    if not ok:
        raise ValueError(message)


def v2(n):
    need(type(n) is int and n != 0, 'nonzero exact integer')
    n = abs(n)
    return (n & -n).bit_length()-1


def step(n):
    need(type(n) is int and n % 2, 'odd exact source')
    a = v2(3*n+1)
    return (3*n+1)//(1 << a), a


def joint(depth, i, j):
    need(type(depth) is int and depth >= 1 and type(i) is int and i >= 1
         and type(j) is int and j >= 1, 'positive exact indices')
    if i < depth:
        return F(1, 1 << i) if j == i else F(0)
    if i > depth:
        return F(1, 1 << i) if j == depth else F(0)
    return F(1, 1 << j) if j > depth else F(0)


def kernel(depth):
    need(type(depth) is int and depth >= 1, 'positive exact depth')
    return 2-F(6, 1 << depth)


def main():
    checks = 0
    def check(ok, message):
        nonlocal checks
        checks += 1
        if not ok: raise RuntimeError(message)

    # Exact finite heads plus analytic geometric tails: no tail is discarded.
    for m in range(1, 41):
        cross = sum((F(i*i, 1 << i) for i in range(1,m)), F(0))
        cross += 2*m*F(m+2, 1 << m)
        info = sum((F(i, 1 << i) for i in range(1,m)), F(0)) + F(2*m, 1 << m)
        agree = sum((F(1, 1 << i) for i in range(1,m)), F(0))
        check(cross-4 == kernel(m), 'covariance kernel')
        check(info == 2*(1-F(1,1 << m)), 'conditional information')
        check(agree == 1-F(2,1 << m), 'agreement probability')
        for i in range(1, 16):
            row = sum((joint(m,i,j) for j in range(1,m+1)), F(0))
            if i == m: row += F(1,1 << m)
            check(row == F(1,1 << i), 'exact geometric marginal')

    # Independent literal affine pair y, y+2^m, all odd residues at precision10.
    bits = 10
    for m in range(1,9):
        counts = Counter()
        for y in range(1,1 << (bits+1),2):
            _,a=step(y); _,b=step(y+(1 << m))
            counts[min(a,bits+1),min(b,bits+1)] += 1
        for i in range(1,bits+2):
            for j in range(1,bits+2):
                if i<=bits and j<=bits: expected=joint(m,i,j)
                elif i==bits+1 and j==m: expected=F(1,1 << bits)
                elif j==bits+1 and i==m: expected=F(1,1 << bits)
                else: expected=F(0)
                check(F(counts[i,j],1 << bits)==expected,'literal two-valuation law')

    weights={1:F(1,3),2:F(2,3)}
    moment=sum((w*F(1,1 << m) for m,w in weights.items()),F(0))
    check(moment==F(1,3),'non-Haar depth moment')
    cov=sum((w*kernel(m) for m,w in weights.items()),F(0))
    agree=sum((w*(1-F(2,1 << m)) for m,w in weights.items()),F(0))
    info=sum((w*2*(1-F(1,1 << m)) for m,w in weights.items()),F(0))
    check((cov,agree,info)==(0,F(1,3),F(4,3)),'three redundant scalar summaries')
    p11=sum((w*joint(m,1,1) for m,w in weights.items()),F(0))
    p22=sum((w*joint(m,2,2) for m,w in weights.items()),F(0))
    check(p11==F(1,3) and p11!=F(1,4) and p22==0,'dependence remains')
    for j in range(1,10):
        diagonal=sum((w*joint(m,j,j) for m,w in weights.items()),F(0))
        check(diagonal==F(1,1 << j)*sum(w for m,w in weights.items() if m>j),
              'full diagonal profile retains the depth tail')
    for i in range(1,10):
        for j in range(1,10):
            if i==j: mixture=F(1,1 << (2*i))
            else:
                m=min(i,j)
                mixture=F(1,1 << m)*joint(m,i,j)
            check(mixture==F(1,1 << (i+j)),'geometric depth yields independent pair')

    raw=19; v=2
    y,a0=step(raw); x=3*(1 << v)*raw+1
    ys=[y]; xs=[x]; ay=[]; bx=[]
    for _ in range(4):
        y,a=step(y); x,b=step(x)
        ys.append(y); xs.append(x); ay.append(a); bx.append(b)
    s,t=2,3
    check(sum(bx[:t])==v+a0+sum(ay[:s]),'lawful Mersenne alignment')
    k=t-s; anchor=xs[t]-3**k*ys[s]
    e=3*anchor+1-3**k; delta=v2(e); nu=v2(3**k-1); a=ay[s]
    check(a<delta,'lawful lockstep guard')
    e_next=F(3*e,1 << a)-(3**k-1)
    check(e_next.denominator==1 and v2(e_next.numerator)>nu,'tie exceeds alleged cap')
    check(e_next == 3*(xs[t+1]-3**k*ys[s+1])+1-3**k,
          'independent literal next alignment')
    cap_hostile={'raw_source':raw,'v':v,'Y':ys,'X':xs,'s':s,'t':t,
                 'E':e,'a':a,'delta':delta,'nu':nu,'E_next':str(e_next)}
    # Algebraic tie depths can be arbitrarily large with fixed k=1,a=1:
    # E=2*(2+2^h)/3 is a dyadic integer for even h; E'=2^h.
    for h in range(2,22,2):
        e=F(2*(2+(1 << h)),3)
        check(e.denominator==1 and v2(e.numerator)==2,'tie source depth')
        check(F(3*e,2)-2 == 1 << h,'unbounded cancellation in tie')

    for k in range(1,17):
        for raw_y in (3,7,11,15):
            paired_y, first=step(raw_y)
            x_start=3*(1 << (4*k))*raw_y+1
            cur=x_start; vals=[]
            for _ in range(2*k):
                cur, exponent=step(cur); vals.append(exponent)
            check(first==1 and vals==[2]*(2*k-1)+[3] and sum(vals)==4*k+1,
                  'unbounded aligned gaps on a positive-probability source cylinder')
            kappa=cur-3**(2*k)*paired_y
            depth=v2(3*kappa+1-3**(2*k))
            check(kappa==F(1-3**(2*k),2) and depth==2+v2(k),
                  'past-determined depth on unbounded aligned gaps')

    for k in range(1,17):
        nu=3+v2(k); depths=Counter()
        for raw_y in range(1,1 << (nu+1),2):
            paired_y,a0=step(raw_y)
            cur=3*(1 << (4*k))*raw_y+1; total=0
            for _ in range(2*k):
                cur,a=step(cur); total+=a
            aligned=total==4*k+a0
            check(aligned==(a0<nu),'exact coarse initial-alignment event')
            if aligned:
                kappa=cur-3**(2*k)*paired_y
                depth=v2(3*kappa+1-3**(2*k))
                check(depth==nu-a0,'actual non-geometric alignment depth')
                depths[depth]+=1
        check(depths=={d:1 << d for d in range(1,nu)},'complete coarse depth census')

    print(json.dumps({'status':'PASS; depth law and return control remain OPEN',
        'checks':checks,'literal_precision_bits':bits,'depth_mixture':{'1':'1/3','2':'2/3'},
        'covariance':str(cov),'agreement':str(agree),'conditional_information_bits':str(info),
        'dependent_p11':str(p11),'dependent_p22':str(p22),
        'cap_hostile':cap_hostile},indent=2,sort_keys=True))


if __name__=='__main__': main()
