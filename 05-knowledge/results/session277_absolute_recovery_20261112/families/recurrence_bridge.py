#!/usr/bin/env python3
"""Exact run/Stirling compiler and certified pole recovery, stdlib only."""
import json, math, itertools, sys
if hasattr(sys,"set_int_max_str_digits"):sys.set_int_max_str_digits(0)
from fractions import Fraction as F
from decimal import Decimal as D, localcontext, ROUND_FLOOR, ROUND_CEILING
from pathlib import Path
PREC=100

def op(a,b,kind,rounding):
    with localcontext() as c:
        c.prec=PREC;c.rounding=rounding
        return {'+':lambda:a+b,'-':lambda:a-b,'*':lambda:a*b,'/':lambda:a/b}[kind]()
class I:
    def __init__(self,lo,hi=None): self.lo=D(lo);self.hi=D(lo if hi is None else hi)
    @staticmethod
    def frac(q): return I(op(D(q.numerator),D(q.denominator),'/',ROUND_FLOOR),op(D(q.numerator),D(q.denominator),'/',ROUND_CEILING))
    def __add__(a,b):
        if not isinstance(b,I):b=I(b)
        return I(op(a.lo,b.lo,'+',ROUND_FLOOR),op(a.hi,b.hi,'+',ROUND_CEILING))
    __radd__=__add__
    def __neg__(a):return I(a.hi.copy_negate(),a.lo.copy_negate())
    def __sub__(a,b):return a+-b
    def __mul__(a,b):
        if not isinstance(b,I):b=I(b)
        return I(min(op(x,y,'*',ROUND_FLOOR) for x in [a.lo,a.hi] for y in [b.lo,b.hi]),max(op(x,y,'*',ROUND_CEILING) for x in [a.lo,a.hi] for y in [b.lo,b.hi]))
    __rmul__=__mul__
    def inv(a):
        assert not a.lo<=0<=a.hi
        return I(op(D(1),a.hi,'/',ROUND_FLOOR),op(D(1),a.lo,'/',ROUND_CEILING))
    def __truediv__(a,b):return a*b.inv()
    def __pow__(a,n):
        b=I(1)
        while n:
            if n&1:b=b*a
            a=a*a;n//=2
        return b

def atan_bounds(q,terms):
    s=sum(((-1)**j*q**(2*j+1)/F(2*j+1) for j in range(terms)),F(0))
    nxt=(-1)**terms*q**(2*terms+1)/F(2*terms+1)
    return min(s,s+nxt),max(s,s+nxt)
def constants(d,digits):
    z=F(1,2*d+1);terms=max(10,math.ceil(digits/math.log10((2*d+1)**2))+3)
    t=2*sum((z**(2*j+1)/F(2*j+1) for j in range(terms)),F(0))
    te=2*z**(2*terms+1)/((2*terms+1)*(1-z*z))
    tau=I(I.frac(t).lo,I.frac(t+te).hi)
    a,b=atan_bounds(F(1,5),math.ceil(digits/math.log10(25))+3)
    c,e=atan_bounds(F(1,239),math.ceil(digits/math.log10(239**2))+3)
    pi=I(I.frac(16*a-4*e).lo,I.frac(16*b-4*c).hi)
    return tau,pi

def cmul(z,w):return (z[0]*w[0]-z[1]*w[1],z[0]*w[1]+z[1]*w[0])
def cpow(z,n):
    r=(I(1),I(0))
    while n:
        if n&1:r=cmul(r,z)
        z=cmul(z,z);n//=2
    return r

def regular_exact(q,n):
    d=(q-1)//2;s=[1]+[0]*n
    for m in range(1,n+1):
        s=[0]+[k*s[k]+s[k-1] for k in range(1,n+1)]
    return sum(q*d**(k-1)*math.factorial(k)*s[k] for k in range(1,n+1)) if n else 1

def recover(q,n,extra_precision=50):
    """Certified interval whose only integer is B_n. Not A000568."""
    global PREC
    assert q>=3 and q%2 and n>=1
    d=(q-1)//2
    # A priori B_n <= q^n n!: no exact-value oracle chooses precision.
    PREC=len(str(math.factorial(n)))+n*len(str(q))+extra_precision
    tau,pi=constants(d,PREC+10)
    factor=I.frac(F(q*math.factorial(n),d*(d+1)))
    head=tau.inv()**(n+1)
    dominant=factor*head
    k=0
    while True:
        k+=1;y=pi*(2*k);den=tau*tau+y*y
        term=cpow((tau/den,-y/den),n+1)[0]
        head=head+term*2
        # Integral bound on the omitted symmetric translates.
        tail=factor*2/(I(n)*(pi*2)**(n+1)*I(k)**n)
        val=factor*head
        lo=op(val.lo,tail.hi,'-',ROUND_FLOOR);hi=op(val.hi,tail.hi,'+',ROUND_CEILING)
        low=int(lo.to_integral_value(rounding=ROUND_CEILING));high=int(hi.to_integral_value(rounding=ROUND_FLOOR))
        if low==high:
            exact=regular_exact(q,n)  # Independent control after certification.
            assert low==exact
            with localcontext() as c:
                c.prec=PREC
                mid=(dominant.lo+dominant.hi)/2
            with localcontext() as c:
                c.prec=16
                err=+(mid-D(exact))
            return {'q':q,'n':n,'exact':str(exact),'pole_pairs':k,'decimal_precision':PREC,'interval_lo':str(lo),'interval_hi':str(hi),'dominant_error_estimate':str(err),'error_definition':'dominant minus exact; display estimate only; interval recovery is directed and certified','tail_bound':str(tail.hi)}
        if tail.hi < D("0.125"):
            # No fixed arithmetic error floor or arbitrary pole-count cutoff.
            return recover(q,n,extra_precision*2)

def H(rows):
    n=len(rows);dp=[[0]*n for _ in range(1<<n)]
    for v in range(n):dp[1<<v][v]=1
    for mask in range(1,1<<n):
        for v in range(n):
            if not dp[mask][v]:continue
            opts=rows[v]&~mask
            while opts:
                bit=opts&-opts;opts-=bit;w=bit.bit_length()-1
                dp[mask|bit][w]+=dp[mask][v]
    return sum(dp[-1]) if n else 1

def inflate(A,sizes):
    labels=[i for i,s in enumerate(sizes) for _ in range(s)]
    return [sum(1<<j for j,b in enumerate(labels) if (a==b and i<j) or (a!=b and A[a][b])) for i,a in enumerate(labels)]
def stirling(n,k):
    s=[1]+[0]*k
    for _ in range(n):s=[0]+[j*s[j]+s[j-1] for j in range(1,k+1)]
    return s[k]
def compiled_H(A,sizes):
    """Run counts bounded by sizes; exact coefficient in rational walk series."""
    q=len(A);counts={};target=tuple(sizes)
    states={((0,)*q,-1):1}
    total=1 if not any(sizes) else 0
    for length in range(1,sum(sizes)+1):
        nxt={}
        for (ks,last),v in states.items():
            for j in range(q):
                if ks[j]>=sizes[j] or (last>=0 and not A[last][j]):continue
                kk=list(ks);kk[j]+=1;kk=tuple(kk)
                nxt[kk,j]=nxt.get((kk,j),0)+v
        for (ks,last),v in nxt.items():
            total+=v*math.prod(math.factorial(k)*stirling(n,k) for n,k in zip(sizes,ks))
        states=nxt
    return total

def compositions(n,q):
    if q==1:yield (n,);return
    for j in range(n+1):
        for rest in compositions(n-j,q-1):yield (j,)+rest

def inverse_target(A,c):
    """Complete positive transitive-inflation inverse for nontrivial strong A."""
    q=len(A)
    assert q>=3
    # Verify strong connectivity rather than silently accepting a bad bound.
    for root in range(q):
        seen={root};todo=[root]
        while todo:
            v=todo.pop()
            for j in range(q):
                if A[v][j] and j not in seen:seen.add(j);todo.append(j)
        assert len(seen)==q
    if c<=1:return []
    bound=(c-1).bit_length()-1
    hits=[]
    for sizes in itertools.product(range(1,bound+1),repeat=q):
        if compiled_H(A,sizes)==c:
            rows=inflate(A,sizes)
            assert H(rows)==c
            hits.append({'sizes':sizes,'rows':rows})
    return hits

def test():
    tested=0
    for q in range(1,5):
        pairs=list(itertools.combinations(range(q),2))
        for bits in range(1<<len(pairs)):
            A=[[0]*q for _ in range(q)]
            for b,(i,j) in enumerate(pairs):A[i][j]=(bits>>b)&1;A[j][i]=1-A[i][j]
            for sizes in itertools.product(range(3),repeat=q):
                if sum(sizes)>7:continue
                assert compiled_H(A,sizes)==H(inflate(A,sizes)),(A,sizes)
                tested+=1
    C=[[0,1,0],[0,0,1],[1,0,0]]
    for n in range(8):
        b=sum(math.factorial(n)//math.prod(math.factorial(x) for x in sizes)*H(inflate(C,sizes)) for sizes in compositions(n,3))
        assert b==regular_exact(3,n)
    targets={str(c):inverse_target(C,c) for c in [7,21,23,25,27,29,31,33,45,47]}
    Path(__file__).with_name('TARGETS.json').write_text(json.dumps(targets,indent=2)+'\n')
    for a in range(1,5):
        for b in range(1,5):assert compiled_H(C,(a,b,1))==2+(2**a-1)*(2**b-1)
    certificates=[recover(q,n) for q in [3,5,7] for n in [10,20,30,50,100,250]]
    return {'direct_inflation_checks':tested,'aggregate_checks':8,'certificates':certificates}
if __name__=='__main__':
    result=test();out=Path(__file__).with_name('RESULTS.json');out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'direct_inflation_checks':result['direct_inflation_checks'],'aggregate_checks':8,'certified_recoveries':len(result['certificates']),'output':str(out)}))
