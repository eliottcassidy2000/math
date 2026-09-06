"""Independent exact polynomial audit; standard library only, no interpolation.

Derive all factorial contents and all three characteristic coefficients in
Q[g], then compare the producer's published rational coefficient arrays.
No producer/root-audit code imports; no numerical roots or parameter census.
"""
from fractions import Fraction as F
from math import factorial, comb
from pathlib import Path
from itertools import permutations
import hashlib
import json
import sys

sys.stdout.reconfigure(newline="\n")
GATES = 0


def need(ok, label):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(label)


def trim(a):
    a=list(map(F,a))
    while len(a)>1 and a[-1]==0:
        a.pop()
    return tuple(a)


ZERO=(F(0),)
ONE=(F(1),)


def add(a,b):
    c=[F(0)]*max(len(a),len(b))
    for i,x in enumerate(a): c[i]+=x
    for i,x in enumerate(b): c[i]+=x
    return trim(c)


def scale(a,x):
    return trim([F(x)*v for v in a])


def sub(a,b):
    return add(a,scale(b,-1))


def mul(a,b):
    c=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): c[i+j]+=x*y
    return trim(c)


def prods(items):
    ans=ONE
    for item in items: ans=mul(ans,item)
    return ans


def div_exact(a,b):
    a=list(a)
    q=[F(0)]*max(1,len(a)-len(b)+1)
    while len(a)>=len(b) and trim(a)!=ZERO:
        k=len(a)-len(b); value=a[-1]/b[-1];q[k]+=value
        for j,v in enumerate(b): a[j+k]-=value*v
        a=list(trim(a))
    need(trim(a)==ZERO,"exact polynomial content division")
    return trim(q)


def fall(lead,constant,n):
    return prods(((F(constant-j),F(lead)) for j in range(n)))


def at(a,g):
    value=F(0)
    for x in reversed(a): value=value*g+x
    return value


def shift(a,c):
    result=[F(0)]*len(a)
    for j,v in enumerate(a):
        for i in range(j+1): result[i]+=v*comb(j,i)*c**(j-i)
    return trim(result)


def rem_t(a,p):
    a=list(a)
    while len(a)>3:
        h=a.pop(); k=len(a)-3
        for j in range(3): a[k+j]=sub(a[k+j],mul(h,p[j]))
    return a+[ZERO]*(3-len(a))


def raw_row(g,m):
    ans=[]
    for x in range(m+1):
        for y in range(m-x+1):
            z=m-x-y
            if -21*x+(2*g-21)*y+(3*g-21)*z==0:
                ans.append(((x,y,z),factorial(m)//(factorial(x)*factorial(y)*factorial(z))))
    return sorted(ans,key=lambda entry:entry[0][2])


def main():
    base=Path(__file__).parent
    record=base/"overnight6_20260906_laurent_cubic_carry.out"
    if not record.exists():
        record=base.parent/"05-knowledge/results/overnight6_20260906_laurent_cubic_carry.out"
    line=next(line for line in record.read_text().splitlines()
              if line.startswith("characteristic_positive_shift_certificates "))
    records=json.loads(line.split(" ",1)[1])
    first=[scale(fall(1,0,10-j),F(1,factorial(9-3*j)*factorial(1+2*j))) for j in range(4)]
    content=scale(fall(1,0,7),F(1,factorial(9)))
    integer_p=[div_exact(a,content) for a in first]
    expected=[prods(((-7,1),(-8,1),(-9,1))),
              scale(prods(((-7,1),(-8,1))),84),scale((-7,1),504),(F(72),)]
    need(integer_p==expected,"complete first factorial content")
    p=[div_exact(a,first[-1]) for a in first]
    need(p[-1]==ONE,"independent monic first row")
    need(all(len(p[j])-1==3-j for j in range(4)),"monic row weighted-degree premises")
    rawq=[scale(fall(2,0,21-j),F(1,factorial(21-3*j)*factorial(2*j))) for j in range(8)]
    kfactor=fall(2,0,14)
    q=[div_exact(a,kfactor) for a in rawq]
    for j in range(8):
        need(q[j]==scale(fall(2,-14,7-j),F(1,factorial(21-3*j)*factorial(2*j))),
             "complete second row positive content")
    ratio=div_exact(q[0],p[0])
    need(ratio==scale(prods(((-10,1),(-15,2),(-17,2),(-19,2))),F(1152,factorial(21))),
         "inverse denominator cancels before reduction")
    # q0/t = -(q0/p0)*(t^2+p2*t+p1); remaining powers have no inverse.
    carried=list(q[1:])
    for j,coef in enumerate((p[1],p[2],ONE)):
        carried[j]=sub(carried[j],mul(ratio,coef))
    r=rem_t(carried,p)
    need(all(len(r[j])-1<=6-j for j in range(3)),"three exact response degree bounds")
    cols=[rem_t([ZERO]*j+r,p) for j in range(3)]
    v=[[cols[j][i] for j in range(3)] for i in range(3)]
    b1=scale(add(add(v[0][0],v[1][1]),v[2][2]),-1)
    b2=ZERO
    for i in range(3):
        for j in range(i+1,3):
            b2=add(b2,sub(mul(v[i][i],v[j][j]),mul(v[i][j],v[j][i])))
    determinant=ZERO
    for pp in permutations(range(3)):
        sign=(-1)**sum(pp[i]>pp[j] for i in range(3) for j in range(i+1,3))
        determinant=add(determinant,scale(prods(v[i][pp[i]] for i in range(3)),sign))
    bs=(b1,b2,scale(determinant,-1))
    for i,b in enumerate(bs,1):
        need(len(b)-1==6*i,"exact characteristic degree")
        shifted=shift(b,11)
        need(all(c>0 for c in shifted),"strict positive full shifted coefficient polynomial")
        rec=records[str(i)]
        expected=scale(tuple(F(c) for c in reversed(rec["coefficients_descending"])),F(1,int(rec["denominator"])))
        need(rec["shift_g"]==11 and shifted==expected,"every published rational coefficient independently derived")
    # Cubic discriminant from the elementary five-term formula.
    b,c,d=p[2],p[1],p[0]
    disc=add(sub(sub(mul(mul(b,b),mul(c,c)),scale(prods((c,c,c)),4)),
                 scale(prods((b,b,b,d)),4)),
             sub(scale(prods((b,c,d)),18),scale(mul(d,d),27)))
    expected_disc=scale(prods(((-8,1),(-7,1),(-7,1),(-27511000,11532575,-1610102,74863))),F(1,1728))
    need(disc==expected_disc,"independent five-term cubic discriminant")
    need(all(c>0 for c in shift(disc,11)),"strict discriminant on the complete real parameter ray")
    for g in (11,13,16,17,19):
        for mass,offset,formula in ((g,1,first),(2*g,0,rawq)):
            row=raw_row(g,mass)
            target=([(g-10+j,9-3*j,1+2*j) for j in range(4)] if mass==g
                    else [(2*g-21+j,21-3*j,2*j) for j in range(8)])
            need([v for v,_ in row]==target,"literal two-loop complete multiplicity fiber")
            need([F(w) for _,w in row]==[at(a,g) for a in formula],"literal multinomial coefficient normalization")
            anchor=target[0]
            for j,(vector,_) in enumerate(row):
                need(tuple(anchor[i]+j*(1,-3,2)[i] for i in range(3))==vector,"every phase-ratio monomial")
        first_anchor=(g-10,9,1)
        second_anchor=(2*g-21,21,0)
        need(second_anchor==tuple(2*first_anchor[i]-(1,-3,2)[i] for i in range(3)),
             "actual doubled anchor has lower carry minus one")
    need(next(m for m in range(1,13) if raw_row(12,m))==4,"gcd-dropped first-return hostile")
    need(at((1,42,84,3),-28)==-1175 and at((1,42,84,3),-27)==1054,
         "real cancellation phase outside ordinary-core real-rooted sector")
    # z^3+z^2+z+2 has one real root in (-2,-1); its remaining root sum is positive.
    need(at((2,1,1,1),-2)<0<at((2,1,1,1),-1),"positive coefficients without real spectrum hostile")
    need(2**2-4*3*1<0,"hostile cubic derivative strictly positive")
    print("Standard-library exact Q[g] reduction: factorial contents and inverse carry PASS")
    print("Independent symbolic B1/B2/B3 via principal minors/determinant: degrees 6/12/18 PASS")
    print("All 39 shifted rational characteristic coefficients equal the published records and are positive")
    print("Independent discriminant formula, literal source fibers, factorial weights and carry exponents PASS")
    print("No interpolation, numerical roots, producer imports, or additional parameter census")
    print("published_record_sha256="+hashlib.sha256(json.dumps(records,sort_keys=True,separators=(",",":")).encode()).hexdigest())
    print(f"PASS: {GATES} optimization-live exact gates")


if __name__=="__main__":
    main()
