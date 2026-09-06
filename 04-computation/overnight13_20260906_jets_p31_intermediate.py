"""Universal odd-minor packets and the p31 intermediate ideal; exact stdlib."""
from math import comb,gcd,factorial
from fractions import Fraction
from functools import reduce

gates=0
def check(ok,message):
    global gates
    gates+=1
    if not ok: raise RuntimeError(message)

def packet(s):
    coeff=[(-1)**(s-i)*sum(comb(s+j,j)*comb(2*s-j,s-j)*comb(j,i-s+j)
                         for j in range(s-i,s+1)) for i in range(s+1)]
    content=reduce(gcd,coeff)
    return content,tuple(x//content for x in coeff)

def evaluate(coeff,a):
    value=0
    for c in reversed(coeff):value=value*a+c
    return value

def valuation(x,p):
    if x==0:return 10**9
    x=abs(x);out=0
    while x%p==0:x//=p;out+=1
    return out

def constant(m,s):
    content,_=packet(s);h=m-s-1
    check(h>=0,'shift domain')
    ratio=Fraction(content)
    for c in range(s+1,3*s+2):ratio*=Fraction(factorial(c+h),factorial(c))
    ratio/=factorial(h)
    for j in range(1,s+1):ratio/=Fraction(factorial(j+h),factorial(j))**2
    check(ratio.denominator==1,'integral constant')
    return ratio.numerator

def determinant(matrix):
    A=[row[:] for row in matrix];n=len(A);previous=1;sign=1
    if n==0:return 1
    for k in range(n-1):
        r=next((r for r in range(k,n) if A[r][k]),None)
        if r is None:return 0
        if r!=k:A[r],A[k]=A[k],A[r];sign=-sign
        pivot=A[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                numerator=A[i][j]*pivot-A[i][k]*A[k][j]
                check(numerator%previous==0,'Bareiss exact division')
                A[i][j]=numerator//previous
        for i in range(k+1,n):A[i][k]=0
        previous=pivot
    return sign*A[-1][-1]

def minimal_rows(m,s):
    h=m-s-1
    high=[(x,j) for x in (0,1) for j in range(h+1,m)]
    return [sorted(high+[(x,h)]) for x in (0,1)]

def normalized_minor(rows,columns,a):
    return [[comb(c,j)*(1 if x==0 else a**(c-j)) for c in columns] for x,j in rows]

def smith_layers(p,a,e=1):
    m=(p+1)//2;precision=max(1,(3*m-1)*e+1);modulus=p**precision
    A=[[comb(c,j)*pow(x,c-j,modulus)%modulus if c>=j else 0
        for c in range(3*m)] for x in (0,p**e,p**e*a) for j in range(m)]
    exponents=[];level=0
    while A:
        n=len(A)
        hit=next(((r,c) for r in range(n) for c in range(n) if A[r][c]%p),None)
        if hit is None:
            check(all(x%p==0 for row in A for x in row),'common layer division')
            A=[[x//p for x in row] for row in A]
            modulus//=p;level+=1
            check(modulus>1,'sufficient modulus')
            continue
        r,c=hit;A[0],A[r]=A[r],A[0]
        for row in A:row[0],row[c]=row[c],row[0]
        inv=pow(A[0][0],-1,modulus)
        pivot=[x*inv%modulus for x in A[0]]
        check(pivot[0]==1,'unit pivot')
        A=[[(row[j]-row[0]*pivot[j])%modulus for j in range(1,n)] for row in A[1:]]
        exponents.append(level)
    check(sum(exponents)==3*m*m*e,'confluent determinant')
    return tuple(exponents)

def deuring(p,a):
    k=(p-1)//2
    return sum(comb(k,j)**2*pow(a,j,p) for j in range(k+1))%p==0

def main():
    print('Universal odd-minor packet and p31 residual E29=631e+1+kappa')
    for s in range(7):
        content,q=packet(s)
        for m in (s+1,s+2,s+4):
            K=constant(m,s);columns=list(range(m,m+2*s+1))
            for a in (2,3,5):
                for which,rows in enumerate(minimal_rows(m,s)):
                    expected=K*a**(s*s if which==0 else (s+1)**2)*(a-1)**(s*s)*evaluate(q if which==0 else tuple(reversed(q)),a)
                    check(abs(determinant(normalized_minor(rows,columns,a)))==abs(expected),'universal minor packet')
    s=14;p=31;m=16;content,q=packet(s);q2=tuple(reversed(q))
    check(content==19380,'primitive packet content')
    K=constant(m,s)
    check(K==23038504627568008043520 and valuation(K,p)==1,'shifted p31 content')
    companion=[]
    for i in range(s+3):
        num=(q2[i] if i<len(q2) else 0)-(q[i-2] if 2<=i<len(q)+2 else 0)
        check(num%p==0,'divided companion coefficient')
        companion.append(num//p)
    roots=[]
    for a in range(2,p):
        left=evaluate(q,a)%p;right=evaluate(q2,a)%p
        check(right==a*a*left%p,'residue proportionality')
        if left==right==0:
            roots.append(a)
            check(evaluate(companion,a)%p!=0,'all-lift common loss one')
        else:check(left!=0 or right!=0,'ordinary unit packet')
    check(roots==[3,11,15,17,21,29],'complete exceptional alphabet')
    for a in range(2,p):
        orbit={a,pow(a,-1,p),(1-a)%p,pow(1-a,-1,p),a*pow(a-1,-1,p)%p,(a-1)*pow(a,-1,p)%p}
        check(all((b in roots)==(a in roots) for b in orbit),'intrinsic cross-ratio orbit')
    check(all(not deuring(p,a) and a not in (2,16,30) for a in roots),'exceptional alphabet Deuring and AP ordinary')
    print('q14 coefficients increasing',q)
    print('content',content,'shifted K',K)
    print('exceptional residues/divided units',[(a,evaluate(companion,a)%p) for a in roots])
    lifts=0
    for a0 in range(2,p):
        for digit in range(p):
            a=a0+p*digit
            check(min(valuation(evaluate(q,a),p),valuation(evaluate(q2,a),p))==int(a0 in roots),'second-digit common loss')
            lifts+=1
    print('all admissible second-digit packet lifts',lifts)
    for a in roots:
        lift=next(a+p*k for k in range(p) if evaluate(q,a+p*k)%(p*p)==0)
        check(valuation(evaluate(q2,lift),p)==1,'individual cancellation hostile')
        print('individual q1 lift',lift,'v31(q1,q2)',valuation(evaluate(q,lift),p),valuation(evaluate(q2,lift),p))
    # Independently enumerate every row complement of size three and column complement of size three.
    from itertools import combinations
    allrows=[(x,j) for x in (0,1) for j in range(m)]
    allcols=list(range(m,3*m));best_rows={};best_cols={}
    maxsum=239;minsum=sum(range(16,45))
    for omitted in combinations(range(32),3):
        rem=set(omitted)
        rows=[row for i,row in enumerate(allrows) if i not in rem]
        loss=maxsum-sum(j for x,j in rows)
        if loss<=1:best_rows.setdefault(loss,[]).append(rows)
        cols=[c for i,c in enumerate(allcols) if i not in rem]
        loss=sum(cols)-minsum
        if loss<=1:best_cols.setdefault(loss,[]).append(cols)
    check({k:len(v) for k,v in best_rows.items()}=={0:2,1:4},'complete row band')
    check({k:len(v) for k,v in best_cols.items()}=={0:1,1:1},'complete column band')
    next_bands=[]
    for rowloss,rowsets in best_rows.items():
        for colloss,colsets in best_cols.items():
            if rowloss+colloss!=1:continue
            for rows in rowsets:
                for cols in colsets:
                    check(sum(cols)-sum(j for x,j in rows)==632,'next weight')
                    if all(j>=1 for x,j in rows):
                        check(31 in cols and all(comb(31,j)%31==0 for x,j in rows),'zero column31')
                    else:
                        check(sum(j==0 for x,j in rows)==1 and all(j==0 or j>=2 for x,j in rows),'single zero derivative')
                        check(31 in cols and 32 in cols,'two sparse columns retained')
                        check(all(comb(c,j)%31==0 for x,j in rows if j>=2 for c in (31,32)),'two columns one row')
                    next_bands.append((rows,cols))
    check(len(next_bands)==6,'complete next-band content proof')
    for a in (2,3,4):
        for which,rows in enumerate(minimal_rows(m,s)):
            expected=K*a**(s*s if which==0 else (s+1)**2)*(a-1)**(s*s)*evaluate(q if which==0 else q2,a)
            check(abs(determinant(normalized_minor(rows,list(range(16,45)),a)))==abs(expected),'literal p31 minimal determinant')
    matrices=0
    for e in (0,1,2,4):
        for a in range(2,p):
            lifted=a if e==1 else a+p*p*(a-17)
            parts=smith_layers(p,lifted,e)
            check(sum(parts[:45])==(0 if e==0 else 631*e+1+int(a in roots)),'full E45')
            check(parts[-1]==(0 if e==0 else 47*e-int(deuring(p,a))),'inherited largest factor')
            matrices+=1
    for a in (3,4):
        parts=smith_layers(p,a)
        check(not deuring(p,a) and a not in (2,16,30),'Deuring and AP blind')
        print('a',a,'full e1 partition',parts)
        print('kernel exponent at precision43',sum(min(43,v) for v in parts))
        matrices+=1
    print('original full Hasse matrix controls',matrices)
    print('PASS',gates,'always-active gates')

if __name__=='__main__':main()
