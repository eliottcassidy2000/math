"""Independent p11/p13 proof certificate and modular Smith referee.

No producer imports. Integer determinant values plus exact degree bounds
certify all polynomial coefficients. Full Hasse matrices are peeled by
unit pivots modulo p^B and repeated division by p, not rational pivots.
"""
from fractions import Fraction as F
from itertools import combinations, permutations
from math import comb, gcd
from pathlib import Path
from collections import Counter
import hashlib
import json
import sys

sys.stdout.reconfigure(newline="\n")
gates=0
CERT_SHA="05f79d3dba323be53a1308ea174f61cc2cb79dff053c61cd4c93ec80a15751ec"


def need(ok,label):
    global gates
    gates+=1
    if not ok:
        raise ArithmeticError(label)


def val(x,p):
    if x==0:
        return 10**9
    x=F(x)
    def integer(n):
        v=0
        while n%p==0:
            n//=p
            v+=1
        return v
    return integer(abs(x.numerator))-integer(x.denominator)


def peval(P,x):
    result=0
    for a in reversed(P):
        result=result*x+a
    return result


def mul(A,B):
    C=[0]*(len(A)+len(B)-1)
    for i,a in enumerate(A):
        for j,b in enumerate(B):
            C[i+j]+=a*b
    return C


def determinant(A):
    # Fraction-free full pivoting: row and column choices differ from
    # the producer's fixed-column Bareiss samples.
    A=[r[:] for r in A]
    n=len(A)
    sign=1
    old=1
    for k in range(n-1):
        options=[(abs(A[i][j]),i,j) for i in range(k,n) for j in range(k,n) if A[i][j]]
        if not options:
            return 0
        _,i,j=min(options)
        if i!=k:
            A[i],A[k]=A[k],A[i]
            sign=-sign
        if j!=k:
            for row in A:
                row[j],row[k]=row[k],row[j]
            sign=-sign
        pivot=A[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                num=pivot*A[i][j]-A[i][k]*A[k][j]
                need(num%old==0,"exact referee fraction-free determinant division")
                A[i][j]=num//old
        for i in range(k+1,n):
            A[i][k]=0
        old=pivot
    return sign*A[-1][-1]


def key(p,r,rows,cols):
    return (p,r,tuple(tuple(v) for v in rows),tuple(cols))


def band_keys(p,r,depth):
    m=(p+1)//2
    R=[(node,j) for node in (1,2) for j in range(m)]
    C=list(range(m,3*m))
    top=sum(sorted((j for node,j in R),reverse=True)[:r])
    base=sum(C[:r])
    # Enumerate row and column levels separately, then take the permitted
    # product. No assumption about a last-column shift or boundary row.
    row_levels={d:[] for d in range(depth+1)}
    col_levels={d:[] for d in range(depth+1)}
    for rows in combinations(R,r):
        d=top-sum(j for node,j in rows)
        if d in row_levels:
            row_levels[d].append(rows)
    for cols in combinations(C,r):
        d=sum(cols)-base
        if d in col_levels:
            col_levels[d].append(cols)
    result={}
    for rd in row_levels:
        for cd in col_levels:
            if rd+cd<=depth:
                for rows in row_levels[rd]:
                    for cols in col_levels[cd]:
                        result[key(p,r,rows,cols)]=(base-top+rd+cd,rd+cd)
    return result


def packet_certificate(records):
    extra={11:{7,8},13:{8,9,10,11}}
    expected={}
    for p in (11,13):
        m=(p+1)//2
        for r in range(1,2*m-1):
            expected.update(band_keys(p,r,int(r in extra[p])))
    actual={key(z['p'],z['rank'],z['rows'],z['columns']):z for z in records}
    need(set(expected)==set(actual) and len(actual)==len(records)==66,
         "complete independent minimal/next-weight index universe")
    values=0
    packets={0:[1],1:[-1,2],2:[2,-7,7],3:[-1,5,-9,6],
             4:[14,-91,234,-286,143],5:[-3,24,-80,140,-130,52]}
    constants={11:[6,126,2940,185220,740880,40748400,32016600,124864740,252252,756756],
               13:[7,196,8232,1037232,11409552,1882576080,6409723320,116656964424,
                   1767529764,42420714336,80047968,17153136]}
    for k,z in actual.items():
        p,r,rows,cols=k
        need((z['weight'],z['loss'])==expected[k],"exact independently derived minor weight")
        second=[j for node,j in rows if node==2]
        count=len(second)
        lo=sum(cols[:count])-sum(second)
        hi=sum(cols[-count:])-sum(second) if count else 0
        if count==0:
            lo=0
        P=z['coefficients']
        need(all(c==0 for c in P[:lo]) and len(P)-1<=hi,"structural monomial factor and degree bound")
        # Each determinant term assigns exactly count distinct columns to
        # the a-node rows. Hence P/a^lo has degree at most hi-lo.
        for a in range(1,hi-lo+2):
            A=[[comb(c,j)*(a if node==2 else 1)**(c-j) for c in cols] for node,j in rows]
            need(determinant(A)==peval(P,a),"complete degree-bounded polynomial identity")
            values+=1
        content=0
        for c in P:
            content=gcd(content,abs(c))
        need(content==z['content'] and val(content,p)==z['content_valuation'],"all exact integer polynomial contents")
        if z['loss']==0:
            s=r//2
            shape=[1]
            for _ in range(s*s):
                shape=mul(shape,[-1,1])
            exponent=s*s
            if r%2:
                if sum(node==1 for node,j in rows)==s+1:
                    shape=mul(shape,packets[s])
                else:
                    shape=mul(shape,[(-1)**s*v for v in reversed(packets[s])])
                    exponent=(s+1)**2
            need(content==constants[p][r-1]
                 and P==[0]*exponent+[content*v for v in shape],"entire displayed factored packet and constant table")
    shifts=0
    for k,z in actual.items():
        p,r,rows,cols=k
        if p!=11:
            continue
        shifted=key(13,r,[(node,j+1) for node,j in rows],[c+1 for c in cols])
        if shifted not in actual:
            continue
        scale=F(1)
        for c in cols:
            scale*=c+1
        for node,j in rows:
            scale/=j+1
        need(actual[shifted]['coefficients']==[scale*c for c in z['coefficients']],
             "all available derivative/degree shifts as exact scalar identities")
        shifts+=1
    # Derive ideal corrections directly from complete minimal packets.
    corrections={}
    for p in (11,13):
        m=(p+1)//2
        corrections[p]=[]
        for r in range(1,2*m-1):
            group=[z for z in records if z['p']==p and z['rank']==r]
            minimal=[z for z in group if z['loss']==0]
            need(len(minimal)==1+(r%2),"number of minimum-weight polynomials")
            base=min(z['content_valuation'] for z in minimal)
            corrections[p].append(base)
            for a in range(2,p):
                exceptional=(p==13 and r==11 and a in (2,7,12))
                target=base+int(exceptional)
                need(min(val(peval(z['coefficients'],a),p) for z in minimal)==target,
                     "all finite-field minimum-packet attainments")
                if not exceptional:
                    need(any(peval(z['coefficients'],a)//(p**base)%p for z in minimal),
                         "universal unit packet after integral content removal")
            if base==2 or (p==13 and r==11):
                near=[z for z in group if z['loss']==1]
                need(len(near)==5+(r%2) and all(z['content_valuation']>=1 for z in near),
                     "all next-band obligations at e=1")
    need(corrections=={11:[0,0,0,0,0,1,2,2,1,1],
                       13:[0,0,0,0,0,0,1,2,2,2,1,1]},"independently reconstructed base corrections")
    return values,shifts


def predicted(p,e,a):
    m=(p+1)//2
    if e==0:
        return [0]*(3*m)
    # Construct cumulative ideals first, independently of the producer's
    # hard-coded factor list, and use the inherited exact largest loss.
    correction=([0,0,0,0,0,1,2,2,1,1] if p==11 else [0,0,0,0,0,0,1,2,2,2,1,1])
    E=[0]
    for r,shift in enumerate(correction,1):
        s=r//2
        weight=3*s*s if r%2==0 else 3*s*s+3*s+1
        if p==13 and r==11 and a%p in (2,7,12):
            shift+=1
        E.append(weight*e+shift)
    k=(p-1)//2
    sigma=int(sum(comb(k,j)**2*pow(a,j,p) for j in range(k+1))%p==0)
    total=3*m*m*e
    E.extend((total-((3*m-1)*e-sigma),total))
    return [0]*m+[v-u for u,v in zip(E,E[1:])]


def modular_smith(nodes,m,p):
    total=m*m*sum(val(a-b,p) for a,b in combinations(nodes,2))
    modulus=p**(total+1)
    N=3*m
    A=[[comb(c,j)*pow(x,c-j,modulus)%modulus if c>=j else 0 for c in range(N)]
       for x in nodes for j in range(m)]
    answer=[]
    layer=0
    while A:
        n=len(A)
        unit=next(((i,j) for i in range(n) for j in range(n) if A[i][j]%p),None)
        if unit is None:
            need(modulus>p and all(x%p==0 for row in A for x in row),"exact p-layer deflation")
            A=[[x//p for x in row] for row in A]
            modulus//=p
            layer+=1
            continue
        i,j=unit
        A[0],A[i]=A[i],A[0]
        for row in A:
            row[0],row[j]=row[j],row[0]
        inverse=pow(A[0][0],-1,modulus)
        A=[[(A[i][j]-A[i][0]*inverse*A[0][j])%modulus for j in range(1,n)] for i in range(1,n)]
        answer.append(layer)
    need(len(answer)==N and sum(answer)==total,"full original Hasse matrix determinant and all factor count")
    return answer


def divided_and_symmetry():
    q1=mul([-1,2],[3,-18,44,-52,26])
    q2=mul([-2,1],[26,-52,44,-18,3])
    lhs=[0]*8
    for i,c in enumerate(q2):
        lhs[i]+=c
    for i,c in enumerate(q1):
        lhs[i+2]+=c
    rhs=mul(mul([-1,1],[1,-1,1]),[4,-2,-1,-2,4])
    need(lhs==[13*c for c in rhs],"integer divided rank-eleven companion identity")
    need([peval(rhs,a)%13 for a in (2,7,12)]==[2,1,12],"all companion exceptional residues are units")
    need([a for a in range(2,13) if peval(q1,a)%13==peval(q2,a)%13==0]==[2,7,12],
         "exact AP exceptional packet")
    deep=[]
    for r in (2,7,12):
        for lift in range(13*13):
            a=r+13*lift
            need(min(val(peval(q1,a),13),val(peval(q2,a),13))==1,"every p3 lift of common residue has joint loss one")
            if val(peval(q1,a),13)>=3:
                deep.append((r,a,val(peval(q1,a),13),val(peval(q2,a),13)))
    need(len(deep)==3,"individual deep cancellation is not an ideal jump")
    for p in (3,5,7,11,13):
        k=(p-1)//2
        roots={a for a in range(2,p) if sum(comb(k,j)**2*a**j for j in range(k+1))%p==0}
        need(roots=={3:{2},5:set(),7:{2,4,6},11:{2,6,10},13:set()}[p],"precise smaller-prime Deuring comparison")
        for a in range(2,p):
            transforms={a%p,(1-a)%p,pow(a,-1,p),pow(1-a,-1,p),
                        (a-1)*pow(a,-1,p)%p,a*pow(a-1,-1,p)%p}
            need(all((b in roots)==(a in roots) for b in transforms),"S3 intrinsic Deuring bit")
            if p==13:
                need(all((b in (2,7,12))==(a in (2,7,12)) for b in transforms),"S3 intrinsic midpoint bit")
    return deep


def full_bank():
    count=0
    for p in (11,13):
        m=(p+1)//2
        for a in range(p*p):
            if a%p in (0,1):
                continue
            need(modular_smith((0,p,p*a),m,p)==predicted(p,1,a),"independent complete p2 parameter bank")
            count+=1
        for e in (0,3):
            for r in range(2,p):
                a=r+2*p**4
                need(modular_smith((0,p**e,p**e*a),m,p)==predicted(p,e,a),"outer depths and unrelated fourth-digit lifts")
                count+=1
        # All six label orders plus a different signed affine unit map.
        for order in permutations((0,p,p*(2-p**3))):
            nodes=tuple(5-3*x for x in order)
            need(modular_smith(nodes,m,p)==predicted(p,1,2),"complete node permutations and affine-unit normalization")
            count+=1
        for r in range(2,p):
            first,second=predicted(p,1,r),predicted(p,2,r)
            diffs1=[v-u for u,v in zip(first,first[1:])]
            diffs2=[v-u for u,v in zip(second,second[1:])]
            need(all(a>=0 and b-a>=0 for a,b in zip(diffs1,diffs2)),
                 "all-height ordering from nonnegative affine gap slopes and e1 endpoints")
    a=modular_smith((0,13,26),7,13)
    b=modular_smith((0,13,39),7,13)
    need(a[-1]==b[-1]==20 and a!=b,"literal full observer Deuring-blind intermediate hostile")
    need(sum(min(s,16) for s in a)==141 and sum(min(s,16) for s in b)==140,"actual finite-precision kernel orders")
    return count+2,a,b


def main():
    path=Path(__file__).with_name("overnight11_20260906_smith_prime_banks_certificates.json")
    raw=path.read_bytes()
    need(hashlib.sha256(raw).hexdigest()==CERT_SHA,"frozen standalone polynomial data")
    records=json.loads(raw)['records']
    values,shifts=packet_certificate(records)
    deep=divided_and_symmetry()
    matrices,a,b=full_bank()
    print("STATUS: PASS; full p11/p13 partitions, all-lift ideal proof, and precision hostile")
    print("POLYNOMIALS:",len(records),"; COMPLETE DEGREE-BOUNDED DETERMINANT VALUES:",values)
    print("EXACT SHARED DERIVATIVE/DEGREE SHIFTS:",shifts)
    print("INDEPENDENT BAND CONTENT TABLE:")
    for p,r in ((11,7),(11,8),(13,8),(13,9),(13,10),(13,11)):
        print(p,r,dict(Counter(z['content_valuation'] for z in records if z['p']==p and z['rank']==r and z['loss']==1)))
    print("INDIVIDUAL DEEP q1 LIFTS (residue,a,vq1,vq2):",deep)
    print("FULL ORIGINAL HASSE MODULAR SMITH MATRICES:",matrices)
    print("p13 HOSTILE:",a,b)
    print("KERNEL ORDERS MOD13^16: 13^141 and13^140, with common largest factor20")
    print("SCOPE: equilateral complete prime-order three-node banks only; no general-prime full law")
    print("ACTIVE GATES:",gates)


if __name__=="__main__":
    main()
