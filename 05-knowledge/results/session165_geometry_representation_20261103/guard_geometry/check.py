#!/usr/bin/env python3
"""Self-contained exact checks for the S165 guard-state geometry note."""
from fractions import Fraction
from itertools import product
from math import comb
import json


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


C=((0,1,0),(1,1,1),(1,0,0))
D=((1,0,1),(0,2,0),(0,0,0))
I=((1,0,0),(0,1,0),(0,0,1))
v0=(1,2,0)


def mm(a,b):
    return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(len(b)))
                       for j in range(len(b[0]))) for i in range(len(a)))


def mv(a,v):
    return tuple(sum(x*y for x,y in zip(row,v)) for row in a)


def power(a,n):
    b=I
    for _ in range(n):
        b=mm(a,b)
    return b


def word_matrix(word):
    m=I
    for letter in word:
        m=mm(C if letter=='C' else D,m)
    return m


def det(a):
    a=[list(map(Fraction,row)) for row in a]
    answer=Fraction(1)
    for k in range(len(a)):
        pivot=next((i for i in range(k,len(a)) if a[i][k]),None)
        if pivot is None:
            return 0
        if pivot!=k:
            a[k],a[pivot]=a[pivot],a[k]
            answer=-answer
        p=a[k][k]
        answer*=p
        for j in range(k,len(a)):
            a[k][j]/=p
        for i in range(k+1,len(a)):
            p=a[i][k]
            if p:
                for j in range(k,len(a)):
                    a[i][j]-=p*a[k][j]
    require(answer.denominator==1,'integral determinant')
    return int(answer)


def rank(rows):
    """Exact incremental rational row rank, stopping once all columns pivot."""
    if not rows:
        return 0
    basis={}
    width=len(rows[0])
    for raw in rows:
        row=list(map(Fraction,raw))
        for j in sorted(basis):
            if row[j]:
                q=row[j]
                row=[a-q*b for a,b in zip(row,basis[j])]
        pivot=next((j for j,x in enumerate(row) if x),None)
        if pivot is not None:
            q=row[pivot]
            basis[pivot]=[a/q for a in row]
            if len(basis)==width:
                return width
    return len(basis)


def transpose(a):
    return tuple(zip(*a))


def flat(a):
    return tuple(x for row in a for x in row)


def tournament(n,bits):
    rows=[0]*n
    k=0
    for i in range(n):
        for j in range(i+1,n):
            if (bits>>k)&1:
                rows[i]|=1<<j
            else:
                rows[j]|=1<<i
            k+=1
    return tuple(rows)


def starts(rows):
    """Independent subset/first-endpoint DP. No transfer formulas."""
    n=len(rows)
    dp=[[0]*n for _ in range(1<<n)]
    for i in range(n):
        dp[1<<i][i]=1
    for mask in range(1,1<<n):
        remaining=mask
        while remaining:
            low=remaining&-remaining
            i=low.bit_length()-1
            rest=mask^low
            if rest:
                choices=rows[i]&rest
                total=0
                while choices:
                    nxt=choices&-choices
                    total+=dp[rest][nxt.bit_length()-1]
                    choices^=nxt
                dp[mask][i]=total
            remaining^=low
    return tuple(dp[-1])


def induced(rows,vertices):
    return tuple(sum(1<<j for j,v in enumerate(vertices) if (rows[u]>>v)&1)
                 for u in vertices)


def state(rows,q):
    incoming=[r for r in range(len(rows)) if (rows[r]>>q)&1]
    require(len(incoming)==1,'guard must have indegree one')
    hp=starts(rows)
    vertices=[r for r in range(len(rows)) if r!=q]
    smaller=starts(induced(rows,vertices))
    return hp[q],sum(hp)-hp[q],smaller[vertices.index(incoming[0])]


def strong(rows):
    full=(1<<len(rows))-1
    for q in range(len(rows)):
        reached=1<<q
        while True:
            nxt=reached
            for i,r in enumerate(rows):
                if (reached>>i)&1:
                    nxt|=r
            if nxt==reached:
                break
            reached=nxt
        if reached!=full:
            return False
    return True


def attach(rows,q,letter):
    new=len(rows)
    out=list(rows)+[0]
    out[q]|=1<<new
    out[new]=((1<<new)-1)^(1<<q)
    return tuple(out),new if letter=='C' else q


def literal_word(word):
    rows=(2,4,1)
    q=0
    for letter in word:
        rows,q=attach(rows,q,letter)
    return rows,q


def gap_ok(v):
    s,u,d=v
    if d:
        return True
    m=1
    while s*(1<<(m+1))<=u:
        m+=1
    return 3*u<=s*(1<<(m+2)) and s*(1<<m)<=u


def family(n):
    rows=[sum(1<<j for j in range(i+1,n)) for i in range(n)]
    for q in (0,2):
        rows[q]^=1<<(n-1)
        rows[n-1]|=1<<q
    return tuple(rows)


def parity_conditions(m):
    a,b,c,e,f,g,h,i,j=flat(m)
    return all(x%2==0 for x in (e-h,f-h-j,g-h-i,a-b-c-h-i-j))


def main():
    report={}
    O=((1,1,0),(1,2,1),(3,3,2))
    V=transpose((v0,mv(C,v0),mv(power(C,2),v0)))
    L=((1,2,1),(3,3,2),(5,6,3))
    R=mm(mm(D,C),D)
    require(R==((1,2,1),(2,4,2),(0,0,0)),'rank-one reset')
    require(det(O)==2 and det(V)==-2 and det(L)==2,'observation/control determinants')
    require(det(mm(O,V))==-4,'Hankel determinant')
    rank_one=[flat(mm(mm(power(C,i),R),power(C,j))) for i in range(3) for j in range(3)]
    require(rank(rank_one)==9,'rank-one matrix algebra basis')
    basis_words=('','C','D','CC','CD','DC','DD','CCD','CDD')
    integral_basis=[flat(word_matrix(w)) for w in basis_words]
    require(det(transpose(integral_basis))==16,'integral matrix basis determinant')
    report['matrix_certificates']={'det_O':2,'det_V':-2,'det_L':2,'hankel_det':-4,
                                   'integral_index':16,'basis_temporal_words':basis_words}

    transfer_cases=0
    for n in range(2,5):
        for bits in range(1<<(n*(n-1)//2)):
            rows=tournament(n,bits)
            for q in range(n):
                if sum((r>>q)&1 for r in rows)!=1:
                    continue
                v=state(rows,q)
                for letter,mat in (('C',C),('D',D)):
                    out,p=attach(rows,q,letter)
                    require(state(out,p)==mv(mat,v),'literal transfer on all small guards')
                    if strong(rows):
                        require(strong(out),'strongness transfer')
                out,p=attach(rows,q,'C')
                out2,p2=attach(out,p,'C')
                h0=sum(starts(rows));h1=sum(starts(out));h2=sum(starts(out2))
                require((h0-2*h1+h2,h0+2*h1-h2,h2-3*h0)==tuple(2*x for x in v),
                        'literal suffix tomography')
                transfer_cases+=1
    require(transfer_cases==110,'small guarded census')
    report['small_guarded_sources']=transfer_cases

    word_count=0;boundary_count=0;literal_count=0
    for n in range(13):
        for letters in product('CD',repeat=n):
            word=''.join(letters)
            m=word_matrix(word);v=mv(m,v0)
            s,u,d=v
            require(s>0 and u>0 and d>=0 and u>=s+d,'state cone')
            require((s+u)%2==1 and (u-d)%2==0,'physical parity')
            require(gap_ok(v),'dyadic projective gap')
            require(parity_conditions(m),'integral algebra congruences')
            boundary_count+=(d==0)
            if n<=6:
                rows,q=literal_word(word)
                require(state(rows,q)==v and strong(rows),'literal word state')
                literal_count+=1
            word_count+=1
    report['finite_word_controls']={'maximum_depth':12,'words':word_count,
                                    'defect_zero_states_with_multiplicity':boundary_count,
                                    'literal_words_through_depth_6':literal_count}

    polynomial_ranks=[]
    for delta in range(5):
        exponents=[(a,b,c) for a in range(delta+1) for b in range(delta+1-a)
                   for c in range(delta+1-a-b)]
        rows=[]
        for m in range(delta+1):
            for n in range(1,delta+2):
                for k in range(delta+1):
                    v=mv(word_matrix('D'*m+'C'+'D'*n+'C'*k),v0)
                    A=1+(1<<(m+1));B=1<<n
                    require(v==mv(power(C,k),(A,A*B,0)),'polynomial bank parameterization')
                    rows.append([v[0]**a*v[1]**b*v[2]**c for a,b,c in exponents])
        r=rank(rows)
        require(r==comb(delta+3,3),'full polynomial bank rank')
        polynomial_ranks.append({'degree':delta,'samples':len(rows),'rank':r})
    report['inhomogeneous_polynomial_bank_ranks']=polynomial_ranks

    homogeneous_ranks=[]
    for delta in range(7):
        exponents=[(a,b,delta-a-b) for a in range(delta+1) for b in range(delta+1-a)]
        rows=[]
        for n in range(delta+1):
            for k in range(delta+1):
                v=mv(word_matrix('D'*n+'C'*k),v0)
                rows.append([v[0]**a*v[1]**b*v[2]**c for a,b,c in exponents])
        r=rank(rows)
        require(r==comb(delta+2,2),'homogeneous polynomial bank rank')
        homogeneous_ranks.append({'degree':delta,'samples':len(rows),'rank':r})
    report['homogeneous_polynomial_bank_ranks']=homogeneous_ranks

    # The inverse follows from predecessor (s,u,d) -> (d,s,u-s-d).
    Cinv=((0,0,1),(1,0,0),(-1,1,-1))
    require(mm(C,Cinv)==I,'inverse matrix')
    normals=[]
    for k in range(31):
        normal=power(Cinv,k)[2]
        require(normal not in normals and tuple(-a for a in normal) not in normals,
                'distinct primitive plane normals')
        normals.append(normal)
    report['distinct_reset_planes_checked']=len(normals)

    family_data=[]
    for n in range(6,16):
        rows=family(n);v=state(rows,0);t=1<<(n-4)
        require(v==(1+t,4*t-2,0),'literal two-backedge family state')
        require(strong(rows),'family strongness')
        require(4*v[0]>v[1] and 8*v[0]<3*v[1],'strict forbidden ratio')
        require(sum(v[:2])==5*t-1,'family Hamilton count')
        reset_rows,reset_q=attach(rows,0,'C')
        reset_rows,reset_q=attach(reset_rows,reset_q,'D')
        # Direct endpoint enumeration of the reset is kept to <=13 vertices.
        if n<=11:
            require(state(reset_rows,reset_q)==tuple((5*t-1)*x for x in v0),
                    'literal family canonical reset')
        family_data.append({'order':n,'state':v,'H':sum(v[:2])})
    report['literal_escape_family']=family_data
    report['six_vertex_rows']=family(6)
    report['status']='all checks passed'
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
