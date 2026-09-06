"""Exact p11/p13 three-bank Smith laws; bounded weighted-minor certificate.

Standard library only. Symbolic Laplace coefficients, literal Bareiss minors,
and full p-local matrix elimination are separate computation paths.
"""
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations
from math import comb, gcd
from pathlib import Path
import hashlib
import json
import sys

sys.stdout.reconfigure(newline='\n')
GATES=0

def need(ok,label):
    global GATES
    GATES+=1
    if not ok:raise RuntimeError(label)

def vp(value,p):
    if not value:return 10**9
    value=Q(value);a,b=abs(value.numerator),value.denominator;v=0
    while a%p==0:a//=p;v+=1
    while b%p==0:b//=p;v-=1
    return v

def trim(a):
    while len(a)>1 and a[-1]==0:a.pop()
    return a

def mul(a,b):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return trim(c)

def power(a,n):
    b=[1]
    for _ in range(n):b=mul(b,a)
    return b

def evaluate(a,x):
    y=0
    for c in a[::-1]:y=y*x+c
    return y

def det(M):
    if not M:return 1
    M=[row[:] for row in M];previous=1;sign=1
    for k in range(len(M)-1):
        if M[k][k]==0:
            i=next((i for i in range(k+1,len(M)) if M[i][k]),None)
            if i is None:return 0
            M[k],M[i]=M[i],M[k];sign=-sign
        pivot=M[k][k]
        for i in range(k+1,len(M)):
            for j in range(k+1,len(M)):
                value=M[i][j]*pivot-M[i][k]*M[k][j]
                if value%previous:raise RuntimeError('nonexact integer Bareiss division')
                M[i][j]=value//previous
        for i in range(k+1,len(M)):M[i][k]=0
        previous=pivot
    return sign*M[-1][-1]

def laplace(rows,cols):
    first=[j for x,j in rows if x==1]
    second=[j for x,j in rows if x==2]
    r=len(rows);result={}
    for chosen in combinations(range(r),len(first)):
        rest=[i for i in range(r) if i not in chosen]
        sign=(-1)**(sum(chosen)-len(first)*(len(first)-1)//2)
        coefficient=sign*det([[comb(cols[c],j) for c in chosen] for j in first])
        coefficient*=det([[comb(cols[c],j) for c in rest] for j in second])
        degree=sum(cols[c] for c in rest)-sum(second)
        result[degree]=result.get(degree,0)+coefficient
    maximum=max(result,default=0)
    return trim([result.get(i,0) for i in range(maximum+1)])

QPOLY={0:[1],1:[-1,2],2:[2,-7,7],3:[-1,5,-9,6],
       4:[14,-91,234,-286,143],5:[-3,24,-80,140,-130,52]}
CONTENTS={
  11:[6,126,2940,185220,740880,40748400,32016600,124864740,252252,756756],
  13:[7,196,8232,1037232,11409552,1882576080,6409723320,116656964424,
      1767529764,42420714336,80047968,17153136]}
CORRECTIONS={11:[0,0,0,0,0,1,2,2,1,1],13:[0,0,0,0,0,0,1,2,2,2,1,1]}
SECOND={11:{7,8},13:{8,9,10,11}}

def predicted(p,e,a):
    m=(p+1)//2
    if e==0:return [0]*(3*m)
    if p==11:
        sigma=int(a%11 in (2,6,10))
        tail=[e,2*e,4*e,5*e,7*e,8*e+1,10*e+1,11*e,13*e-1,14*e,
              16*e-1+sigma,17*e-sigma]
    elif p==13:
        tau=int(a%13 in (2,7,12))
        tail=[e,2*e,4*e,5*e,7*e,8*e,10*e+1,11*e+1,13*e,14*e,
              16*e-1+tau,17*e-tau,19*e-1,20*e]
    else:raise ValueError('proved scopes are p11 and p13')
    return [0]*m+tail

def local_smith(nodes,m,p):
    n=len(nodes)*m
    A=[[Q(comb(c,j)*x**(c-j) if c>=j else 0) for c in range(n)]
       for x in nodes for j in range(m)]
    result=[]
    for k in range(n):
        value,i,j=min((vp(A[i][j],p),i,j) for i in range(k,n) for j in range(k,n))
        need(0<=value<10**9,'finite p-integral pivot')
        A[k],A[i]=A[i],A[k]
        for row in A:row[k],row[j]=row[j],row[k]
        pivot=A[k][k];result.append(value)
        for i in range(k+1,n):
            factor=A[i][k]/pivot
            need(vp(factor,p)>=0,'p-integral unimodular elimination multiplier')
            for j in range(k+1,n):A[i][j]-=factor*A[k][j]
            A[i][k]=Q(0)
        for j in range(k+1,n):A[k][j]=Q(0)
    need(result==sorted(result),'ordered full Smith factors')
    need(sum(result)==m*m*sum(vp(x-y,p) for x,y in combinations(nodes,2)),
         'literal full observer determinant valuation')
    return result

def certificate():
    records=[]
    for p in (11,13):
        m=(p+1)//2
        allrows=[(x,j) for x in (1,2) for j in range(m)]
        for r in range(1,2*m-1):
            cols=list(range(m,m+r))
            top=sum(sorted([j for x,j in allrows],reverse=True)[:r])
            weight=sum(cols)-top
            cases=[]
            for rows in combinations(allrows,r):
                loss=top-sum(j for x,j in rows)
                if loss==0 or (loss==1 and r in SECOND[p]):cases.append((rows,cols,loss))
                if loss==0 and r in SECOND[p]:cases.append((rows,cols[:-1]+[cols[-1]+1],1))
            need(sum(loss==0 for _,_,loss in cases)==(2 if r%2 else 1),
                 'all equality-weight row and column sets')
            if r in SECOND[p]:
                need(sum(loss==1 for _,_,loss in cases)==(6 if r%2 else 5),
                     'complete next-weight band')
            for rows,columns,loss in cases:
                poly=laplace(rows,columns)
                content=0
                for c in poly:content=gcd(content,abs(c))
                need(content>0,'nonzero retained minor')
                if loss==0:
                    s=r//2
                    exponent=s*s
                    shape=power([-1,1],s*s)
                    if r%2:
                        if sum(x==1 for x,j in rows)==s+1:
                            shape=mul(shape,QPOLY[s])
                        else:
                            shape=mul(shape,[(-1)**s*c for c in QPOLY[s][::-1]])
                            exponent=(s+1)**2
                    expected=[0]*exponent+[CONTENTS[p][r-1]*c for c in shape]
                    need(poly==trim(expected),'complete factored minimal-weight identity')
                    need(vp(content,p)==CORRECTIONS[p][r-1], 'minimal content valuation')
                else:
                    need(vp(content,p)>=1,'next-weight band divisibility needed at e1')
                for a in (-3,0,1,2,17):
                    literal=[[comb(c,j)*(1 if x==1 else a)**(c-j) for c in columns]
                             for x,j in rows]
                    need(det(literal)==evaluate(poly,a),'independent literal Bareiss minor')
                records.append(dict(p=p,rank=r,rows=rows,columns=columns,weight=weight+loss,
                                    loss=loss,coefficients=poly,content=content,
                                    content_valuation=vp(content,p)))
            minimal=[laplace(rows,columns) for rows,columns,loss in cases if loss==0]
            for a in range(2,p):
                delta=CORRECTIONS[p][r-1]+int(p==13 and r==11 and a in (2,7,12))
                need(min(vp(evaluate(poly,a),p) for poly in minimal)==delta,
                     'attaining minimal-weight companion at every residue')
                # Reduced packet checks, distinct from testing full valuations at one lift.
                if not (p==13 and r==11 and a in (2,7,12)):
                    need(any(evaluate(poly,a)//p**CORRECTIONS[p][r-1]%p for poly in minimal),
                         'universal all-lift unit in minimal packet')
    need(len(records)==66,'bounded certificate:26 p11 and40 p13 polynomials')
    return records

def divided_companion():
    first=QPOLY[5]
    second=[-c for c in first[::-1]]
    combined=[0]*8
    for i,c in enumerate(second):combined[i]+=c
    for i,c in enumerate(first):combined[i+2]+=c
    need(all(c%13==0 for c in combined),'integral divided rank11 companion')
    expected=mul(mul([-1,1],[1,-1,1]),[4,-2,-1,-2,4])
    need(trim([c//13 for c in combined])==expected,'exact divided companion factorization')
    for a,value in ((2,2),(7,1),(12,12)):
        need(evaluate(expected,a)%13==value,'companion unit on each AP residue')
    for p in (11,13):
        k=(p-1)//2
        roots=[a for a in range(2,p) if sum(comb(k,j)**2*a**j for j in range(k+1))%p==0]
        need(roots==([2,6,10] if p==11 else []),'inherited exact Deuring roots')
    # The same lowest-weight normalized determinant is independent of m,
    # because binom(c+1,j+1)=(c+1)/(j+1) binom(c,j).
    for r in range(1,11):
        s=r//2
        rows=[(1,j) for j in range(6-(s+1 if r%2 else s),6)]
        rows += [(2,j) for j in range(6-s,6)]
        cols=list(range(6,6+r))
        old=laplace(rows,cols)
        new=laplace([(x,j+1) for x,j in rows],[c+1 for c in cols])
        ratio=Q(1)
        for c in cols:ratio*=c+1
        for x,j in rows:ratio/=j+1
        need([Q(c) for c in new]==[ratio*c for c in old], 'exact derivative/degree shift map')

def full_controls():
    count=0; examples=[]
    for p in (11,13):
        m=(p+1)//2
        for a in range(p*p):
            if a%p in (0,1):continue
            need(local_smith((0,p,p*a),m,p)==predicted(p,1,a),'complete admissible p2 lift bank e1')
            count+=1
        for e in (0,2,5):
            for residue in range(2,p):
                a=residue-p**3
                need(local_smith((0,p**e,p**e*a),m,p)==predicted(p,e,a),
                     'negative higher unit lifts and independent outer depths')
                count+=1
        for a in (2,3):
            base=(0,p,p*a)
            for transformed in ((base[2],base[0],base[1]),tuple(7-2*x for x in base)):
                need(local_smith(transformed,m,p)==predicted(p,1,a),'permutation and affine-unit controls')
                count+=1
            examples.append(dict(p=p,a=a,partition=predicted(p,1,a)))
    ap=predicted(13,1,2);ordinary=predicted(13,1,3)
    need(sum(min(x,16) for x in ap)==141 and sum(min(x,16) for x in ordinary)==140,
         'actual finite-precision kernel-order separation')
    for p in (11,13):
        for e in range(1,11):
            for a in range(2,p):
                values=predicted(p,e,a)
                need(values==sorted(values) and sum(values)==3*((p+1)//2)**2*e,
                     'formula ordering and determinant at all branch boundaries')
    return count,examples

def main():
    records=certificate()
    divided_companion()
    matrices,examples=full_controls()
    data=dict(scope='PROVED p11 and p13 equilateral three-bank prime-order laws',
              records=records,full_matrix_controls=matrices,examples=examples)
    encoded=(json.dumps(data,sort_keys=True,indent=2)+'\n').encode()
    destination=Path(__file__).with_name(Path(__file__).stem+'_certificates.json')
    destination.write_bytes(encoded)
    print('STATUS: PASS; two exact full partitions and an intermediate-ideal Deuring hostile')
    print('MINOR CERTIFICATE: 26 p11 plus 40 p13; all other minors controlled by weight gap')
    print('NEXT-BAND CONTENT VALUATIONS:')
    for p in (11,13):
        for r in sorted(SECOND[p]):
            print(p,r,dict(sorted(Counter(z['content_valuation'] for z in records
                                          if z['p']==p and z['rank']==r and z['loss']==1).items())))
    for example in examples:print(json.dumps(example,sort_keys=True))
    print('FULL MATRIX CONTROLS:',matrices)
    print('HOSTILE: p13 a2/a3 both Deuring ordinary, largest 20; mod13^16 kernel orders 13^141/13^140')
    print('ALL-LIFT SIDECAR: (q2+a^2 q1)/13=(a-1)(a^2-a+1)(4a^4-2a^3-a^2-2a+4)')
    print('CERTIFICATE SHA256:',hashlib.sha256(encoded).hexdigest())
    print('ACTIVE GATES:',GATES)

if __name__=='__main__':main()
