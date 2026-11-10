"""Exact recurrence extraction followed by finite-dimensional proof certificate."""
from apex_words import *
from fractions import Fraction as F

def bm(s):
    C=[F(1)]; B=[F(1)]; L=0; m=1; b=F(1)
    for n in range(len(s)):
        d=sum(C[j]*s[n-j] for j in range(len(C)))
        if not d: m+=1;continue
        old=C[:]; coef=d/b
        C += [F(0)]*max(0,len(B)+m-len(C))
        for j in range(len(B)): C[j+m]-=coef*B[j]
        if 2*L<=n: L=n+1-L; B=old; b=d;m=1
        else:m+=1
    return C

def main():
    # Pair states, with multiplicities, are exactly the tensor-square carrier.
    trans={}
    starts=[(ZERO,ZERO),(ZERO,18+ZERO),(18+ZERO,ZERO),(18+ZERO,18+ZERO)]
    seen=set(starts); todo=list(starts)
    while todo:
        i,j=todo.pop(); row=defaultdict(int)
        for a in range(4):
            for u in moves(ORDERS[i//18],STATES[i%18],a):
                for v in moves(ORDERS[j//18],STATES[j%18],a):
                    dest=(18*(i//18)+u,18*(j//18)+v)
                    row[dest]+=1
        trans[i,j]=row
        for d in row:
            if d not in seen:seen.add(d);todo.append(d)
    accept={p:int(accepting(ORDERS[p[0]//18],STATES[p[0]%18]) and
                  accepting(ORDERS[p[1]//18],STATES[p[1]%18])) for p in seen}
    v=dict.fromkeys(starts,1); seq=[]
    # Find a candidate, but do not trust a prefix fit as proof.
    for n in range(120):
        seq.append(sum(x*accept[p] for p,x in v.items()))
        new=defaultdict(int)
        for p,x in v.items():
            for q,m in trans[p].items():new[q]+=m*x
        v=new
    C=bm(seq); degree=len(C)-1
    assert all(x.denominator==1 for x in C)
    C=list(map(int,C))
    # Certificate: residuals vanish for dim(B) consecutive indices.
    # Each residual is f B^n p(B)c; Cayley-Hamilton propagates all later ones.
    needed=degree+len(seen)
    while len(seq)<needed:
        seq.append(sum(x*accept[p] for p,x in v.items()))
        new=defaultdict(int)
        for p,x in v.items():
            for q,m in trans[p].items():new[q]+=m*x
        v=new
    for n in range(degree,len(seq)):
        assert sum(C[j]*seq[n-j] for j in range(degree+1))==0
    assert seq[:6]==[1,12,192,4040,101568,2883072]
    out={'status':'EXACT-RECURRENCE-CERTIFICATE','reachable_tensor_states':len(seen),
         'recurrence_degree':degree,'backward_coefficients':C,
         'residuals_checked':len(seq)-degree,'initial_terms':seq[:degree],
         'proof':'Residual is a matrix coefficient of powers of B; dim(B) zero residuals imply every later one by Cayley-Hamilton.'}
    roots=[36]+[24]*2+[16]*2+[12]*2+[8]*4+[4]*3
    q=[1]
    for root in roots:
        row=[0]*(len(q)+1)
        for i,a in enumerate(q):row[i]+=a;row[i+1]-=root*a
        q=row
    assert C==q+[0,0,0]
    numerator=[sum(q[j]*seq[n-j] for j in range(min(14,n)+1)) for n in range(17)]
    at=F(1,36)
    leading=sum(a*at**i for i,a in enumerate(numerator))
    for root in roots[1:]:leading/=1-root*at
    assert leading==F(841,25600)
    out.update(characteristic_factorization='X^3*(X-36)*(X-24)^2*(X-16)^2*(X-12)^2*(X-8)^4*(X-4)^3',
               ogf_numerator_coefficients=numerator, leading_36_coefficient=str(leading),
               limiting_relative_variance=str(64*leading-1))
    Path(__file__).with_name('moment_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
