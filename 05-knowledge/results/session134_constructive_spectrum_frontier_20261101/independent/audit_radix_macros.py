#!/usr/bin/env python3
"""Independent all-product, all-integer-eigenvalue audit through four operations.
No producer import; local matrices are retained independently audited inputs.
"""
import itertools,json,pathlib
from fractions import Fraction
from audit_migration_closure import conjugate,S0,S1
from audit_migration_compiler import build,state as direct_state,literal

def need(c,msg):
    if not c:raise RuntimeError(msg)
def det3(a,b,c,d,e,f,g,h,i):return a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g)
def determinant(m):
    a,b,c,d,e,f,g,h,i,j,k,l,o,p,q,r=m
    return a*det3(f,g,h,j,k,l,p,q,r)-b*det3(e,g,h,i,k,l,o,q,r)+c*det3(e,f,h,i,j,l,o,p,r)-d*det3(e,f,g,i,j,k,o,p,q)
def det_shift(m,B):
    x=list(m)
    for i in (0,5,10,15):x[i]-=B
    return determinant(x)
def matrix_mul(A,B):
    return tuple(sum(A[4*i+k]*B[4*k+j] for k in range(4)) for i in range(4) for j in range(4))
def charpoly(m):
    trace=sum(m[i] for i in (0,5,10,15))
    e2=sum(m[4*i+i]*m[4*j+j]-m[4*i+j]*m[4*j+i] for i,j in itertools.combinations(range(4),2))
    e3=0
    for ids in itertools.combinations(range(4),3):e3+=det3(*(m[4*i+j] for i in ids for j in ids))
    return [Fraction(determinant(m)),Fraction(-e3),Fraction(e2),Fraction(-trace),Fraction(1)]
def trim(p):
    p=list(p)
    while p and p[-1]==0:p.pop()
    return p

def remainder(a,b):
    a=trim(a);b=trim(b)
    while len(a)>=len(b):
        c=a[-1]/b[-1];shift=len(a)-len(b)
        for j,v in enumerate(b):a[j+shift]-=c*v
        a=trim(a)
    return a

def evaluate(p,x):
    ans=Fraction(0)
    for c in p[::-1]:ans=ans*x+c
    return ans

def integer_roots_ge2(p,bound):
    seq=[trim(p),[i*p[i] for i in range(1,len(p))]]
    while seq[-1]:
        r=remainder(seq[-2],seq[-1])
        if not r:break
        seq.append([-v for v in r])
    cache={}
    def variation(twox):
        if twox not in cache:
            vals=[evaluate(s,Fraction(twox,2)) for s in seq]
            need(all(v!=0 for v in vals[:1]),'monic polynomial has half-integer root')
            signs=[1 if v>0 else -1 for v in vals if v]
            cache[twox]=sum(a!=b for a,b in zip(signs,signs[1:]))
        return cache[twox]
    answer=[];stack=[(3,2*bound+1)] if bound>=2 else []
    while stack:
        lo,hi=stack.pop()
        if variation(lo)==variation(hi):continue
        if hi-lo==2:
            B=(lo+hi)//4
            if evaluate(p,B)==0:answer.append(B)
        else:
            mid=2*((lo+hi)//4)+1
            need(lo<mid<hi,'bad half-integer split');stack.extend(((lo,mid),(mid,hi)))
    return sorted(answer)

def nullspace(m,B):
    a=[[Fraction(m[4*i+j]-(B if i==j else 0)) for j in range(4)] for i in range(4)]
    piv=[];row=0
    for col in range(4):
        pivot=next((r for r in range(row,4) if a[r][col]),None)
        if pivot is None:continue
        a[row],a[pivot]=a[pivot],a[row];scale=a[row][col];a[row]=[v/scale for v in a[row]]
        for r in range(4):
            if r!=row and a[r][col]:
                c=a[r][col];a[r]=[u-c*v for u,v in zip(a[r],a[row])]
        piv.append(col);row+=1
    free=[i for i in range(4) if i not in piv];basis=[]
    for f in free:
        v=[Fraction(0)]*4;v[f]=1
        for r,c in enumerate(piv):v[c]=-a[r][f]
        basis.append(tuple(v))
    return basis

def main():
    root=pathlib.Path(__file__).resolve().parent
    raw=json.loads((root/'migration_matrix_bank.json').read_text())['records'];bank={S0,S1}
    for rec in raw:
        a=rec['local'];fresh=range(1,len(a));cycle=any((a[i]>>j&1 and a[j]>>k&1 and a[k]>>i&1) or (a[i]>>k&1 and a[k]>>j&1 and a[j]>>i&1) for i,j,k in itertools.combinations(fresh,3))
        if not cycle:bank.add(conjugate(rec['matrix']))
    need(len(bank)==32,'bank size');bank=tuple(sorted(tuple(v for row in M for v in row) for M in bank))
    # Independent Sturm controls, including repeated and noninteger roots.
    for diag,expect in (((1,2,3,4),[2,3,4]),((0,1,2,2),[2])):
        M=tuple(diag[i] if i==j else 0 for i in range(4) for j in range(4))
        need(integer_roots_ge2(charpoly(M),max(diag))==expect,'Sturm diagonal control')
    golden=(1,1,0,0,1,0,0,0,0,0,1,0,0,0,0,0)
    need(integer_roots_ge2(charpoly(golden),2)==[],'Sturm algebraic-noninteger control')
    current={tuple(int(i==j) for i in range(4) for j in range(4))};union=set();depths=[]
    # Sparse left-action implementation, independently deduplicating every depth.
    sparse=[tuple(tuple((k,M[4*i+k]) for k in range(4) if M[4*i+k]) for i in range(4)) for M in bank]
    for depth in range(1,5):
        nxt=set()
        for B in current:
            for A in sparse:
                nxt.add(tuple(sum(v*B[4*k+j] for k,v in row) for row in A for j in range(4)))
        current=nxt;union.update(nxt)
        depths.append({'depth':depth,'exact_length_distinct_products':len(nxt),'eigen1_products':sum(det_shift(M,1)==0 for M in nxt)})
    eligible={};eigen1=0;kernel_dims={};all_integer_eigen_pairs=0;zero_output=0;sign_changing=0
    for M in sorted(union):
        if det_shift(M,1):continue
        eigen1+=1;poly=charpoly(M)
        need(evaluate(poly,1)==0,'characteristic determinant discrepancy')
        roots=integer_roots_ge2(poly,max(sum(M[4*i:4*i+4]) for i in range(4)))
        for B in roots:
            need(det_shift(M,B)==0,'integer root determinant discrepancy')
            basis=nullspace(M,B);kernel_dims[len(basis)]=kernel_dims.get(len(basis),0)+1;all_integer_eigen_pairs+=1
            need(len(basis)==1,'higher-dimensional B-eigenspace requires separate cone audit')
            a=basis[0];out=2*a[0]+a[2]+a[3]
            if not out:zero_output+=1;continue
            a=tuple(v/out for v in a)
            if min(a)<0:sign_changing+=1;continue
            key=(B,a)
            need(key not in eligible or eligible[key]==M,'two distinct matrices share normalized B-ray')
            eligible[key]=M
    need(len(eligible)==84,'eligible ray count changed')
    need(len(set(eligible.values()))==84,'eligible matrices not distinct')
    # Actual affine survivor, with signed intercept allowed.
    def recipe_matrix(n,code,z):
        rows=[0]*n;bit=0
        for i in range(n):
            for j in range(i+1,n):
                a,b=(i,j) if code>>bit&1 else (j,i);rows[a]|=1<<b;bit+=1
        rec=next(r for r in raw if r['local']==rows and r['new_root']==z)
        return tuple(v for row in conjugate(rec['matrix']) for v in row)
    W=tuple(int(i==j) for i in range(4) for j in range(4))
    parent=[0];cut=0
    for p in ((3,0,2),(4,29,3),(3,5,2)):
        W=matrix_mul(recipe_matrix(*p),W)
        n,code,z=p;rows=[0]*n;bit=0
        for i in range(n):
            for j in range(i+1,n):
                u,v=(i,j) if code>>bit&1 else (j,i);rows[u]|=1<<v;bit+=1
        parent,cut=build(parent,cut,rows,z)
    actual=direct_state(parent,cut)
    need(actual==(17,53,36,0) and literal(parent,cut)==actual,'literal first-W witness')
    need((W[0],W[3],W[12],W[15])==(64,18,56,17),'positive survivor restriction')
    need(all(W[4*i+j]==0 for i in (1,2) for j in (0,3)),'survivor plane not invariant')
    a=(Fraction(9,26),0,0,Fraction(4,13));b=(Fraction(-9,26),0,0,Fraction(9,13))
    need(2*a[0]+a[2]+a[3]==1 and 2*b[0]+b[2]+b[3]==0,'signed affine normalization')
    def mv(M,v):return tuple(sum(M[4*i+j]*v[j] for j in range(4)) for i in range(4))
    need(mv(W,a)==tuple(80*x for x in a),'survivor slope')
    need(mv(W,b)==tuple(b[i]-27*a[i] for i in range(4)),'signed nonzero digit')
    state=(0,0,0,1);family=[]
    for k in range(6):
        F=2*state[0]+state[2]+state[3]
        need(F==(52*80**k+27)//79 and (52*80**k+27)%79==0,'survivor family')
        family.append(F);state=mv(W,state)
    print(json.dumps({'status':'PASS','generator_count':32,'depths':depths,'union_distinct_products_lengths_1_to_4':len(union),'union_eigen1_products':eigen1,'all_integer_B_eigen_pairs':all_integer_eigen_pairs,'B_eigenspace_dimension_histogram':kernel_dims,'zero_output_rays':zero_output,'sign_changing_output_normalized_rays':sign_changing,'eligible_positive_output_rays':len(eligible),'eligible_radices':sorted(set(B for B,a in eligible)),'normalized_B_ray_injective':True,'survivor_literal_first_native_state':actual,'survivor_parent_masks':parent,'survivor_cut':cut,'survivor_signed_intercept':[str(v) for v in b],'survivor_digit':-27,'survivor_radix':80,'survivor_first_values':family,'scope':'Every product of1..4 operations in the exact32-generator bank; all integer B>=2 tested by exact Sturm isolation. No two distinct NONZERO digits on a common output-normalized affine ray in this finite macro universe. Zero-digit macros and longer words remain outside this obstruction.'},indent=2,sort_keys=True))
if __name__=='__main__':main()
