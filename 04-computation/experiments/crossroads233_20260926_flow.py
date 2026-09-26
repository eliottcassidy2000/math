"""Two-cutoff incompatibility and finite cutoff kernels in global P2 pairings.

Exact integers/Fractions only. Run normally or with -O; no dependencies.
The proof is in crossroads233_20260926_flow.md.
"""
from array import array
from fractions import Fraction
from functools import lru_cache
from itertools import product
import json


def require(test,message):
    if not test:raise AssertionError(message)


def constraints(X):
    equal=[];less=[]
    for i in range(2,X+1):
        if i%2==0:
            j=3*i//2
            if j<=X:equal.append((i,j))
        else:
            j,k=(3*i-1)//2,(3*i+1)//2
            if j<=X:less.append((j,i))
            if k<=X:less.append((i,k))
    return equal,less


def optimize(X,small=0):
    zero=[0]*(X+1);one=[0]*(X+1)
    for i in range(X,1,-1):
        a,b=0,1+int(i<=small)
        if i%2==0:
            j=3*i//2
            if j<=X:a+=one[j];b+=zero[j]
        else:
            j,k=(3*i-1)//2,(3*i+1)//2
            if j<=X:a+=zero[j];b+=min(zero[j],one[j])
            if k<=X:a+=min(zero[k],one[k]);b+=one[k]
        zero[i]=a;one[i]=b
    bits=[0]*(X+1)
    if X>=2:
        stack=[(2,int(one[2]<zero[2]))]
        while stack:
            i,b=stack.pop();bits[i]=b
            if i%2==0:
                j=3*i//2
                if j<=X:stack.append((j,1-b))
            else:
                j,k=(3*i-1)//2,(3*i+1)//2
                if j<=X:stack.append((j,0 if not b else int(one[j]<zero[j])))
                if k<=X:stack.append((k,1 if b else int(one[k]<zero[k])))
    eq,le=constraints(X)
    require(all(bits[i]+bits[j]==1 for i,j in eq),'complement clause')
    require(all(bits[i]<=bits[j] for i,j in le),'monotone clause')
    value=min(zero[2],one[2]) if X>=2 else 0
    require(sum(bits)+sum(bits[:small+1])==value,'objective reconstruction')
    return value,bits


@lru_cache(None)
def direct_tree(h,i,weighted):
    w=2 if weighted and h>=2 else 1
    if h==0:return 0,w,0
    if i%2==0:
        x,y,_=direct_tree(h-1,3*i//2,weighted)
        a,b=y,w+x
        return a,b,min(a,b)-min(x,y)
    a,b,_=direct_tree(h-1,(3*i-1)//2,weighted)
    x,y,_=direct_tree(h-1,(3*i+1)//2,weighted)
    zero,one=a+min(x,y),w+min(a,b)+y
    return zero,one,min(zero,one)-min(a,b)-min(x,y)


def coefficients(weighted,depth=20):
    differences=array('i',[1]);s=Fraction(0);rows=[];audits=0
    for h in range(1,depth+1):
        modulus=len(differences);new=array('i');total=0
        w=2 if weighted and h>=2 else 1
        bound=2 if weighted else 1
        for i in range(2*modulus):
            if i%2==0:
                x=differences[(3*i//2)%modulus]
                delta,toll=w-x,min(max(x,0),w)
            else:
                x=differences[((3*i-1)//2)%modulus]
                y=differences[((3*i+1)//2)%modulus]
                delta=w+min(x,0)+max(y,0)
                toll=min(max(-x,0),w+max(y,0))
            require(abs(delta)<=bound*(h+1) and 0<=toll<=bound*h,'depth bound')
            if h<=7:
                a,b,t=direct_tree(h,i,weighted)
                require((b-a,t)==(delta,toll),'independent full-tree audit')
                audits+=1
            new.append(delta);total+=toll
        differences=new
        require(sum(differences)==w*len(differences),'mean-difference identity')
        s+=Fraction(total,3**(h+1))
        rows.append((h,total,s))
    return rows,audits


@lru_cache(None)
def kernel_direct_tree(h,i,root_weights):
    """Independent fixed-root optimization, without difference recurrences."""
    w=root_weights[h]
    if h==0:return 0,w,0
    if i%2==0:
        x,y,_=kernel_direct_tree(h-1,3*i//2,root_weights)
        a,b=y,w+x
        return a,b,min(a,b)-min(x,y)
    a,b,_=kernel_direct_tree(h-1,(3*i-1)//2,root_weights)
    x,y,_=kernel_direct_tree(h-1,(3*i+1)//2,root_weights)
    zero,one=a+min(x,y),w+min(a,b)+y
    return zero,one,min(zero,one)-min(a,b)-min(x,y)


def kernel_coefficients(kernel,depth):
    require(all(a>=0 for a in kernel) and sum(kernel)>0,'nonnegative nonzero kernel')
    W=sum(kernel);normalizer=sum((a*Fraction(2,3)**k for k,a in enumerate(kernel)),Fraction(0))
    root_weights=tuple(sum(kernel[:h+1]) for h in range(depth+1))
    differences=array('q',[root_weights[0]]);partial=Fraction(0);rows=[];audits=0
    for h in range(1,depth+1):
        modulus=len(differences);new=array('q');total=0;w=root_weights[h]
        for i in range(2*modulus):
            if i%2==0:
                x=differences[(3*i//2)%modulus]
                delta,toll=w-x,min(max(x,0),w)
            else:
                x=differences[((3*i-1)//2)%modulus]
                y=differences[((3*i+1)//2)%modulus]
                delta=w+min(x,0)+max(y,0)
                toll=min(max(-x,0),w+max(y,0))
            require(abs(delta)<=W*(h+1) and 0<=toll<=W*h,'kernel depth bound')
            if h<=7:
                a,b,t=kernel_direct_tree(h,i,root_weights)
                require((b-a,t)==(delta,toll),'kernel independent full-tree audit')
                audits+=1
            new.append(delta);total+=toll
        differences=new
        require(sum(differences)==w*len(differences),'kernel mean-difference identity')
        partial+=Fraction(total,3**(h+1));rows.append((h,total,partial))
    return rows,normalizer,audits


def kernel_prefix_optimum(X,kernel):
    cutoffs=[2**k*X//3**k for k in range(len(kernel))]
    weights=[sum(a for a,cut in zip(kernel,cutoffs) if i<=cut) for i in range(X+1)]
    zero=[0]*(X+1);one=[0]*(X+1)
    for i in range(X,1,-1):
        a,b=0,weights[i]
        if i%2==0:
            j=3*i//2
            if j<=X:a+=one[j];b+=zero[j]
        else:
            j,k=(3*i-1)//2,(3*i+1)//2
            if j<=X:a+=zero[j];b+=min(zero[j],one[j])
            if k<=X:a+=min(zero[k],one[k]);b+=one[k]
        zero[i]=a;one[i]=b
    return min(zero[2],one[2]),weights


def audit_finite_kernels():
    audits=0
    for m in range(2,9):
        kernel=[3**j*2**(m-1-j) for j in range(m)]
        rows,normalizer,count=kernel_coefficients(kernel,16);audits+=count
        require(normalizer==m*2**(m-1),'adjacent equal-mass normalization')
        bound=rows[-1][2]/normalizer
        print(json.dumps({'kernel_spacing':1,'cuts':m,'depth':16,
                          'density_lower_exact':str(bound),'density_lower_decimal':float(bound)}))
    spaced=[0]*11
    for j in range(6):spaced[2*j]=9**j*4**(5-j)
    rows,normalizer,count=kernel_coefficients(spaced,20);audits+=count
    require(normalizer==6*4**5,'spacing-two equal-mass normalization')
    bound=rows[-1][2]/normalizer
    require(bound==Fraction(4915373767757,16067102519808),'spacing-two comparison')
    print(json.dumps({'kernel_spacing':2,'cuts':6,'depth':20,
                      'density_lower_exact':str(bound),'density_lower_decimal':float(bound)}))
    kernel=[3**j*2**(7-j) for j in range(8)]
    rows,normalizer,count=kernel_coefficients(kernel,20);audits+=count
    require(kernel==[128,192,288,432,648,972,1458,2187],'winning kernel')
    bound=rows[-1][2]/normalizer
    require(sum(kernel)==6305 and normalizer==1024,'winning kernel normalization')
    require(rows[-1][2]*3**21==3286041555236,'winning aggregate numerator')
    require(rows[-1][1]==3259399538,'winning final toll coefficient')
    require(bound==Fraction(821510388809,2677850419968),'winning rational certificate')
    cases=0
    for X in range(2,13):
        value,weights=kernel_prefix_optimum(X,kernel);eq,le=constraints(X);best=sum(weights)
        for choice in product([0,1],repeat=X-1):
            bits=[0,0]+list(choice);cases+=1
            if all(bits[i]+bits[j]==1 for i,j in eq) and all(bits[i]<=bits[j] for i,j in le):
                best=min(best,sum(w*b for w,b in zip(weights,bits)))
        require(value==best,'kernel exhaustive-prefix audit')
    print(json.dumps({'kernel_spacing':1,'cuts':8,'depth':20,'kernel':kernel,
                      'root_weight_bound':sum(kernel),'normalizer':str(normalizer),
                      'toll_coefficients':[row[1] for row in rows],
                      'aggregate_numerator':3286041555236,
                      'density_lower_exact':str(bound),'density_lower_decimal':float(bound),
                      'independent_kernel_full_tree_cases':audits,
                      'independent_kernel_binary_assignments':cases}))


def main():
    cases=0
    for X in range(2,16):
        small=4*X//9;eq,le=constraints(X);best=X+small
        for choice in product([0,1],repeat=X-1):
            bits=[0,0]+list(choice);cases+=1
            if all(bits[i]+bits[j]==1 for i,j in eq) and all(bits[i]<=bits[j] for i,j in le):
                best=min(best,sum(bits)+sum(bits[:small+1]))
        require(best==optimize(X,small)[0],'exhaustive/DP discrepancy')
    print(json.dumps({'independent_binary_assignments':cases,'cutoffs_checked':[2,15]}))
    first=None
    for X in range(1,234):
        small=4*X//9
        if optimize(X,small)[0]>optimize(X)[0]+optimize(small)[0]:first=(X,small);break
    require(first==(5,2),'smallest two-cutoff conflict')
    print(json.dumps({'first_conflict_large_cutoff':5,'small_cutoff':2,'separate_optima':[0,1],'joint_optimum':2}))
    for X in [233,2250,22500,225000]:
        small=4*X//9;val,_=optimize(X,small)
        print(json.dumps({'large_cutoff':X,'small_cutoff':small,'exact_joint_optimum':val,
                          'separate_optima':[optimize(small)[0],optimize(X)[0]],'joint_ratio':str(Fraction(val,X+small))}))
    old,old_audits=coefficients(False);new,new_audits=coefficients(True)
    old_upper=old[-1][2]+23*Fraction(2,3)**21
    require(old_upper==Fraction(3089623223,10460353203),'inherited alpha upper bound')
    for h,A,s in new:
        beta=Fraction(9,13)*s
        print(json.dumps({'depth':h,'weighted_toll_sum':A,'beta_partial_exact':str(beta),'beta_partial_decimal':float(beta)}))
    beta12=Fraction(9,13)*new[11][2];beta20=Fraction(9,13)*new[-1][2]
    require(beta12>old_upper,'twelve terms fail to separate densities')
    require(beta20==Fraction(4538724347,15109399071),'weighted twenty-term coefficient sum')
    forced_upper=Fraction(13,4)*beta20-Fraction(9,4)*old_upper
    oscillation=forced_upper-old_upper
    require(forced_upper==Fraction(13417603,43046721),'limsup certificate')
    require(oscillation==Fraction(170854306,10460353203),'oscillation certificate')
    print(json.dumps({'independent_full_tree_cases':old_audits+new_audits,
                      'alpha_upper':str(old_upper),'natural_density_lower':str(beta20),
                      'forced_limsup_for_alpha_attainers':str(forced_upper),'forced_oscillation':str(oscillation),
                      'forced_limsup_decimal':float(forced_upper),'forced_oscillation_decimal':float(oscillation)}))
    audit_finite_kernels()
    print('ALL CHECKS PASSED')


if __name__=='__main__':main()
