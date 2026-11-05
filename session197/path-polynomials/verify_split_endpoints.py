#!/usr/bin/env python3
from verify_marked_contexts import *
from fractions import Fraction

@lru_cache(None)
def from_root(a,u):
    n=len(a);d=[[0]*n for _ in range(1<<n)];d[1<<u][u]=1
    for mask in range(1,1<<n):
        for i,c in enumerate(d[mask]):
            if not c:continue
            todo=a[i]&~mask
            while todo:
                bit=todo&-todo;todo-=bit;d[mask|bit][bit.bit_length()-1]+=c
    return d

def linked(a,u,v):
    if u==v:
        mask=((1<<len(a))-1)^(1<<u);return all_bins(a)[mask]
    bins=all_bins(a);d=from_root(a,u);full=(1<<len(a))-1;out=[0]
    for mask in range(1,1<<len(a)):
        if d[mask][v]:out=add(out,scale(bins[full^mask],d[mask][v]))
    return out

def on_columns(cols,p):return add(*[scale(cols[m],v) for m,v in enumerate(p) if v])

def top_shift_from_columns(cols,b):
    # Centered operator moments sum_s b_s(x)s^m, then Lagrange isolate s=b.
    moments=[]
    for m in range(2*b+1):
        moments.append(add(*[scale(mul(mon(m-j),cols[j]),comb(m,j)*(-1)**(m-j)) for j in range(m+1)]))
    lag=[Fraction(1)]
    for s in range(-b,b):lag=mul(lag,[Fraction(-s,b-s),Fraction(1,b-s)])
    out=add(*[scale(moments[j],v) for j,v in enumerate(lag)])
    require(all(Fraction(v).denominator==1 for v in out),'integral top shift')
    return [int(v) for v in out]

def run():
    contexts=[(tournament(n,m),r) for n in range(1,4) for m in range(1<<(n*(n-1)//2)) for r in range(n)]
    inputs=[(tournament(n,m),u,v) for n in range(2,4) for m in range(1<<(n*(n-1)//2)) for u,v in permutations(range(n),2)]
    checks=0;virtual=0;singletons=0
    for q,r in contexts:
        A=[];B=[]
        for m in range(5):
            A.append(bin_polynomial(clone(q,r,m)) if len(q)-1+m else [1])
            a=transitive(m+2);out=substitute(q,r,a);offset=len(q)-1
            B.append(linked(out,offset+m+1,offset))
            require(linked(a,m+1,0)==[0] and marked(a,(m+1,0),('S','E'))==mon(m),'split test module')
        for par,k in product((0,1),range(2)):
            require(positive_extended(on_columns(B,basis(k,par))),'positive connecting block');virtual+=1
        for a,u,v in inputs:
            out=substitute(q,r,a);offset=len(q)-1;j=linked(a,u,v);k=marked(a,(u,v),('S','E'))
            require(linked(out,offset+u,offset+v)==add(on_columns(A,j) if j!=[0] else [0],on_columns(B,k)),('triangular full profile',q,r,a,u,v));checks+=1
        for n in range(1,4):
            for mask in range(1<<(n*(n-1)//2)):
                a=tournament(n,mask);out=substitute(q,r,a);offset=len(q)-1
                for u in range(n):
                    f=linked(a,u,u);require(linked(out,offset+u,offset+u)==on_columns(A,f),'unsplittable singleton');singletons+=1
    wchecks=0
    for n in range(2,5):
        for mask in range(1<<(n*(n-1)//2)):
            a=tournament(n,mask);w=wrapper(a)
            for u,v in permutations(range(n),2):
                j=linked(a,u,v);k=marked(a,(u,v),('S','E'))
                require(linked(w,u,v)==add(D(j),k),'exact W upper block')
                require(marked(w,(u,v),('S','E'))==markedD(k,1,1),'exact W lower block');wchecks+=1
    hostile=[]
    for n in (2,3):
        a=transitive(n);u=n-1;v=0;j=linked(a,u,v);k=marked(a,(u,v),('S','E'));out=linked(wrapper(a),u,v)
        require(j==[0] and out==k,'j-only hostile');hostile.append({'input_order':n,'j':j,'k':k,'j_after_W':out})
    from verify_matching_amplitude import balanced,cross_matching
    high=[]
    for q in (balanced(1,1,0),balanced(2,0b1001,0),balanced(2,0b1111,0),balanced(3,0b100010001,0)):
        nu,mu=cross_matching(q,0);b=nu-1;cols=[];offset=len(q)-1
        for m in range(2*b+1):
            out=substitute(q,0,transitive(m+2));cols.append(linked(out,offset+m+1,offset))
        top=top_shift_from_columns(cols,b);want=[mu*nu]
        for z in range(1,nu):want=mul(want,[z,1])
        require(top==want,('balanced open coefficient',q,nu,mu,top,want));high.append({'rows':q,'nu':nu,'mu':mu,'forward_coefficient':top})
    print(json.dumps({'rooted_contexts':len(contexts),'full_split_context_checks':checks,'formal_connecting_basis_checks':virtual,'unsplittable_singleton_checks':singletons,'explicit_W_pair_checks':wchecks,'same_j_different_outputs':hostile,'higher_matching_controls':high},sort_keys=True,indent=2))
if __name__=='__main__':run()
