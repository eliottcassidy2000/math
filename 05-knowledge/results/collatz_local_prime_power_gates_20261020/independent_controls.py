#!/usr/bin/env python3
"""Independent exact controls: no producer or repository mathematics imports."""
from fractions import Fraction
from itertools import combinations, product
from math import gcd
from hashlib import sha256
import json

counts = {}
def check(ok, tag):
    counts[tag] = counts.get(tag, 0) + 1
    if not ok:
        raise RuntimeError(tag)

def affine(w):
    a,b=Fraction(1),Fraction(0)
    for bit in w:
        a = a*(3 if bit else 1)/2
        b = (b*(3 if bit else 1)+bit)/2
    return a,b

def B(w):
    positions=[i for i,b in enumerate(w) if b]
    return sum((1<<i)*3**(len(positions)-j-1) for j,i in enumerate(positions))

def rot(w,r):return w[r:]+w[:r]

def core(u,v):
    a=0;b=len(u)
    while a<b and u[a]==v[a]:a+=1
    while b>a and u[b-1]==v[b-1]:b-=1
    return u[a:b],v[a:b]

for q in range(2,11):
    for w in product((0,1),repeat=q):
        p=sum(w)
        if p in (0,q):continue
        D=2**q-3**p; carry=B(w);g=gcd(abs(D),carry)
        a,b=affine(w)
        check(a==Fraction(3**p,2**q) and b==Fraction(carry,2**q),'affine_composition')
        check(b/(1-a)==Fraction(carry,D),'rational_fixed_point')
        for r in range(1,q):
            v=rot(w,r);a=sum(w[:r]);delta=B(v)-carry
            check(2**r*B(v)==3**a*carry+D*B(w[:r]),'rotation_affine_identity')
            check(gcd(abs(D),delta)==g*gcd(abs(D)//g,3**a-2**r),'arbitrary_rotation_content')
        u,v=core(w,rot(w,1))
        check(gcd(abs(D),B(u)-B(v))==g,'single_core_exact')
        ell=len(u);h=sum(u)
        check(0<abs(B(u)-B(v))<=(2**(ell-h)-1)*(3**h-2**h),'core_width_bound')

for q in range(3,501):
    for d in range(1,q//2+1):
        D=2**q-9;carry=3+2**d
        g=gcd(abs(D),carry)
        check(g==gcd(carry,2**(q-2*d)-1),'all_period_two_ones')
        if q>=5:check(g<D,'two_one_positive_exclusion')

# All two-successive-unequal-swap histories through q=12, retaining reversals.
for q in range(3,13):
    for p in range(1,q):
        if gcd(p,q)>1:continue
        c=tuple((i+1)*p//q-i*p//q for i in range(q))
        D=2**q-3**p
        for i in range(q):
            if c[i]==c[(i+1)%q]:continue
            v=list(c);v[i],v[(i+1)%q]=v[(i+1)%q],v[i];v=tuple(v)
            for j in range(q):
                if v[j]==v[(j+1)%q]:continue
                w=list(v);w[j],w[(j+1)%q]=w[(j+1)%q],w[j];w=tuple(w)
                exception={(i+1)%q,(j+1)%q};g=gcd(abs(D),B(w));clean=0
                inverse=pow(p,-1,q)
                for t in range(q-1):
                    r=t*inverse%q;s=(t+1)*inverse%q
                    if r in exception or s in exception:continue
                    clean+=1
                    perturb=[]
                    for cut in [r,s]:
                        total=0
                        for before,at in [(c,i),(v,j)]:
                            loc=(at-cut)%q;z=rot(before,cut)
                            check(loc<q-1,'exception_cut_removed')
                            total+=(1 if z[loc] else -1)*2**loc*3**sum(z[loc+2:])
                        check(B(rot(w,cut))-B(rot(c,cut))==total,'two_swap_perturbation')
                        perturb.append(total)
                    smooth=B(rot(c,r))-B(rot(c,s))
                    z=abs(smooth)
                    while z and z%2==0:z//=2
                    while z and z%3==0:z//=3
                    check(z==1,'balanced_edge_smooth')
                    delta=B(rot(w,r))-B(rot(w,s))
                    check(delta==smooth+perturb[0]-perturb[1] and delta!=0,'five_term_edge')
                    check(delta%g==0,'bounded_swap_gcd_gate')
                check(clean>=max(0,q-5),'clean_edge_bound')

# Independent literal positional sums, with no optimized swap update.
hist={};digest=sha256();occ=nonunit=integral=0
for q in range(3,121):
    for p in range(1,q):
        if gcd(q,p)!=1:continue
        c=tuple((i+1)*p//q-i*p//q for i in range(q));D=2**q-3**p
        for i in range(q):
            j=(i+1)%q
            if c[i]==c[j]:continue
            w=list(c);w[i],w[j]=w[j],w[i]
            g=gcd(abs(D),B(w));occ+=1;nonunit+=g>1;integral+=abs(D)>1 and g==abs(D)
            digest.update(f'{q}:{p}:{i}:{g}\n'.encode())
check(integral==0,'direct_swap_no_integral_abs_D_gt_1')
print(json.dumps({'status':'PASS','checks':counts,'total_checks':sum(counts.values()),
 'universes':{'all_words':'nonconstant binary words 2<=q<=10, all rotation shifts',
 'two_ones':'3<=q<=500, every 1<=d<=floor(q/2)',
 'two_swaps':'all two-successive unequal cyclic swaps of coprime mechanical words, 3<=q<=12',
 'direct_swaps':'all directed unequal cyclic swaps of coprime mechanical words, 3<=q<=120, repeated outputs retained'},
 'direct_swap_occurrences':occ,'direct_swap_nonunit':nonunit,'direct_swap_integral_abs_D_gt_1':integral,
 'direct_swap_digest_sha256':digest.hexdigest()},indent=2,sort_keys=True))
