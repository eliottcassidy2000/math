#!/usr/bin/env python3
"""Standalone exact controls for the credited concurrent extensions."""
import json, math, itertools
from fractions import Fraction
from collections import Counter
from pathlib import Path

def require(p,why):
    if not p: raise RuntimeError(why)

def jacobi(a,n):
    require(n>0 and n%2,'odd positive modulus'); a%=n; s=1
    while a:
        while a%2==0:
            a//=2
            if n%8 in (3,5): s=-s
        a,n=n,a
        if a%4==n%4==3: s=-s
        a%=n
    return s if n==1 else 0

def factor(n):
    f={}; p=2
    while p*p<=n:
        while n%p==0: f[p]=f.get(p,0)+1; n//=p
        p+=1
    if n>1:f[n]=f.get(n,0)+1
    return f

def sigma(n):
    return math.prod((p**(e+1)-1)//(p-1) for p,e in factor(n).items())

def masks(ns,basis):
    return [sum((factor(n).get(p,0)%2)<<i for i,p in enumerate(basis)) for n in ns]

def rowcode(columns,r):
    return {sum(((z&c).bit_count()%2)<<i for i,c in enumerate(columns)) for z in range(1<<r)}

def main():
    out={}; bank=[4,2,5,3,7,13]; sums=list(map(sigma,bank))
    require(sums==[7,3,6,4,8,14],'sigma bank')
    columns=masks(sums,[2,3,7]); code=rowcode(columns,3)
    require(len(code)==8,'sigma rank three')
    for word in code:
        c={b:(word>>i)&1 for i,b in enumerate(bank)}
        require(c[3]==0 and c[2]^c[5]^c[7]==0 and c[4]^c[7]^c[13]==0,'complete sigma relations')
    counts=Counter(tuple(jacobi(a,q) for a in sums) for q in range(3,168,4) if math.gcd(q,168)==1)
    require(len(counts)==8 and set(counts.values())=={3},'CRT sigma classes')
    require(masks([sigma(11),sigma(2)],[2,3,7])==[2,2],'eleven duplicate coordinate')
    out['sigma_bank']={'bases':bank,'sigma':sums,'columns':columns,'colors':8,'CRT_modulus':168,'classes_per_color':3}
    seven=[2,5,7,19,29,59,89]; cols=masks([p+1 for p in seven],[2,3,5])
    require(set(cols)==set(range(1,8)),'all seven nonzero columns')
    dual=rowcode(cols,3); kernel=set()
    for w in range(128):
        total=0
        for i,c in enumerate(cols):
            if w>>i&1:total^=c
        if total==0:kernel.add(w)
    require(len(kernel)==16 and Counter(w.bit_count() for w in kernel)=={0:1,3:7,4:7,7:1},'Hamming weights')
    require(Counter(w.bit_count() for w in dual)=={0:1,4:7},'simplex weights')
    require(all((a&b).bit_count()%2==0 for a in kernel for b in dual),'orthogonality')
    for i in range(7):
        for j in range(i):
            require(set(Counter(((w>>i)&1,(w>>j)&1) for w in dual).values())=={2},'pair independence')
    out['Hamming_simplex']={'bank':seven,'columns':cols,'kernel_weight_counts':dict(Counter(w.bit_count() for w in kernel)),'rowspace_weight_counts':dict(Counter(w.bit_count() for w in dual))}
    s=[2]
    for _ in range(8):s.append(s[-1]**2-s[-1]+1)
    for j in range(1,len(s)):
        require(s[j]==1+math.prod(s[:j]),'Sylvester product')
        require(s[j]%4==3,'odd Sylvester mod four')
        for i in range(1,j):
            require(math.gcd(s[i],s[j])==1 and jacobi(s[j],s[i])==1 and jacobi(s[i],s[j])==-1,'all odd Jacobi orientation')
    out['Sylvester']={'terms_checked':len(s),'odd_terms_checked':len(s)-1,'includes_composite_1807':s[4]==1807,'orientation':'later beats earlier'}
    exps={2:5,3:7,5:3,7:4,11:1,13:2,29:1,41:1,59:1,71:1}
    def exponent(b,es):return min(es.get(p,0)//e for p,e in factor(b).items())
    f={b:exponent(b,exps) for b in bank+[11]}
    require(f=={4:2,2:5,5:3,3:7,7:4,13:2,11:1},'composite-base cycle')
    adjusted=dict(exps); adjusted[2]-=1
    g={b:exponent(b,adjusted) for b in bank+[11]}
    require(g=={4:2,2:4,5:3,3:7,7:4,13:2,11:1},'central-two adjustment')
    out['valuation_maps']={'R':f,'R_divided_by_2':g,'terminal':1}
    # Fermat identity is exact for all m; finite prime controls do not assert later primality.
    out['Fermat_prime_checks']=[]
    for m in range(5):
        p=2**(2**m)+1; d=2**(m+1)
        require(2**d-1==p*(p-2),'Fermat exact factorization')
        require(pow(2,d,p)==1 and pow(2,d//2,p)==p-1,'Fermat order certificate')
        require(pow(2,p-1,p*p)!=1,'non-Wieferich finite guard')
        out['Fermat_prime_checks'].append({'p':p,'order':d})
    for p,d in [(1093,364),(3511,1755)]:
        require(pow(2,d,p*p)==1 and pow(2,d,p**3)!=1,'exact valuation two')
        require(all(pow(2,d//r,p)!=1 for r in factor(d)),'minimal order')
    require(pow(2,182,1093**2)==1093**2-1,'negative one in doubling subgroup')
    out['Paley_Gram_moduli']=[]
    for q in [3,7,11,19,43]:
        K=[[jacobi(y-x,q) for y in range(q)] for x in range(q)]
        require(all(sum(row)==0 for row in K),'Paley row sum')
        for x in range(q):
            for z in range(q):
                require(sum(K[x][y]*K[z][y] for y in range(q))==(q-1 if x==z else -1),'Paley exact Gram')
        out['Paley_Gram_moduli'].append(q)
    colored_counts={}
    for n in range(1,5):
        pairs=list(itertools.combinations(range(n),2)); permutations=list(itertools.permutations(range(n)))
        for q in (1,2):
            canonical=set()
            for bits in range(1<<len(pairs)):
                edges={(i,j) if bits>>k&1 else (j,i) for k,(i,j) in enumerate(pairs)}
                for colors in itertools.product(range(q),repeat=n):
                    forms=[]
                    for perm in permutations:
                        renamed={(perm[i],perm[j]) for i,j in edges}
                        renamed_colors=[0]*n
                        for i in range(n):renamed_colors[perm[i]]=colors[i]
                        mask=sum(((i,j) in renamed)<<k for k,(i,j) in enumerate(pairs))
                        forms.append((tuple(renamed_colors),mask))
                    canonical.add(min(forms))
            fixed_sum=0
            for perm in permutations:
                seen=set(); lengths=[]
                for i in range(n):
                    if i in seen:continue
                    j=i; length=0
                    while j not in seen:seen.add(j); length+=1; j=perm[j]
                    lengths.append(length)
                if any(l%2==0 for l in lengths):continue
                e=sum((l-1)//2 for l in lengths)+sum(math.gcd(a,b) for a,b in itertools.combinations(lengths,2))
                fixed_sum+=(2**e)*(q**len(lengths))
            formula=Fraction(fixed_sum,math.factorial(n))
            require(formula==len(canonical),'colored Burnside enumeration')
            if n==3:require(formula==Fraction(4*q**3+2*q,3),'A3 polynomial')
            if n==4:require(formula==Fraction(8*q**4+4*q**2,3),'A4 polynomial')
            colored_counts[f'n={n},q={q}']=len(canonical)
    out['colored_tournament_counts']=colored_counts
    out['status']='PASS: exact finite controls; general claims proved in companion note'
    Path(__file__).with_name('concurrent_extensions_receipt.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
