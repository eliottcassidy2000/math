"""Ternary rigid palettes: Mahler sections, discrepancy, exact circulant norms."""
from pathlib import Path
from math import prod
from collections import Counter
import importlib.util
import json

CHECKS=0
def need(ok,msg):
    global CHECKS
    CHECKS+=1
    if not ok: raise RuntimeError(msg)

def lucas(n):
    a,b=2,1
    for _ in range(n): a,b=b,a+b
    return a

def digitsign(n):
    sign=1
    while n:
        if n%3==2: sign=-sign
        n//=3
    return sign

def bareiss(matrix):
    a=[row[:] for row in matrix]; n=len(a); prev=1; sign=1
    for k in range(n-1):
        if not a[k][k]:
            j=next((j for j in range(k+1,n) if a[j][k]),None)
            if j is None: return 0
            a[k],a[j]=a[j],a[k]; sign=-sign
        pivot=a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                z=a[i][j]*pivot-a[i][k]*a[k][j]
                need(z%prev==0,'Bareiss exact division'); a[i][j]=z//prev
            a[i][k]=0
        prev=pivot
    return sign*a[-1][-1]

def trim(a):
    while a and a[-1]==0: a.pop()
    return a

def gcd_degree(a,b,p):
    a=trim([x%p for x in a]); b=trim([x%p for x in b])
    while b:
        r=a[:]; inv=pow(b[-1],-1,p)
        while len(r)>=len(b):
            c=r[-1]*inv%p; shift=len(r)-len(b)
            for i,x in enumerate(b): r[i+shift]=(r[i+shift]-c*x)%p
            trim(r)
        a,b=b,r
    return len(a)-1

def gaussian_eval(coeff):
    re=im=0
    for n,x in enumerate(coeff):
        if n%4==0: re+=x
        elif n%4==1: im+=x
        elif n%4==2: re-=x
        else: im-=x
    return re,im

def main():
    rows=[]; d=[1]; determinant=1
    for k in range(11):
        n=len(d)
        need(all(v==digitsign(i) for i,v in enumerate(d)),'digit-rule equals substitution')
        prefixes=[0]
        for x in d: prefixes.append(prefixes[-1]+x)
        need(prefixes[-1]==1,'imbalance one')
        need(max(prefixes)==k+1 and min(prefixes)==min(0,2-k),'exact prefix extrema')
        re,im=gaussian_eval(d)
        need(re*re+im*im==5**k,'quarter-turn phase resonance')
        if k: determinant=determinant**3*lucas(n)
        row={'depth':k,'N':n,'prefix_min':min(prefixes),'prefix_max':max(prefixes),
             'weight_minus_signs':d.count(-1),'value_at_i':[re,im]}
        if k<=5:
            other=prod(lucas(3**r)**(3**(k-r)) for r in range(1,k+1))
            need(determinant==other,'determinant recurrence/product')
            row['determinant']=str(determinant)
            if k<=4:
                actual=bareiss([[d[(j-i)%n] for j in range(n)] for i in range(n)])
                need(actual==determinant,'literal integer circulant determinant')
            # Rank via polynomial gcd holds even in inseparable characteristic3.
            ranks={str(p):n-gcd_degree(d,[-1]+[0]*(n-1)+[1],p) for p in [2,3,19]}
            need(ranks['2']==1 and ranks['3']==n,'Frobenius rank contrast')
            row['finite_field_ranks']=ranks
        if k<=6:
            # In F3, (1+z-z²)D(z)^2=1 through degree N-1.
            square=[sum(d[j]*d[i-j] for j in range(i+1))%3 for i in range(n)]
            lhs=[(square[i]+(square[i-1] if i>=1 else 0)-(square[i-2] if i>=2 else 0))%3 for i in range(n)]
            need(lhs==[1]+[0]*(n-1),'Mahler/Frobenius quadratic equation')
        rows.append(row)
        d=d+d+[-x for x in d]
    # Clock-sensitive hostile: length9 primitive necklaces versus full W2 freedom.
    primitive_necklaces=(2**9-2**3)//9
    need(primitive_necklaces==56,'global primitive count')
    spec=importlib.util.spec_from_file_location('inherited_canonical',Path(__file__).resolve().parent.parent/'weighted_crt_clocks_20261028/check_canonical.py')
    old=importlib.util.module_from_spec(spec); spec.loader.exec_module(old)
    d9=tuple(digitsign(i) for i in range(9)); group=old.w.permutations_a(2); t=old.w.negation(2)
    spectrum=Counter(); central=[]
    for g in group:
        v=old.moved(d9,g); determinant=bareiss([[v[(j-i)%9] for j in range(9)] for i in range(9)])
        spectrum[determinant]+=1
        if all(g[t[i]]==t[g[i]] for i in range(9)): central.append((g,v,determinant))
    need(dict(spectrum)=={4864:30,1792:18,256:15,1024:18},'fixed-clock hostile spectrum')
    profiles=[]
    for _,v,_ in central:
        profile=[]
        for g,_,_ in central:
            moved=old.moved(v,g)
            profile.append(bareiss([[moved[(j-i)%9] for j in range(9)] for i in range(9)]))
        profiles.append(tuple(profile))
    need(len(set(profiles))==3,'future re-markings recover all three relative frames')
    result={'status':'FINITE-EXACT controls for supplied all-depth proofs',
            'rows':rows,'global_primitive_length9':primitive_necklaces,'hierarchically_rigid_orbits':2,
            'fixed_cyclic_clock_hostile':{'A2_orbit_determinants':dict(spectrum),
                   'C2_determinants':[x[2] for x in central],'C2_future_profiles':profiles},
            'checks':CHECKS}
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__': main()
