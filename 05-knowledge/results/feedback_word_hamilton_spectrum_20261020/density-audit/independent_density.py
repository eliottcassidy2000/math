#!/usr/bin/env python3
"""Clean-room controls with rational logarithm enclosures, no producer imports."""
from bisect import bisect_right
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product, permutations
import json


def check(p,why):
    if not p:
        raise RuntimeError(why)


def H(word):
    # Direct closed combinatorial sum, not the producer's transfer recurrence.
    return 1+sum(2 if j==i+1 else 1 << (j-i-1)
                 for i in range(len(word)) for j in range(i+1,len(word))
                 if word[i] and not word[j])


def literal(word):
    # Apex is vertex zero. Chain labels are 1,...,m.
    def arc(a,b):
        return bool(word[b-1]) if a==0 else not word[a-1] if b==0 else a<b
    return sum(all(arc(p[i],p[i+1]) for i in range(len(word)))
               for p in permutations(range(len(word)+1)))


def integer(word):
    return sum(b << (len(word)-1-i) for i,b in enumerate(word))


def dual(word):
    return tuple(1-b for b in word[::-1])


def identity_check(word):
    m=len(word); scale=1 << m
    n=integer(word);r=integer(dual(word))
    adjacent=sum(word[i]==1 and word[i+1]==0 for i in range(m-1))
    backward=sum(1 << (m+j-i-1) for i in range(m) for j in range(i)
                 if word[i] and not word[j])
    error=scale*H(word)-n*r
    check(error==scale*(1+adjacent)-backward,('scaled identity',word))
    check(2*abs(error)<=(m+2)*scale,('uniform product error',word))
    check(2*adjacent<=m and 2*backward<=m*scale,('ingredients including empty word',word))


@lru_cache(None)
def log_interval(z):
    """Rational enclosure of ln(z), for 1<=z<=2, via atanh series."""
    z=Q(z);check(1<=z<=2,('log domain',z))
    a=(z-1)/(z+1)
    terms=28
    lo=2*sum((a**(2*j+1)/Q(2*j+1) for j in range(terms)),Q(0))
    remainder=2*a**(2*terms+1)/(Q(2*terms+1)*(1-a*a))
    return lo,lo+remainder


def F_interval(z):
    z=Q(z)
    if z<=Q(1,4):
        return Q(0),Q(0)
    if z<=Q(1,2):
        lo,hi=log_interval(4*z)
        return 4*z*lo-4*z+1,4*z*hi-4*z+1
    if z<1:
        lo,hi=log_interval(1/z)
        return 4*z-3+4*z*lo,4*z-3+4*z*hi
    return Q(1),Q(1)


def main():
    report={}
    all_controls=0
    for m in range(11):
        for word in product((0,1),repeat=m):
            identity_check(word);all_controls+=1
    report['all_binary_words_identity_controls_including_empty']=all_controls
    ln2lo,ln2hi=log_interval(Q(2))
    by_length={};moments=[];prefix_controls=literal_controls=0
    cdf_controls=0
    for m in range(2,15):
        s=(m-2)//2;scale=1 << m
        counts=Counter();prefixes=Counter();fixed=0
        for mid in product((0,1),repeat=m-2):
            word=(1,)+mid+(0,);dw=dual(word);h=H(word)
            identity_check(word)
            check(h==H(dw),('duality',word))
            check(Q(scale,4)<h<scale,('length cutoff',word))
            n,r=integer(word),integer(dw)
            cell_n=n >> (m-s-1); cell_r=r >> (m-s-1)
            check(1 << s<=cell_n<1 << (s+1) and 1 << s<=cell_r<1 << (s+1),'first digit')
            prefixes[cell_n,cell_r]+=1
            # The common-prefix interval is exactly the deterministic coupling cell.
            width=1 << (m-s-1)
            check(cell_n*width<=n<(cell_n+1)*width and cell_r*width<=r<(cell_r+1)*width,'coupling cells')
            fixed+=word==dw;counts[h]+=1;prefix_controls+=1
            if m<=6:
                check(h==literal(word),('literal apex permutations',word));literal_controls+=1
        total=1 << (m-2)
        check(len(prefixes)==1 << (2*s) and set(prefixes.values())=={1 << (m-2-2*s)},('independent prefixes',m))
        check(fixed==(0 if m%2 else 1 << (m//2-1)),('fixed points',m))
        delta=Q(1,1 << s)+Q(m+2,2*scale)
        mean=Q(sum(h*c for h,c in counts.items()),total)
        if m>=3:
            check(mean==Q(9*scale,16)+Q(1,2),('exact mean',m))
        mean_normalized=mean/scale
        second=Q(sum(h*h*c for h,c in counts.items()),total*scale*scale)
        variance=second-mean_normalized**2
        check(abs(mean_normalized-Q(9,16))<=delta,('coupled mean',m))
        check(abs(second-Q(49,144))<=2*delta,('coupled second moment',m))
        check(abs(variance-Q(55,2304))<=4*delta,('coupled variance',m))
        seen=0;worst=Q(0)
        for h,c in sorted(counts.items()):
            flo,fhi=F_interval(Q(h,scale))
            # Left/right limits at every jump certify the supremum distance.
            for empirical in [Q(seen,total),Q(seen+c,total)]:
                distance=max(abs(empirical-flo),abs(empirical-fhi))
                worst=max(worst,distance)
                check(distance<=4*ln2lo*delta,('rationally certified uniform CDF bound',m,h))
                cdf_controls+=1
            seen+=c
        moments.append({'m':m,'exact_mean':str(mean),'normalized_variance':str(variance),
                        'cdf_distance_rational_upper_approx':float(worst)})
        by_length[m]=sorted(counts.elements())
    report['reduced_word_prefix_and_identity_controls']=prefix_controls
    report['literal_permutation_controls']=literal_controls
    report['rigorous_rational_CDF_jump_checks']=cdf_controls
    report['moments_and_CDF']=moments

    # Exact coefficient cancellation in basis 1,t,t*ln(t),t*ln(2).
    constant=[1,0,0,0]
    Ft2=[-3,2,-2,2]
    twice_Ft4=[2,-2,2,0]
    check([a+b+c for a,b,c in zip(constant,Ft2,twice_Ft4)]==[0,0,0,2], 'symbolic dyadic phase cancellation')
    report['phase_identity']='exact coefficient cancellation'

    bank=sorted(h for values in by_length.values() for h in values)
    spectrum=sorted(set(bank)|{1})
    cutoffs=[Q(n) for n in range(4,8192)]
    cutoffs += [Q((1 << L)*r,16) for L in range(2,13) for r in range(17,32)]
    for X in cutoffs:
        L=(X.numerator//X.denominator).bit_length()-1
        W=bisect_right(bank,X)
        check(W==(1 << (L-1))-1+bisect_right(by_length[L+1],X)+bisect_right(by_length[L+2],X),('exact two levels',X))
        QL=Q((1 << (L//2))+(1 << ((L+1)//2)))+Q(2*L+7,8)
        E_upper=max(abs(Q(W+1)-X*ln2lo),abs(Q(W+1)-X*ln2hi))
        check(E_upper<=4*ln2lo*QL,('rationally certified all-cutoff error',X))
        spectrum_upper=X*ln2lo/2+2*ln2lo*QL+(1 << ((L+2)//2-1))
        check(bisect_right(spectrum,X)<=spectrum_upper,('distinct values upper bound',X))
    report['exact_or_rationally_certified_cutoffs']=len(cutoffs)
    report['word_ratios']=[{'X':1 << L,'W':bisect_right(bank,1 << L),
                            'W_over_X':str(Q(bisect_right(bank,1 << L),1 << L))} for L in range(4,13)]
    report['log2_enclosure_width']=str(ln2hi-ln2lo)
    report['status']='PASS'
    print(json.dumps(report,indent=2,sort_keys=True))


if __name__=='__main__':
    main()
