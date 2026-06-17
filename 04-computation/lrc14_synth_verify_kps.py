from fractions import Fraction as Fr
from math import gcd
from functools import reduce
from itertools import combinations

def danger_arcs(v):
    w=Fr(1,14*v); A=[]
    for k in range(v+1):
        c=Fr(k,v); lo=c-w; hi=c+w
        if lo<0: A+=[(Fr(0),hi),(1+lo,Fr(1))]
        elif hi>1: A+=[(lo,Fr(1)),(Fr(0),hi-1)]
        else: A.append((lo,hi))
    return A

def L_exact(S):
    A=[]; [A.extend(danger_arcs(v)) for v in S]; A=sorted((a,b) for a,b in A if b>a)
    tot=Fr(0); cl=ch=None
    for a,b in A:
        if ch is None: cl,ch=a,b
        elif a<=ch: ch=max(ch,b)
        else: tot+=ch-cl; cl,ch=a,b
    if ch is not None: tot+=ch-cl
    return 1-tot

def lcm(a,b): return a*b//gcd(a,b)
def lcmset(S): return reduce(lcm,S,1)

# 1. The two claimed tight configs
AP=list(range(1,14))
SP=[1,2,3,4,5,6,7,8,9,10,11,13,24]
print("L(AP {1..13})      =", L_exact(AP), " lcm=", lcmset(AP))
print("L(sporadic 12->24) =", L_exact(SP), " lcm=", lcmset(SP))

# 2. The 1/1260 extremizer
CH=[1,2,3,4,5,6,7,8,9,10,11,13,36]
print("L({1..11,13,36})   =", L_exact(CH), " == 1/1260 ?", L_exact(CH)==Fr(1,1260))


# 3. max-min := max_tau min_v ||v tau||  via exact candidate enumeration.
# The maximizer of min_v||v tau|| occurs at a breakpoint where two danger boundaries cross
# or at a danger-arc boundary. Candidates: tau = (2k+-1)/(2*14*?) ... simpler: the optimum
# of min_v ||v tau|| is attained at tau where some ||v tau|| boundary meets. Use the standard:
# critical tau are rationals p/q with q <= max(S) where derivative sign changes. We enumerate
# tau = j/(2*L) for fine grid won't be exact. Instead: candidate optima are at points where
# two functions ||a tau||, ||b tau|| are equal & both maximal => tau = (i+-... ). 
# Robust exact method: the function f(tau)=min_v||v tau|| is piecewise linear with breakpoints
# at tau = k/v (peaks of individual sawtooth) and crossings. Local maxima of f occur at crossings
# of two sawtooths: a*tau - i = -(b*tau - j) type. Enumerate all candidate crossing taus.

def maxmin_exact(S):
    cands=set([Fr(0)])
    M=max(S)
    # individual sawtooth peaks ||v tau|| local maxima at tau=(2k+1)/(2v)
    for v in S:
        for k in range(2*v+1):
            t=Fr(2*k+1,2*v)
            if 0<=t<1: cands.add(t)
    # crossings of two sawtooths: v*t - kv = +-(w*t - kw)
    Sl=list(S)
    for ai in range(len(Sl)):
        for bi in range(len(Sl)):
            a=Sl[ai]; b=Sl[bi]
            for ka in range(-1,a+1):
                for kb in range(-1,b+1):
                    # a t - ka = w t - kb  => (a-b) t = ka-kb
                    if a!=b:
                        t=Fr(ka-kb,a-b)
                        if 0<=t<1: cands.add(t)
                    # a t - ka = -(b t - kb) => (a+b) t = ka+kb
                    t2=Fr(ka+kb,a+b)
                    if 0<=t2<1: cands.add(t2)
    best=Fr(0)
    for t in cands:
        m=min(min(frac:=(v*t - int(v*t)), 1-frac) for v in S)
        # careful with Fraction floor
        vals=[]
        for v in S:
            x=v*t
            fl=x - Fr(int(x)) if x>=0 else x-Fr(int(x))+ (0 if x==int(x) else 1)
            # use proper nearest-integer distance
            n=round(float(x))
            # exact nearest: 
            lower=x.numerator//x.denominator
            for cand_int in (lower,lower+1):
                pass
            d1=x-Fr(lower); d2=Fr(lower+1)-x
            vals.append(min(d1,d2))
        m=min(vals)
        if m>best: best=m
    return best

print()
print("max-min AP        =", maxmin_exact(AP), "== 1/14 ?", maxmin_exact(AP)==Fr(1,14))
print("max-min sporadic  =", maxmin_exact(SP), "== 1/14 ?", maxmin_exact(SP)==Fr(1,14))
print("max-min 12->36    =", maxmin_exact(CH), "(should be > 1/14, loose config)")
