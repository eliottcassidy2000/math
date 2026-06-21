from fractions import Fraction as Fr
# Carefully understand the trapezoid h_j and confirm TV <= 2/7 always (so Koksma OK).
# h_j(a) = overlap length of [a, a+L) with 1-periodic strips I_j=[j/7,(j+1)/7) (mod 1).
# Width of strip = 1/7. Window width = L = p/(7q).
# As a slides over period 1, the window [a,a+L) overlaps the (periodic) strip.
# In ONE period the strip appears... actually strips of width 1/7 spaced 1/7 apart = the
# whole line is covered by 7 strips per unit, sector j strip is one of every 7.
# So in [0,1) there is exactly ONE copy of strip I_j (length 1/7) PLUS wrap.
# h_j(a): trapezoid rising from 0 to min(L,1/7) and back, over a within one period.
# Wait: there are MULTIPLE strips within distance L if L>1/7? No: strips of sector j are
# spaced 1 apart (one per unit interval, at [j/7,(j+1)/7)+n). Distance between consecutive
# sector-j strips = 1. Since L = p/(7q) <= 43/140 < 1, the window hits at most ONE j-strip
# (plus its periodic translate when wrapping). 
# Over a in [0,1): window [a,a+L). The function h_j(a) is the standard "tent/trapezoid"
# convolution of two intervals of widths 1/7 and L, which has TV = 2*min(1/7,L).
# Since L>1/7 (because p/q>1 => L=p/(7q)>1/(7)... wait L=p/(7q), p/q>1 => p>q => L>1/7? 
#   L>1/7 iff p/(7q)>1/7 iff p/q>1. YES. So min(1/7,L)=1/7, TV=2/7.
# BUT my audit found one sector with TV<2/7. WHY? Because of the PERIODIC WRAP: when the
# strip and window both wrap mod 1, on a FULL period the rises/falls can MERGE if L+1/7>1,
# but here L+1/7 < 43/140+20/140=63/140<1, no merge. So why TV<2/7 for one j?
# Let me recompute very carefully, sampling h_j on a dense rational grid and computing TV.

def hj(a, L, j):
    lo=a; hi=a+L; tot=Fr(0)
    for n in range(-2,4):
        s0=Fr(j,7)+n; s1=Fr(j+1,7)+n
        left=max(lo,s0); right=min(hi,s1)
        if right>left: tot+=right-left
    return tot

def TV_exact(L,j,steps=2*3*5*7*4):
    # exact piecewise-linear TV via breakpoints
    bps=set()
    for n in range(-2,4):
        for edge in (Fr(j,7)+n, Fr(j+1,7)+n):
            for shift in (Fr(0), -L):
                bps.add((edge+shift)%1)
    bps=sorted(set(b for b in bps if 0<=b<1))+[Fr(1)]
    tv=Fr(0)
    prev=hj(Fr(0),L,j)
    for b in bps[1:]:
        cur=hj(b%1 if b<1 else Fr(0),L,j)
        # but at a=1 -> wrap to 0; use limit. Evaluate at b for b<1, and close loop to a=0 value.
        val = hj(b if b<1 else Fr(0), L, j)
        tv+=abs(val-prev); prev=val
    return tv

for (p,q) in [(3,2),(2,1),(5,3),(43,20)]:
    L=Fr(p,7*q)
    print(f"p/q={p}/{q}, L={L}, L>1/7={L>Fr(1,7)}, L+1/7={L+Fr(1,7)}")
    for j in range(7):
        # sample plateau max and TV
        # max of h_j = min(L,1/7) normally; check
        import fractions
        grid=[Fr(t,420) for t in range(420)]
        vals=[hj(g,L,j) for g in grid]
        mx=max(vals)
        tv=TV_exact(L,j)
        print(f"   j={j}: max h_j={mx}={float(mx):.4f}  TV={tv}={float(tv):.4f}  (2/7={float(Fr(2,7)):.4f})")
