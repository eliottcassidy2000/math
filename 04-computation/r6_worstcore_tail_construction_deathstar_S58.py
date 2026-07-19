from fractions import Fraction as F
P=[1,2,4,7,9,11,12]
def nd(x): x=x%1; return min(x,1-x)
al,ah=F(71,154),F(13,28)

print("=== PART 1: tail construction rigor (b>=400) ===")
# (a) candidate-exists formula: (b+8)*(1/308 - 18/(77b)) >= 1 for all b>=400 ?
def lb_len(b): return (b+8)*(F(1,308)-F(18,77*b))
print("  (b+8)|A_good| lower bound at b=400:",float(lb_len(400)),"  b=380:",float(lb_len(380)),"  b=375:",float(lb_len(375)))
print("  monotone increasing? diff 400->100000:",float(lb_len(400)),"->",float(lb_len(100000)))
# verify >=1 for all b>=400 by checking it's increasing and >=1 at 400 (derivative>0)
ok_formula=all(lb_len(b)>=1 for b in range(400,2000))  # spot; it's (b+8)/308 - 18(b+8)/(77b), increasing
print("  (b+8)|A_good|>=1 for b in [400,2000]:",ok_formula)
# (b) verify a valid t* actually exists AND margins hold, for b in [400, 4000] exactly
def check_construct(b):
    K=b+8; rho_needed=F(9,77*b)
    # candidate half-integers in A_good = (al+9/(77b), ah-9/(77b))
    glo=al+F(9,77*b); ghi=ah-F(9,77*b)
    nlo=int(glo*K-F(1,2)); nhi=int(ghi*K-F(1,2))+2
    for n in range(nlo,nhi+1):
        t=(F(n)+F(1,2))/K
        if not (glo<=t<=ghi): continue
        # verify killer b margin and all: min over speeds of (nd(vt)-1/14)/v >= rho_needed
        good=True; rho=None
        for v in P+[b,b+2,b+4,b+6,b+8]:
            d=nd(v*t)-F(1,14)
            if d<=0: good=False;break
            r=d/v
            if rho is None or r<rho: rho=r
        if good and rho>=rho_needed and 2*rho>F(1,7*(b+8)):
            return True
    return False
fails=[b for b in range(400,4001) if not check_construct(b)]
print("  construction failures in [400,4000]:",len(fails),fails[:10])
# (c) confirm killer b is binding with margin exactly 9/(77b) at a sample
b=1000; K=b+8; n=int((al+ah)/2*K-F(1,2)); t=(F(n)+F(1,2))/K
if al<t<ah:
    dists={v:nd(v*t) for v in [b,b+2,b+4,b+6,b+8]}
    print("  sample b=1000: killer dists from Z:",{k:float(v) for k,v in dists.items()},"(min killer b =%.4f, >=29/154=%.4f)"%(float(min(dists.values())),29/154))
