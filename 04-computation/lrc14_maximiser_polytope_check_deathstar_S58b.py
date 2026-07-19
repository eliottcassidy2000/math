from fractions import Fraction as F
import itertools
# their exact bad_from_g
def bad_kps(g):
    cuts=[]
    for x in g:
        h=F(7,6)*x; a=max(h-F(1,6),F(0)); b=min(h,F(1))
        if b>a: cuts.append((a,b))
    cuts.sort(); cur=F(0); L=F(0)
    for a,b in cuts:
        if a>cur and a-cur>L: L=a-cur
        if b>cur: cur=b
    if 1-cur>L: L=1-cur
    return L<=F(1,6)
# my derived polytope (sorted right-edges e_i=(7/6)g_i): e1<=1/3, e2-e1<=1/3, e3-e2<=1/3, e3>=5/6
def bad_poly(g):
    e=sorted(F(7,6)*x for x in g)
    return e[0]<=F(1,3) and e[1]-e[0]<=F(1,3) and e[2]-e[1]<=F(1,3) and e[2]>=F(5,6)
# verify they agree on a grid (away from clamping the two should match; check ALL and report mismatches)
mism=0; tot=0
for a in range(0,42):
    for b in range(0,42):
        for c in range(0,42):
            g=(F(a,42),F(b,42),F(c,42)); tot+=1
            if bad_kps(g)!=bad_poly(g): mism+=1
print("bad_kps vs bad_poly on (1/42)^3 grid: %d mismatches / %d"%(mism,tot))
# exact sojourn via fine cells: for direction d, scan u, exact bad, measure via interval merge on fine grid
def sojourn_exact(DS, D):
    # D = common denominator; find bad u-subintervals by checking bad at cell midpoints on a partition
    # partition at all j/d_i and at fine grid; use N=lcm-ish
    N=D
    inb=[bad_poly(tuple((-d*F(t,N))%1 for d in DS)) for t in range(N)]
    # measure = (#bad)/N is approximate; refine by using N large. Return fraction.
    return F(sum(inb),N)
for DS in [(1,2,3),(1,3,2),(3,9,6),(1,2,7),(2,4,7)]:
    s=sojourn_exact(DS, 2*3*4*5*7*9)  # 7560
    print("d=%s sojourn(N=7560)=%s=%.6f  (2/21=%.6f)"%(DS,s,float(s),2/21))
