# Is the stacked-jump edge EMPTY for the GMC radial D=(1-b t)^2-4 gamma r t^2?
# Fold t-value at pinch root r (root of g(r):=r b'(r)^2 = gamma):  t_*(r)=1/phi(r), phi(r)=b(r)-2 r b'(r).
# "Stacked" = two DISTINCT pinch roots r1!=r2 at the SAME gamma with the SAME fold t-value => phi(r1)=phi(r2).
# So stacking <=> the planar curve r|->(g(r),phi(r)) SELF-INTERSECTS (same g AND same phi at r1!=r2).
import numpy as np
def analyze(bcoef,label,rmax=8.0,N=400000):
    b=np.polynomial.Polynomial(bcoef); bp=b.deriv()
    r=np.linspace(0.02,rmax,N)
    g=r*bp(r)**2; phi=b(r)-2*r*bp(r)
    # detect self-intersections: find r1<r2 with (g,phi) close. Coarse bucket then refine.
    found=[]
    # scan pairs via sorting by g then checking phi matches at equal g (two branches crossing)
    idx=np.argsort(g)
    gs=g[idx]; ps=phi[idx]; rs=r[idx]
    for i in range(len(gs)-1):
        # look ahead a small window in g-sorted order for a phi crossing at nearly equal g with distant r
        for j in range(i+1,min(i+60,len(gs))):
            if abs(rs[i]-rs[j])<0.05: continue
            if abs(gs[i]-gs[j])<3e-3 and abs(ps[i]-ps[j])<3e-3:
                found.append((round(rs[i],3),round(rs[j],3),round(gs[i],3),round(ps[i],4)))
    # dedup by (r1,r2) rounded
    uniq={}
    for f in found:
        key=(round(f[0],1),round(f[1],1))
        uniq[key]=f
    U=list(uniq.values())
    print(f"[{label}] b={bcoef}: self-intersections of (g,phi) [=(gamma, 1/t_*)] with r1!=r2:")
    if U:
        for f in U[:6]:
            print(f"    r1={f[0]}, r2={f[1]}: gamma={f[2]}, phi={f[3]}  => STACKED (same gamma, same fold t_*=1/phi)")
        print(f"   => STACKING OCCURS: {len(U)} distinct real self-crossing(s). Edge NOT empty for this b.")
    else:
        print("   => no real self-intersection found in scan.")
    return U
analyze([0,0,1.0],"b=r^2")
analyze([0,1.0,1.0],"b=r+r^2")
analyze([0,0,0,1.0],"b=r^3")
analyze([0,1.0,0,-0.3],"b=r-0.3r^3")
analyze([0,2.0,-1.0,0.1],"b=2r-r^2+0.1r^3")
print("\nNote: g'=b'(b'+2rb''), phi'=-(b'+2rb''), so g'=-b'*phi' => the curve (g,phi) is generically immersed")
print("and self-intersects (transversal, distinct slopes -b'(r1)!=-b'(r2)). So stacking is generic, not empty.")
