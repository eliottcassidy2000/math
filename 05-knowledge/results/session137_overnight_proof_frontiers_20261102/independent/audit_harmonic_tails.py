"""Independent exact harmonic maxima, full ceilings, and terminal mass audit."""
from audit_harmonic import safe,require,F
from math import gcd
import json

def run():
    maxima=ceilings=masses=0;isolated_extra=0
    for q in range(3,7):
        for mask in range(1<<(q-1)):
            A={i for i in range(1,q) if mask>>(i-1)&1}
            C=list(range(1,q))+[q+i for i in sorted(A)]
            d=next((j for j in range(1,q) if j not in A),q);n=q+d
            E=[u+q*(u in A) for u in range(1,q) if gcd(u,q)==1];s=min(E)
            for m in range(1,5):
                if q<=5 and m<=3:
                    for k in range(1,n):
                        V=sorted(set(C+[j*q*m for j in range(1,k+1)]));H=max(V)
                        if d*m>=s:
                            r=min(F(m,q*m+s),F(1,k+1))
                            require(bool(safe(V,r)),'max lower witness')
                            require(not safe(V,r+F(1,8*H*H)),'max upper rational certificate')
                        else:
                            N=n//gcd(q*m,n)
                            require(bool(safe(V,F(1,n)))==(k<N),'below-floor decision')
                        maxima+=1
                for h in range(1,d+1):
                    nh=q+h;k=nh-1;r=F(1,nh)
                    V=sorted(set(C+[j*q*m for j in range(1,k+1)]));H=max(V)
                    actual=safe(V,r)
                    points=set();weighted=F(0)
                    for e in E:
                        a=pow(e,-1,q)
                        for j in range(1,m*h//e+1):
                            if gcd(j,nh)!=1:continue
                            t=F(a,q)-F(j,q*m*nh);points.update((t,1-t))
                            ell=pow(j,-1,nh)
                            weighted+=F(2,q*m)*(F(1,ell)+F(1,nh-ell))
                    if h==d and gcd(q*m,nh)==1:
                        grid={F(b,nh) for b in range(1,nh) if gcd(b,nh)==1}
                        isolated_extra+=len(grid-points);points|=grid
                        harmonic=sum((F(1,u) for u in range(1,nh) if gcd(u,nh)==1),F(0))
                        ell=pow(m,-1,nh)
                        weighted+=F(2,q*m)*(harmonic-F(1,ell)-F(1,nh-ell))
                    require(actual==[(t,t) for t in sorted(points)],'complete ceiling maximizers')
                    require(bool(points)==(m*h>=s),'ceiling criterion')
                    ceilings+=1
                    if points:
                        eta=F(1,16*H*H)
                        growth=sum((b-a for a,b in safe(V,r-eta)),F(0))
                        require(growth/eta==weighted,'terminal mass coefficient');masses+=1
    print(json.dumps({'scope':'FINITE-EXACT all promotion subsets q3..6,m1..4, every ceiling; q3..5,m1..3 all k<n for maximum certificates',
                      'maximum_checks':maxima,'ceiling_checks':ceilings,'terminal_mass_checks':masses,
                      'primitive_grid_extra_points_counted':isolated_extra},sort_keys=True,indent=2))
if __name__=='__main__':run()
