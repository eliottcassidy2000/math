import sympy as sp
u,a=sp.symbols('u a')
def CT(R,N,mm):
    return sp.Poly(sp.expand(R**mm),u).coeff_monomial(u**(N*mm))
print("STRUCTURE: in gauge r0=r_d=1, CT(Lambda^m) = sum_y c_{m,y} a^y with c_{m,y} = POSITIVE")
print("multinomials.  A nullcone point a* is a common root of these POSITIVE-coeff polynomials.")
print("Roots of a positive-coeff poly lie in an annulus (Enestrom-Kakeya): min(c_y/c_{y+1})")
print("<= |a*| <= max(c_y/c_{y+1}) on the support.  If the annuli of CT(m0) and CT(2m0) are")
print("DISJOINT, no common root -> coprime -> TNC.  (Roots of unity = the ANGULAR positions;")
print("the radius comes from multinomial magnitudes = the RECURRENCE growth.)")
print()
def support_coeffs(c):
    p=sp.Poly(sp.expand(c),a)
    d={e[0]:co for e,co in zip(p.monoms(),p.coeffs())}
    ys=sorted(d)
    return ys,[d[y] for y in ys]
def annulus(c):
    ys,cs=support_coeffs(c)
    if len(cs)<2: return None
    # consecutive ratios on the actual (possibly gapped) support
    rs=[]
    for i in range(len(cs)-1):
        gap=ys[i+1]-ys[i]
        rs.append((abs(cs[i]/cs[i+1]))**sp.Rational(1,gap))
    return min(rs), max(rs)
print("   pattern         CT(m0) annulus |a|       CT(2m0) annulus |a|      disjoint?")
allsep=True; cnt=0
for N in [2,3,4]:
    for j in range(1,9):
        for d in range(j+1,12):
            if j==N or d==N: continue
            R=1+a*u**j+u**d
            m0=next((m for m in range(1,16) if sp.expand(CT(R,N,m))!=0 and not sp.expand(CT(R,N,m)).is_number),None)
            if m0 is None: continue
            c0=sp.expand(CT(R,N,m0))
            if len([x for x in sp.Poly(c0,a).all_coeffs() if x!=0])<2: continue
            cnt+=1
            an0=annulus(c0); an1=annulus(sp.expand(CT(R,N,2*m0)))
            if an0 is None or an1 is None: continue
            lo0,hi0=float(an0[0]),float(an0[1]); lo1,hi1=float(an1[0]),float(an1[1])
            disjoint = hi0<lo1 or hi1<lo0
            allsep=allsep and disjoint
            print(f"   N={N}(-{N},{j-N},{d-N}): [{lo0:.3f},{hi0:.3f}]        "
                  f"[{lo1:.3f},{hi1:.3f}]      {disjoint}")
print(f"\n   tunable trinomials: {cnt};  Enestrom-Kakeya annuli DISJOINT (proves coprime): {allsep}")
