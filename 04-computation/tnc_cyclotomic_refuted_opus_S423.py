import sympy as sp
u,a=sp.symbols('u a')
def CT(R,N,mm):
    return sp.Poly(sp.expand(R**mm),u).coeff_monomial(u**(N*mm))
print("BROAD SEARCH: trinomials whose CT(m0) is a GENUINE polynomial in a (>=2 terms),")
print("and whether ALL roots have |root|=1 (roots of unity).  Hunting a counterexample.")
print()
tunable=[]
for N in [2,3,4]:
    for j in range(1,9):
        for d in range(j+1,12):
            if j==N or d==N: continue
            R=1+a*u**j+u**d
            m0=next((m for m in range(1,16) if sp.expand(CT(R,N,m))!=0 and not sp.expand(CT(R,N,m)).is_number),None)
            if m0 is None: continue
            c=sp.expand(CT(R,N,m0)); p=sp.Poly(c,a)
            # strip a^k
            co=p.all_coeffs()
            nz=[x for x in co if x!=0]
            if len(nz)>=2:   # genuine polynomial -> tunable
                # roots (nonzero)
                pr=p
                while pr.eval(0)==0 and pr.total_degree()>0: pr=sp.Poly(sp.quo(pr.as_expr(),a),a)
                rts=[complex(r) for r in sp.nroots(pr)]
                mags=[abs(r) for r in rts]
                allunit=all(abs(m-1)<1e-9 for m in mags)
                tunable.append((N,j,d,m0,sp.factor(c),allunit,[round(m,4) for m in mags]))
print(f"tunable trinomial patterns found: {len(tunable)}")
for N,j,d,m0,c,unit,mags in tunable:
    tag = "ALL |root|=1 (cyclotomic)" if unit else "*** NON-UNIT ROOT -- CYCLOTOMIC ROUTE REFUTED ***"
    print(f"   N={N} (-{N},{j-N},{d-N}) m0={m0}: {str(c)[:30]:30s} |roots|={mags}  {tag}")
print()
counter=[x for x in tunable if not x[5]]
if counter:
    print(f"COUNTEREXAMPLES to cyclotomic (non-root-of-unity tuned points): {len(counter)}")
    print("=> the single-shot cyclotomic route as stated is FALSE. Tuned points need NOT be")
    print("   roots of unity; |a| is a ratio of multinomial coefficients.")
else:
    print("NO counterexample: every tunable trinomial has root-of-unity tuned points.")
    print("=> cyclotomic route SURVIVES for trinomials.")
