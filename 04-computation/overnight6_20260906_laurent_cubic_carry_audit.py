"""Independent source/resultant interpolation certificate for the cubic family.

No producer code imports or copied characteristic coefficients. Reads only
the published coefficient record and checks it using independent raw fibers
and resultants at 19 points, with the degree bound proved in the audit note.
"""
from pathlib import Path
from math import factorial,gcd,prod
import json,hashlib,sys
import sympy as S
T,X,G=S.symbols('t x g')
CHECKS=0
def need(ok,label):
    global CHECKS
    CHECKS+=1
    if not ok:raise RuntimeError(label)

def raw(g,m):
    # Independently solve in the gamma multiplicity; use charge 2y+3z.
    ans=[]
    for z in range(m+1):
        rem=21*m-3*g*z
        if rem<0 or rem%(2*g):continue
        y=rem//(2*g);x=m-y-z
        if x>=0:ans.append(((x,y,z),factorial(m)//(factorial(x)*factorial(y)*factorial(z))))
    return ans

def main():
    sys.stdout.reconfigure(newline='\n')
    base=Path(__file__).parent
    transcript=base/'overnight6_20260906_laurent_cubic_carry.out'
    if not transcript.exists():
        transcript=base.parent/'05-knowledge/results/overnight6_20260906_laurent_cubic_carry.out'
    line=next(z for z in transcript.read_text(encoding='utf-8').splitlines()
              if z.startswith('characteristic_positive_shift_certificates '))
    records=json.loads(line.split(' ',1)[1])
    wanted={}
    for k in (1,2,3):
        rec=records[str(k)];cs=list(map(S.Rational,rec['coefficients_descending']))
        d=S.Integer(rec['denominator'])
        need(rec['shift_g']==11 and len(cs)==6*k+1 and d>0 and all(x>0 for x in cs),'positive full coefficient record')
        wanted[k]=S.Poly.from_list(cs,G).as_expr().subs(G,G-11)/d
    # Normalize raw factorial coefficients by their leading coefficient.
    psym=S.Poly(sum(S.prod(G-7-z for z in range(3-j))*S.Rational(factorial(7),
                       factorial(9-3*j)*factorial(1+2*j))*T**j for j in range(4)),T)
    for j in range(4):
        coeff=S.Poly(S.expand(psym.nth(j).subs(G,G+11)),G)
        need(all(c>0 for c in coeff.all_coeffs()),'positive full-domain first coefficients')
    disc=S.Poly(S.expand(S.discriminant(psym.as_expr(),T).subs(G,G+11)),G)
    need(all(c>0 for c in disc.all_coeffs()),'positive full-domain cubic discriminant')
    sampled={k:[] for k in (1,2,3)}
    for g in range(11,30):
        one,two=raw(g,g),raw(g,2*g)
        need([v for v,c in one]==[(g-10+j,9-3*j,1+2*j) for j in range(4)],'first source fiber')
        need([v for v,c in two]==[(2*g-21+j,21-3*j,2*j) for j in range(8)],'doubled source fiber')
        first=S.Poly(sum(c*T**j for j,(v,c) in enumerate(one)),T)
        p=first.monic();K=prod(range(2*g-13,2*g+1))
        need(S.expand(p.as_expr()-psym.as_expr().subs(G,g))==0,'independent monic source formula')
        q=S.Poly(sum(S.Rational(c,K)*T**j for j,(v,c) in enumerate(two)),T)
        # The resultant has roots Qbar(t_i)/t_i. Normalize its nonzero
        # leading coefficient explicitly, independent of Sylvester sign.
        ch=S.Poly(S.resultant(p.as_expr(),T*X-q.as_expr(),T),X).monic()
        need(ch.degree()==3 and ch.LC()==1,'monic resultant characteristic')
        for k in (1,2,3):
            value=ch.nth(3-k);sampled[k].append((g,value))
            need(value==wanted[k].subs(G,g),'published polynomial exact specialization')
        # Independent Sturm root count, rather than producer discriminant formula.
        need(p.count_roots(-S.oo,0)==3 and S.gcd(p,p.diff()).degree()==0,'simple negative first roots')
        need(ch.count_roots(0,S.oo)==0 and ch.nth(0)>0,'no nonnegative response root')
        # Sign reversal at the three first roots is encoded by tau<0.
        need(S.gcd(p,q).degree()==0,'no common first/doubled zero')
        if gcd(g,21)==1:
            need(all(not raw(g,m) for m in range(1,2*g) if m!=g),'no earlier or intermediate source return')
    for k in (1,2,3):
        interpolated=S.interpolate(sampled[k],G)
        need(S.Poly(interpolated,G).degree()==6*k,'interpolated exact degree')
        need(S.expand(interpolated-wanted[k])==0,'complete degree-bounded coefficient identity')
    # Independent symbolic denominator cancellation and weighted-degree premise.
    p0=psym.nth(0)
    q0=S.prod(2*G-14-j for j in range(7))/factorial(21)
    ratio=S.cancel(q0/p0)
    need(S.denom(ratio).is_Integer and S.Poly(ratio,G).degree()==4,'sole inverse denominator cancels')
    # Consumer hostiles: gcd is required; trace and norm alone are insufficient.
    need(next(m for m in range(1,13) if raw(12,m))==4,'gcd-dropped first mass')
    hostile=S.Poly((X+5)*(X-1)**2,X)
    need(hostile.nth(2)>0 and hostile.nth(0)>0 and hostile.nth(1)<0,'trace and norm miss middle coefficient')
    print('Independent 19-parameter raw-fiber/resultant interpolation: degrees6/12/18 PASS')
    print('All3 full shifted characteristic polynomials strictly positive on g>=11 PASS')
    print('Source fibers, Sturm roots, exact denominator cancellation, gcd and missing-middle hostiles PASS')
    print('Polynomial identities use analytic weighted-degree bounds; no bounded extrapolation')
    print('PASS_OPTIMIZATION_LIVE_GATES',CHECKS)
    print('published_record_sha256',hashlib.sha256(json.dumps(records,sort_keys=True,separators=(',',':')).encode()).hexdigest())

if __name__=='__main__':main()
