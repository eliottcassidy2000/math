#!/usr/bin/env python3
"""Exact identities and independent finite-field boundary-fibre audits."""
import json
import sympy as s

w,v,t,x,y,z,C,A,B = s.symbols('w v t x y z C A B')


def check(ok, message):
    if not ok:
        raise RuntimeError(message)


def polynomial(expr, *gens):
    return s.Poly(s.cancel(expr),*gens,domain=s.QQ)


def make_map(R):
    p = s.diff(R,w)
    c = -p.subs(w,1)
    kappa = s.diff(p,w).subs(w,1)/c
    check(c!=0 and kappa!=-2 and R.subs(w,0)==0 and R.subs(w,1)==0
          and p.subs(w,0)==0, "seed hypotheses")
    a = -(1+kappa)/(2+kappa)
    gamma = 1+a*v+t
    u = 1+v
    q = s.integrate(w*s.diff(p,w)/c,w)
    q -= q.subs(w,0)
    alpha = s.cancel(u+q.subs(w,u*gamma)/gamma**2)
    beta = s.cancel(c+p.subs(w,u*gamma)/gamma)
    sub = {v:x*y,t:x*x*z}
    F = [polynomial(alpha.subs(sub)/x**2,x,y,z).as_expr(),
         polynomial(beta.subs(sub)/x,x,y,z).as_expr(),
         x*(1+a*x*y+x*x*z)]
    det = s.expand(s.det(s.Matrix(F).jacobian([x,y,z])))
    check(det==c, "constant Jacobian")
    check(s.expand(F[1].subs(x,0))==-c*y/(kappa+2), "boundary B")
    check(s.diff(F[0].subs(x,0),z)==kappa+2, "boundary A coefficient")
    return F,c,a,kappa,p


def finite_audit(F,c,p1,prime):
    # Enumerate actual source triples, independently of the inverse equation.
    terms = []
    for f in F:
        row = []
        for powers,coef in s.Poly(f,x,y,z).terms():
            num,den = map(int,s.fraction(coef))
            row.append((powers,num*pow(den,-1,prime)%prime))
        terms.append(row)
    hist = {}
    for xx in range(prime):
        for yy in range(prime):
            for zz in range(prime):
                values = [sum(coef*pow(xx,i,prime)*pow(yy,j,prime)*pow(zz,k,prime)
                              for (i,j,k),coef in row)%prime for row in terms]
                if values[2]==0:
                    key = tuple(values[:2])
                    hist[key] = hist.get(key,0)+1
    def mod(r):
        num,den=map(int,s.fraction(r))
        return num*pow(den,-1,prime)%prime
    for aa in range(prime):
        for bb in range(prime):
            if mod(p1)==0:
                predicted = 1+int(bb!=0)
            else:
                delta=(bb*bb-2*mod(c*p1)*aa)%prime
                predicted=1+sum(delta*xx*xx%prime==mod(c*c)
                                for xx in range(1,prime))
            check(hist.get((aa,bb),0)==predicted,
                  f"finite field boundary: p={prime}, A={aa}, B={bb}")
    return prime**3


def main():
    # Last two seeds have repeated nonzero roots, excluded in THM-3448 §§1–5.
    seeds = [('cubic',w*w-w**3),
             ('quartic_zero_linear',w**4-w**3),
             ('quintic_nonzero_double',-w*w*(w-1)*(w-2)**2),
             ('degree_seven_two_three',w**3*(w-1)*(w-2)**3)]
    rows=[]
    for name,R in seeds:
        R=s.expand(R)
        F,c,a,kappa,p=make_map(R)
        mu=min(i[0] for i,coef in s.Poly(R,w).terms() if coef)
        _,fs=s.factor_list(R,w)
        nonzero=[multiplicity for f,multiplicity in fs if f.subs(w,0)!=0]
        expected_disc=mu+sum(e-1 for e in nonzero)
        # Full two-parameter discriminant, not a potentially special slice.
        T=R-B*C*w+c*A*C*C
        disc=s.Poly(s.discriminant(T,w),C,A,B)
        valuation=min(powers[0] for powers,coef in disc.terms() if coef)
        check(valuation==expected_disc,"raw discriminant valuation")
        cycles=sorted(([mu-1] if mu>2 else [])+[e for e in nonzero if e>1])
        defect=sum(e-1 for e in cycles)
        check(valuation-defect==2,"universal index one")
        finite_checks=sum(finite_audit(F,c,s.diff(p,w).subs(w,0),pp)
                          for pp in (5,7,11)
                          if all(s.denom(coef)%pp for f in F
                                 for coef in s.Poly(f,x,y,z).coeffs()) and c%pp)
        rows.append(dict(name=name,R=str(R),degree=s.degree(R,w),c=str(c),
                         a=str(a),kappa=str(kappa),mu=mu,cycles=cycles,
                         raw_discriminant_valuation=valuation,
                         inertia_defect=defect,index=1,
                         finite_source_triples=finite_checks))
    # The known quartic's complete chart boundary, independently expanded.
    u=1+3*x*y; gamma=1-4*x*y-x*x*z
    G=[s.cancel((2*u+u*u-3*u**4*gamma**2)/x**2),
       s.cancel((1+u-2*u**3*gamma**2)/x),x*gamma]
    check(s.expand(s.det(s.Matrix(G).jacobian([x,y,z])))==-6,"quartic det")
    check([s.expand(f.subs(x,0)) for f in G]==[87*y*y+6*z,y,0],"quartic boundary")
    check([f.subs({x:1,y:0,z:0}) for f in G]==[0,0,1],"quartic collision 1")
    check([f.subs({x:-1,y:0,z:2}) for f in G]==[0,0,1],"quartic collision 2")
    # A universal simple-root count identity, checked independently by factor multiplicities.
    tests=0
    for i in range(1,5):
        for j in range(1,5):
            for k in range(1,5):
                T=s.Poly(w**i*(w-1)**j*(w-2)**k,w)
                D=s.gcd(T,T.diff())
                formula=T.degree()-D.degree()-D.sqf_part().degree()
                check(formula==sum(e==1 for e in (i,j,k)),"simple-root count")
                check(formula==T.degree()-2*D.degree()+s.gcd(D,D.diff()).degree(),
                      'two-gcd implementation')
                tests+=1
    # MRB's finite-prefix and all-orders-asymptotic blind perturbation.
    q=s.symbols('q')
    mrb=[]
    for K in range(1,9):
        gen=s.factor(s.factorial(K)*(-q)**(K+1)/(1+q)**(K+1))
        check(s.expand(s.series(gen,q,0,K+5).removeO())==
              sum((-1)**n*s.prod(n-j for j in range(1,K+1))*q**n
                  for n in range(K+1,K+5)),"MRB shift series")
        mrb.append(str(gen))
    def partitions(n,upper=None):
        if n==0:
            yield ()
            return
        for first in range(min(n,upper or n),0,-1):
            for tail in partitions(n-first,first):
                yield (first,)+tail
    # Different exceptional-component ledgers, not raw coordinate-polynomial counts.
    partition_counts=[1]+[0]*14
    for part in range(1,15):
        for total in range(part,15):
            partition_counts[total]+=partition_counts[total-part]
    lower_bounds=[]
    degree_five=[]
    fixed_scale_checks=0
    for n in range(4,15):
        ledgers=set()
        for finite,partitions_of in ((3,n-3),(2,n-2)):
            for lam in partitions(partitions_of):
                if finite==2 and max(lam)==1:
                    continue
                mu,remaining=(2,lam) if finite==3 else (lam[0]+1,lam[1:])
                H=w**mu
                for j,e in enumerate(remaining,1):
                    H*=((w+j)/s.Integer(1+j))**e
                R=-H*(w-1)
                pp=s.diff(R,w)
                check(s.degree(R,w)==n and pp.subs(w,1)==-1,
                      'partition seed normalization')
                kappa=s.diff(pp,w).subs(w,1)
                check(kappa < -2,'partition seed admissibility')
                cpair=(finite,tuple(sorted(e for e in lam if e>1)))
                lpair=(n-2,(2,))
                ledger=tuple(sorted((cpair,lpair)))
                check(ledger not in ledgers,'inequivalence invariant distinct')
                ledgers.add(ledger)
                if n==5:
                    degree_five.append(dict(finite=finite,partition=lam,
                                            R=str(s.factor(R)),kappa=str(kappa),
                                            component_ledger=ledger))
                if finite==3:
                    r=len(lam)
                    roots=[1+s.Rational(e*(2**r-1),2**i)
                           for i,e in enumerate(lam,1)]
                    check(len(set(roots))==r and all(rr>1 for rr in roots),
                          'fixed-scale root separation')
                    RR=-w*w*(w-1)
                    for rr,e in zip(roots,lam):
                        RR*=((w-rr)/(1-rr))**e
                    check(s.diff(RR,w).subs(w,1)==-1 and
                          s.diff(RR,w,2).subs(w,1)==0,
                          'identical positive component')
                    fixed_scale_checks+=1
        lower=partition_counts[n-3]+partition_counts[n-2]-1
        check(len(ledgers)==lower,'independent partition recurrence')
        lower_bounds.append([n,lower])
    quartic_counts=[]
    for PP,QQ in [(s.Rational(1,8),-s.Rational(1,64)),(0,0),(-1,-1),(0,1)]:
        TT=s.Poly(w**3-w**4-PP*w+QQ,w)
        DD=s.gcd(TT,TT.diff())
        quartic_counts.append(4-DD.degree()-DD.sqf_part().degree())
    check(quartic_counts==[0,1,2,4],'quartic no-three-point spectrum')
    print(json.dumps(dict(seeds=rows,quartic_boundary='(87*y**2+6*z,y,0)',
                         root_multiplicity_controls=tests,mrb_series=mrb,
                         inequivalent_lower_bounds=lower_bounds,degree_five=degree_five,
                         fixed_scale_partition_checks=fixed_scale_checks,
                         second_quartic_witness_counts=quartic_counts,
                         checks='PASS'),indent=2,sort_keys=True,default=int))


if __name__=='__main__':
    main()
