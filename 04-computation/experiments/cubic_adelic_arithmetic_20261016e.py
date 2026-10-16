#!/usr/bin/env python3
"""Exact local arithmetic and an elliptic torsion carrier for the cubic field.
PARI is an independent check of field discriminant and elliptic conductor.
No elliptic-curve database, modular-form table, or conjectural rank is used.
"""
import json
import subprocess
import sympy as s
import cubic_adelic_20261016e as A
import cubic_crossedge_packets_20261016d as C

checks=0
def need(ok,label):
    global checks
    checks+=1
    if not ok:raise RuntimeError(label)

def run():
    x,j=s.symbols('x j');f=x**3+2*x*x-4*x-4;g=x**3-x*x-3*x+1
    need(all(int(f.subs(x,a))%3 for a in range(3)),'irreducibility of cubic modulo three')
    need(s.discriminant(f,x)==592 and s.discriminant(g,x)==148,'two integral polynomial discriminants')
    r=x*x/2+x-1
    need(s.rem(g.subs(x,r),f,x)==0,'maximal-order generator polynomial')
    need(s.rem(r*r-3-x,f,x)==0,'inverse integral-generator transformation')
    need(s.rem((x*x/2)**3-6*(x*x/2)**2+8*(x*x/2)-2,f,x)==0,
         'half-square is integral')
    correction=s.expand((g-(x+1)**3)/2)
    need(s.gcd(s.Poly(x+1,x,modulus=2),s.Poly(correction,x,modulus=2)).degree()==0,
         'Dedekind index criterion at two')
    need(s.gcd(148,5**30*11**36)==1 and s.gcd(148,5**6*11**4)==1,
         'coprime cubic and both old-field discriminants')
    for sign in [-1,1]:
        need(s.rem(s.expand(f.subs(x,-3+sign*j)),j*j+11,j)==76,'both local denominator evaluations')
    u=A.new_direction();one=C.I;w=C.W
    a=C.sc(C.add(C.sc(one,20),C.sc(C.lift(C.G.b),11)),C.F(1,19))
    need(C.add(C.add(C.power(u,3),C.mul(a,C.power(u,2))),
          C.add(C.sc(C.mul(C.conj(a),u),-1),C.sc(one,-1)))==C.O,
         'self-inversive cubic polynomial of the new unit')
    need(C.sc(C.add(C.add(u,C.conj(u)),C.sc(one,-2)),C.F(1,11))==
         C.sc(C.sub(w,C.sc(one,4)),C.F(1,19)),'mixed-prime displacement identity')

    # Exact metric in the basis 1,u,u^2, over the original field L.
    G=C.G;aa=G.sc(G.add(G.sc(G.I,20),G.sc(G.b,11)),C.F(1,19))
    def pmul(p,q):
        out=[G.O for _ in range(5)]
        for i in range(3):
            for h in range(3):out[i+h]=G.add(out[i+h],G.mul(p[i],q[h]))
        for d in [4,3]:
            for h,t in enumerate([G.I,G.conj(aa),G.sc(aa,-1)]):
                out[d-3+h]=G.add(out[d-3+h],G.mul(out[d],t))
        return out[:3]
    W=[G.sc(G.add(G.sc(G.I,-14),G.sc(G.b,11)),C.F(1,11)),
       G.sc(G.add(G.sc(G.I,39),G.sc(G.b,11)),C.F(1,11)),G.sc(G.I,C.F(19,11))]
    Wpowers=[[G.I,G.O,G.O],W,pmul(W,W)]
    def from_coords(cs):return A.sum_vectors([C.mul(C.lift(c),C.power(u,i)) for i,c in enumerate(cs)])
    need(from_coords(W)==w,'explicit cubic change of basis')
    def coords(z):
        out=[]
        for h in range(3):
            c=G.O
            for i in range(3):c=G.add(c,G.mul(z[40*i:40*i+40],Wpowers[i][h]))
            out.append(c)
        return out
    def metric(cs):
        p,q,r=cs;d1=G.add(G.mul(q,G.conj(p)),G.mul(r,G.conj(q)));d2=G.mul(r,G.conj(p))
        cross=G.add(G.add(d1,G.sc(G.mul(aa,d2),-1)),G.conj(d2))
        scalar=G.add(G.add(G.mul(p,G.conj(p)),G.mul(q,G.conj(q))),G.mul(r,G.conj(r)))
        scalar=G.add(scalar,G.add(G.mul(G.conj(aa),d2),G.mul(aa,G.conj(d2))))
        return cross==G.O and scalar==G.I
    for p in range(-2,3):
        for q in range(-2,3):
            for r in range(-2,3):
                cs=[G.sc(G.I,x) for x in [p,q,r]]
                need(metric(cs)==C.unit(from_coords(cs)),'all 125 rational coefficient triples: metric versus independent tensor norm')
    v=C.generator();vc=coords(v)
    need(from_coords(vc)==v and metric(vc) and all(any(c) for c in vc),
         'prior generator is a genuinely three-coordinate unit')
    need(G.mul(aa,G.conj(aa))==G.sc(G.I,C.F(1731,361)), 'middle-zero support obstruction')

    # F_(19^2), and its unramified modulus-19^2 lift, represented by j^2=-11.
    def add(a,b,m):return ((a[0]+b[0])%m,(a[1]+b[1])%m)
    def neg(a,m):return (-a[0]%m,-a[1]%m)
    def mul(a,b,m):return ((a[0]*b[0]-11*a[1]*b[1])%m,(a[0]*b[1]+a[1]*b[0])%m)
    def scale(a,c,m):return (a[0]*c%m,a[1]*c%m)
    def inv(a,m):return scale((a[0],-a[1]),pow((a[0]**2+11*a[1]**2)%m,-1,m),m)
    def evalf(a,m):return add(add(mul(mul(a,a,m),a,m),scale(mul(a,a,m),2,m),m),add(scale(a,-4,m),(-4,0),m),m)
    roots=[]
    for r0 in [(4,0),(-3,-1),(-3,1)]:
        need(evalf(r0,19)==(0,0),'all three residue roots')
        m=19**2;der=add(add(scale(mul(r0,r0,m),3,m),scale(r0,4,m),m),(-4,0),m)
        r1=add(r0,neg(mul(evalf(r0,m),inv(der,m),m),m),m)
        need(evalf(r1,m)==(0,0),'one Hensel step independently checked')
        num=add(r1,(3,-1),m);den=add(r1,(3,1),m)
        vn=1 if all(c%19==0 for c in num) else 0
        vd=1 if all(c%19==0 for c in den) else 0
        need((not vn or num!=(0,0)) and (not vd or den!=(0,0)),
             'zero valuations are exactly one, not merely at least one')
        roots.append(dict(residue=[c%19 for c in r0],hensel=list(r1),unit_valuation=vn-vd))
    need([r['unit_valuation'] for r in roots]==[0,-1,1],'fixed branch and conjugate pole-zero branches')
    primes=list(s.primerange(3,44));traces=[]
    for p in primes:
        if p==37:continue
        count=1+sum(1 if (v:=int(f.subs(x,a))%p)==0 else 2 if pow(v,(p-1)//2,p)==1 else 0 for a in range(p))
        ap=p+1-count;traces.append([int(p),ap])
    need(dict(traces)[13]==0 and pow((-4*13)%5,2,5)==4,'irreducible Frobenius polynomial on E[5] at thirteen')
    b2,b4,b6=8,-8,-16;b8=-48
    c4=b2*b2-24*b4;c6=-b2**3+36*b2*b4-216*b6
    disc=-(b2*b2*b8+8*b4**3+27*b6**2-9*b2*b4*b6)
    need((c4,c6,disc)==(256,640,9472),'direct Weierstrass invariants')
    gp='f=x^3+2*x^2-4*x-4; E=ellinit([0,2,0,-4,-4]); print([nfdisc(f),ellglobalred(E)[1],E.c4,E.c6,E.disc]); print(vector(12,i,ellap(E,[3,5,7,11,13,17,19,23,29,31,41,43][i])));\n'
    result=subprocess.run(['gp','-fq'],input=gp,text=True,capture_output=True,check=True)
    need(not result.stderr.strip(),'PARI emitted no errors')
    lines=result.stdout.strip().splitlines();pari=json.loads(lines[0]);pa=json.loads(lines[1])
    need(pari==[148,148,256,640,9472],'independent PARI arithmetic')
    need(pa==[a for p,a in traces],'independent exact Frobenius point counts')
    return dict(status='PROVED algebra identities; FINITE-EXACT independent conductor computation; modularity is CITED',
      cubic=str(f),integer_generator=str(g),field_discriminant=148,polynomial_discriminant=592,
      local_nineteen=roots,elliptic_coefficients=[0,2,0,-4,-4],elliptic_conductor=148,
      elliptic_invariants=dict(c4=c4,c6=c6,discriminant=disc,j='65536/37'),frobenius_traces=traces,
      mod_five_irreducibility_prime=13,pari_commands=gp,pari_result=pari,checks=checks,
      old_generator_in_unit_basis=[[str(x) for x in c] for c in vc],metric_rational_controls=125,
      unit_checks=A.checks,engine_checks=C.CHECKS)

if __name__=='__main__':print(json.dumps(run(),separators=(',',':'),sort_keys=True))
