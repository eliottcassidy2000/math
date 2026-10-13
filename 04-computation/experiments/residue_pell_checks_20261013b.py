#!/usr/bin/env python3
"""Exact algebra/finite universes, plus explicitly numerical spectra.

Dependencies: numpy and sympy. No chromatic-number or Collatz-convergence
claim is inferred from floating-point output. The general proofs are in
tree_resonance_composite_pell_20261013b.md.
"""
import argparse
import json
import math
import numpy as np
import sympy as s

CHECKS = 0


def require(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(label)


def graph_spectrum(modulus, nonsquare=-1):
    v = np.arange(modulus, dtype=np.int64)
    points = ((v[:, None]**2-nonsquare*v[None, :]**2) % modulus == 1)
    z = np.fft.fft2(points)
    require(float(np.max(abs(z.imag))) < 1e-8, 'real spectrum')
    lam = np.sort(z.real.ravel())
    return int(points.sum()), lam


def local_checks():
    rows = []
    specs = {}
    for p, maxr in [(3,4),(7,3),(11,2)]:
        for r in range(1,maxr+1):
            degree, spectrum = graph_spectrum(p**r)
            require(degree == (p+1)*p**(r-1), 'norm-one count')
            require(spectrum[0] >= -degree*2*math.sqrt(p)/(p+1)-1e-8,
                    'all-character Weil/stationary bound')
            specs[p**r] = (degree,spectrum)
            # Audit conductor-specific bounds on every character.
            M=p**r
            a=np.arange(M,dtype=np.int64)
            mask=((a[:,None]**2+a[None,:]**2)%M == 1)
            eig=np.fft.fft2(mask).real
            for i in range(M):
                for j in range(M):
                    if i==j==0:
                        continue
                    u=math.gcd(math.gcd(i,j),M)
                    ss=r
                    while u>1:
                        u//=p;ss-=1
                    bound=p**(r-ss)*2*p**(ss/2)
                    require(abs(eig[i,j]) <= bound+1e-8,'conductor bound')
            rows.append({'p':p,'r':r,'vertices':M*M,'degree':degree,
                         'least_eigenvalue_numerical':float(spectrum[0]),
                         'hoffman_numerical':1-degree/float(spectrum[0])})
    composite=[]
    for m,n in [(3,7),(9,7),(9,49),(27,11)]:
        d,sp=graph_spectrum(m*n)
        dm,sm=specs[m];dn,sn=specs[n]
        target=np.sort(np.multiply.outer(sm,sn).ravel())
        require(d==dm*dn,'CRT degree')
        require(np.max(abs(sp-target))<1e-7,'complete CRT spectrum')
        require(abs(sp[0]/d-min(sm[0]/dm,sn[0]/dn))<1e-9,
                'tensor least normalized eigenvalue')
        composite.append({'modulus':m*n,'degree':d,'vertices':(m*n)**2,
                          'hoffman_numerical':1-d/float(sp[0])})
    # The relevant first fixed prime: u^2=2 is unramified quadratic at 109.
    p=109;M=p*p;D=2
    base=[(a,b) for a in range(p) for b in range(p) if (a*a-D*b*b-1)%p==0]
    points=[]
    for a,b in base:
        carry=(a*a-D*b*b-1)//p
        if b:
            for i in range(p):
                j=((carry+2*a*i)*pow(2*D*b,-1,p))%p
                points.append((a+p*i,b+p*j))
        else:
            for j in range(p):
                i=((-carry+2*D*b*j)*pow(2*a,-1,p))%p
                points.append((a+p*i,b+p*j))
    require(len(set(points)) == p*(p+1),'109^2 torus count')
    for a,b in points:
        require((a*a-D*b*b-1)%M==0,'109^2 exact norm')
    points=np.array(points,dtype=np.int64)
    roots={}
    for z in range(M):
        roots.setdefault(z*z%M,[]).append(z)
    rng=np.random.default_rng(20261013)
    maxerr=0.
    for i,j in [(1,0),(0,1),((M-1)//2,0)]+[tuple(map(int,rng.integers(0,M,2))) for _ in range(128)]:
        if i%p==j%p==0:
            continue
        direct=np.exp(2j*np.pi*((i*points[:,0]+j*points[:,1])%M)/M).sum()
        norm=(i*i-j*j*pow(D,-1,M))%M
        predicted=p*sum(np.exp(2j*np.pi*z/M) for z in roots.get(norm,[]))
        maxerr=max(maxerr,float(abs(direct-predicted)))
        require(abs(direct-predicted)<1e-8,'primitive stationary-phase formula')
    db,sb=graph_spectrum(p,D)
    least=max(-float(sb[0]),2.)
    p109={'p':p,'r':2,'vertices':M*M,'degree':p*(p+1),
          'exact_norm_points':len(points),'sample_stationary_max_error':maxerr,
          'prime_least_eigenvalue_numerical':float(sb[0]),
          'prime_power_least_eigenvalue_numerical':p*float(sb[0]),
          'hoffman_numerical':1+(p+1)/least,
          'weil_hoffman_proved_floor':1+(p+1)/(2*math.sqrt(p))}
    return rows,composite,p109


def pell_checks():
    z,L=s.symbols('z L');t=2*z-1
    rows=[]
    for n in range(1,21):
        a=s.Rational(1,2*(n+1));b=s.Rational(1,2*n)
        P=s.Poly(a*s.chebyshevt(n+1,t)-b*s.chebyshevt(n,t),z)
        S=s.Poly(a*s.chebyshevu(n,t)-b*s.chebyshevu(n-1,t),z)
        c=P.eval(0)
        chi=P-s.Poly(c,z)
        G=s.Poly(s.chebyshevt(n,t)+(t-1)*s.chebyshevu(n-1,t),z)
        require((P.diff()-G).is_zero,'Chebyshev primitive')
        rem=P*P-s.Poly(t*t-1,z)*S*S
        require((rem-s.Poly(a*a+b*b-2*a*b*t,z)).is_zero,'two-term Pell norm')
        delta=s.Poly((L-chi.as_expr())**2-(t*t-1)*S.as_expr()**2,z)
        formula=s.Poly((L+c)**2-2*(L+c)*P.as_expr()+a*a+b*b-2*a*b*t,z)
        require((delta-formula).is_zero,'generic fibre exact formula')
        require(delta.degree()==n+1,'exact degree')
        rows.append({'depth':n,'generic_degree':delta.degree(),'genus':n//2})
    # Universal binary-refinement genus formula for the marked fish cover.
    branch=[{0,1},{0,4},{1,-3},{s.Rational(-1,2),s.Rational(3,2)}]
    genera=[]
    for rank in [3,4]:
        union=set().union(*branch[:rank]);d=2**rank
        g=1+(-2*d+len(union)*d//2)//2
        quotient=[]
        for mask in range(1,2**rank):
            support=set()
            for i in range(rank):
                if mask>>i&1:support^=branch[i]
            require(bool(support),'independent square classes')
            quotient.append((len(support)-2)//2)
        require(sum(quotient)==g,'character genus decomposition')
        genera.append({'rank':rank,'genus':g,'character_quotient_genera':quotient})
    return rows,genera


def golden_checks():
    phi=(1+s.sqrt(5))/2
    bar=(1-s.sqrt(5))/2
    vertices=[(a,b) for a in range(-3,4) for b in range(-3,4)
              if abs(a+b*phi)<=7 and abs(a+b*bar)<=1]
    def G(v):
        a,b=v
        return (a//2,b//2) if a%2==b%2==0 else (
            (b+b%2)//2,(a+b+(a+b)%2)//2)
    roots={(0,0),(1,1),(1,2),(2,3)}
    edges=[]
    maximum=0
    for v in vertices:
        w=G(v)
        require(w in vertices,'invariant golden box')
        path=[v]
        while path[-1] not in roots:
            nxt=G(path[-1]);require(nxt not in path,'no other golden cycle')
            path.append(nxt)
        maximum=max(maximum,len(path)-1)
        edges.append({'source':v,'target':w,'remaining_steps':len(path)-1})
    for v in roots:require(G(v)==v,'golden fixed point')
    require(bool((8*phi+1)/14<1),'outside-box strict contraction')
    return {'exact_vertices':len(vertices),'edges':edges,
            'fixed_points':sorted(roots),'maximum_steps_inside_box':maximum}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',required=True);args=ap.parse_args()
    local,composite,p109=local_checks()
    pell,fish=pell_checks()
    golden=golden_checks()
    data={'checks':CHECKS,'local_spectra_NUMERICAL':local,
          'composite_spectra_NUMERICAL':composite,'fixed_prime_power':p109,
          'pell_exact':pell,'fish_exact':fish,'golden_exact':golden}
    with open(args.output,'w') as f:json.dump(data,f,indent=2);f.write('\n')
    print(json.dumps(data,indent=2))


if __name__=='__main__':main()
