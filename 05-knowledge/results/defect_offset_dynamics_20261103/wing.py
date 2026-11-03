"""A converse-closed four-module family and its exact two-coordinate wings."""
from run_polynomial import *


def wing(m,l):
    M,L=2**m,2**l
    return M*(L-1)-L//2+2,L//2-1


def hwing(m,l,k,n):
    x,b=wing(m,l);y,c=wing(n,k)
    return 2*(x*y+b*c)-1


def source4(m,l,k,n):
    sizes=(m,l,1,1,k,n)
    ad=[(v,j) for v,c in enumerate(sizes) for j in range(c)]
    T=tuple(sum(1<<j for j,(w,s) in enumerate(ad) if
                (v==w and r<s) or (v!=w and BASE[v]>>w&1)) for v,r in ad)
    return T,ad


def clock():
    states=list(product(range(1,4),repeat=4))
    def grow(s,a):return tuple((u+v-1)%3+1 for u,v in zip(s,a))
    moves=list(product(range(3),repeat=4))
    def obs(s):return tuple(z%7 for pair in (wing(s[0],s[1]),wing(s[3],s[2])) for z in pair)
    observed={obs(s) for s in states}
    signatures={tuple(hwing(*grow(s,a))%7 for a in moves) for s in states}
    require(len(observed)==49 and len(signatures)==81,'49 current wing pairs, 81 future scalar classes')
    # One wing loses its first clock exactly at l=0 mod3.
    wings={(m,l):tuple(x%7 for x in wing(m,l)) for m,l in product(range(1,4),repeat=2)}
    require(len(set(wings.values()))==7,'seven current wing observations')
    require(len({(wings[m,l],wings[m,l%3+1]) for m,l in wings})==9,'one next-l query restores nine wing phases')
    return dict(wing_states=9,wing_observations=7,joint_states=81,joint_observations=49,
                future_scalar_classes=len(signatures),scope='modulo-seven Hamilton response under four independent block growths')


def main():
    poly=compile_runs(BASE,(0,1,4,5))
    # The sixteen positive coefficients are a complete symbolic certificate.
    import sympy as sp
    M,L,K,N=sp.symbols('m2 l2 k2 n2')
    x=M*(L-1)-L/2+2;b=L/2-1;y=N*(K-1)-K/2+2;c=K/2-1
    require(sp.expand(symbolic(poly,('m','l','k','n'))-(2*(x*y+b*c)-1))==0,'all-size wing identity')
    checks=0
    for sizes in product((1,2),repeat=4):
        T,_=source4(*sizes)
        require(tables(T)[0][-1]==hwing(*sizes)==evaluate_runs(poly,sizes),'independent four-block DP')
        require(strong(T),'strong substitution');checks+=1;tables.cache_clear()
    require(induced(opposite(BASE),list(reversed(range(6))))==BASE,'literal converse core transport')
    # Prime quotient and rigidity, by every proper nontrivial subset and degree-preserving permutation.
    from branch_search import canon
    for mask in range(1,63):
        vs=[i for i in range(6) if mask>>i&1]
        if len(vs)<2:continue
        require(any(len({BASE[u]>>v&1 for v in vs})>1 for u in range(6) if u not in vs),'prime six-vertex quotient')
    from itertools import permutations
    require(sum(induced(BASE,p)==BASE for p in permutations(range(6)))==1,'rigid quotient')
    z=sp.symbols('z');dyadic=2*z/(1-2*z);ordinary=z/(1-z)
    X=dyadic**2-sp.Rational(3,2)*dyadic*ordinary+2*ordinary**2
    Y=dyadic*ordinary/2-ordinary**2
    series=sp.factor(z**2*(2*(X**2+Y**2)-ordinary**4))
    expected=z**6*(17-76*z+148*z**2-144*z**3+64*z**4)/((1-z)**4*(1-2*z)**4)
    require(sp.cancel(series-expected)==0,'marked-class generating function')
    out=dict(coefficients=[[list(p),w] for p,w in sorted(poly.items())],literal_checks=checks,clock=clock(),
             marked_series=str(series),marked_coefficients=[int(sp.expand(sp.series(series,z,0,13).removeO()).coeff(z,j)) for j in range(6,13)])
    (HERE/'wing.json').write_text(json.dumps(out,indent=2)+'\n')
    print('wing identity and class parameters proved; checks',checks,'clock',out['clock'],flush=True)


if __name__=='__main__':main()
