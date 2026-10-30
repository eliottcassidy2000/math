"""Independent Session114 support/channel and cyclic-parent audit.

Python 3 + SymPy. No producer code is imported. Exception gates survive -O.
"""
import itertools
import json
import math
import sympy as sp


def require(test, message):
    if not test:
        raise RuntimeError(message)


def phi(n):
    return sum(math.gcd(a, n) == 1 for a in range(n))


def weight_support(w, m, degree, target):
    return {(i, j) for i in range(degree + 1)
            for j in range(degree + 1 - i)
            if (i - w*j - target) % m == 0}


def jacobian_audit():
    cases = 0
    for p,q in itertools.product(range(1,8),repeat=2):
        if math.gcd(p,q) != 1:
            continue
        for d in range(1,16):
            for target,cutoff in ((-q,max(q*(d-1),p*d+q)),
                                  (p,max(q*d+p,p*(d-1)))):
                m = cutoff+1
                require(all(p*i-q*j == target for i in range(d+1)
                            for j in range(d+1-i)
                            if (p*i-q*j-target)%m == 0),
                        "Mixed-weight finite certificate")
                cases += 1
    for w in range(1, 13):
        for d in range(1, 31):
            cutoff = max(w*(d-1), d+w)
            for m in (cutoff+1, cutoff+2, 2*cutoff+1):
                require(all(i-w*j == -w for i, j in
                            weight_support(w, m, d, -w)), "Q alias")
                cases += 1
    for w in range(1, 11):
        for d in range(2, 31):
            m = w*d+1
            for target in (1, -w):
                support = weight_support(w, m, d, target)
                aliases = {(i, j) for i, j in support if i-w*j != target}
                expected = ({(0, d)} if target == 1 else
                            ({(d, 0)} if w == 1 else set()))
                require(aliases == expected, "Boundary support mismatch")
                cases += 1
    for m in range(3, 65):
        for target in (1, -1):
            aliases = {(i, j) for i, j in weight_support(1, m, m, target)
                       if i-j != target}
            require(aliases == ({(0, m-1)} if target == 1 else {(m-1, 0)}),
                    "Balanced degree-m support mismatch")
            cases += 1
    x, y, t = sp.symbols('x y t')
    for p,q in itertools.product(range(1,6),repeat=2):
        if math.gcd(p,q) != 1:
            continue
        A, B = 2-3*t+5*t**2, 3+7*t-2*t**2+t**3
        P,Q = x*A.subs(t,x**q*y**p),y*B.subs(t,x**q*y**p)
        J = sp.Poly(sp.diff(P,x)*sp.diff(Q,y)-sp.diff(P,y)*sp.diff(Q,x),x,y)
        require(J.coeff_monomial(x**(5*q)*y**(5*p)) == 5*(1+2*q+3*p),
                "Mixed-weight top coefficient")
        extra = sum((i+j+1)*x**i*y**j for i in range(6)
                    for j in range(6-i) if p*i-q*j != p)
        JJ = sp.Poly(sp.diff(P+extra,x)*sp.diff(Q,y)-
                     sp.diff(P+extra,y)*sp.diff(Q,x),x,y)
        zero = sum(c*x**i*y**j for (i,j),c in JJ.terms() if p*i-q*j == 0)
        require(sp.expand(zero-J.as_expr()) == 0,"Mixed-weight projection")
        cases += 2
    for w in range(1, 8):
        for da, db in itertools.product(range(4), repeat=2):
            # All coefficients are nonzero, with varying signs; top degree
            # extraction is checked directly in x,y, not by the claimed ODE.
            A = sum((-1)**i*(i+1)*t**i for i in range(da+1))
            B = sum((i+2)*t**i for i in range(db+1))
            P, Q = x*A.subs(t, x**w*y), y*B.subs(t, x**w*y)
            J = sp.Poly(sp.diff(P,x)*sp.diff(Q,y)-sp.diff(P,y)*sp.diff(Q,x),x,y)
            top = J.coeff_monomial(x**(w*(da+db))*y**(da+db))
            require(top == (1+w*da+db)*(-1)**da*(da+1)*(db+2),
                    "Top coefficient mismatch")
            # Add an arbitrary collection of other integer-weight components;
            # direct monomial filtering must recover the same weight-zero J.
            extra = sum((i+j+1)*x**i*y**j for i in range(5)
                        for j in range(5-i) if i-w*j != 1)
            JJ = sp.Poly(sp.diff(P+extra,x)*sp.diff(Q,y)-
                         sp.diff(P+extra,y)*sp.diff(Q,x),x,y)
            zero = sum(c*x**i*y**j for (i,j),c in JJ.terms() if i-w*j == 0)
            require(sp.expand(zero-J.as_expr()) == 0, "Projection identity")
            cases += 2
    for d in range(2, 10):
        alpha, beta = sp.symbols('alpha beta')
        A, B = 1+2*t+3*t**2, 2-3*t+5*t**3
        P = x*A.subs(t,x*y)+alpha*y**d
        Q = y*B.subs(t,x*y)+beta*x**d
        J = sp.Poly(sp.diff(P,x)*sp.diff(Q,y)-sp.diff(P,y)*sp.diff(Q,x),x,y)
        components = {}
        for (i,j),c in J.terms():
            components[i-j] = components.get(i-j,0)+c*x**i*y**j
        require(sp.expand(components.get(d+1,0)+
                          beta*d*x**(d+1)*sp.diff(A,t).subs(t,x*y)) == 0,
                "Positive alias channel")
        require(sp.expand(components.get(-d-1,0)+
                          alpha*d*y**(d+1)*sp.diff(B,t).subs(t,x*y)) == 0,
                "Negative alias channel")
        cases += 2
    # Important excluded boundary: degree1, m=w+1 has a scalar source action.
    for w in range(1, 10):
        for target in (1,-w):
            require(weight_support(w,w+1,1,target) == {(1,0),(0,1)},
                    "Degree-one scalar action")
            cases += 1
    return cases


def compose(p, q):
    return tuple(p[q[i]] for i in range(len(p)))


def powers(p):
    n = len(p)
    out = [tuple(range(n))]
    for k in range(1,n):
        out.append(compose(p,out[-1]))
    return out


def is_full_cycle(p):
    seen, at = set(), 0
    while at not in seen:
        seen.add(at)
        at = p[at]
    return at == 0 and len(seen) == len(p)


def parent_audit(r,m):
    n = r*m
    count = 0
    parents = {}
    parent_slopes = {}
    descriptions = {}
    g = math.gcd(r,m)
    for shift in range(r):
        for z in itertools.product(range(m),repeat=r):
            p = tuple((i+shift)%r+r*((j+z[i])%m)
                      for j in range(m) for i in range(r))
            descriptions[p] = (shift,sum(z)%m)
            full = is_full_cycle(p)
            require(full == (math.gcd(shift,r)==1 and math.gcd(sum(z),m)==1),
                    "Cycle criterion")
            if full:
                count += 1
                pp = powers(p)
                key = min(pp[k] for k in range(n) if math.gcd(k,n)==1)
                parents.setdefault(key,0)
                parents[key] += 1
                slope = (sum(z)*pow(shift,-1,g))%g if g>1 else 0
                parent_slopes.setdefault(key,set()).add(slope)
    expected = phi(r)*phi(m)*m**(r-1)//phi(n)
    require(count == phi(r)*phi(m)*m**(r-1), "Generator count")
    require(len(parents) == expected, "Parent count")
    require(set(parents.values()) == {phi(n)}, "Generator multiplicity")
    require(all(len(v)==1 for v in parent_slopes.values()),"Parent slope changes")
    slopes = sorted({next(iter(v)) for v in parent_slopes.values()})
    require(len(slopes)==phi(g),"Number of slopes")
    for rho in slopes:
        require(sum(next(iter(v))==rho for v in parent_slopes.values()) ==
                m**(r-1)//g,"Slope roster size")
        nu = next(k for k in range(1,m) if math.gcd(k,m)==1 and k%g==rho)
        h = tuple(i+r*((nu*j)%m) for j in range(m) for i in range(r))
        hi = tuple(h.index(k) for k in range(n))
        a = tuple((k+1)%n for k in range(n))
        v = tuple((k+r)%n if k%r==0 else k for k in range(n))
        vi = tuple(v.index(k) for k in range(n))
        b = compose(v,compose(a,vi))
        gens = [compose(h,compose(p,hi)) for p in (a,b)]
        identity = tuple(range(n))
        generated,stack = {identity},[identity]
        while stack:
            cur = stack.pop()
            for p in gens:
                new = compose(p,cur)
                if new not in generated:
                    generated.add(new)
                    stack.append(new)
        intended = {p for p,(tt,ss) in descriptions.items() if (ss-rho*tt)%g==0}
        require(generated == intended,"Slope-group generation")
    # Directly enumerate the common W-invariant tournament bank and verify
    # every parent acts on every member; no DSU/compiler from producers.
    natural = tuple(range(n))
    alt = lambda q: tuple((2*(k//2 if k%2==0 else q-1-k//2))%q for k in range(q))
    selector = tuple(i+r*j for i in alt(r) for j in alt(m))
    require(selector[0]==0 and selector[-1]==n-1 and len(set(selector))==n,
            "Selector ports")
    counts = [0,0]
    for bits in itertools.product((0,1),repeat=(r+m-2)//2):
        def edge(x,y):
            i,j = x%r,x//r
            k,l = y%r,y//r
            q, delta, off = (m,(l-j)%m,(r-1)//2) if i==k else (r,(k-i)%r,0)
            return bool(bits[off+min(delta,q-delta)-1] ^ (delta>q//2))
        for p in parents:
            require(all(edge(x,y)==edge(p[x],p[y]) for x in range(n)
                        for y in range(x+1,n)), "Parent not invariant")
        for k,path in enumerate((natural,selector)):
            if all(edge(x,y) for x,y in zip(path,path[1:])):
                counts[k] += 1
                # For a fixed path, native complementation is injective and
                # each of these distinct named parents owns the same output.
                output = tuple(edge(x,y) ^ ((x,y) not in set(zip(path,path[1:]))
                    and (y,x) not in set(zip(path,path[1:])))
                    for x in range(n) for y in range(x+1,n))
                require(len(output)==n*(n-1)//2,"Output size")
    require(counts == [2**((r+m-4)//2),1],"Common path banks")
    return {'r':r,'m':m,'full_cycles':count,'parents':expected,
            'phase_slopes':slopes,'parents_per_slope':m**(r-1)//g,
            'natural_common_sources':counts[0],'selector_common_sources':counts[1]}


if __name__ == '__main__':
    print(json.dumps({'jacobian_exact_gates':jacobian_audit(),
                      'parent_cases':[parent_audit(r,m) for r,m in
                                      ((3,3),(3,5),(5,3),(3,7))]},sort_keys=True,indent=2))
