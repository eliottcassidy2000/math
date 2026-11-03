"""Exact exploratory path-bin calculus; no third-party dependencies."""
from functools import lru_cache
from itertools import permutations, product
from math import comb


def require(ok, detail):
    if not ok:
        raise RuntimeError(detail)


def opposite(T):
    n = len(T)
    return tuple(((1 << n) - 1) ^ (1 << i) ^ row for i, row in enumerate(T))


def induced(T, vertices):
    return tuple(sum(((T[u] >> v) & 1) << j for j, v in enumerate(vertices)) for u in vertices)


def transitive(n):
    return tuple(((1 << n) - 1) ^ ((1 << (i + 1)) - 1) for i in range(n))


def tiling(n, bits):
    T = list(transitive(n))
    k = 0
    for i in range(n):
        for j in range(i + 2, n):
            if bits >> k & 1:
                T[i] ^= 1 << j
                T[j] ^= 1 << i
            k += 1
    return tuple(T)


@lru_cache(None)
def tables(T):
    n = len(T)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if mask >> v & 1:
                rest = mask ^ (1 << v)
                if rest:
                    dp[mask][v] = sum(dp[rest][u] for u in range(n) if rest >> u & 1 and T[u] >> v & 1)
    H = tuple(sum(row) for row in dp)
    return (1,) + H[1:], tuple(tuple(row) for row in dp)


def convolution(a, b):
    out = []
    for mask in range(len(a)):
        v, sub = 0, mask
        while True:
            v += a[sub] * b[mask ^ sub]
            if not sub:
                break
            sub = (sub - 1) & mask
        out.append(v)
    return tuple(out)


def bin_values(T, root, maximum):
    H, ends = tables(T)
    _, starts = tables(opposite(T))
    E = tuple(row[root] for row in ends)
    S = tuple(row[root] for row in starts)
    current = (1,) + (0,) * (len(H) - 1)
    p, e, s = [0], [0], [0]
    for m in range(1, maximum + 1):
        e.append(convolution(E, current)[-1])
        s.append(convolution(S, current)[-1])
        current = convolution(H, current)
        p.append(current[-1])
    return p, e, s


def seam(X, a, Y, b):
    """Left X, right Y^op, central c; returns literal rows and c."""
    nx, ny = len(X), len(Y)
    c = nx + ny
    Z = list(X) + [row << nx for row in opposite(Y)] + [0]
    for y in range(ny):
        Z[nx + y] |= (1 << nx) - 1
    for x in range(nx):
        if x == a:
            Z[x] |= 1 << c
        else:
            Z[c] |= 1 << x
    for y in range(ny):
        if y == b:
            Z[c] |= 1 << (nx + y)
        else:
            Z[nx + y] |= 1 << c
    return tuple(Z), c


def out_birth(T, root):
    n = len(T)
    return tuple(row | (0 if i == root else 1 << n) for i, row in enumerate(T)) + (1 << root,), n


def strong(T):
    n = len(T)
    for start in range(n):
        reached, old = 1 << start, 0
        while old != reached:
            old = reached
            for v in range(n):
                if reached >> v & 1:
                    reached |= T[v]
        if reached != (1 << n) - 1:
            return False
    return True


def polynomial(values):
    """Integer Newton coefficients, f(m)=sum c_j binom(m,j)."""
    row, answer = list(values), []
    while row:
        answer.append(row[0])
        row = [b-a for a, b in zip(row, row[1:])]
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return tuple(answer)


def evaluate(coeff, m):
    if m < 0:
        return sum(c * (-1)**j * comb(j-m-1, j) for j, c in enumerate(coeff))
    return sum(c * comb(m, j) for j, c in enumerate(coeff) if j <= m)


def as_polynomials(values):
    p,e,s=values
    n=len(p)-1
    def rooted(v):
        # Rooted polynomials have degree <= n-1 and are combinatorial at m>=1.
        # Their formal value at zero need not be zero.
        zero=sum((-1)**(j+1)*comb(n,j)*v[j] for j in range(1,n+1))
        return polynomial([zero]+v[1:])
    return polynomial(p),rooted(e),rooted(s)


def seam_poly(left, right):
    n = len(left[0]) + len(right[0]) - 1
    p, e, s = [], [], []
    for m in range(n+1):
        px, ex, sx = [evaluate(poly, m) for poly in left]
        py, ey, sy = [evaluate(poly, m) for poly in right]
        ax, ay = evaluate(left[1], m+1), evaluate(right[1], m+1)
        p.append(m * (ax*ay + ax*py + px*ay + (px-sx)*(py-sy)))
        e.append(ex*py + evaluate(left[0], m-1)*(py-sy))
        s.append(ey*px + evaluate(right[0], m-1)*(px-sx))
    return tuple(polynomial(v) for v in (p,e,s))


def birth_poly(source):
    n = len(source[0])
    p, e, s = [], [], []
    for m in range(n+1):
        pm, em, sm = [evaluate(poly, m) for poly in source]
        p.append(m*(evaluate(source[2],m+1)+pm-em))
        e.append(pm-em)
        s.append(sm+evaluate(source[0],m-1))
    return tuple(polynomial(v) for v in (p,e,s))


def tail_coefficients(word, m):
    """1 means tail vertex -> new root. Six-state weighted word automaton."""
    t=len(word)
    L=(m-1)**t+sum(m**i*(m-1)**(t-i-1) for i,b in enumerate(word) if b)
    D=(m-1)**t+sum((m-1)**i*m**(t-i-1) for i,b in enumerate(word) if not b)
    states={(0,0):1}
    for b in word:
        nxt={}
        def add(key,value):
            nxt[key]=nxt.get(key,0)+value
        for (last,started),weight in states.items():
            add((last,started),(m-1)*weight)
            add((1 if b else 2,started),weight)
            if started or not b:
                add((last,1),weight)
        states=nxt
    B=sum(v for (last,started),v in states.items() if last!=2)
    return L,D,B


def funnel(T, ports, word):
    n,t=len(T),len(word)
    root=n+t
    tailmask=((1 << t)-1) << n
    rows=[row|tailmask|(0 if i in ports else 1 << root) for i,row in enumerate(T)]
    rows.extend((transitive(t)[i] << n) | ((1 << root) if b else 0) for i,b in enumerate(word))
    rows.append(sum(1 << i for i in ports)+sum(1 << (n+i) for i,b in enumerate(word) if not b))
    return tuple(rows),root


def funnel_poly(source, word):
    """One core port, arbitrary binary tail; output rooted at new vertex."""
    n=len(source[0])+len(word)
    p,e,s=[],[],[]
    for m in range(n+1):
        pm,em,sm=[evaluate(poly,m) for poly in source]
        L,D,B=tail_coefficients(word,m)
        Lplus=tail_coefficients(word,m+1)[0]
        p.append(m*(Lplus*evaluate(source[2],m+1)+B*pm-D*em))
        e.append(L*pm-(m-1)**len(word)*em)
        s.append(m**len(word)*sm+D*evaluate(source[0],m-1))
    return tuple(polynomial(v) for v in (p,e,s))


def joint_values(T, ports, m):
    """C_J(m+1): first bin ends in J, second starts in J, m-1 other bins."""
    H,ends=tables(T)
    _,starts=tables(opposite(T))
    E=tuple(sum(row[r] for r in ports) for row in ends)
    S=tuple(sum(row[r] for r in ports) for row in starts)
    current=convolution(E,S)
    for _ in range(m-1):
        current=convolution(current,H)
    return current[-1]


def components(T):
    n=len(T)
    reach=[]
    for v in range(n):
        mask,old=1 << v,0
        while old!=mask:
            old=mask
            for u in range(n):
                if mask >> u & 1:
                    mask|=T[u]
        reach.append(mask)
    unseen=set(range(n)); result=[]
    while unseen:
        v=min(unseen)
        block=tuple(u for u in sorted(unseen) if reach[v] >> u & 1 and reach[u] >> v & 1)
        result.append(block);unseen.difference_update(block)
    blocks=tuple(result)
    result.sort(key=lambda block: -sum(1 for other in blocks if other!=block and T[block[0]] >> other[0] & 1))
    return result


def main():
    roots = [(tiling(n,bits),r) for n in range(1,4) for bits in range(1 << ((n-1)*(n-2)//2)) for r in range(n)]
    cases = 0
    for X,a in roots:
        for Y,b in roots:
            Z,c = seam(X,a,Y,b)
            x = bin_values(X,a,4)
            y = bin_values(Y,b,4)
            z = bin_values(Z,c,3)
            for m in range(1,4):
                px,ex,sx = [v[m] for v in x]
                py,ey,sy = [v[m] for v in y]
                expected = (m*(x[1][m+1]*y[1][m+1]+x[1][m+1]*py+px*y[1][m+1]+(px-sx)*(py-sy)), ex*py+x[0][m-1]*(py-sy), ey*px+y[0][m-1]*(px-sx))
                require(tuple(v[m] for v in z)==expected, (X,a,Y,b,m,z,expected))
            if x[1][1] and y[1][1]:
                require(strong(Z), ('strong seam',X,a,Y,b))
            ordinal=tuple(row|(((1 << len(Y))-1) << len(X)) for row in X)+tuple(row << len(X) for row in Y)
            ox=bin_values(ordinal,a,3);oy=bin_values(ordinal,len(X)+b,3)
            for m in range(1,4):
                require(tuple(v[m] for v in ox)==(x[0][m]*y[0][m],x[1][m]*y[0][m-1],x[2][m]*y[0][m]),'ordinal root X')
                require(tuple(v[m] for v in oy)==(x[0][m]*y[0][m],x[0][m]*y[1][m],x[0][m-1]*y[2][m]),'ordinal root Y')
            cases += 1
    print('seam root pairs',cases,'m=1,2,3 OK')
    print('ordinal root choices',2*cases,'m=1,2,3 OK')
    for T,r in roots:
        Z,c=out_birth(T,r)
        p=as_polynomials(bin_values(T,r,len(T)))
        q=birth_poly(p)
        observed=bin_values(Z,c,len(Z)+1)
        require(all(evaluate(q[i],m)==observed[i][m] for i in range(3) for m in range(1,len(Z)+2)), ('birth',T,r))
    print('outdegree-one births',len(roots),'OK')
    funnel_cases=0
    for T,r in roots:
        for t in range(4):
            for word in product((0,1),repeat=t):
                Z,c=funnel(T,{r},word)
                source=as_polynomials(bin_values(T,r,len(T)))
                predicted=funnel_poly(source,word)
                actual=bin_values(Z,c,3)
                require(all(evaluate(predicted[i],m)==actual[i][m] for i in range(3) for m in range(1,4)), ('funnel',T,r,word))
                funnel_cases+=1
    print('one-port funnel cases',funnel_cases,'m=1,2,3 OK')
    joint_cases=0
    for T in (tiling(3,0),tiling(3,1),tiling(4,1),tiling(4,3)):
        n=len(T)
        for jmask in range(1,1 << n):
            ports={r for r in range(n) if jmask >> r & 1}
            for word in product((0,1),repeat=2):
                Z,c=funnel(T,ports,word)
                observed=bin_values(Z,c,2)
                v=[bin_values(T,r,3) for r in ports]
                for m in (1,2):
                    p=v[0][0][m]; pm1=v[0][0][m-1]
                    e=sum(q[1][m] for q in v);s=sum(q[2][m] for q in v)
                    sp=sum(q[2][m+1] for q in v)
                    L,D,B=tail_coefficients(word,m)
                    lp=tail_coefficients(word,m+1)[0]
                    prediction=(m*(lp*sp-m**len(word)*joint_values(T,ports,m)+B*p-D*e),L*p-(m-1)**len(word)*e,m**len(word)*s+D*pm1)
                    require(prediction==tuple(q[m] for q in observed),('joint funnel',T,ports,word,m,prediction,observed))
                joint_cases+=1
    print('multiport funnels',joint_cases,'m=1,2 OK')
    decomposition_cases=0
    for n in range(3,7):
        totals={1:set(),2:set()}
        ready={1:set(),2:set()}
        for bits in range(1 << ((n-1)*(n-2)//2)):
            T=tiling(n,bits)
            if not strong(T):
                continue
            h=tables(T)[0][-1]
            _,starts=tables(opposite(T))
            for r in range(n):
                s=starts[-1][r]
                if s in totals:
                    totals[s].add(h)
                    keep=[v for v in range(n) if v!=r]
                    U=induced(T,keep)
                    blocks=components(U)
                    require(all(len(block)==1 for block in blocks[1:]),('tail',T,r,s,blocks))
                    W=induced(U,blocks[0])
                    J=[i for i,v in enumerate(blocks[0]) if T[r] >> keep[v] & 1]
                    root_starts=tables(opposite(W))[1][-1]
                    require(sum(root_starts[i] for i in J)==s,('budget',T,r))
                    require(len(J)<=s,('ports',T,r))
                    if s==1:
                        require(len(J)==1 and root_starts[J[0]]==1,('one-port',T,r))
                    decomposition_cases+=1
                v=bin_values(T,r,2)
                d=v[1][2]-h
                if d in ready:
                    ready[d].add(h)
        print('order',n,'start budgets', {k:sorted(v) for k,v in totals.items()}, 'F-h', {k:sorted(v) for k,v in ready.items()})
    print('small-start decompositions',decomposition_cases,'OK')


if __name__=='__main__':
    main()
