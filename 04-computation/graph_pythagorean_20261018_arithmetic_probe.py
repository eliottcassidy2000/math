"""Exact graph/Pythagorean arithmetic controls, 2026-10-18.

Universe: all positive shore pairs x >= y with x+y <= 200; primitive
Euclid pairs 1 <= n < m <= 100; selected complete, bipartite and cubic
controls. Standard library only, integer arithmetic throughout.
"""
from itertools import combinations
from math import gcd, isqrt


def check(condition, message):
    if not condition:
        raise AssertionError(message)


def complete(n):
    return set(combinations(range(n), 2))


def bipartite(m, n):
    return {(i, m+j) for i in range(m) for j in range(n)}


def disjointness(edges):
    items = sorted(edges)
    return items, {(i, j) for i, j in combinations(range(len(items)), 2)
                   if set(items[i]).isdisjoint(items[j])}


def seidel(n, edges):
    return [[0 if i == j else -1 if tuple(sorted((i,j))) in edges else 1
             for j in range(n)] for i in range(n)]


def multiply(a, b):
    return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]


def switch(n, edges, shore):
    return {e for e in combinations(range(n),2)
            if (e in edges) != ((e[0] in shore) != (e[1] in shore))}


def degree(n, edges):
    return [sum(v in e for e in edges) for v in range(n)]


def paley_square_prime(p):
    """F_(p^2)=F_p[sqrt(d)], with d a quadratic nonresidue."""
    d = next(x for x in range(2,p) if pow(x,(p-1)//2,p)==p-1)
    points = [(a,b) for a in range(p) for b in range(p)]
    def chi(a,b):
        norm = (a*a-d*b*b)%p
        if not norm:
            return 0
        return 1 if pow(norm,(p-1)//2,p)==1 else -1
    edges = {(i+1,j+1) for i,j in combinations(range(p*p),2)
             if chi(points[i][0]-points[j][0],points[i][1]-points[j][1])==1}
    alpha = next(z for z in points if chi(*z)==-1)
    k = (p-1)//2
    selected = {1+points.index(((b+u*alpha[0])%p,(u*alpha[1])%p))
                for b in range(k) for u in range(p)}
    check(len(selected)==k*p,"independent affine-line switch size")
    return edges, selected


def main():
    count = primitive = 0
    for c in range(2, 201):
        for y in range(1, c//2+1):
            x = c-y
            z = isqrt(x*y)
            if z*z != x*y:
                continue
            a, b = x-y, 2*z
            check(a*a+b*b == c*c, "square-shore light cone")
            count += 1
            if a > 0 and gcd(gcd(a,b),c) == 1:
                m, n = isqrt(x), isqrt(y)
                check(m*m == x and n*n == y, "primitive shores are squares")
                check(gcd(m,n) == 1 and (m-n)%2 == 1, "primitive Euclid gate")
                primitive += 1
    print(f"square-shore universe c<=200: integral={count}, primitive_positive={primitive}")

    pairs = gap = 0
    for m in range(2,101):
        for n in range(1,m):
            if gcd(m,n) != 1 or (m-n)%2 == 0:
                continue
            a,b,c = m*m-n*n,2*m*n,m*m+n*n
            check(c-b == (m-n)**2, "hypotenuse/even-leg square gap")
            check(c-a == 2*n*n, "hypotenuse/odd-leg double square gap")
            check((c-b == 1) == (m-n == 1), "unit-gap characterization")
            pairs += 1
            gap += c-b == 1
    print(f"primitive Euclid universe m<=100: pairs={pairs}, unit_gap={gap}")

    print("central branch: m n a b c a-b")
    m,n = 2,1
    for _ in range(8):
        a,b,c = m*m-n*n,2*m*n,m*m+n*n
        check(abs(a-b) == 1, "Pell consecutive legs")
        print(m,n,a,b,c,a-b)
        m,n = 2*m+n,m

    print("parabolic unit-gap branch: m n a b c")
    for n in range(1,7):
        m=n+1
        print(m,n,m*m-n*n,2*m*n,m*m+n*n)

    items, pet = disjointness(complete(5))
    check(len(items)==10 and len(pet)==15 and set(degree(10,pet))=={3}, "Petersen size")
    standard_map = {}
    for i in range(5):
        standard_map[("o",i)] = items.index(tuple(sorted((i,(i+2)%5))))
        standard_map[("i",i)] = items.index(tuple(sorted(((i+3)%5,(i+4)%5))))
    gp_edges = set()
    for i in range(5):
        for p,q in [(("o",i),("o",(i+1)%5)),(("i",i),("i",(i+2)%5)),(("o",i),("i",i))]:
            gp_edges.add(tuple(sorted((standard_map[p],standard_map[q]))))
    check(gp_edges==pet, "explicit standard Petersen isomorphism")

    ps = seidel(10,pet)
    check(multiply(ps,ps)==[[9 if i==j else 0 for j in range(10)] for i in range(10)], "Petersen Seidel conference")
    print("Petersen: Seidel S^2=9 I; standard GP(5,2) map verified")

    v = items.index((0,1))
    neighbors = {i for i in range(10) if tuple(sorted((v,i))) in pet}
    switched = switch(10,pet,neighbors)
    check(degree(10,switched)[v]==0, "switched isolated root")
    remap = {}
    t = [2,3,4]
    for i in range(3):
        remap[(0,i)] = items.index((0,t[i]))
        remap[(1,i)] = items.index((1,t[i]))
        remap[(2,i)] = items.index(tuple(x for x in t if x!=t[i]))
    expected = {tuple(sorted((remap[(a,i)],remap[(b,j)])))
                for (a,i),(b,j) in combinations(remap,2) if a!=b and i!=j}
    check(switched==expected, "Petersen switched to K1 plus D(K3,3)")
    for (a,i),(b,j) in combinations(remap,2):
        p=((a+i)%3,(a-i)%3)
        q=((b+j)%3,(b-j)%3)
        check((a!=b and i!=j)==((p[0]==q[0]) != (p[1]==q[1])), "D(K3,3) rook map")
    print("Petersen switch N(v): K1 + D(K3,3); F3 map to L(K3,3) verified")

    print("Seidel root radii: Petersen=3, K5=4, K3,3=5; 3^2+4^2=5^2")
    for n in range(2,21):
        sk=seidel(n,complete(n))
        # Characteristic equation: eigenvalues 1 and -(n-1).
        lhs=multiply(sk,sk)
        check(all(lhs[i][j]+(n-2)*sk[i][j]-(n-1)*(i==j)==0
                  for i in range(n) for j in range(n)), "complete Seidel minimal polynomial")
        for p in range(1,n):
            sb=seidel(n,bipartite(p,n-p))
            # Switching the p-shore transforms complete bipartite to empty.
            check(switch(n,bipartite(p,n-p),set(range(p)))==set(), "bipartite switching erases cut")
            lhs=multiply(sb,sb)
            check(all(lhs[i][j]-(n-2)*sb[i][j]-(n-1)*(i==j)==0
                      for i in range(n) for j in range(n)), "bipartite Seidel polynomial")
    print("Seidel complete/bipartite controls: all orders 2..20 and every shore split")

    print("Paley-square conference continuation: t order degree Seidel-radius triangle clique-order bipartite-shore")
    for t in [3,5,7,11]:
        ed, selected = paley_square_prime(t)
        order=t*t+1
        h=switch(order,ed,selected)
        check(set(degree(order,h))=={t*(t-1)//2}, "regular switched Paley conference")
        s=seidel(order,h)
        check(multiply(s,s)==[[t*t if i==j else 0 for j in range(order)]
                             for i in range(order)],"Paley conference square")
        a,b,c=t,(t*t-1)//2,(t*t+1)//2
        check(a*a+b*b==c*c and c-b==1,"conference Pythagorean family")
        print(t,order,t*(t-1)//2,t,(a,b,c),c,(t*t+3)//4)

    # Independently count the disjointness degree; no spectral approximation.
    print("clique edge-disjointness continuation: order_c, vertex_count, degree")
    for c in [5,6,13,29]:
        vertices, ed = disjointness(complete(c))
        dg=degree(len(vertices),ed)
        check(set(dg)=={(c-2)*(c-3)//2}, "Kneser degree")
        print(c,len(vertices),dg[0])

    # Small graph polarization checks, including one missing edge.
    print("clique/cut polarization controls: m n Pythagorean_defect")
    for m,n in [(2,1),(3,2),(5,2),(4,3)]:
        fg=m+2*len(complete(m))
        fh=n+2*len(complete(n))
        a,b,c=fg-fh,2*len(bipartite(m,n)),fg+fh
        check(c*c-a*a-b*b==0, "clique/cut exact light cone")
        print(m,n,c*c-a*a-b*b)
        if m>=2:
            fg-=2
            a,c=fg-fh,fg+fh
            check(c*c-a*a-b*b==-8*n*n, "one missing edge shell defect")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
