#!/usr/bin/env python3
"""Fresh exact verification for the reconstructed Session225 graft packet.
Standard library only. Explicit exceptions remain active under python -O.
Every tournament is encoded by outgoing-neighbor bit masks.
"""
import itertools
import json
from fractions import Fraction

K = [[5,-5,27,-22], [1,-1,6,-6], [5,-4,29,-20], [1,-1,7,-6]]
B4 = [8,1,3,6]
O = [[1,0,0,0], [0,0,0,0], [0,0,2,0], [0,0,0,1]]
Z = [[1,-1,1,0], [0,0,0,1], [0,0,2,-1], [0,0,0,1]]
J = [[1,-1,1,0], [1,-1,0,0], [1,0,1,0], [1,0,0,0]]

def require(test, message):
    if not test:
        raise RuntimeError(message)

def valid(rows, root):
    n = len(rows)
    require(n > 0 and 0 <= root < n, "root domain")
    for i, row in enumerate(rows):
        require(0 <= row < (1 << n) and not (row >> i & 1), "row domain")
        for j in range(i):
            require(((row >> j) & 1) + ((rows[j] >> i) & 1) == 1, "orientation")

def counts(rows, root):
    """Literal subset DP, independent of the proposed graft formulas."""
    valid(rows, root)
    n = len(rows)
    size = 1 << n
    dp = [[0] * n for _ in range(size)]
    start = [[0] * n for _ in range(size)]
    for v in range(n):
        dp[1 << v][v] = 1
    start[1 << root][root] = 1
    for mask in range(1, size):
        for last in range(n):
            a, b = dp[mask][last], start[mask][last]
            if a == 0 and b == 0:
                continue
            available = rows[last] & ~mask
            while available:
                bit = available & -available
                available -= bit
                v = bit.bit_length() - 1
                dp[mask | bit][v] += a
                start[mask | bit][v] += b
    hp = [sum(row) for row in dp]
    hp[0] = 1
    starts = [sum(row) for row in start]
    ends = [row[root] for row in dp]
    full = size - 1
    alpha = sum(starts[m] * hp[full ^ m] for m in range(size))
    return hp, starts, ends, [hp[full], ends[full], alpha, hp[full ^ (1 << root)]]

def profile(rows, root):
    return counts(rows, root)[3]

def subsets(mask):
    s = mask
    while True:
        yield s
        if s == 0:
            return
        s = (s - 1) & mask

def stats(rows, root):
    hp, starts, ends, x = counts(rows, root)
    full = (1 << len(rows)) - 1
    without = full ^ (1 << root)
    Q = sum(hp[m] * hp[full ^ m] for m in range(full + 1))
    delta = sum(ends[m] * hp[full ^ m] for m in range(full + 1))
    T = sum(starts[m] * sum(hp[v] * hp[(full ^ m) ^ v]
                            for v in subsets(full ^ m))
            for m in range(full + 1))
    b = sum(hp[v] * hp[without ^ v] for v in subsets(without))
    return dict(H=x[0], Q=Q, s=starts[full], e=x[1], alpha=x[2],
                delta=delta, T=T, q=x[3], b=b)

def operator(rows, root):
    d = stats(rows, root)
    H,Q,s,e,a,de,T,q,b = [d[k] for k in
                         ("H","Q","s","e","alpha","delta","T","q","b")]
    return [[H,-H,Q-H,2*H-Q],
            [e,-e,de-e,e-de],
            [a,s-a,2*T-a,3*a-s-2*T],
            [q,-q,b-q,2*q-b] if len(rows)>1 else [1,0,0,0]]

def graft(a, root, b, gadget_root=0):
    """Keep both induced tournaments; B -> root, A-root -> B."""
    n = len(a)
    new_bits = ((1 << len(b)) - 1) << n
    rows = [row | (new_bits if i != root else 0)
            for i, row in enumerate(a)]
    rows.extend((row << n) | (1 << root) for row in b)
    return rows, n + gadget_root

def mv(a, x):
    return [sum(u*v for u,v in zip(row,x)) for row in a]

def mm(a, b):
    return [[sum(a[i][k]*b[k][j] for k in range(4))
             for j in range(4)] for i in range(4)]

def determinant(matrix):
    a = [[Fraction(v) for v in row] for row in matrix]
    out = Fraction(1)
    for i in range(len(a)):
        pivot_row = next((j for j in range(i,len(a)) if a[j][i]), None)
        if pivot_row is None:
            return 0
        if pivot_row != i:
            a[i],a[pivot_row] = a[pivot_row],a[i]
            out = -out
        pivot = a[i][i]
        out *= pivot
        for j in range(i+1,len(a)):
            factor = a[j][i]/pivot
            for k in range(i+1,len(a)):
                a[j][k] -= factor*a[i][k]
    require(out.denominator == 1, "integer determinant")
    return int(out)

def reaches_all(rows, v):
    seen = todo = 1 << v
    while todo:
        bit = todo & -todo
        todo -= bit
        i = bit.bit_length() - 1
        nxt = rows[i] & ~seen
        seen |= nxt
        todo |= nxt
    return seen == (1 << len(rows)) - 1

def strong(rows):
    return all(reaches_all(rows,v) for v in range(len(rows)))

def tournaments(n):
    pairs = list(itertools.combinations(range(n),2))
    for code in range(1 << len(pairs)):
        rows = [0]*n
        for k,(i,j) in enumerate(pairs):
            if code >> k & 1:
                rows[i] |= 1 << j
            else:
                rows[j] |= 1 << i
        yield rows


def permutation_profile(rows, root):
    """Separate direct permutation/split enumeration for small controls."""
    vertices = list(range(len(rows)))
    def good(path):
        return all(rows[u] >> v & 1 for u,v in zip(path,path[1:]))
    D = E = U = 0
    for order in itertools.permutations(vertices):
        if good(order):
            D += 1
            E += int(order[-1] == root)
        for cut in range(len(order)):
            if order[cut] == root and good(order[:cut]) and good(order[cut:]):
                U += 1
    d = sum(good(order) for order in itertools.permutations(
        [v for v in vertices if v != root]))
    return [D,E,U,d]

def permutation_stats(rows, root):
    def good(path):
        return all(rows[u] >> v & 1 for u,v in zip(path,path[1:]))
    n = len(rows)
    H = Q = s = e = alpha = delta = T = 0
    for order in itertools.permutations(range(n)):
        if good(order):
            H += 1
            s += int(order[0] == root)
            e += int(order[-1] == root)
        for i in range(n+1):
            if good(order[:i]) and good(order[i:]):
                Q += 1
                alpha += int(i > 0 and order[0] == root)
                delta += int(i > 0 and order[i-1] == root)
        if order[0] == root:
            for i in range(1,n+1):
                for j in range(i,n+1):
                    if good(order[:i]) and good(order[i:j]) and good(order[j:]):
                        T += 1
    q = b = 0
    for order in itertools.permutations([v for v in range(n) if v != root]):
        q += int(good(order))
        for i in range(len(order)+1):
            b += int(good(order[:i]) and good(order[i:]))
    return dict(H=H,Q=Q,s=s,e=e,alpha=alpha,delta=delta,T=T,q=q,b=b)

def main():
    bank = [rows for n in range(1,5) for rows in tournaments(n)]
    require(len(bank)==75,"bank size")
    xs = [profile(a,0) for a in bank]
    matrices = [operator(b,0) for b in bank]
    for a,x in zip(bank,xs):
        require(permutation_profile(a,0)==x,"permutation profile control")
        require(permutation_stats(a,0)==stats(a,0),"permutation gadget-stat control")
    first,first_root = graft([2,0],0,B4)
    require(permutation_profile(first,first_root)==profile(first,first_root),
            "six-vertex permutation control")
    pairs = strong_checks = 0
    for a,x in zip(bank,xs):
        require(x[2] >= x[3] + x[0] - x[1], "source positive-split inequality")
        for b,m in zip(bank,matrices):
            rows,r = graft(a,0,b)
            require(profile(rows,r)==mv(m,x),"general graft closure")
            pairs += 1
            if len(a)>=2 and reaches_all(a,0) and strong(b):
                require(strong(rows),"strongness lemma")
                strong_checks += 1
    require(operator(B4,0)==K,"literal gadget matrix")
    require(stats(B4,0)==dict(H=5,Q=32,s=1,e=1,alpha=5,delta=7,T=17,q=1,b=8),
            "gadget statistics")
    require(operator([0],0)==J,"singleton exception")
    generator_checks = 0
    for a,x in zip(bank,xs):
        n = len(a)
        ao = [row | (1 << n) for row in a]+[0]
        az,r = graft(a,0,[0])
        require(profile(ao,0)==mv(O,x),"O")
        require(profile(az,0)==mv(Z,x),"Z")
        require(profile(az,r)==mv(J,x),"J")
        generator_checks += 3
    basis_graphs = [([0],0),([2,0],0),([2,4,1],0),([6,4,0],0)]
    basis = [profile(a,r) for a,r in basis_graphs]
    require(basis==[[1,1,1,1],[1,0,2,1],[3,1,3,1],[1,0,4,1]],"basis profiles")
    require(determinant(basis)==-4,"physical span")
    require(determinant(K)==-8,"K determinant")
    require([determinant(g) for g in (O,Z,J)]==[0,0,0],"singular generators")
    k2,k3 = mm(K,K),mm(mm(K,K),K)
    k4 = mm(k3,K)
    identity = [[k4[i][j]-27*k3[i][j]-37*k2[i][j]+31*K[i][j]
                 -8*int(i==j) for j in range(4)] for i in range(4)]
    require(all(v==0 for row in identity for v in row),"recurrence identity")
    a,root,x = [2,0],0,[1,0,2,1]
    literal_iterates = []
    for step in range(1,4):
        a,root = graft(a,root,B4)
        x = mv(K,x)
        require(profile(a,root)==x,"literal iterate")
        require(strong(a),"literal iterate strong")
        literal_iterates.append(dict(step=step,vertices=len(a),profile=x[:],strong=True))
    x,sequence = [1,0,2,1],[]
    for _ in range(8):
        sequence.append(x[0])
        x = mv(K,x)
    require(sequence[:5]==[1,37,1113,31349,886465],"initial recurrence values")
    for i in range(4,len(sequence)):
        require(sequence[i]==27*sequence[i-1]+37*sequence[i-2]
                -31*sequence[i-3]+8*sequence[i-4],"scalar recurrence")
    print(json.dumps(dict(status="PASS",sources=75,gadgets=75,
          universe="All labeled tournaments of orders 1..4, root fixed at label 0",
          general_closure_pairs=pairs,strongness_controls=strong_checks,
          permutation_profile_controls=76,permutation_stat_controls=75,
          generator_controls=generator_checks,gadget_rows=B4,gadget_root=0,
          gadget_statistics=stats(B4,0),K=K,K_determinant=-8,
          physical_basis=basis,physical_basis_determinant=-4,
          literal_iterates=literal_iterates,sequence=sequence,
          recurrence=[27,37,-31,8],
          scope="Fresh finite controls plus exact matrix arithmetic; all-order proof is separate",
          original_artifacts_recovered=False),indent=2,sort_keys=True))

if __name__ == "__main__":
    main()
