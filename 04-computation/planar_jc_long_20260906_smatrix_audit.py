#!/usr/bin/env python3
"""Independent exact audit: raw dual testers and actual Student lattices.

Does not import the producer.  Integer invariants are computed by elementary
Euclidean row/column operations, not a Smith library or the two-path DP.
Small determinantal ideals use direct matrix-support matching enumeration.
"""
from fractions import Fraction
from functools import reduce
from itertools import permutations
from math import gcd, lcm, prod
from pathlib import Path
import hashlib
import json
import sympy as sp

GATES = 0


def check(ok, label):
    global GATES
    if not ok:
        raise RuntimeError(label)
    GATES += 1


def emit(tag, *items):
    print(tag, json.dumps(items, separators=(",", ":")), flush=True)


def matrix(m):
    # Direct coefficient action on the monomial basis; no parity decomposition.
    A = [[0]*m for _ in range(m+1)]
    for j in range(m):
        A[j+1][j] -= 2*m
        A[j+1][j] += j
        if j:
            A[j-1][j] += 6*j
    return A


def integer_diagonal(A):
    """Unimodular Euclidean diagonalization with the divisibility repair.

    Once pivot row/column are clear, a trailing entry not divisible by the
    pivot is moved into its row. The next Euclidean step strictly reduces
    the nonzero pivot absolute value. This ensures Smith divisibility.
    """
    M = [row[:] for row in A]
    height, width = len(M), len(M[0])
    out = []
    for k in range(min(height,width)):
        choices = [(abs(M[i][j]),i,j) for i in range(k,height)
                   for j in range(k,width) if M[i][j]]
        if not choices:
            break
        _,i,j = min(choices)
        M[k],M[i] = M[i],M[k]
        for row in M:
            row[k],row[j] = row[j],row[k]
        while True:
            restart = False
            for i in range(k+1,height):
                if M[i][k]:
                    quotient = M[i][k]//M[k][k]
                    M[i] = [a-quotient*b for a,b in zip(M[i],M[k])]
                    if M[i][k]:
                        M[k],M[i] = M[i],M[k]
                    restart = True
                    break
            if restart:
                continue
            for j in range(k+1,width):
                if M[k][j]:
                    quotient = M[k][j]//M[k][k]
                    for i in range(height):
                        M[i][j] -= quotient*M[i][k]
                    if M[k][j]:
                        for row in M:
                            row[k],row[j] = row[j],row[k]
                    restart = True
                    break
            if restart:
                continue
            offender = next(((i,j) for i in range(k+1,height)
                             for j in range(k+1,width)
                             if M[i][j] % M[k][k]),None)
            if offender is None:
                break
            i,_ = offender
            M[k] = [a+b for a,b in zip(M[k],M[i])]
        if M[k][k] < 0:
            M[k] = [-a for a in M[k]]
        out.append(M[k][k])
    check(all(M[i][j] == (out[i] if i == j and i < len(out) else 0)
              for i in range(height) for j in range(width)), "complete integer diagonal")
    check(all(b%a == 0 for a,b in zip(out,out[1:])), "diagonal divisibility")
    return out


def matching_ideals(A):
    # Enumerate actual matrix support by columns, with a used-row bitmask.
    # This is not nonadjacent-edge recursion and does not use parity paths.
    m = len(A[0])
    edges = [[(i,A[i][j]) for i in range(len(A)) if A[i][j]] for j in range(m)]
    ideals = [0]*(m+1)
    count = [0]*(m+1)

    def visit(j,used,size,value):
        if j == m:
            ideals[size] = gcd(ideals[size],abs(value))
            count[size] += 1
            return
        visit(j+1,used,size,value)
        for i,a in edges[j]:
            if not used & (1<<i):
                visit(j+1,used | (1<<i),size+1,value*a)
    visit(0,0,0,1)
    return ideals,count


def component_diagonal(A):
    # Extract connected support components from the actual matrix, without
    # assuming their parity description or path ordering. This prevents
    # unnecessary intermediate coefficient growth in large dense reductions.
    height,width = len(A),len(A[0])
    rowcols = [[j for j in range(width) if A[i][j]] for i in range(height)]
    colrows = [[i for i in range(height) if A[i][j]] for j in range(width)]
    remaining = set(range(width))
    diagonals = []
    while remaining:
        seed = min(remaining)
        columns,rows = {seed},set()
        pending = [seed]
        while pending:
            j = pending.pop()
            for i in colrows[j]:
                if i not in rows:
                    rows.add(i)
                    for jj in rowcols[i]:
                        if jj not in columns:
                            columns.add(jj)
                            pending.append(jj)
        remaining -= columns
        block = [[A[i][j] for j in sorted(columns)] for i in sorted(rows)]
        diagonals.extend(integer_diagonal(block))
    # Smith-reduce the direct sum of already diagonal matrices. On two
    # coordinates diag(a,b) is equivalent to diag(gcd(a,b),lcm(a,b)).
    for i in range(len(diagonals)):
        for j in range(i+1,len(diagonals)):
            a,b = diagonals[i],diagonals[j]
            d = gcd(a,b)
            diagonals[i],diagonals[j] = d,a//d*b
    check(all(b%a == 0 for a,b in zip(diagonals,diagonals[1:])),"merged Smith divisibility")
    return diagonals


def declared_path_formula(m):
    # Independent valuation recursion on the abstract edge paths.  Gcds of
    # products become minima of valuations; this uses no large product DP.
    # Every possible prime is bounded by max(m,3), proved by a maximal minor.
    paths = [[],[]]
    for j in range(0,m,2):
        if j:
            paths[0].append(6*j)
        paths[0].append(2*m-j)
    for j in range(1,m,2):
        paths[1].extend((6*j,2*m-j))
    divisors = [1]*(m+1)
    for p in list(sp.primerange(2,max(m,3)+1)):
        minima = []
        for weights in paths:
            capacity = (len(weights)+1)//2
            old2 = old1 = [0]+[10**9]*capacity
            for w in weights:
                valuation = 0
                while w%p == 0:
                    valuation += 1
                    w //= p
                new = [0]+[min(old1[k],valuation+old2[k-1])
                            for k in range(1,capacity+1)]
                old2,old1 = old1,new
            minima.append(old1)
        for k in range(m+1):
            v = min(minima[0][a]+minima[1][k-a]
                    for a in range(len(minima[0]))
                    if 0 <= k-a < len(minima[1]))
            check(v < 10**9, "full column matching exists")
            divisors[k] *= p**v
    return divisors


def generic_testers():
    for n,r,s,seed in [(3,1,2,0),(4,2,2,1),(4,2,3,2),(4,3,2,3)]:
        A = sp.Matrix(n,r,lambda i,j: (i+2)**(j+1)+seed*(i==j))
        B = sp.Matrix(n,s,lambda i,j: (i+1)**(j+r+1)-seed)
        variables = sp.symbols('u:'+str(r*n))
        U = sp.Matrix(r,n,variables)
        G = U*A
        residual = sp.expand(G.det())*B-A*G.adjugate()*U*B
        expected_support = set()
        minors = []
        coefficients = []
        for ordered in permutations(range(n),r):
            exponents = tuple(int(any(a*n+ordered[a] == k for a in range(r)))
                              for k in range(r*n))
            expected_support.add(exponents)
            for j in range(n):
                for alpha in range(s):
                    rows = list(ordered)+[j]
                    bordered = A[rows,:].row_join(B[rows,alpha]).det()
                    coefficient = sp.Poly(sp.expand(residual[j,alpha]),*variables).coeff_monomial(exponents)
                    check(coefficient == bordered, "ordered bordered coefficient")
                    coefficients.append(int(coefficient))
                    minors.append(int(bordered))
        for value in residual:
            poly = sp.Poly(sp.expand(value),*variables)
            check(all(monom in expected_support or coefficient == 0
                      for monom,coefficient in poly.terms()), "no extra tester monomials")
        check(reduce(gcd,coefficients,0) == reduce(gcd,minors,0), "coefficient ideal over Z")
        emit("DUAL",n,r,s,"coefficient-gcd",abs(reduce(gcd,coefficients,0)))
    # Explicit module and rank-drop boundaries; arbitrary-ring identity does
    # not supply either a full-rank chart or global integral membership.
    check(not any((2*a)%6 == 1 for a in range(6)), "zero-divisor ring membership hostile")
    check(sp.Matrix([[2,1],[0,0]]).det() == 0, "all mixed minors zero at integral hostile")
    A = sp.zeros(3,2)
    B = sp.Matrix([1,0,0])
    U = sp.Matrix(2,3,sp.symbols('v:6'))
    check((U*A).det()*B-A*(U*A).adjugate()*U*B == sp.zeros(3,1), "rank-drop tester blind")
    check(A.rank() < A.row_join(B).rank(), "rank-drop nonmembership")
    A = sp.Matrix([1,sp.I,0])
    B = sp.Matrix([0,0,1])
    check(A.T*A == sp.zeros(1,1) and A.T*B == sp.zeros(1,1), "complex Gram blind")
    check(A.row_join(B).rank() == 2, "raw complex measurement survives")


generic_testers()
rows = {}
for m in range(2,49):
    A = matrix(m)
    diagonal = component_diagonal(A)
    if m <= 14:
        check(diagonal == integer_diagonal(A),"undecomposed Euclidean control")
    check(len(diagonal) == m, "actual Student full rank")
    direct_divisors = [1]
    for d in diagonal:
        direct_divisors.append(direct_divisors[-1]*d)
    formula = declared_path_formula(m)
    check(direct_divisors == formula, "Euclidean diagonal versus valuation-path formula m="+str(m))
    if m <= 10:
        exhaustive,counts = matching_ideals(A)
        check(exhaustive == direct_divisors, "all support matching ideals m="+str(m))
    rows[m] = diagonal
    if m in (2,5,8,14,24,48):
        emit("EUCLIDEAN_SMITH",m,diagonal,"torsion-order",prod(diagonal))

inherited = {5:[21,0,14,0,36,0],6:[77,0,42,0,84,0,360],
             7:[143,0,66,0,108,0,360,0],8:[715,0,286,0,396,0,1080,0,5040]}
for m,expected in inherited.items():
    null = sp.Matrix(matrix(m)).T.nullspace()
    check(len(null) == 1, "dense rational nullity")
    den = lcm(*(v.q for v in null[0]))
    row = [int(v*den) for v in null[0]]
    content = reduce(gcd,row,0)
    row = [v//content for v in row]
    if row[0] < 0:
        row = [-v for v in row]
    check(row == expected, "inherited Student vector from dense nullspace")

A2 = sp.Matrix(matrix(2))
B2 = sp.Matrix([[0,2],[1,0],[0,-1]])
X2 = A2.gauss_jordan_solve(B2)[0]
check(X2 == sp.Matrix([[sp.Rational(-1,4),0],[0,sp.Rational(1,3)]]), "two-column dense inverse")
check(lcm(*(q.q for q in X2)) == 12 == rows[2][-1], "two-column sharp common denominator")

A14 = sp.Matrix(matrix(14))
B14 = sp.Matrix(15,7,lambda i,j: int(i == 2*j+1))
X14 = A14.gauss_jordan_solve(B14)[0]
check(A14*X14 == B14, "all raw row14 equations")
denominators = [lcm(*(q.q for q in X14[:,j])) for j in range(7)]
check(denominators == [28,182,2184,4004,20020,18018,16016], "row14 individual denominators")
check(lcm(*denominators) == 720720 == rows[14][-1], "row14 sharp complete-bank exponent")
check(rows[14] == [1,1,1,1,1,1,1,2,2,6,12,12,840,720720], "row14 entire Smith form")
check(prod(rows[14]) == 2092278988800, "row14 torsion order")
check(720720 == 2**4*3**2*7*715 and 715 == 5*11*13,"localized exponent")
emit("DENSE_BANK",2,12,"row14",denominators,"joint",720720,"localized",715)
emit("SCOPE","field iff needs full column rank","ring iff requires a unit minor chart",
     "integral saturation quotient is retained","torsion vanishes over characteristic-zero fields")
emit("UNIVERSE","actual integer matrices m=2..48","all support matchings m=2..10",
     "four generic-dual size/seed controls","four inherited dense nullspaces","two sharp dense banks")
emit("PASS",GATES)
emit("SOURCE_SHA256",hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
