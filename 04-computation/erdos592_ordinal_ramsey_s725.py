#!/usr/bin/env python3
"""
S725 — Erdos problem 592 (OPEN, $1000, ordinal partition calculus), creatively bridged to the repo.

592: determine the countable ordinals beta with  omega^beta -> (omega^beta, 3)^2  : every red/blue
coloring of edges of K_{omega^beta} has a RED clique of order type omega^beta OR a BLUE triangle K_3.
Known: Specker omega^2 YES; omega^n NO for 3<=n<omega; Chang omega^omega YES; beta=omega^2 OPEN.

This is INFINITE combinatorics, outside our finite repo's direct reach. The creative contribution is the
set of GENUINE bridges + an honest finite exploration:

(A) REFRAME via independence (THM-209, H=I(Omega,2)): "no blue K_3" <=> blue graph triangle-free <=>
    RED graph has independence number <= 2. A red clique = a blue-independent set. So
        592 <=>  every TRIANGLE-FREE graph on omega^beta has an independent set of order type omega^beta.
    Triangle-free graphs are the repo's chromatic-plane / Moser-spindle world; the "3" (blue triangle) is
    the repo's pervasive triangle (the C_3 that starts tournament theory).

(B) THE FINITE SEED IS IN THE REPO: the blue-triangle obstruction's finite floor is R(3,3)=6, whose
    UNIQUE extremal coloring of K_5 (no monochromatic triangle) is C_5 = the Paley graph on 5 vertices
    (THM-130/309 C_5-Paley, THM-436 C_5 threshold). 592 is the transfinite tower over this finite seed.

(C) THE TOWER: omega^beta in Cantor normal form is a SHELL TOWER (a tower of omega's) -- the same shape
    as the LRC cyclotomic witness tower (THM-439) and the covering-system smooth-modulus tower (S724).
    "Stepping up" beta -> beta+1 is a TRANSFER/recursion (the S720/721 temperature ladder). The
    non-monotone pattern YES(1,2) NO(3..finite) YES(omega) is a TOWER/renormalization phenomenon: the
    obstruction lives at finite heights, the limit omega is a fixed point that washes it out.

(D) honest finite exploration: the ULTRAMETRIC tree coloring (color a pair by its first-difference level)
    -- shows a blue level with branching >=3 forces a blue triangle, so ultrametric colorings CANNOT
    exploit omega-branching to avoid blue triangles => the omega^n (n>=3) counterexamples must be
    NON-ultrametric (use within-level order). That is a correct insight into why n>=3 is subtle.

No numpy/sympy.
"""
from itertools import combinations, product

def is_triangle_free(n, edges):
    E=set(map(frozenset,edges))
    for a,b,c in combinations(range(n),3):
        if frozenset((a,b)) in E and frozenset((b,c)) in E and frozenset((a,c)) in E: return False
    return True
def independence_number(n, edges):
    E=set(map(frozenset,edges)); best=0
    # brute max independent set (n small)
    for r in range(n,0,-1):
        for S in combinations(range(n),r):
            if all(frozenset((a,b)) not in E for a,b in combinations(S,2)):
                return r
    return 0

if __name__=="__main__":
    print("="*86)
    print("S725 — Erdos 592 (omega^beta -> (omega^beta,3)^2), bridged to the repo")
    print("="*86)

    # (B) the finite seed: R(3,3)=6, extremal C_5 = Paley graph on 5 vertices
    print("\n(B) FINITE SEED (in the repo): R(3,3)=6; extremal coloring of K_5 = C_5 = Paley graph P_5")
    # Paley graph on Z/5: i~j iff i-j is a QR mod 5 (QRs={1,4})
    V=range(5); blue=[(i,j) for i,j in combinations(V,2) if (i-j)%5 in (1,4) or (j-i)%5 in (1,4)]
    red=[(i,j) for i,j in combinations(V,2) if (i,j) not in set(blue)]
    print(f"   blue=P_5 edges {blue}")
    print(f"   blue triangle-free: {is_triangle_free(5,blue)}; red (complement) triangle-free: {is_triangle_free(5,red)}")
    print(f"   => K_5 has a 2-coloring with NO monochromatic triangle (C_5); K_6 does NOT (R(3,3)=6).")
    print(f"   This C_5/Paley extremal is THM-130/309/436 in the repo: the finite floor of 592's blue-K_3.")

    # (A) reframe: blue triangle-free <=> red independence <= 2
    print("\n(A) REFRAME: blue triangle-free  <=>  red graph independence number <= 2")
    print(f"   for C_5: red independence number = {independence_number(5,red)} (=2, since red is also C_5)")
    print("   => 592 = 'every triangle-free graph on omega^beta has an independent set of order type omega^beta'.")
    print("   (independence = the repo's H = I(Omega,2), THM-209.)")

    # (D) ultrametric tree coloring: blue level branching >=3 forces a blue triangle
    print("\n(D) ULTRAMETRIC TREE coloring on b-ary depth-n tree (model of omega^n truncated to branching b):")
    print("    d(x,y)=first differing level in {1..n}; color a pair by chi(level). KEY ultrametric fact:")
    print("    for x<y<z lex, d(x,z)=min(d(x,y),d(y,z)) => every triangle's 3 levels are {p,p,q}, p<=q.")
    for b in (2,3,4):
        n=3
        pts=list(product(range(b),repeat=n))
        def dlev(x,y):
            for i in range(n):
                if x[i]!=y[i]: return i+1
            return n+1
        # make level 2 BLUE, others RED; count blue triangles
        chi={1:'R',2:'B',3:'R'}
        blue_pairs=[(i,j) for i,j in combinations(range(len(pts)),2) if chi[dlev(pts[i],pts[j])]=='B']
        bt = not is_triangle_free(len(pts), blue_pairs)
        print(f"   branching b={b}, level-2 BLUE: blue triangle present? {bt}  "
              f"({'b>=3 at a blue level => blue triangle (cannot use omega-branching at a blue level)' if b>=3 else 'b=2: no blue triangle, but branching truncated to 2 (not omega)'})")
    print("   CONCLUSION: any BLUE level with omega-branching forces a blue triangle; so to avoid blue K_3")
    print("   an ultrametric coloring must keep blue levels finite-branching -- it cannot realize the")
    print("   omega^n counterexamples (n>=3). Hence Specker's NO-constructions are NON-ultrametric: they")
    print("   exploit the ORDER WITHIN levels, not just the tree. (Correct structural insight into 592.)")

    # (C) the tower / renormalization heuristic
    print("\n(C) THE TOWER (heuristic): omega^beta = a shell tower (CNF) = LRC cyclotomic tower (THM-439) /")
    print("    covering smooth-modulus tower (S724). beta-pattern YES(1,2) NO(3..finite) YES(omega):")
    print("    - finite beta: 'hot/generic' -- a (non-ultrametric) obstruction coloring exists for beta>=3;")
    print("    - limit omega: 'cold crystalline fixed point' -- the obstruction cannot be sustained across")
    print("      ALL finite levels at once (a compactness flip), so the relation is restored (Chang).")
    print("    Conjectural reading of OPEN beta=omega^2: it is a SECOND limit (limit of limits); if the")
    print("    'fixed point restores at limits' principle iterates, omega^2 should be YES (like omega), with")
    print("    the FALSE band sitting at the successor/finite levels between limits. (Matches Schipperus-type")
    print("    positive results at limit ordinals; a tower-renormalization heuristic, NOT a proof.)")
    print("="*86)
