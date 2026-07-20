#!/usr/bin/env python3
"""
death-star-2026-07-20-S61b (HYP claim below) -- the ODD/EVEN duality the owner asked for:
odd-valued functions (Redei: H(T)=#Ham paths is ODD) on tournaments, vs EVEN concepts
(even graphs E_n = cycle space; even functions = complement-invariant). Central object:
H is an ODD-VALUED EVEN FUNCTION. Plus the bridge Sum_T H(T) = n! * 2^C(n-1,2) tying the
odd invariant to the even-graph count 2^C(n-1,2) via S_n. Verified exactly n=3..6.
"""
from itertools import permutations, product
from math import comb, factorial

def ham_paths(adj, n):
    """# directed Hamiltonian paths in tournament adj (adj[i][j]=1 iff i->j)."""
    cnt = 0
    for perm in permutations(range(n)):
        if all(adj[perm[k]][perm[k+1]] for k in range(n-1)):
            cnt += 1
    return cnt

def tournaments(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for bits in range(1<<m):
        adj = [[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(edges):
            if (bits>>k)&1: adj[i][j]=1
            else: adj[j][i]=1
        yield adj, bits, edges

def complement(adj, n):
    return [[adj[j][i] for j in range(n)] for i in range(n)]

print("=== H is an ODD-VALUED EVEN FUNCTION, and the odd/even bridge ===\n")
for n in range(3,7):
    m = comb(n,2)
    allH = []
    all_odd = True
    complement_invariant = True
    for adj, bits, edges in tournaments(n):
        h = ham_paths(adj, n)
        allH.append(h)
        if h % 2 == 0: all_odd = False
        hc = ham_paths(complement(adj,n), n)
        if hc != h: complement_invariant = False
    S = sum(allH)
    bridge = factorial(n) * (2**comb(n-1,2))
    print(f"n={n}: #tournaments=2^{m}={1<<m}")
    print(f"  (1) ODD-VALUED (Redei): every H(T) odd = {all_odd}")
    print(f"  (2) EVEN FUNCTION (complement-invariant H(T)=H(T^op)): {complement_invariant}")
    print(f"  (3) BRIDGE  Sum_T H(T) = {S}   n!*2^C(n-1,2) = {bridge}   match={S==bridge}")
    print(f"      avg H = {S}/{1<<m} = n!/2^(n-1) = {factorial(n)}/{2**(n-1)}")
    print(f"  (4) EVEN-GRAPH count |cycle space of K_n| = 2^C(n-1,2) = {2**comb(n-1,2)}"
          f"  (= #tilings = #switching classes)")
    print()

# Walsh/parity: H is R-even => only even-degree Walsh coefficients (HYP-534, reconfirmed n=3,4)
def walsh_even_only(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    # f on {0,1}^m ; Walsh coeff hat f(S) = (1/2^m) sum_x f(x)(-1)^{<S,x>}
    vals = {}
    for adj,bits,ed in tournaments(n):
        vals[bits] = ham_paths(adj,n)
    maxodd = 0.0
    for S in range(1<<m):
        c = 0
        for x in range(1<<m):
            sign = -1 if bin(S & x).count("1")%2 else 1
            c += sign*vals[x]
        c /= (1<<m)
        if bin(S).count("1")%2==1 and abs(c)>1e-9:
            maxodd = max(maxodd, abs(c))
    return maxodd
print("=== H is R-EVEN: largest |odd-degree Walsh coefficient| (should be 0) ===")
for n in [3,4]:
    print(f"  n={n}: max |odd-degree Walsh coeff of H| = {walsh_even_only(n):.2e}  (0 => H is an even function)")

print("\n=== JC bridge: the counterexample's ODD fiber = H(3-cycle) ===")
# the JC counterexample F has a triple collision (fiber size 3). 3 = H(C_3), smallest odd H>1.
C3 = [[0,1,0],[0,0,1],[1,0,0]]   # 3-cycle 0->1->2->0
print(f"  H(3-cycle C_3) = {ham_paths(C3,3)}  = JC-counterexample fiber size (3) = smallest odd H>1")
print(f"  det JF = -2 (EVEN); fiber = 3 (ODD): the JC witness is 'odd fiber over even determinant',")
print(f"  the same odd(value)/even(structure) tension as H (odd values, even function).")
