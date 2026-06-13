#!/usr/bin/env python3
"""
FREE POSITION UNIVERSALITY v2 — corrected approach.

The average W over ALL tournaments (complete graph model) gives only
the leading r^f term. F_f(r) is NOT the average W.

But F_f(r) IS the W-polynomial of ALL (f+1)! orderings of f+1 vertices
where each step has sign ±1/2 (the "sign-sum"):

F_f(r) = sum_{sigma in S_{f+1}} prod_{i=0}^{f-1} (r + sign(sigma(i+1)-sigma(i))/2)

This is already proved in THM-059. The question is: why does this appear
as the coefficient in the W-decomposition?

Actually, let me re-read opus's THM-059 statement more carefully.
The coefficient C_I(r) for invariant I is 2^{parts(I)} * F_f(r).
And F_f(r) = sum_k A(f+1,k) * (r+1/2)^{f-k} * (r-1/2)^k.

The universality question is: why does the SAME F_f appear regardless
of which specific cycle vertices are pinned?

NEW APPROACH: Think about it in terms of the BACKBONE MODEL.

In the backbone model: vertices 0,1,...,n-1 with backbone 0->1->...->n-1.
Non-backbone edges are variables.

A 3-cycle on {a,b,c} involves backbone edges among {a,b,c} and non-backbone
edges. The "pinned" positions are where the cycle sits in the path.

The free positions involve the remaining vertices, which form an arbitrary
subpath with known backbone structure.

The key insight might be that the BACKBONE MODEL treats all positions
symmetrically — each non-backbone edge has ±1/2 independently.

Let me verify: compute F_f from the backbone model on f+1 CONSECUTIVE
vertices vs f+1 NON-CONSECUTIVE vertices.

kind-pasteur-2026-03-07-S27
"""
from itertools import permutations
from fractions import Fraction
from collections import defaultdict
from math import comb

def eulerian_number(n, k):
    return sum((-1)**j * comb(n+1, j) * (k+1-j)**n for j in range(k+2))

def master_poly_F(f):
    result = defaultdict(lambda: Fraction(0))
    for k in range(f+1):
        Ank = eulerian_number(f+1, k)
        poly_term = {0: Fraction(Ank)}
        for _ in range(f-k):
            new = {}
            for pw, c in poly_term.items():
                new[pw+1] = new.get(pw+1, Fraction(0)) + c
                new[pw] = new.get(pw, Fraction(0)) + c * Fraction(1, 2)
            poly_term = new
        for _ in range(k):
            new = {}
            for pw, c in poly_term.items():
                new[pw+1] = new.get(pw+1, Fraction(0)) + c
                new[pw] = new.get(pw, Fraction(0)) + c * Fraction(-1, 2)
            poly_term = new
        for pw, c in poly_term.items():
            result[pw] += c
    return dict(result)

def poly_str(poly, var='r'):
    terms = []
    for power in sorted(poly.keys(), reverse=True):
        c = poly[power]
        if c == 0: continue
        terms.append(f"{c}*{var}^{power}" if power > 0 else str(c))
    return " + ".join(terms) if terms else "0"

# The SIGN-SUM definition of F_f:
# F_f(r) = sum_{sigma in S_{f+1}} prod_{i=0}^{f-1} (r + sign(sigma(i+1)-sigma(i))/2)
# where sign(x) = +1 if x > 0, -1 if x < 0.

print("F_f via sign-sum:")
for f in range(5):
    result = defaultdict(lambda: Fraction(0))
    for perm in permutations(range(f+1)):
        poly = {0: Fraction(1)}
        for i in range(f):
            s = Fraction(1, 2) if perm[i+1] > perm[i] else Fraction(-1, 2)
            new = {}
            for pw, c in poly.items():
                new[pw+1] = new.get(pw+1, Fraction(0)) + c
                new[pw] = new.get(pw, Fraction(0)) + c * s
            poly = new
        for pw, c in poly.items():
            result[pw] += c

    Ff = master_poly_F(f)
    match = all(result.get(p, Fraction(0)) == Ff.get(p, Fraction(0))
               for p in set(list(result.keys()) + list(Ff.keys())))
    print(f"  F_{f}(r) via sign-sum = {poly_str(dict(result))}, match={match}")

# So F_f(r) uses sign(sigma(i+1)-sigma(i)) which is +1 for ascent, -1 for descent.
# This is a property of the ORDERING, not of any specific tournament.
# The fact that W's decomposition uses this is because:
#
# In the W-polynomial, each step's contribution is (r + s_e) where s_e = ±1/2
# depends on the edge orientation. When we extract the coefficient of a
# specific Fourier monomial (cycle invariant), the free edges are NOT
# constrained. Their contribution sums over all permutations of the free
# vertices, and at each step, the sign depends on whether the edge goes
# "forward" (ascending vertex label) or "backward" (descending label) in
# the path... NO WAIT. The sign depends on the TOURNAMENT, not the vertex labels.
#
# But in the sign-sum formula, the sign depends on the vertex label ordering.
# How do these two things match?
#
# The answer: for the "free" contribution in the Fourier decomposition,
# each free edge has BOTH orientations summed over (since we're extracting
# the degree-0 Fourier component of the free edges). The sum over both
# orientations of edge (u,v) gives:
# (r + 1/2) + (r - 1/2) = 2r ... but that would give just r^f.
#
# WAIT. That's not right. The sum is over PRODUCTS, not individual factors.
# If we have f free steps with independent ±1/2 signs, the total is:
# sum over all 2^{free_edges} orientations of prod_{step i} (r + s_i)
# But the signs at different steps are NOT independent — they share vertices!
# Edge (u,v) appears in the product whenever u and v are consecutive in the path.
#
# OK let me think about this completely differently.
#
# F_f(r) = sum over sigma in S_{f+1} of prod (r + sgn(sigma(i+1)-sigma(i))/2)
#
# This can be rewritten: sigma(i+1) > sigma(i) iff the pair (sigma(i), sigma(i+1))
# is an "ascent" in the NATURAL ORDER. The sign is NOT about tournament edges
# but about vertex label comparisons.
#
# In the TOURNAMENT context, the free vertices are labeled, and for each
# ordered pair (u,v), the tournament determines the sign. The sign-sum
# formula replaces tournament-dependent signs with label-order-dependent signs.
#
# This is possible because the Fourier decomposition treats the cycle edges
# separately. The free edges' contribution, after summing over all tournaments,
# produces the sign-sum polynomial.

# PROOF IDEA: Consider the contribution of free vertex set V_free = {v_1,...,v_{f+1}}
# to the W-coefficient of a cycle invariant. For each pair (v_i, v_j),
# the edge orientation s_{v_i,v_j} is independent of the cycle edges.
# When we sum over both orientations of each free edge:
#
# sum_{orientations} prod_{steps involving free edges} (r + s) =
# sum_{sigma in S_{f+1}} sum_{orientations consistent with sigma}
#   prod_{i=0}^{f-1} (r + s_{sigma(i),sigma(i+1)})
#
# For a FIXED permutation sigma of the free vertices, the edges
# (sigma(i), sigma(i+1)) are determined. Each such edge has a fixed
# orientation in the sum. But we're summing over ALL orientations of
# ALL free edges.
#
# Hmm, this isn't quite factoring correctly because different permutations
# use different edges, and edges shared between permutations are NOT
# independently summed.
#
# Let me try yet another approach: direct computation.

print("\n" + "="*70)
print("DIRECT TEST: Does the free-edge sum give F_f(r)?")
print("="*70)

# Take n=5. Fix a 3-cycle on vertices {0,1,2} (using backbone edges 0->1, 1->2
# and non-backbone edge 2->0 or 0->2).
# The "pinned" contribution from the cycle is: the path visits {0,1,2}
# in some order using 2 steps. The "free" vertices are {3,4} which use 2 more steps.
#
# Actually the decomposition is more subtle. Let me compute W(r) for ALL
# n=5 tournaments and extract the coefficient of s_{02} (the centered variable
# for edge (0,2)). This should relate to F_f somehow.

n = 5

def compute_W_poly(A, n):
    result = defaultdict(lambda: Fraction(0))
    for perm in permutations(range(n)):
        poly = {0: Fraction(1)}
        for i in range(n-1):
            s = Fraction(A[perm[i]][perm[i+1]]) - Fraction(1, 2)
            new_poly = {}
            for power, coeff in poly.items():
                new_poly[power+1] = new_poly.get(power+1, Fraction(0)) + coeff
                new_poly[power] = new_poly.get(power, Fraction(0)) + coeff * s
            poly = new_poly
        for power, coeff in poly.items():
            result[power] += coeff
    return dict(result)

def tournament_from_tiling(n, tiling_bits):
    A = [[0]*n for _ in range(n)]
    for i in range(n-1):
        A[i][i+1] = 1
    idx = 0
    for i in range(n):
        for j in range(i+2, n):
            if (tiling_bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def num_tiling_bits(n):
    return n*(n-1)//2 - (n-1)

m = num_tiling_bits(n)

# The tiling edges at n=5:
edges = []
for i in range(n):
    for j in range(i+2, n):
        edges.append((i,j))
print(f"  Non-backbone edges at n={n}: {edges}")
# edges: [(0,2), (0,3), (0,4), (1,3), (1,4), (2,4)]

# For each tiling bit pattern, compute W(r).
# Then: compute the "Fourier component" for a specific s-variable pair.
# The Fourier component for s_{02} * s_{12}:
# W = sum_bits W(bits) = sum over all 2^m tournaments of W_T(r)
# Each W_T(r) = sum_perm prod(r + s_e)
# s_e = A_e - 1/2 = (bit - 1/2) for non-backbone, = +1/2 for backbone.

# Fourier: W_T(r) = sum of r-polynomials * products of s-variables
# The coefficient of s_{02} (say, the first tiling variable) in W_T is:
# [partial W / partial s_{02}]
# which in practice means: sum over perms of [derivative of product w.r.t. s_{02}]

# Instead of derivatives, let me use the indicator:
# W_T = W_even(r) + s_{02} * W_odd(r) + higher order
# where W_even is the part with s_{02}^0 and W_odd is with s_{02}^1.

# Since s_{02} = A_{02} - 1/2, and A_{02} is a specific tiling bit:
# If bit=1: s_{02} = +1/2
# If bit=0: s_{02} = -1/2
# So: W(bit=1) - W(bit=0) = 2 * s_{02} * W_odd(r) (at s_{02} from +1/2 to -1/2)
# Actually: W(s=+1/2) - W(s=-1/2) = 2 * W_odd coefficient (since all other
# variables have definite values that may contribute)

# Hmm, this is getting complicated. Let me take the cleanest approach:
# compute (1/2) * [W(s_{02}=+1/2) - W(s_{02}=-1/2)] for ALL combinations
# of other s-variables.

# Actually, the simplest Fourier extraction is:
# coefficient of s_{e_1} * s_{e_2} * ... in W(r) =
# (1/2^m) * sum over all 2^m orientations of W_T(r) * prod sign(s_{e_k})

# This is the Walsh-Hadamard coefficient. For the monomial proportional
# to s_{02} alone (degree 1):
# coeff = (1/2^m) * sum_T W_T(r) * sign(s_{02}(T))
# where sign(s_{02}) = +1 if A_{02}=1, -1 if A_{02}=0.

# For a 3-cycle invariant, we need degree-2 Fourier components.
# But let me check: what IS the Fourier decomposition of t3?

# t3 = number of directed 3-cycles.
# t3 = sum over triples {a,b,c} of [a->b->c->a or a<-b<-c<-a]
# In s-variables: for triple {a,b,c}, the cycle indicator is:
# (1/2 + s_{ab})(1/2 + s_{bc})(1/2 + s_{ca}) + (1/2 - s_{ab})(1/2 - s_{bc})(1/2 - s_{ca})
# = 2 * [1/8 + s_{ab}*s_{bc}/2 + s_{bc}*s_{ca}/2 + s_{ca}*s_{ab}/2]... hmm

# Actually, for directed 3-cycle a->b->c->a:
# I(a->b->c->a) = A_{ab} * A_{bc} * A_{ca}
# = (1/2 + s_{ab})(1/2 + s_{bc})(1/2 + s_{ca})
# = 1/8 + (s_{ab} + s_{bc} + s_{ca})/4 + (s_{ab}*s_{bc} + s_{bc}*s_{ca} + s_{ca}*s_{ab})/2 + s_{ab}*s_{bc}*s_{ca}

# t3 counts BOTH directions for each cyclic triple:
# I(a->b->c->a) + I(a->c->b->a)
# = (1/2+s_{ab})(1/2+s_{bc})(1/2+s_{ca}) + (1/2-s_{ab})(1/2-s_{bc})(1/2-s_{ca})
# = 2*(1/8 + (s_{ab}*s_{bc} + s_{bc}*s_{ca} + s_{ca}*s_{ab})/2)
# = 1/4 + s_{ab}*s_{bc} + s_{bc}*s_{ca} + s_{ca}*s_{ab}

# So the degree-2 Fourier part of t3 involves sum of pairwise products
# of s-variables around each triangle.

# The W-coefficient of t3 at n=5 is 2*F_2(r) = 12r^2 - 1.
# This should equal the coefficient of each degree-2 monomial (from cycles)
# in W, multiplied appropriately.

# This is getting quite involved. Let me just verify the RESULT:
# At n=5, W(r) = F_4(r) + 2*F_2(r)*t3 + 2*F_0(r)*t5.
# I've already verified this holds for all 64 tournaments (earlier scripts).
# The universality is an ALGEBRAIC FACT.

# The PROOF of universality comes from the permutation formula:
# F_f(r) = sum_{sigma in S_{f+1}} prod(r + sgn(sigma(i+1)-sigma(i))/2)
# This is independent of vertex labels by construction.
# The question is only: why does THIS function appear as the W-coefficient?

# The answer (which I now understand) is:
# 1. W(r) = sum_P prod_{steps} (r + s_e)
# 2. Each edge e has s_e = A_e - 1/2
# 3. In the Fourier decomposition, term proportional to cycle monomial M involves:
#    - Cycle edges contributing their s_e values (fixed by M)
#    - Free edges contributing (r + s_e) where s_e is summed over ±1/2
# 4. The sum over permutations of free vertices with free edge signs:
#    sum_{sigma} prod_{free step i} (r + s_e)
#    where s_e = +1/2 if sigma(i+1) > sigma(i) and -1/2 otherwise
#    (by convention, since each free edge has both orientations)
# 5. This gives exactly F_f(r) by the sign-sum formula.

# BUT WAIT: step 4 says s_e depends on the LABEL ORDER, not the tournament.
# Why? Because we're summing over BOTH orientations of each free edge.
# For a free edge (u,v) with u < v in the natural ordering:
# - Orientation A_{uv} = 1: s_e = +1/2 (this happens when v > u in the path, an "ascent")
# - Orientation A_{uv} = 0: s_e = -1/2 (v < u in the path, a "descent")
# When we sum over both: the factor (r + 1/2) + (r - 1/2) = 2r appears...
# but that's for the sum, not the product.

# The product structure is more subtle. Each permutation sigma of f+1 vertices
# defines a sequence of edges. The product prod(r + s_e(sigma)) depends on
# which edges are ascents vs descents in sigma. Summing over all 2^{f free edges}
# orientations of the free edges, each permutation sigma's contribution is:
# prod_{i: ascent} (r + 1/2) * prod_{i: descent} (r - 1/2)
# because the edge orientation that makes step i an ascent in sigma is the
# ONLY one that "matches" (the other orientation would make it a descent).

# NO WAIT. The free edges are summed over both orientations, but different
# permutations share edges. So we can't sum each permutation independently.

# I think the correct argument is:
# In the TILING MODEL, each non-backbone edge has an independent bit.
# The Fourier decomposition extracts degree-k terms.
# For the degree-0 term (no cycle variables), we get the "background" C_0(r).
# C_0(r) involves summing W(r) over ALL 2^m edge orientations and dividing by 2^m.
# The result is F_{n-1}(r).

# For the degree-2 term proportional to a specific cycle's Fourier monomial,
# we FIX the cycle edges and sum over all free-edge orientations.
# The free edges contribute independently of the cycle, giving F_f(r).

# The key: when we FIX cycle edges and sum over free edges, the result
# factors as (cycle contribution) * (free contribution).
# The free contribution sums over all 2^{free_edges} orientations of free edges
# AND all orderings of free vertices in the path.
# This is equivalent to summing over all tournaments on the free vertices,
# which gives... ah but earlier that gave only r^f.

# Let me re-examine. The free edges are between free vertices and ALSO
# between free vertices and cycle vertices. Only the "internal" free edges
# (between free vertices) are independent of the cycle.

# So the free contribution involves:
# - Edges between free vertices (summed over both orientations)
# - Edges between free and cycle vertices (also summed? or fixed by cycle?)

# If the cycle invariant only pins the CYCLE EDGES (between cycle vertices),
# then edges between free and cycle vertices are also free.
# In the backbone model: backbone edges are always forward.
# Non-backbone edges between free and cycle vertices are tiling bits,
# independent of cycle bits.

# So the "free edges" include:
# - All non-backbone edges between free vertices: C(f+1,2) - f edges
# - All non-backbone edges between free and cycle vertices

# The total free contribution sums over all orientations of these edges.
# This should give F_f(r) * (some factor from free-cycle edges).

# Hmm, this is getting complicated. Let me just VERIFY the factorization.

print("\n" + "="*70)
print("TEST: F_f as sum over ALL tournaments on f+1 vertices in BACKBONE model")
print("="*70)

for nv in [2, 3, 4, 5]:
    f = nv - 1
    m_bb = nv - 1  # backbone edges
    m_nb = nv*(nv-1)//2 - m_bb  # non-backbone edges
    total_poly = defaultdict(lambda: Fraction(0))

    for bits in range(2**m_nb):
        A = [[0]*nv for _ in range(nv)]
        for i in range(nv-1):
            A[i][i+1] = 1  # backbone
        idx = 0
        for i in range(nv):
            for j in range(i+2, nv):
                if (bits >> idx) & 1:
                    A[i][j] = 1
                else:
                    A[j][i] = 1
                idx += 1

        for perm in permutations(range(nv)):
            poly = {0: Fraction(1)}
            for i in range(nv-1):
                s = Fraction(A[perm[i]][perm[i+1]]) - Fraction(1, 2)
                new = {}
                for pw, c in poly.items():
                    new[pw+1] = new.get(pw+1, Fraction(0)) + c
                    new[pw] = new.get(pw, Fraction(0)) + c * s
                poly = new
            for pw, c in poly.items():
                total_poly[pw] += c

    # Divide by 2^{m_nb} (average over non-backbone orientations)
    avg_poly = {p: c / 2**m_nb for p, c in total_poly.items()}

    Ff = master_poly_F(f)
    match = all(abs(float(avg_poly.get(p, Fraction(0)) - Ff.get(p, Fraction(0)))) < 1e-10
               for p in set(list(avg_poly.keys()) + list(Ff.keys())))

    print(f"  nv={nv} (f={f}): m_nb={m_nb}, avg over backbone model:")
    print(f"    avg W = {poly_str(avg_poly)}")
    print(f"    F_{f}  = {poly_str(Ff)}")
    print(f"    MATCH: {match}")
