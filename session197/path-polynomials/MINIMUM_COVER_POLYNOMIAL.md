# Minimum-cover factorization of the highest clone shift

Session 197. **PROVED polynomial strengthening of the minimum-cover scalar factorization.** This packet is separate from the previously frozen matching-amplitude theorem.

## 1. Data and the exact identity

Let Q be a tournament with distinguished root r. Put O=N^+(r), I=N^-(r), and let G be the bipartite graph whose edges are the arcs u->v with u in O and v in I. Let nu be its maximum matching size. Choose any minimum vertex cover C of G, so |C|=nu. Write

    C_O=C intersect O, C_I=C intersect I,
    O'=O minus C, I'=I minus C.

Every arc between I' and O' points from I' to O', since C covers every O-to-I edge.

The full labelled-bin context operator from `UNIVERSAL_CONTEXTS.md` has a finite shift expression

    A_Q f(x)=sum_j p_j(x) f(x+j).

`CONSEQUENCES.md` proves that its greatest forward shift is nu. We factor that entire polynomial p_nu(x), not merely its value at x=1.

Define a **left arrangement** to be a collection of pairwise vertex-disjoint directed paths P_c in Q[O'], one labelled path for each c in C, satisfying:

- If c is in C_I, P_c must be nonempty.
- If c is in C_O, P_c may be empty.
- A nonempty P_c ends at a vertex pointing to c.

Let U_L be the unused vertices of O'. Put

    L_C(x)=sum_left_arrangements B_(Q[U_L])(x).            (1)

Define a **right arrangement** similarly using paths S_c in Q[I']:

- If c is in C_O, S_c must be nonempty.
- If c is in C_I, S_c may be empty.
- A nonempty S_c starts at a vertex to which c points.

With U_R the unused vertices of I', put

    R_C(x)=sum_right_arrangements B_(Q[U_R])(x).           (2)

Every directed path is an ordered vertex word. Two different path orders are counted separately. Empty paths consume no vertices and have one choice. Empty residual graphs have B_empty(x)=1.

Then

    p_nu(x)=x(x+1)...(x+nu-1) L_C(x) R_C(x).              (3)

The empty product at nu=0 is 1. Both L_C and R_C have nonnegative integer coefficients in ordinary powers of x. The factorization holds for every minimum cover, though the two factors may change with the cover. No assertion is made that arbitrary positive polynomial pairs are physically realizable.

## 2. Why maximum-shift skeletons use each cover vertex once

In the exact odd-cycle context expansion, a crossing cycle visits the substituted module in one or more nonempty runs. Each exterior fragment between consecutive runs is a directed path starting in O and ending in I. If there are nu runs in total, there are nu disjoint exterior fragments.

Each such fragment contains an O-to-I arc: along its vertex sequence, consider the first transition from the O side to the I side. That arc has an endpoint in C. Since the fragments are disjoint and |C|=nu, all nu cover vertices are used, exactly one per fragment. No cover vertex is left for a pure exterior cycle or as an unused vertex.

For the fragment using c, every vertex before c lies in O'. Otherwise its first O-to-I transition would already use a cover vertex different from c. Similarly every vertex after c lies in I': an O vertex after c would require a later O-to-I transition using another cover vertex. Thus the fragment is precisely

    P_c, c, S_c.

Its first vertex must be in O and its last in I. This gives exactly the mandatory/nonmandatory conditions in (1) and (2). Conversely any two arrangements in (1) and (2) make nu disjoint O-to-I fragments of this form, each containing its distinct cover vertex.

This is a labelled-fragment bijection. It does not discard endpoint incidence conditions, introduce a virtual graph or assume positivity after a general change of coordinates.

## 3. Permutations, crossing-cycle parity and the rising factorial

Once the nu fragments are fixed, their organization into crossing cycles is encoded by a permutation sigma of their cover labels: each cycle of sigma gives one cyclically ordered crossing cycle, with one module run between adjacent exterior fragments. The cycle notation already identifies cyclic rotations and treats different crossing cycles as unordered, so there is no additional factorial quotient.

The odd-cycle context expansion has two important features established in `UNIVERSAL_CONTEXTS.md`:

1. A crossing cycle with k module runs and exterior parity epsilon contributes the shift factor

       (S-1)^k + (-1)^epsilon (S^(-1)-1)^k,

   where S f(x)=f(x+1). Its coefficient of the extreme forward shift S^k is exactly 1, for every k>=1 and every exterior parity. The factor 2 for an odd cycle and the factor 1/2 from parity projection have already cancelled in this formula.

2. Its exterior weight includes one factor x for each crossing cycle. Pure exterior odd cycles and unused exterior vertices contribute the ordinary odd-cycle expansion of B on the residual exterior graph.

Therefore every sigma contributes x^(number of cycles of sigma) to p_nu(x), regardless of the parity of its exterior fragment lengths. This does not assert that every unexpanded clone assignment is an odd cycle. It extracts the extreme shift coefficient from the exact parity-projected operator.

The sum over all permutations is

    sum_sigma x^(number of cycles of sigma)
       = x(x+1)...(x+nu-1).

One elementary proof inserts the last label: it either forms its own cycle, contributing x, or is inserted after one of the previous nu-1 labels, contributing nu-1. This proves the product recursively, including the empty-permutation value 1.

All remaining vertices lie in U_L or U_R. A pure exterior directed cycle cannot meet both sets, because every arc between I' and O' points from I' to O'. Its cycle collection therefore splits independently into the two sides. Summing their odd-cycle weights gives B_(Q[U_L])(x)B_(Q[U_R])(x). The left and right arrangements are independent as well. Multiplying these contributions proves (3).

## 4. Zero rank and scalar interpretation

When nu=0, C is empty and G has no edges, so I->O. There is one empty path arrangement per side, giving L_C=B_(Q[O]) and R_C=B_(Q[I]). The actual substituted tournament is the ordered join Q[I]->T->Q[O], and hence

    A_Q f=B_(Q[I]) B_(Q[O]) f.

Thus (3) holds directly even with empty O or I, and no crossing-cycle argument is needed at this boundary.

For any nu, evaluating (1) at x=1 counts the labelled prefixes together with one directed Hamilton path on their unused O' vertices, allowed to be empty. Likewise R_C(1) counts the labelled suffixes together with one initial path on unused I'. Therefore

    a_(nu+1)=p_nu(1)=nu! L_C(1)R_C(1),                   (4)

where a_(nu+1) is the maximum-run coefficient. This is the exact scalar minimum-cover factorization. The polynomial identity strengthens it by keeping all residual labelled-bin information.

## 5. Degree, matching multiplicity and inherited lower bound

Each of the |C_I| mandatory left paths consumes at least one O' vertex. Therefore

    deg L_C <= |O'|-|C_I|.

The top coefficient counts the choices in which those paths are singletons and every optional path is empty. These are precisely the matchings from distinct vertices of O' into all vertices of C_I. Let this number be ell. Similarly the top coefficient of R_C is the number rho of matchings from all C_O into distinct I', at degree |I'|-|C_O|.

A maximum matching of G has nu edges, each touching C. Because |C|=nu and matching edges are disjoint, every such edge touches exactly one cover vertex, and every cover vertex is used. In particular no edge between C_O and C_I can occur in a maximum matching. Thus a maximum matching decomposes uniquely into the two independent matchings counted by ell and rho. Their product is the number mu of maximum matchings:

    mu=ell*rho>0.

Both degree bounds are consequently attained. If N=|O|+|I|, identity (3) now yields

    deg p_nu=N-nu, leading coefficient mu,
    p_nu(x) >= mu*x^(N-2nu)*x(x+1)...(x+nu-1)

coefficientwise. If N=2nu, both side polynomials are constant and equality holds. These recover the previously audited matching-amplitude results, while (3) retains the extra exact factors.

Nonnegative ordinary coefficients in (1) and (2) follow from the physical odd-cycle expansion of every residual B polynomial. This is a particular physical factorization; no general real-rootedness, total-positivity or factor-realization conclusion follows.

## 6. Reproduction and dependencies

Run `python3 verify_minimum_cover_polynomial.py` and `python3 -O verify_minimum_cover_polynomial.py`. The checker enumerates directed side paths and unused residual bin polynomials directly. For every minimum cover it compares (3) with the previously frozen literal odd-cycle skeleton operator, checking all rooted labelled contexts through order four, every balanced rank-at-most-two crossing pattern, and a rank-three control. It also checks the degree, leading coefficient and scalar amplitude identities. The finite controls supplement the all-order bijection above.

Dependencies:

- `UNIVERSAL_CONTEXTS.md`, SHA256 `8615845e7bf06a6cdb1925fdbd8198a77b1a8864aa9c22704b0459678ff2ed51`.
- `CONSEQUENCES.md`, SHA256 `5616a1d747edad5828f48fb36efb615cdb01d60eaf770bef9e00029d81b74e1d`.
- `MATCHING_AMPLITUDE.md`, SHA256 `443097af71c4d1c364c4bc8ad5d417453df6aea71d734ccb5700b74b93bf48ee`.
- `verify_matching_amplitude.py`, SHA256 `e062734acf04e628a5fee3ed15d57123a22bcb6f7a2f1058edd1972f0566773c`, and its explicit named dependency chain.

Existence of a cover of size nu is the bipartite matching/cover theorem. Constructively, start from any maximum matching, follow alternating paths from unmatched O vertices, and let the cover be the unreached O vertices together with the reached I vertices. Standard alternating-path counting gives size nu and covers every edge. The identity itself holds for any cover of that size, so it is not tied to one matching or one alternating-path choice.
