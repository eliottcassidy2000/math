# A divisor-breaking joint-entry adapter, with exact dyadic spacing

**PROVED, scoped all-order construction and obstruction.** These results
concern literal tournaments and their named endpoint responses. They do not
establish reusable `3H+2 / 3H+4` return closure or any new global Hamilton
spectrum classification. The packet has independent algorithms, but awaits
the parent session's separate review before publication.

## 1. Inheritance and the changed interface

The incoming [joint-entry packet](https://github.com/eliottcassidy2000/math/blob/8f14b96ab020f7d5ceeb88e7bda413bb6f8bc9b7/05-knowledge/results/joint_entry_return_20261104/ENTRY.md)
proves that a same-bin joint count `J_U(v,r)=1` is the exact old-port
unit-start reentry condition for the three-birth dial `N^3`. Its C3-based
all-size adapter restores three ports, but every wrapper at those ports
has Hamilton count divisible by three. The [48-field analysis](https://github.com/eliottcassidy2000/math/blob/8f14b96ab020f7d5ceeb88e7bda413bb6f8bc9b7/05-knowledge/results/joint_entry_return_20261104/RESPONSE.md)
states that obstruction and asks for a divisor-breaking change of interface.

We replace its C3 core by the incoming four-vertex two-start atom. This
produces an all-size adapter whose very same literal component supports
wrappers in both nonzero residues modulo three. However, an exact spacing
theorem shows that *no two such wrappers of a fixed component*, even at
different prefix clocks, differ by two. Thus breaking the divisor is a
real advance in the interface, not the missing sibling return.

The closest mechanism is the incoming unit-joint transfer. The canonical
hostile is its all-wrapper divisor-three obstruction. The corrected near
miss is the inference from two residue lanes to the two exact desired
children. The least-used sidecar is the four-vertex atom's three different
joint ports. The live concepts are physical names, joint entry, marked
defect, modular observation, dyadic spacing, and target existence.

## 2. General core transfer

For a tournament U and a named r, write

    h=H(U), f=E_r(U,2),

where `E_r(U,2)` counts ordered two-bin covers, the first bin nonempty and
ending at r, the second possibly empty. Let U be nontrivial strong.
For m>=1 create D_m(U,r) by adding TT_m=Z and a root w, with

    w -> Z -> U,       r -> w -> U\{r}.

All arcs of U remain. The inherited zero-ready wrapper formula gives

    H(D_m)=h+2^m f.                                      (1)

Its retained-root data are

    s_w(D_m)=h,  H(D_m-w)=h,  E_w(D_m,2)=h+2^m f.         (2)

For the first two identities, a path starting at w must traverse Z before
U, and deleting w leaves TT_m->U. For the last identity, w has sole
in-neighbor r. The cut-after-w injection gives
`E_w(D_m,2)-H(D_m)=s_r(D_m-w)=0`, since the nonempty TT_m precedes r's
component. This is the converse of the inherited sole-exit failed-join
identity. Equivalently it follows by direct bin partition in the wrapper.

For v in U, every simple v-to-w path is a simple v-to-r path inside U,
followed by r->w. It cannot visit Z, which dominates U. The complement of
this path is TT_m followed by its complement in U. Therefore

    J_(D_m)(v,w)=J_U(v,r).                               (3)

Apply the literal three-birth dial: add A beaten only by w, then B beaten
only by A, then C beaten only by B. The newborns induce A->B->C->A.
Call the final tournament X_m(U,r). The incoming dial identity yields

    H(X_m)=3H(D_m)+3E_w(D_m,2)-s_w(D_m)
          =5h+6*2^m f,                                  (4)
    s_v(X_m)=J_U(v,r)  for v in U.                        (5)

These are all-order identities on the literal named graphs. Each D_m is
strong: w reaches Z and U, while every U vertex reaches r and hence w;
Z reaches U. Each dial birth preserves strongness, since the current
strong tournament has at least two vertices, the newborn is reached
through its designated parent, and it reaches all other old vertices.

For any core with h not divisible by three and at least one J=1 port,
equations (4)-(5) produce an all-size unit-start adapter whose *component*
Hamilton count is not divisible by three. This does not yet assert anything
about the counts of its port wrappers, which require the separate E2 data.

## 3. The explicit four-vertex core and its complete response

Take U on names a=0,b=1,c=2,r=3, with outgoing bit rows

    (2,4,9,3).

Thus a->b->c->a, r->a,b, and c->r. The inherited JE4 atom has

    H(U)=5, E_r(U,2)=6, J_U(a,r)=J_U(b,r)=J_U(c,r)=1.

The three joint equalities also follow directly: each port's only simple
path to r follows C3 to c and then r, and its complement is transitive.
Use the construction above. It has m+8 vertices and

    H_m=25+36*2^m,                                      (6)
    s_a=s_b=s_c=1.

The canonical physical roles are `a,b,c,r,w,Z_1,...,Z_m,A,B,C`, with w=4.
`count.cpp`, `verify.py`, and `response.json` retain every actual adjacency
matrix and the physical index map. In particular the Hamilton count at
m=1 is 97; this does not identify this graph with any older graph having
the same scalar count.

Let F_v=E_v(X_m,2). The three restored root responses are

    F_a=190*3^m+288*2^m-45,
    F_b=190*3^m+224*2^m-27,
    F_c=114*3^m+160*2^m-45.                             (7)

The first m=1 instance has `(F_a,F_b,F_c)=(1101,991,617)`.
Deleting those marks gives respectively

    q_a=15+18*2^m, q_b=5+12*2^m, q_c=5+6*2^m.            (8)

Thus every restored port has nontransitive deletion. The incoming JE3
obstruction applies immediately: another uninterrupted N^3 application
from any restored unit-start port has no unit-start output. This family
does not evade that structural no-two-ready-returns theorem.

For each pair `(v,w)`, v in {a,b,c}, the packet retains all 16 coordinates

    H, P(2), E_v(1), S_v(1), q_v, E_v(2), S_v(2),
    E_w(1), S_w(1), q_w, E_w(2), S_w(2),
    W2_vv, W2_vw, W2_wv, W2_ww.

W2 uses the end of bin0 and start of bin1, not two independent marginals.
The repeated-name coordinates are retained zeros. Every one of the 48
coordinates has an exact integer expansion in `1,2^m,3^m`, recorded in
`response.json`. No coordinate at bases 4 through 10 was discarded.

### All-order coefficient proof, not recurrence fitting

The transitive Z module has a fixed eight-name exterior. Contract every
maximal Z run in a one- or two-bin path cover. Choose an exterior
permutation and, for two bins, a cut. A gap is available exactly when its
adjacent exterior vertices have the required edges to/from Z. Every
backward exterior adjacency within one bin forces its gap to be nonempty.
An endpoint constraint on an exterior name forbids its adjacent boundary
gap. For k available gaps and q required nonempty gaps, the exact count is

    sum_(j=0)^q (-1)^j binom(q,j) (k-j)^m.               (9)

Indeed each of the m ordered Z names chooses a gap; its order inside that
gap is then forced. Inclusion-exclusion forces required occupation.
There are at most eight exterior names plus two bins, hence at most ten
gaps. `verify.py` directly enumerates all 8! exterior permutations, all
nine cuts, and the smaller deleted-name skeletons. It produces every
integer coefficient, including the zero-base boundary, without using
the DP samples to fit a recurrence. Formula (9) proves the results for
all m. The zero-base coefficient explains the m=0 exception: actual
H(X_0)=57, not the value 61 obtained by misapplying (6).

## 4. Both nonzero residues, but no sibling pair

The zero-ready wrapper at any restored port v and prefix clock k>=1 has

    Y_(m,k,v)=H_m+2^k F_v.                              (10)

For m>=1, H_m=1 mod3. Equation (7) gives

    F_a=0, F_b=2*2^m, F_c=2^m  (mod3).

Consequently Y at a is always 1 modulo3; at b and c the residues are 0
and 2 in one order, with the order determined by the parity of m+k.
The same actual component therefore supplies both nonzero residues at
every prefix clock. This removes the incoming C3 adapter's divisor
obstruction rather than merely exchanging its scalar source for another
multiple of three.

**Exact spacing theorem.** For every fixed m>=1, the minimum positive
difference between any two distinct values in

    {Y_(m,k,v): k>=1, v in {a,b,c}}

is exactly

    128*2^m-36 >=220.                                  (11)

In particular no pair, even with unequal clocks, can be `3H+2,3H+4` for
one common H.

Proof. Put M=2^m. The responses satisfy `0<F_c<F_b<F_a<2F_c`, with

    F_a-F_b=64M-18,
    F_b-F_c=76*3^m+64M+18,
    2F_c-F_a=38*3^m+32M-45.

The last two expressions are larger than `F_a-F_b`. For the last one,
subtracting gives `38*3^m-32*2^m-27`, which is 23 at m=1 and strictly
increasing thereafter. At equal clocks the smallest difference is
`2^k(F_a-F_b)`, minimized at k=1. At different clocks, writing the
smaller clock as k, the higher scale is at least `2^(k+1)F_c` and the
lower at most `2^k F_a`; the difference is at least
`2^k(2F_c-F_a)`, strictly larger than the claimed minimum. This also
covers different clocks at the same port. Equality is attained at
ports a,b with k=1. The common H_m cancels. QED.

The residue result is therefore deliberately not called two-branch
closure. Changing the core, changing m between the two branches, or
editing further incidences is outside this fixed-component obstruction.

## 5. The exact observer lattice and physical bulk operation

Write the 48-vector as `y=C(u1,u2,u3)^T`, with the realized family
`(u1,u2,u3)=(1,2^m,3^m)`. Its integer row lattice is exactly

    Z x 2Z x 38Z.                                      (12)

Every row satisfies the indicated column divisibilities. Conversely let
S=S_a(1), E=E_c(1), and use the three F values from (7). Then

    S=u1,
    D=11E+F_b-F_a-73S=2u2,
    2F_c-F_a+45S-16D=38u3.                            (13)

These explicit row identities prove (12), including zero coordinates.
The incoming C3 adapter had second-coordinate divisor six; this one has
divisor two. Its modulo-three loss of that direction is genuinely removed.
No novelty is claimed for Smith-normal-form or observation-lattice methods.

To pass literally from X_m to X_(m+1), birth a vertex beating precisely
the four old core names {a,b,c,r} and losing to every other old name.
It is the last vertex of Z. Every old name and arc survives, and both
endpoints are strong. The canonical indices of A,B,C shift; the receipt
contains the actual role-to-physical-index permutation, not a relabelling
permission. On the three linear modes the update is `diag(1,2,3)`.

The lattice in (12) is stable under that update. For arbitrary modular
boundary weights modulo N, the number of distinct instantaneous responses,
and hence of states distinguished by these observations and all future
bulk insertions, is exactly

    N [N/gcd(N,2)] [N/gcd(N,38)].                       (14)

The linear dimension over Q is three by (13). These minimality claims
are confined to the displayed linear weights, these 48 observations,
and bulk insertion. They say nothing about arbitrary complements,
nonlinear encodings, or reconstructing unrelated tournaments. On the
actual realized family, the integer m already determines the graph.

## 6. Verification, boundaries, and next useful question

Run, with no dependencies installed:

    g++ -std=c++17 -O3 count.cpp -o count
    for m in $(seq 0 10); do ./count "$m"; done > samples.jsonl
    python verify.py
    python -O verify.py

`response.json` retains all DP samples, all rows, all starting counts,
the fixed field order and all ten positive-base coefficient columns.
The DP producer restricts its range to m=0..10 (n<=18); every count
is at most `(n+1)n!<2^64`, so its unsigned 64-bit integers are exact.
`receipt.json` adds the independently derived complete symbolic table
including base zero. The independent algorithm checks all 528 sample
coordinates at m=0..10, every literal graph's strongness, and old-name
bulk transport. A third method enumerates all permutations and cuts at
orders eight and nine. Normal and optimized runs give identical stdout.
The m=20 vector is an evaluation of the proved symbolic identity, not an
additional independent direct-graph holdout and is labelled accordingly.

The [incoming nonzero-preparation reduction](https://github.com/eliottcassidy2000/math/blob/8f14b96ab020f7d5ceeb88e7bda413bb6f8bc9b7/05-knowledge/results/nonzero_preparation_20261103/algebra/PREPARATION_REDUCTION.md)
already separates actual target existence from native reachability.
Our construction does not select ready targets for the ternary descendants
of H1429, H170807 or the latest named frontier. No transport from one of
those saved graphs is claimed merely because a scalar or field phase
matches. The new result is a reusable *bulk entry adapter*, not a reusable
ternary digit grammar.

The next decisive task is an actual source rewrite selecting two targets
with difference exactly two, while preserving a future entry mechanism.
Within the present fixed X_m, neither changing the restored mark nor an
arbitrary dyadic prefix clock can do that. The full joint data, the changed
observer lattice and the exact spacing obstruction identify what a further
incidence edit must change.
