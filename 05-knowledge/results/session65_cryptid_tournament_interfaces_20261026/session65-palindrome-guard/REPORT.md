# Palindrome obstruction and a guarded tape carrier for native black K

2026-10-26. **Verdict: PROVED, with bounded independent controls.**
This is one native-interface gate for Session65's existing tape carrier.
No repository or previously sealed artifacts were edited; no commit or push
was made. The only inherited mathematical inputs are the two reports pinned
in `SOURCES.json`. The proofs below are independent elementary derivations
from their definitions. No external or historical search was performed.

## 1. Definitions and exact color theorem

Let G0 be the tournament on 0,1,2 with arcs 0->1, 1->2, 0->2;
let G1 reverse only 0->2, so G1 is a directed triangle C3. For a finite
binary word w define the order join

    E(w) = G_(w_0) join ... join G_(w_(m-1)),

all inter-block arcs directed from earlier to later blocks. A Hamiltonian
word Q=(q_0,...,q_(n-1)) is directed in its actual tournament. Its positional
reversal rho_Q sends q_i to q_(n-1-i). The word is **blue** exactly when rho_Q
is an anti-automorphism, meaning x->y iff rho_Q(y)->rho_Q(x); otherwise it is
**black**. A tournament is *pure blue* if every Hamiltonian word is blue.
The empty tape may be included with its one empty Hamiltonian word, which
is blue vacuously; no positive-order obstruction depends on that convention.

**Theorem.** For every binary word w, the following are equivalent:

1. w is a palindrome;
2. E(w) is self-converse, up to tournament isomorphism;
3. every Hamiltonian word of E(w) is blue;
4. at least one Hamiltonian word of E(w) is blue.

If w is not a palindrome, every Hamiltonian word of E(w) is black.

### Proof, explicitly handling zero-block mergers

The strong components of E(w), in their canonical source-to-sink order,
are one C3 for each symbol 1 and three consecutive singleton components for
each symbol 0. Adjacent zero blocks may merge into a larger transitive
interval, but they do **not** merge into a nontrivial strong component.
A run of r zeros becomes exactly 3r successive singleton components.
Thus the ordered SCC signature, using sizes 1 and 3, has a unique decoder:
each size-3 component becomes 1, and each maximal singleton run of length
3r becomes r zeros. This includes initial, terminal, and all-zero runs.

An anti-isomorphism reverses the order of strong components. Both singleton
and C3 are individually self-converse. Consequently E(w) can be
self-converse only if its SCC signature equals its reversal, and the decoder
then gives w=w^rev. Conversely, for a palindromic w the positional reversal
of the literal block-order word is already an explicit anti-automorphism.

More strongly, a directed Hamiltonian word must respect the SCC order:
once a directed path enters a later component, no arc permits it to return
to an earlier component. The singleton order is forced. Within each C3 the
three choices are precisely the three cyclic shifts of its directed
three-vertex word. Hence all Hamiltonian words are obtained by choosing
these shifts independently and concatenating; their number is 3^#ones(w).
This statement does not pretend that each zero block is a canonical maximal
module: it follows from the forced order of all its singleton vertices.

For palindromic w, rho_Q pairs corresponding three-vertex blocks of equal
type. On a G0 block, it reverses its unique linear order. On two C3 blocks,
it takes the first, middle, last vertices of one directed three-vertex word
to the last, middle, first vertices of the other; this reverses all three
arcs, for every choice of the two cyclic shifts. Between blocks it reverses
their order and therefore their arc direction. Thus every Q is blue.
A blue word always gives an anti-automorphism, so a nonpalindromic w cannot
have even one blue word. This proves all equivalences.

**Scope:** the equivalence between self-converse and pure blue is proved
for this carrier family. It is not asserted for arbitrary tournaments.

## 2. Own-word color is invariant across a legal native K edge

For n>=1 and a directed Hamiltonian word Q in T, the full-order native operation
K_Q retains the n-1 path arcs and reverses every other unordered pair.
The word Q remains directed in K_Q(T), and applying K_Q again returns T.

Let M_Q be the set of unordered non-path pairs. Positional reversal rho_Q
preserves adjacency of positions, hence preserves both the path-pair set
and M_Q. For any ordered pair (x,y), the two arc indicators compared in
the anti-automorphism condition

    T[x,y] = T[rho_Q(y),rho_Q(x)]

are therefore either both flipped by K_Q or both retained. Their equality
is unchanged. Thus

    Q blue in T  iff  Q blue in K_Q(T).                (1)

In particular Q is black at one endpoint iff it is black at the other.
This concerns the **same word on its own edge**; it does not claim that
arbitrary other words, or pure-blue status itself, are K-invariant.

**Corollary.** A pure-blue tournament has no incident black-K edge.
Indeed, if K_Q(S)=T and Q is black in S, (1) makes Q black in T, contradicting
pure blueness. The converse orientation also follows from involutivity.
This survives arbitrary relabelling: transporting T and Q by an isomorphism
conjugates rho_Q and preserves its anti-automorphism property, while K is
equivariant. Hence a pure-blue *isomorphism class* has no incident black
edge in the unmerged isomorphism quotient either.

Consequently no positive-length black-only K chain can reach such a target,
even with relabellings interspersed. If source and target are nonisomorphic,
a zero-length chain cannot identify them either. This is an endpoint-color
obstruction, not the fixed-label parity obstruction of `session65-k-write`.

## 3. A parity-compatible write still blocked even up to isomorphism

Consider the one-bit write

    E(0001) -> E(0000).

The order is n=12, and each full-order K reverses

    C_n = (n-1)(n-2)/2 = 55

unordered pairs. Thus the literal parity law does not prohibit an odd-length
chain for this one-arc target. Nevertheless the source has exactly three
Hamiltonian words, all black, whereas the target is the transitive
12-vertex tournament with its unique Hamiltonian word, which is blue.
The target's explicit word and anti-automorphism are

    Q=(0,1,2,3,4,5,6,7,8,9,10,11),   rho(i)=11-i.

The literal write reverses only {9,11}. The source has a C3 strong component
and the target has none, so they are nonisomorphic. By the corollary, **no
black-only native-K chain of any length realizes this write, even up to
isomorphism and even with arbitrary relabellings.** A single legal K using
the identity word changes 55 pairs, not one; parity compatibility is not
a construction. No impossibility is claimed here for unrestricted K chains
that are allowed to use blue words.

## 4. Guarding removes this endpoint defect, not the simulation gap

For fixed payload length m define

    F(w)=E(0 w 1).

The first and last symbols differ, so 0w1 is never palindromic. Therefore
every Hamiltonian word of F(w) is black, for every payload w. Since the
prefix and suffix are fixed common symbols, Session65's existing exact
fixed-length isomorphism-edit/Hamming formula gives

    d_iso(F(u),F(v)) = d_H(u,v).

No new metric theorem is needed. The guards are fixed tape context, not a
replacement for the separate head, control, origin, or program sidecars.
This removes the particular pure-blue endpoint obstruction. It does not
supply a legal K word taking F(u) to F(v), a sequence of such words, or a
proof that a sequence remains in any encoding family.

The guard alone also does not remove the fixed-label parity obstruction at
every order. If useful, append p additional leading zero cells, choosing
p from {0,1,2} so that L=m+2+p is 0 or 1 modulo 4. Then

    F_p(w)=E(0^(1+p) w 1),   n=3L,

still has differing endpoint symbols and all words black, and C_n is odd.
For fixed m the same p is used for all words, so the same Hamming isometry
holds. This bounded padding removes the **two stated necessary obstructions
only**. It proves neither native-K simulation nor global black connectivity.

## 5. Independently authored finite controls

`check.py` imports no producer or repository code. It uses adjacency matrices,
explicit native-K pair reversals, and explicit exception checks that remain
active under Python optimization. Product enumeration is exhaustive because
of the Hamiltonian-word decomposition proved above. It is **not** described
as an n! search at orders up to 18 (or the guarded orders up to 24).

For every binary tape of length m=0..6 it checks every product Hamiltonian
word for legality, positional-reversal color, native-K own-word color
invariance and involutivity. It independently computes SCCs by reachability,
decodes them, and checks the decoder after reversing all vertex labels.
Every tape's exact count, a representative word, and its reversal map are
retained in `normal.json`; for palindromes these are explicit blue witnesses.

| m | Tapes | Palindromes | Hamiltonian words | Blue | Black |
|---:|---:|---:|---:|---:|---:|
| 0 | 1 | 1 | 1 | 1 | 0 |
| 1 | 2 | 2 | 4 | 4 | 0 |
| 2 | 4 | 2 | 16 | 10 | 6 |
| 3 | 8 | 4 | 64 | 40 | 24 |
| 4 | 16 | 4 | 256 | 100 | 156 |
| 5 | 32 | 8 | 1,024 | 400 | 624 |
| 6 | 64 | 8 | 4,096 | 1,000 | 3,096 |
| **Total** | **127** | **29** | **5,461** | **1,555** | **3,906** |

The independent count identities are 4^m total words and
10^(m/2) blue words for even m, or 4*10^((m-1)/2) for odd m.
Adjacency-guided DFS independently reproduces the product lists on all
15 tapes through m=3 (85 words). All 2,893 candidate vertex permutations
for the seven tapes through m=2 are tested for anti-isomorphism, without
early stopping, agreeing with the palindrome criterion. These are finite
controls of the proof, not substitutes for its all-length argument.

The guarded family 0w1 is checked for every payload through length 6:
127 guarded tapes and **16,383** complete product Hamiltonian words, all
black. The stated bounded padding rule and odd-C arithmetic are also checked
for each such payload. Its isometry is inherited from the already proved
carrier metric formula, not independently re-proved by a large isomorphism
search here.

A separate small native-K control exhausts every labelled tournament and
every vertex word for n=1..5, filters actual directed Hamiltonian words,
and checks legality after K, K involutivity, own-word color invariance,
exact changed-pair count, a transported transposition relabelling, and K
equivariance. It explicitly checks that no black transition lands in any
pure-blue tournament in the exhaustive universe.

| n | Tournaments | Legal K pairs | Blue pairs | Black pairs | Pure-blue tournaments |
|---:|---:|---:|---:|---:|---:|
| 1 | 1 | 1 | 1 | 0 | 1 |
| 2 | 2 | 2 | 2 | 0 | 2 |
| 3 | 8 | 12 | 12 | 0 | 8 |
| 4 | 64 | 192 | 96 | 96 | 24 |
| 5 | 1,024 | 7,680 | 1,920 | 5,760 | 184 |
| **Total** | **1,099** | **7,887** | **2,031** | **5,856** | **219** |

There are 124,469 candidate tournament/vertex-word pairs before the legality
filter. The exact legal-pair identity n!*2^C_n is checked at each order.
The arbitrary-chain obstruction is proved from (1), not established by
searching chains at order 12.

Reproduce from this directory:

    python3 -B check.py > normal.json
    python3 -B -O check.py > optimized.json
    cmp normal.json optimized.json
    sha256sum -c SHA256SUMS

Both executions passed with byte-identical output. The seal pins the report,
source manifest, authored verifier, and both outputs. The bounded interface
gate is resolved; no new simulation or connectivity claim is made.
