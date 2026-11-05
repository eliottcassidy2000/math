# Induced ternary ancestry fails at 13 -> 41 -> 127

Session 194, 2026-11-05. **PROVED using an exact finite certificate and cited
strong-order/classification inputs.** No feedback-vertex restriction is
assumed. The unrestricted `{7,21}` Hamilton-path spectrum remains **OPEN**.

## 1. Statement and inheritance

**Theorem.** There do not exist nested induced tournaments

    A subset B subset C,   H(A)=13, H(B)=41, H(C)=127.

All three tournaments may have arbitrary orders and may be nonstrong.
This does not exclude any of the three scalar counts, an induced 41->127
extension whose H41 input does not contain H13, a construction changing an
old arc, or selecting a different H41 witness for the next step.

The earlier Session 193 induction packet excludes any induced H41 inside a
feedback-vertex-number-at-most-two H127 tournament. This theorem instead
retains the H13 ancestry and removes every feedback-vertex restriction.
The statements have different scopes; neither subsumes the other. Background
and the reviewed predecessor packets are routed from the published
[Session 193 guide](https://github.com/eliottcassidy2000/math/blob/169792fb5eeb483e351112b5a09706e6d7bd89f0/00-navigation/S193.md).

The inherited classification from Session 191's unit-start boundary is:
there is exactly one strong H13 tournament, up to isomorphism, represented
by outgoing bitmasks

    D=(2,12,9,17,7).

This classification was proved there by the complete five-name census and
Busch's order bound. We cite it rather than rerun that blanket census.

The external order bound used here is A. H. Busch, *A Note on the Number of
Hamiltonian Paths in Strong Tournaments*, EJC 13 (2006), N3,
https://doi.org/10.37236/1141:

    H(U)>=5^((|U|-1)/3), for strong U.

Its published proof is cited, not independently re-proved. The integer
comparisons 5^5>13^3, 5^7>41^3 and 5^10>127^3 imply strong-core orders
at most 5, 7 and 10, respectively.

## 2. Reduction to one exact finite universe

Strong components of a tournament are totally ordered and H is their
product. Every nontrivial strong component contains a directed triangle,
so by Hamilton-path insertion monotonicity its count is at least three.
A tournament of prime count therefore has one nontrivial strong component
with that prime count, with only singleton components before and after it.

Assume the displayed nested chain exists. Its H13 core is a copy of D.
Because D is strong and induced, all its vertices lie in one strong component
of B. That component cannot be a singleton and must be B's H41 component,
call it U. Similarly U lies in C's H127 component V. Hence

    D induced in U induced in V,
    |U|<=7,   |V|<=10,   H(U)=41, H(V)=127.

Fix the names of D as 0,...,4. Every U is obtained by at most two births.
Exhaust all 32 masks for vertex5. There is no one-birth H41. Exhaust all 64
masks for vertex6 after each parent with H<=41. There are exactly ten
labelled two-birth H41 records, all strong. These are extension records,
not asserted to be ten isomorphism classes.

For each record U, exhaust every extension by one, two and three vertices.
At each intermediate step retain every tournament with H<=127, whether
strong or nonstrong. At each next step try every outgoing incidence mask.
Hamilton-path insertion monotonicity makes this pruning rigorous: no
intermediate whose count exceeds the target can later return to the target
without changing old arcs. There are no symmetry, degree, apex, or port
filters.

The exhaustive certificate gives:

    first masks tested                  1,280
    first parents retained                272
    second masks tested                69,632
    second parents retained             1,220
    third masks tested                624,640
    H127 hits at each birth depth      0,0,0.

Thus no eligible V exists. This proves the all-order statement, including
arbitrary neutral padding. It does not assume the intermediate tournaments
remain strong. Indeed some useful intermediates below are nonstrong.

## 3. The independent counting algorithms

For a tournament U, let e_A(i),s_A(i) be the Hamilton-path endpoint counts
on the induced subset A. Let

    c_ij=sum_(A contains i, A excludes j) e_A(i) s_(V-A)(j).

A newborn z has outgoing set M among the old vertices. Splitting every
Hamilton word at z gives the exact formula

    H(U+z)=sum_(i in M) s_V(i)
           +sum_(i notin M) e_V(i)
           +sum_(i notin M,j in M)c_ij.                (1)

The first two terms put z first or last. The third term puts a nonempty
old path before z and a nonempty old path after it, with their boundary
vertices i and j. This is a bijection, not a fitted cut polynomial.

`ancestry_verify.cpp` computes (1) from a reverse subset recurrence. In
`--full` mode it independently recomputes **every candidate count, including
all rejected intermediate masks and all 624,640 final masks**, by forward
endpoint subset DP on the actual enlarged adjacency matrix. Exact integer
agreement is required. All counts here are at most 10! and fit int64_t.
The small source inventory is also independently checked before pruning.

The two recurrences have different control flow, and the birth formula adds
an independent decomposition. These controls do not replace the exhaustive
universe or the all-order reduction in section 2.

## 4. Positive physical controls and a surviving sibling branch

All row lists are outgoing bitmasks: bit j of row i equals one iff i->j.
A birth mask M is the newborn's outgoing set among the existing labels.
Old labels and every old incidence are retained.

Starting at D, add vertex5 with mask31 (it beats every old vertex). This
nonstrong intermediate still has H13. Add vertex6 with mask47 or55:

    H41: A=(2,12,9,17,71,31,47),       second mask47;
    H43: B=(2,12,9,81,7,31,55),       second mask55.

Both final tournaments are strong. These are the standard one-clock rooted
wrappers at roots4 and3, respectively. Direct rooted two-path cover counts
on D are (26,18,18,15,14), so H+2F gives 41 and43.

The lower descendant H125 is available while retaining A: add vertex7 with
mask31, then vertex8 with mask31, giving

    (2,12,9,17,71,415,431,287,31),     H=125.

The intermediate has H69 and the final tournament is strong. By contrast,
no induced descendant of this A has H127, by the theorem's finite search
and strong-core reduction.

The sibling branch survives without changing old arcs. From B, one birth
with mask10 gives

    (130,12,137,81,135,159,183,10),    H=131=3*43+2,

and one birth with mask29 gives

    (2,140,9,81,7,159,183,29),        H=133=3*43+4.

Both are strong. The alternative mask79 also gives H133, but is not needed
for this certificate. Thus the obstruction is not a blanket inability to
grow through two generations: the 13->43->{131,133} branch has an explicit
induced physical realization, whereas 13->41->127 does not.

A direct induced H13-to-H127 route also exists if the intermediate H41 is
bypassed. Starting from D, use masks29,7,1. The successive strong counts are

    13 ->29 ->75 ->127,

with final rows (2,172,137,209,199,221,135,1). This is a relabeling of the
S193 H127 positive control, independently checked here. Thus the theorem
isolates the intermediate41, rather than excluding all H13 ancestry of127.

No outgoing return operations from H131 or H133 are asserted. The reusable
induced-only ternary tree must fail, but reconstruction, old-arc resets,
branch-specific products and existential target selection remain available.

## 5. Reproduction and scope

    g++ -O3 -std=c++17 -Wall -Wextra -pedantic ancestry_verify.cpp -o ancestry_verify
    ./ancestry_verify --full > ancestry_full.json
    python3 verify_ancestry.py
    python3 -O verify_ancestry.py

The Python verifier independently enumerates literal Hamilton words for the
five positive final controls and both source branches, checks every birth
mask, strongness, exact integer cutoffs, and the ten seven-name source
counts. The C++ full receipt retains every source mask and each source's
search totals. A `fast` run omits the last-layer independent DP only and is
not described as a full independent replay.

The global theorem depends on the cited strong H13 classification and Busch
bound. The finite negative result needs neither dependency. No new global
forbidden Hamilton-path count is claimed.
