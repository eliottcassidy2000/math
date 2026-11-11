# One-way physical trades: exact additive digits without pair interactions

**Status: PROVED in this note from the inherited fusion law; literal input data
and replay are FINITE-EXACT. Chan–Pak remains OPEN.** No claim of historical
priority for the elementary triangular-matrix or subset-sum facts.

## 1. Inheritance and the question changed by the experiment

The operation is Session264 `fusion/PROOF.md`: identify roots of A and B and
orient every nonroot A vertex towards every nonroot B vertex. The carrier is

    M(T)=[h,U-h,d-h,H-h; 0,Q,0,V-h; 0,0,k,e-h; 0,0,0,h].

Here h=H(T-root), Q is its empty-slot two-cover count, U,V are weighted masses
of paths starting/ending at the root, and d,e their odd-minus-even masses.
For a nontrivial pointed tournament k=0; the rooted singleton identity has k=1.
Matrices multiply under this physical fusion. This law is inherited, not
re-proved or relabelled here.

Earlier equal-profile switches altered both U and V, so separately chosen
switches interact. Their positive ordering and even-difference gcd do not give
independent digits. We first looked for an H-preserving one-way trade and found
one, then changed the question: why require its scalar H to remain fixed?
Dropping that restriction produces a primitive, coprime-diagonal controller.
This implements the user's proposed move from invariant operations to
non-invariant operations that compose in a controlled way.

## 2. A square-zero condition that survives every physical spacer

**Theorem 1 (one-way affine controller).** Let each choice at slot i be between
actual pointed tournaments A_i and B_i with the same h,Q,k,V,e. No equality of
H,U,d is required. Insert arbitrary fixed pointed fusion contexts between the
slots and at either end. Then the final Hamiltonian count is exactly

    H(bits)=H(all B)+sum_i w_i bit_i,

where w_i is the count difference obtained by changing ONLY slot i from B_i to
A_i. Every mixed discrete derivative of order at least two vanishes. The same
statement holds simultaneously for the whole eight-entry fusion carrier.

**Proof.** D_i=M(A_i)-M(B_i) has nonzero entries only in positions12,13,14.
For every upper triangular matrix C and every j,

    D_i C D_j=0,

because D_j has only row1 while (D_i C) has zero column1. Expand the ordered
matrix product, keeping every spacer in its actual place. Every term with two
or more D factors vanishes, regardless of how many baseline matrices/spacers
separate them. The remaining terms are exactly the baseline and single-slot
variations. Scalar H is the sum of entries11 and14. All choices remain actual
positive constructions; this proof implements no matrix subtraction physically.
The singleton/empty word is the identity matrix and is covered. QED.

This is stronger than a vanishing second derivative at one surrounding word:
the zero is an algebraic identity against *every* upper triangular spacer.
The binary controls are reusable at arbitrary length; no finite fuel is used.
Their amplitudes and absolute baseline are nevertheless not freely chosen.

## 3. Two literal six-vertex trades

Root0 and outgoing-neighbor masks specify all arcs. Direct subset path counting
gives the following pairs (see `replay.py`, `RESULTS.json`).

### 3a. Neutral trade, with common-prime debt

    X=(46,8,2,4,15,30), Y=(46,20,0,6,13,30).
    (H,h,Q,V,d,e)=(9,3,48,9,0,3), U_X=36, U_Y=42.

Thus M(Y)-M(X)=6E12. This already gives an exact additive controller in arbitrary
fusion contexts, despite identical scalar inputs. But every entry of both
matrices is divisible3, so every length-m pure word has H divisible3^m.
It cannot approach the missing prime targets. This was a positive result whose
failure indicated exactly which constraint to relax.

### 3b. Non-neutral trade, without common diagonal prime

    A=(6,16,2,7,13,31), B=(4,37,0,7,15,29).
    Common (h,Q,k,V,e)=(5,64,0,32,0).
    (H,U,d)_A=(9,16,2), (H,U,d)_B=(5,10,0).

Consequently

    M(A)-M(B)=6E12+2E13+4E14.

The h and Q diagonals are coprime. This is not an equal-profile trade; the
scalar inputs5 and9 deliberately differ. Four direct compositions and four
literal spacer cubes verify the zero mixed derivative without using the matrix
law to count their Hamiltonian paths.

## 4. All-length explicit selector and exact inverse

Let T(bits) fuse m copies chosen from A,B in order, bit1 meaning A. It has
1+5m vertices. Put b_m=H(B fused m times), computable by the displayed matrix,
and for r>=0 define

    g_0=4,
    g_r=(162*64^r-44*5^r)/59  (r>=1).

**Theorem 2.** For every m>=1,

    H(T(bits))=b_m+sum_(i=0)^(m-1) bit_i 5^i g_(m-1-i).       (1)

All coefficients are positive and, in their displayed order, strictly
superincreasing towards the left. Therefore exactly 2^m different scalar
counts occur. Given a target N, subtract b_m and greedily test these descending
weights; the target is realized by this family iff the final remainder is0.
When it is, the decoder produces the literal fusion word and hence a witness.

**Proof.** A prefix of length i has first diagonal5^i. Only the first row of
D is nonzero, so no other prefix information enters its contribution. For
r>=1 the baseline suffix has

    (M(B)^r)_24=27(64^r-5^r)/59,
    (M(B)^r)_34=-5^r,
    (M(B)^r)_44=5^r.

Pairing its last column with (6,2,4) gives g_r. For r=0 the identity gives4,
not the value2 obtained by incorrectly extending the r>=1 expression. Theorem1
proves (1). We retain that zero-length boundary explicitly.

For r>=2, g_r>64g_(r-1); and g_1/(5g_0)=43/5.
Thus each descending weight exceeds (43/5) times the next. The sum of all
smaller weights is less than one seventh of the next larger weight, in
particular less than that weight. Hence the first bit on which two words
differ determines the sign of their count difference. Greedy inversion follows:
if the remainder is at least the current weight, that bit must be1; otherwise0.
QED.

Initial cases, directly counted where indicated:

    m=1: baseline5, weights(4), counts5,9.
    m=2: baseline185, weights(172,20), counts185,205,357,377.
    m=3: baseline10365, weights(11228,860,100).

This is genuine exact count selection, but for a sparse explicitly recognized
language. Its output partition polynomial is exactly

    sum_bits z^H(T(bits)) = z^b_m product_i (1+z^(5^i g_(m-1-i))).

It encodes the whole scalar fiber distribution, not a heuristic moment or a
finite fit. Every coefficient is0 or1 by injectivity. These formulas remain
valid at arbitrarily large orders without enumerating tournament classes.

## 5. Why even a primitive independent controller does not yet cover an interval

Every g_r is divisible4, and g_r is coprime to5 for r>=1. The first assertion
follows from g_1=172 and g_r=64g_(r-1)+44*5^(r-1) for r>=2.
The second follows from the closed formula modulo5. Since the last weight is
4*5^(m-1), the gcd of all weights is exactly4 for EVERY m. There is no growing
common prime debt, unlike the neutral controller, but there are enormous gaps
inside its allowed congruence class. This is a structural gap proof, not merely
an unsuccessful bounded experiment.

**General interval criterion, CONDITIONAL on physical digit supply.** In any
one-way controller suppose slot i can choose every multiple0,w_i,...,d_i w_i
of a positive response w_i, independently, and the fixed baseline is B.
For a desired lattice spacing s, sort the weights increasingly. If

    w_1=s,  and  w_i <= s+sum_(j<i) d_j w_j,

with every w_i divisible s, then its actual output set contains every point
of the interval B+s{0,...,sum_i d_i w_i/s}. Proof: induction joins the translated
previous intervals, whose neighboring starts differ by at most their length+s.
Choosing the largest admissible digit recursively supplies an inverse witness.
This elementary complete-sequence criterion is not a newly proved numeration
theorem. The contribution here is its exact lawful consumer in Theorem1.

For Chan–Pak one still needs an s=2 physical digit bank, the displayed overlap
inequalities, and baselines whose intervals cover all sufficiently large odd
integers (plus the remaining initial values). None is supplied by gcd alone.
Our binary bank fails the overlap condition exponentially. Its failure names
a resource requirement: digit menu width must compete with transport expansion,
not merely remain nonzero at every generation.

## 6. An exact boundary against the easiest repair

**Theorem 3 (no nontrivial one-way bank with unit deletion count).** At fixed
order, two pointed tournaments with h=1 and the same V are rooted-isomorphic.
Consequently their entire fusion carriers agree. One cannot remove all
baseline amplification by seeking a nontrivial h=1 one-way trade.

**Proof.** H(T-root)=1 implies the deletion is transitive (the inherited
H=1 characterization). Order its vertices uniquely as v1->...->vn. Every
nontrivial path ending at the root has a unique last old vertex v_i which
points to the root. Its earlier old vertices form an arbitrary subset of
v1,...,v_(i-1), visited in order. Every complementary induced deletion has
H=1. Therefore

    V=1+sum_(i: v_i->root) 2^(i-1).

The unique binary expansion of V-1 recovers every root incidence. Thus it
recovers the entire pointed tournament. For an unspecified order the additional
coordinate Q=2^n recovers n, so the same conclusion holds for equal h,Q,V.
The empty deletion gives singleton V=1,Q=1 separately. QED.

The boundary does not prohibit interacting unit-deletion constructions,
multiple marks, re-rooting, changing the operation, or one-way trades with
h>1 (as exhibited). It rules out precisely the tempting no-amplification repair.

## 7. Connections that generated and constrained this route

* Session131 `zeckendorf/PROOFS.md` retains integral carries after reflection
  quotienting. Here the carrier retains rooted parity and *kills* mixed carries
  by a physical one-way condition. A numeration analogy becomes an exact
  subset-sum selector only after Theorem1. The systems are not isomorphic.
* Seeded random archive sampling (seed27521 over sorted tracked reflection
  paths) selected `berggren-positive-cubes-and-amm-causal-horizons-kps-s181.md`.
  Its methodological lesson is to move the frozen coordinate responsible for
  failure. Here freeing H converts a3-divisible neutral letter into the
  coprime-diagonal5/9 pair. No Pell or AMM theorem is imported.
* The same draw selected `fugacity-axis-and-vanishing-theorem.md`, whose
  alternating-sum theme suggested asking for more than one vanishing top
  Walsh coefficient. Theorem1 supplies ALL Walsh interactions of degree>=2
  zero, in every spacer context. This conclusion is proved by our matrix
  identity, independent of any historical reflection's correctness.
* Session216 `interval/PROOF.md` supplies actual strong substitution contexts
  H(G_d[U])=36Q_U(2)+dH(U), d=-7,-5,...,13. These cannot be substituted into
  (1) while ignoring Q_U(2): for our single letters Q_A(2)=168 and Q_B(2)=128,
  though root-deleted Q is64 for both. Their end-path length polynomials also
  differ despite equal V,e. A new richer one-way carrier or new exact response
  law is needed to combine those eleven-offset decoders with this controller.

This is the precise missing compatibility object, not a vague demand for more
numerical coverage. A promising next step is to find a one-way module bank
whose full two-cover response is controlled and whose menu fills the overlap
inequality, or to rotate interfaces between blocks without restoring mixed
interactions. Neither is asserted here.

## Reproduction and scope

Run `python3 -B 05-knowledge/results/session275_exact_selection_20261111/switches/replay.py`.
`carrier.py` is the relevant helper code copied from Session264 with attribution.
Four pure pairs and four spacer cubes are counted as literal tournaments;
levels1..8 exercise the exact formula/inverse;254 transitive-deletion rooted
words check the binary boundary. The all-order statements rest on the proofs,
not those finite tests. Small search scripts are exploratory discovery tools,
not full-class classification claims or dependencies of the stated theorems.
