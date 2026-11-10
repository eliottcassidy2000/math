# Session248: exact observer thresholds and recursive code banks

**Status: PROVED transfer statements and explicit sharp families below,
using Session247's marked-interface theorem as an inherited input.**
No claim of a new coding-theory theorem. The new object is a physically
realized tournament family whose aggregate Hamilton responses have an
exactly quantified interface threshold. Means are not distributions, and
ensemble symmetries are not vertex automorphisms.

## 1. Recover the right object from three old threads

The old `tournaments-as-codes.md` reflection is a historical source of ideas,
not a license to identify arbitrary tournament score classes with binary
linear codes. Its warning is retained. We put an actual binary linear code
on an actual cube of independent internal arc reversals instead.

Session170 supplied equal-path-cover-profile module switches; Session247
showed that q exceptional vertices can expose at most 2q-way interactions
among their bits. The old operation-congruence reflection distinguishes
closed averaged observables from a quotient of the full state transition.
THM-4524's loop-Walsh degree bound gives another physical precedent for
alternating-sum cancellation, but its switching gauge and degree bound are
not the module interface used here.

The bridge is:

    physical invisible module cube -> bounded-degree response
       -> exact linear-code coset average -> interface complexity.

The code lives on module choices, not on an alleged linear score syndrome.

## 2. Module-bit response and its Fourier support

Fix k labelled, disjoint module pairs M_i^0,M_i^1 on identical vertex sets,
with equal full path-cover profiles. For each experiment, choose a quotient
tournament between modules and hold it fixed while b varies. Different
experiments may choose different quotients. The pair may differ by one designated internal arc, as in the
inherited six-vertex pair A/B, masks 18/22.

An observer template O consists of q labelled exceptional vertices, all arcs
among them, and all arcs from them to each individual module vertex. Those
arcs are fixed independently of b in {0,1}^k. Let

    F_O(b) = H(the resulting tournament with modules M_i^(b_i)).

**INHERITED INPUT (Session247 marked/PROOF.md, §7):** F_O is multilinear of
degree at most 2q. The proof cuts each Hamilton path into module runs;
only modules adjacent to the q exceptional vertices can retain bit-dependent
endpoint restrictions, and at most 2q modules are touched.

Put chi_w(b)=(-1)^(w dot b). Multilinear degree and Walsh support degree agree,
so

    F_O(b)=sum_(|w|<=2q) fhat_O(w) chi_w(b).             (O1)

This is a degree bound in MODULE SWITCH bits, not the total tournament arc
space, and not a claim that isomorphism indicators are low degree.

## 3. Exact coset transfer and joint-moment budget

Let C be any binary linear subspace of F_2^k and D=C^perp. For a coset a+C,
let E_a mean uniform average over its bit words, with literal multiplicity.
Elementary character cancellation gives

    E_a chi_w = chi_w(a) if w in D, and 0 otherwise.

**Proposition O1 (PROVED).** Every q-vertex observer satisfies

    E_a F_O = sum_(w in D, |w|<=2q) fhat_O(w) chi_w(a).  (O2)

For distinct cosets define the separating dual weight

    d_C(a,b)=min{|w|: w in D, w dot(a+b)=1}.

Then all their q-observer means agree whenever 2q<d_C(a,b). In particular,
if the minimum nonzero weight of D is d, all cosets have the same mean for
2q<d, equal to the unrestricted cube mean.

Proof: insert (O1) and average each character. Terms distinguishing the
cosets require chi_w(a) != chi_w(b), hence weight at least d_C(a,b).

**Proposition O2 (PROVED).** For templates O_1,...,O_s with q_j vertices
and nonnegative integers r_j, all cosets have equal joint moment

    E_a product_j F_(O_j)^(r_j)

whenever 2 sum_j r_j q_j < d. More generally replace d by d_C(a,b) for two
specified cosets.

Proof: the product has multilinear degree at most 2 sum_j r_j q_j after
reducing b_i^2=b_i. Apply the same character calculation.

Thus one root can have equal first moments but unequal second moments.
Neither equality of means nor finitely many moments implies equality of
laws, equality of individual H values, equal automorphism groups, or a
lumpable state quotient. Observers that read b and then choose different
incidences are outside the fixed-template hypothesis.

## 4. A sharp parity visibility threshold at every size

Take C to be the even-parity subspace. Its dual is {0,1^k}. The even and odd
sectors therefore have equal means for every q<ceil(k/2), independently of
the chosen fixed observer template.

**Proposition O3 (PROVED, sharp physical realization).** For every k>=1,
the least number of exceptional vertices that can separate these two
parity-sector means for the inherited six-vertex switch family is

    q_min = ceil(k/2).                                 (O3)

Construction: Session247's one-root, two-module template has responses

    f(00),f(01),f(10),f(11)=4161,4297,4933,5085.

Its top Walsh coefficient is

    (4161-4297-4933+5085)/4=4.

For k=2q use the ordinal sum of q such templates. It has q exceptional
vertices and count equal to the product of the q response polynomials.
The full-character coefficient is 4^q, so the even-minus-odd mean is

    2*4^q != 0.                                      (O4)

For k=2q-1 use q-1 two-module templates and the inherited one-root single
module template with counts 41,43. Its top coefficient is -1, giving
mean difference -2*4^(q-1). All joins are actual ordinal joins of literal
tournaments. The lower bound was universal over ALL allowed quotient
orientations and observer incidences; this construction attains it.

Recursive parity provides a concrete depth interpretation: r-ary parity
composition of depth h has k=r^h leaves. Its two tournament ensembles
remain mean-blind to fewer than ceil(r^h/2) marked vertices. The threshold
is exponential in recursion depth, linear in the number of leaf modules.
This is not a claim of exponential computational hardness.

## 5. An affine-address bank with 2k sectors and exact threshold k/4

The parity construction hides one bit. A more interesting recursively
addressed bank hides m+1 bits on k=2^m leaf modules, for m>=2.

Index the modules by z in F_2^m. For b=(b_z), define the actual binary syndrome

    sigma(b) = ( sum_z b_z, sum_z b_z z ) in F_2^(m+1). (O5)

Let C_m=ker(sigma). Its dual D_m is the space of affine Boolean truth tables
w(z)=a_0+a dot z. It has dimension m+1; a nonconstant affine function has
weight k/2, while the nonzero constant has weight k. Hence the 2k syndrome
sectors each contain 2^(k-m-1) literal configurations.

**Proposition O4 (PROVED, exact pairwise interface threshold).**

1. All 2k sectors have identical q-observer means for q<k/4.
2. Every pair of distinct sectors can be separated by a fixed physical
   observer template with exactly q=k/4 exceptional vertices.
3. A family of m+1 such templates suffices jointly to distinguish every
   syndrome sector.

Proof of (1): D_m has minimum weight k/2; apply O1.

For (2), distinct sectors define a nonzero functional on D_m. Some
nonconstant affine w has value 1 under this functional: if every
nonconstant w had value 0, the constant also would, since it is the sum of
two nonconstant affine functions (for example z_1 and 1+z_1). That would
annihilate D_m and contradict distinct sectors.

Let W be the halfspace of k/2 indices where this w is 1. Arrange those
modules into k/4 successive two-module rooted templates, with all other
modules unrooted and joined ordinally outside them. Label the cube so that
its quotient order realizes this arrangement. The Fourier support of the
response lies inside W. The only nonzero word of D_m supported in W is w:
all nonzero weights are at least |W|, and equality forces the same support.
Consequently the sector mean is a constant plus

    4^(k/4) 15^(k/2) chi_w(a).

The difference for the selected two sectors is exactly, up to sign,

    2*4^(k/4)*15^(k/2),                              (O6)

which is nonzero. The observer budget reaches the universal lower bound.
The quotient and indexing can be chosen once for this pairwise realization;
the bound in (1) holds for every fixed quotient. No claim that one fixed
quotient realizes all pairwise-optimal templates is needed.

For (3), choose a nonconstant affine basis z_1,...,z_m,1+z_1. Use the above
construction separately for each basis word. These m+1 signs recover the
syndrome functional. These are separate experiments with their own fixed
quotient/observer templates, not a single tournament with only k/4 roots
that simultaneously contains every template.

### The recursive address law

Split the final address coordinate t into two halves. If their partial
syndromes are (p_0,u_0) and (p_1,u_1), then

    sigma_(m+1)(b)=(p_0+p_1, u_0+u_1, p_1).           (O7)

This is an exact recursion for the bank, not a numerology. Doubling the leaf
modules doubles the minimal observer budget and adds one syndrome bit.
For m=2,3,4,5, the bank sizes are respectively 8,16,32,64 sectors and their
exact pairwise thresholds are 1,2,4,8 roots.

## 6. How an even symmetry coexists with odd automorphism groups

Translation b->b+t is an actual action of the binary cube group on the
labelled FAMILY of tournaments: it reverses specified internal arcs.
Modulo C_m it induces the regular action of F_2^(m+1) on the 2k syndrome
sectors. Below k/4 roots, every permitted mean response is constant on
this entire sector set. Thus the response bank has a full even translation
symmetry, although each individual tournament's vertex automorphism group
has odd order.

There is no paradox and no group embedding into Aut(T). The actions are
on different objects. The metagraph move changes tournament arcs; an
individual automorphism only relabels vertices while preserving arcs.
This distinction makes operational, rather than metaphorical, the phrase
“symmetry emerges after coarse observation.” At the threshold, explicit
physical interfaces break the response-bank symmetry by measuring a
nontrivial character.

There is also no conflict with Session247's divisibility law. Our modules
have charge 3. At the affine-bank threshold q=k/4, every individual observed
count still has the factor 3^(k/2), and the displayed difference (O6) indeed
contains 15^(k/2). Becoming observable is not the same as becoming
arithmetically free. Charge-free switches could remove this obstruction,
but the exact sharp response template would have to be supplied anew.

## 7. Connection ledger and next questions

* **Module topology -> coding:** the map is literal independent arc choices
  b. Preserved: physical cube adjacency and multiplicities. Lost upon code
  quotient: particular configurations. Sidecar: full labelled module word.
* **Bounded interfaces -> Fourier support:** preserved exact H response;
  the lost information is every character above degree 2q. This is an
  observation theorem, not an isomorphism classifier.
* **Code dual -> observer cost:** a distinguishing character must touch at
  least d_C modules, while a path meets at most 2q exceptional incidences.
  O3 and O4 turn the bound into an attained physical threshold.
* **Recursive address -> emergent symmetry:** (O7) is the actual recursive
  carrier. The sector means, not individual tournament states, admit the
  low-budget even action. The restored data are syndrome characters.

Open: determine observer thresholds for other code banks and charge-free
invisible module switches; classify when a fixed quotient order attains all
pairwise thresholds; quantify distinguishability using full response laws
rather than bounded moments; relate character visibility to the separate
vertex-support spectral filtration in ../bridges/PROOF.md.

Sources: `05-knowledge/results/session170_holonomy_homotopy_20261103/metagraph_topology/NOTE.md`;
`05-knowledge/results/session247_metagraph_weave_20261110/marked/PROOF.md`, §§4,7;
`07-reflections/tournaments-as-codes.md` (historical impetus only);
`07-reflections/the-self-similar-object-is-the-operation-congruence-not-the-count-codex-S16.md`;
`01-canon/theorems/THM-4524-selfie-tournaments-loop-gauge-arc-parity-and-the-odd-mallows-sloane-partner.md`.

`python3 observers/probe.py` from the session directory checks the new code
identities and sharp mean formulas, using the credited small response tables
as inputs. It does not rerun or claim to audit the inherited seed census.
