# The F1 shadow agrees; the quadratic and symplectic extensions do not

2026-10-29. **PROVED elementary counts and graph criteria / CITED arithmetic
model / FINITE-EXACT controls.** No F1 model of Sha is constructed.

The incoming [symplectic gluing packet](../session96_sha_pairing_growth_20261029/f1/SYMPLECTIC_SKELETON_AND_GLUE.md)
uses a valid finite symplectic geometry. A further distinction matters for
the arithmetic comparison: BKLPR's Selmer model is built from maximal
isotropics of a **split quadratic** space; its alternating matrices appear
as graph charts. Its predicted arithmetic distribution is a conjecture,
not the theorem that the random model exists. Source: [BKLPR, sections
1.2--1.4 and 4](https://math.mit.edu/~poonen/papers/rst_distribution.pdf).
This note proves the finite geometry below independently.

## F1. The same graph chart tests different matrix symmetries

On V=F_q^n + F_q^n, write

    Q(x,y)=x^T y,
    omega((x,y),(x',y'))=x^T y' - x'^T y.

For an n-by-n matrix A, let Gamma_A={(x,Ax):x in F_q^n}. Direct calculation
gives, in every characteristic,

    Gamma_A is omega-isotropic iff A=A^T,
    Q vanishes on Gamma_A iff A^T=-A and every A_ii=0.     (F1)

For the second equivalence, substitute each standard basis vector and then
each sum of two distinct basis vectors in x^TAx. Thus the counts of graphs
are q^[n(n+1)/2] and q^[n(n-1)/2], respectively. At q=1 both polynomial
counts become one, but their first derivatives differ by n.

In characteristic two, minus equals plus, but the zero-diagonal condition
remains. Even the n=1 graph y=x is symplectic-isotropic and not Q-isotropic.
An alternating matrix over an odd field is generally **not** symmetric;
therefore a tournament skew matrix naturally supplies the quadratic graph
chart, not a symplectic graph chart using the same coordinates.

## F2. Full counts and the first lost jet

Let L_n(q) count symplectic Lagrangians in dimension 2n and O_n(q) count
all maximal Q-isotropic n-subspaces for the split form above. Then

    L_0=O_0=1,
    L_n=(q^n+1)L_(n-1),
    O_n=(q^(n-1)+1)O_(n-1),

so for n>=1,

    L_n(q)=product_(i=1)^n(q^i+1),
    O_n(q)=2 product_(i=1)^(n-1)(q^i+1).                  (F2)

**Proof.** Double-count pairs (isotropic line, maximal isotropic subspace).
In the symplectic space every line is isotropic. The number of lines is
(q^(2n)-1)/(q-1), and each maximal isotropic contains (q^n-1)/(q-1)
lines. The line's perpendicular quotient is symplectic dimension 2n-2.

For Q, the number of nonzero isotropic vectors is
(q^n-1)(q^(n-1)+1): for x=0 there are q^n choices of y, and for each
nonzero x there are q^(n-1) choices, then subtract the zero vector.
After division by q-1 and by the lines per n-space, the second recurrence
follows. An isotropic line splits off a hyperbolic plane, so its quadratic
perpendicular quotient is again split. This holds in characteristic two:
choose a polar partner and adjust it by a multiple of the given vector
to make it Q-isotropic.

Both thin values are 2^n, but their first derivatives differ:

    L_n(1)=O_n(1)=2^n,
    L_n'(1)-O_n'(1)=n 2^(n-1).                           (F3)

Thus for this precise pair of counting polynomials, retaining the first
q-1 jet suffices to distinguish what the constant term loses. It is not a
universal minimal-carrier theorem for arithmetic objects.

At n=2 the complete RREF census gives:

| q | symplectic Lagrangians | quadratic maximal isotropics |
|---:|---:|---:|
| 2 | 15 | 6 |
| 3 | 40 | 8 |
| 5 | 156 | 12 |
| 11 | 1464 | 24 |

Each row counts actual subspaces, including all components on the quadratic
side. At q=1 both polynomials give four. No geometry is being evaluated
over a literal field of cardinality one.

## Relation to previous F1 work

Lorscheid's [Tits--Weyl models](https://arxiv.org/abs/1201.1324) provide a
precise context for Weyl-group shadows. The repo's [Session59 signed flags](../session59_f1_geometry_and_information_20261025/signed-flags/THEOREM.md)
already retains whole q-polynomials before specialization. Here the
underlying issue is especially sharp: the same chosen coordinate-pair
skeleton can have different thick geometries. In both spaces, the 2^n
coordinate choices select one axis from each hyperbolic pair.

The polynomial jet separates their aggregate counts, while the actual
quadratic equation Q supplies the object-level predicate. Neither one
recovers arithmetic localization data by itself. This is consistent with
the older [root-algebra/F1 boundary](../normseven_f1_20261025/KUMMER_AND_F1.md):
additive relations cannot be removed from a construction simply because
its multiplicative or combinatorial shadow remains simple.

The quadratic distinction here and the characteristic element of a
Cassels--Tate pairing are different typed objects. No identification of
them follows from their shared sensitivity to characteristic two.
