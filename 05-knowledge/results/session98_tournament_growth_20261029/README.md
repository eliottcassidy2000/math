# Session98: recursive tournament growth and generating functions

**Status:** proved formulas with exact computational controls; classical and
inherited inputs are explicitly credited. No claim to a polynomial-time
unrestricted enumerator or a complete metagraph atlas.

Source snapshot: `5862c94555eafea4e9a0638a17f7c661f8e45d37`.

## Main result: an exact recurrence with a smaller symmetry state

Let a_n=A000568(n), including a_0=1. For odd indices j, keep an integer state
s_j with s_1=0. Define R_0(s)=1 and

\[
nR_n(s)=\sum_{\substack{p\le n\\p\text{ odd}}}
 2^{(p-1)/2+s_p+n-p}
 R_{n-p}\bigl((s_j+\gcd(p,j)-j)_{j\le n-p}\bigr).
\]

Then **a_n=R_n(0)**. No unknown prime-family bank is needed for this equation.
Only coordinates j<=n matter. Its formal generating hierarchy is

\[
\mathcal R_s'(x)=\sum_{p\text{ odd}}2^{(p-1)/2+s_p}x^{p-1}
 \mathcal R_{s+g_p-\mathrm{id}}(2x),\qquad
 \mathcal R_s(0)=1.
\]

This is an infinite-state formal hierarchy, not a finite-order scalar ODE.
The smaller state removes a response proportional to conserved vertex mass.
Proof, multiplicity implementation, inherited-formula comparison and limits:
[projective response recurrence](recurrence/PROJECTIVE_RESPONSE.md).

Against the same raw fixed-target recurrence, the tested N=100 computation
reduces accumulated stage states **16,330 -> 13,544** and peak states
**2,953 -> 2,721**, with identical exact counts. These are state counts,
not a measured runtime speedup. Direct independent partition comparisons
cover every n<=30; both implementations agree at N=40,60,100.

## A convergent generating function with an exact global classification

The raw A000568 OGF and EGF have radius zero. Instead put

\[
\beta_n=\frac{n!a_n}{2^{\binom n2}},\qquad B(z)=\sum_{n\ge0}\beta_nz^n.
\]

B has radius 1 and a meromorphic continuation to the whole complex plane.
Its poles are **exactly 4^k**, with exact order **3k+1**, for k>=0.
The k-th symmetry bank contains p(k) nontrivial-cycle types, indexed by
ordinary partitions of k. The unique largest core is k three-cycles.
This proves the pole order rather than merely fitting its pattern.

For a core mu of odd parts >=3, write m=sum(mu), c=length(mu),
E=sum_i(mu_i-1)/2+sum_(i<j)gcd(mu_i,mu_j), and
z_mu=product_j j^(m_j)m_j!. Then

\[
 B(z)=\sum_\mu\frac{(m!/z_\mu)2^{E-\binom m2}z^m}
 {(1-2^{-(m-c)}z)^{m+1}}.
\]

The empty core gives 1/(1-z). Normal convergence away from poles is proved.
Both raw and normalized sequences are non-P-recursive, for different reasons:
superfactorial growth in the first case and infinitely many poles in the
second. See [global pole theorem](limits/ALL_POLES_CORE_EXPANSION.md) and
[normalization and SCC convolution](limits/LIMITS_AND_NORMALIZATION.md).
This is an analytic deduction from the classical Burnside count and the
inherited unit-tail core; no literature-priority claim is made.

## Partial structures now have explicit growth classifications

1. **Missing comparisons.** The defect polynomial P_n(t) counts partial
   tournaments with d unknown comparisons. It has an exact all-cycle insertion
   recursion; P_n(0)=a_n. At defect d, even-cycle contributions use at most
   2d vertices, and every partial automorphism group satisfies
   v_2(|Aut|)<=d, sharply. The first layer obeys
   2p(n,1)=e_n+v_(n-1), with marked-pair and marked-vertex counts.
   [Defect grammar](defects/MISSING_COMPARISON_GRAMMAR.md).

2. **Growing stars.** Two or more unknown comparisons incident to one added
   vertex force odd automorphism order. One unknown comparison permits even
   symmetry exactly when its endpoints are twins. Exposure flags cannot be
   forgotten during completion: C3 gives four flagged versus three unflagged
   completions. [Partial-star classification](partials/PARTIAL_STAR_CLASSIFICATION.md).

3. **Module surplus.** Prime-root growth by at most h vertices needs only a
   weighted cycle-index Taylor jet. At h=12 the sufficient jet has 70
   coordinates, but its exact ambient linear quotient has rank57 for formal
   palettes and rank33 for actual ordinary tournament palettes. Two actual
   order-12 prime families, using a pinned inherited type census, agree under
   every ordinary palette through surplus3 and separate at surplus4.
   [Base theorem](composition/RESULT.md), [sharper quotient](composition/QUOTIENT_EXTENSION.md).

4. **Changing the comparison alphabet.** With q pair colors and t colors
   fixed by endpoint reversal, the normalized GF has poles q^k of order2k+1
   when t>0. With t=0, only q^(2k) remain, of order3k+1. This unifies ordinary
   graphs, tournaments and three-state partial comparisons under one proved
   counting mechanism. [Palette pole theorem](limits/PALETTE_POLE_CLASSIFICATION.md).

[Growth synthesis](GROWTH_THEORY.md) states the common conservation/observation
principle and the remaining arbitrary-attachment boundary. These counting
quotients do not by themselves preserve all individual metagraph neighbors.

## Verification and reproduction

All code uses Python's standard library, exact integers/Fractions, and active
checks under `python3 -O`. The input `prime_types12.inherited.txt` is a pinned,
unchanged full type census; its existence data is explicitly inherited, not
newly generated. `composition/inherited_recreate.py` is the unchanged
Session86 dependency, also replayed as a standalone check.

Run all scripts from this directory with

    python3 verify_packet.py

This invokes every published mathematical checker in normal and optimized
modes, compares stdout byte-for-byte, and validates its recorded receipt where
applicable. The new recurrence reaches100; the all-poles coefficient comparison
reaches35; the inherited standalone enumerator reaches50. Finite checks do
not substitute for the proofs of analytic continuation or all-order claims.
Independent controls and review are documented in [audit](audit/AUDIT.md).

`MANIFEST.sha256` lists the shipped files (excluding itself). The final
publication checks compare the remote bytes and every relative Markdown link.
No CI status is inferred from the absence of workflow runs.
