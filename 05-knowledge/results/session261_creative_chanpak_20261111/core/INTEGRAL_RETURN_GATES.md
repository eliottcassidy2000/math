# Integral return gates: a rational cylinder can have no integer seed

**PROVED general reduction; exact finite certificates for the two-apex
single-letter case.** This is not a new Chan–Pak exclusion. It independently
explains a restriction on one proposed proof grammar.

## 1. The randomly imported distinction

Random wave2 chose `THM-2974-discriminant-cover-integral-order-smith-and-owner-boundary.md`.
Its useful mechanism is that a rational or normalized isomorphism need not
preserve an integral order or its distinguished owner. We do not transport
its algebraic-geometric conclusion. We transport the distinction to the
spine compiler: rational affine return equations need not contain a single
integer combination of physically reachable states, much less a physical
state itself. `RANDOM_WAVE2.json` retains all six choices and the exact pool.

## 2. A complete integral linear envelope

For any finite integer-matrix grammar, augment the state with a constant1.
Write v0 for the actual initial vector and B_a for the augmented digit
matrices. Define the orbit lattice

    Λ = span_Z {B_w v0 : w any finite word}.

It is effectively computable. Start from Z v0 and repeatedly adjoin each
B_aΛ. Store the resulting column lattice in Hermite normal form. The chain
is increasing inside Z^D, so it stabilizes by the ascending-chain condition
for submodules of a finitely generated Z-module. At stabilization Λ contains
v0 and is invariant under every B_a; induction contains every orbit vector.
Conversely every generator added is an integer combination of orbit vectors.
Thus the result is exactly Λ, not a guessed congruence lattice.

This is an envelope only. Negative integer combinations are allowed in Λ;
most of its vectors need not be reachable or nonnegative.

Let R be the finite rowspace of all binary-return defects from
`../wildcard/PROOF.md`, section7, and c the last-coordinate row. If

    R v=0,  c v=1,  v in Λ                                  (1)

has no integer solution, no physical seed supports that binary cylinder.
If K is a column basis of Λ, this is the integer column-membership test

    (0,...,0,1)^T in image_Z ([R K; c K]).

Hermite normal form decides it exactly. A positive answer is necessary,
not sufficient; physical reachability and nonnegativity remain unproved.

## 3. Two-apex certificate

`lattice_probe.py` uses the explicit integer matrices from the positive
spine compiler, with core a1→a2. Its orbit lattice stabilizes at rank28
inside the37-dimensional augmented representation. The successive HNF
ranks are5,16,27,28,28. The one-apex ranks are3,6,6. In both cases the gcd
of (H−constant) on the computed lattice basis is exactly2, so the familiar
Hamilton oddness appears as an integral, not merely rational, invariant.

For each of the twelve ordered distinct one-letter pairs, construct the
complete invariant defect rowspace R and test (1). Every pair is excluded.
For eight pairs even the rational equations are inconsistent. For the four
ordered pairs

    (0,1), (1,0), (0,2), (2,0),

the rational equations **are consistent**, yet the integral lattice test
rejects them. In the HNF basis of their constraint image, the target's
unique coordinates include −1/2. All remaining coordinates are integers.
Thus solving the rational affine cylinder equations can produce a false
positive even after restricting to the rational span of all physical words.
The integral order is a genuine additional gate.

This finite calculation is an all-word obstruction because both closure
steps are certified invariant: Λ contains every physical word, and R
contains every future defect row. It does not extrapolate from a finite
scan of seed lengths. The separate combinatorial proof in
`../wildcard/EXTENSIONS.md` explains the same single-letter failure through
three live path slots, inert sinks and a residual one-apex register.

## 4. Why this does not reject the next route

Multiletter blocks, state-dependent digits and module replacements change
the operators or the control graph. Recompute the corresponding integral
envelope and pullback constraints; do not reuse the twelve-pair rejection
outside its scope. Even a successful lattice intersection still needs a
literal seed word. The promising next use is to prune proposed finite-control
return systems before expensive target-state searches.

**Reproduce:** `python3 core/lattice_probe.py` from the packet root, or run
that script by its repository path. Uses SymPy for exact HNF and rational
linear algebra. `LATTICE_RECEIPT.json` records the rank closures and outcomes.

## 5. Other wave2 pulls and stopping decisions

- THM-3253 positive owner-mass cyclicity: positivity plus a spanning orbit
  suggests examining cyclic observable modules, but no map from its Singer
  gauges to tournament states was supplied. Kept as an analogy only.
- THM-851 factorization-deck closure: operation-context refinement supports
  retaining a full continuation representation. This is consonant with the
  marked-tensor trade, but does not identify the two carriers.
- THM-3757 Pell/Chebyshev no-critical-point tower: a locally successful
  condition need not supply its global mate. No useful Hamilton constructor
  emerged from its Chebyshev norm, so that route stopped.
- THM-420 clock/shell witnesses: a finite family of residues can filter
  candidates, but an exact target fibre is not a residue witness. No coverage
  claim was inferred from congruences.
- THM-3944 conductor collapse: changing to a simpler carrier can lose the
  very rank or character needed for the target. The concrete replacement is
  the integer-lattice gate above, not a claimed geometric identification.
