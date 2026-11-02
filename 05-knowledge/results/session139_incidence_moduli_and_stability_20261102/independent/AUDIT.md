# Session139 independent mathematical audit

Status: VERIFIED independently; no Git mutation. Audited against source HEAD
`291b4d9757842c02aa451beb65acfade2848609b` and explicitly named local drafts.
The computations below are checks, not extrapolations to unenumerated orders.

## 1. Weighted five-vertex spectral loss

The custom verifier `audit_bridges.py` uses only Python standard-library
`Fraction` and independently implemented bivariate polynomial arithmetic.
It does not import the discovering agent's code or SymPy.

For `K14=a,K25=b` (opposite entries skew, all other entries zero),
`A=(J-I+K)/2`, `w=(1,1,-2,0,0)`, and `Q=I-ww^T/3`, it verifies coefficientwise:

* `QQ^T=I`, `Q1=1`, and `B=QAQ^T` has diagonal zero and complementary entries;
* every entry of A and B is in `[0,1]` for all `|a|,|b|<=1`;
* the complete polynomial `det(I-zA)=det(I-zB)` and scalar V, tau, u, v agree;
* `H(A)=15/2`, `H(B)-H(A)=a^2 b^2/54`;
* `M23(A)=(10-a^2-b^2)/16`, `M23(B)-M23(A)=a^2 b^2/108`.

At `(a,b)=(1/4,1/2)`, H values are `15/2` and `25921/3456`, V is
`155/64`, tau is `5/2`, u is `233/1024`, v is `-37/256`; hence the
H and Delta differences both equal `1/3456`. The characteristic coefficients
in ascending z order are `1,0,-155/64,-5/2,-1019/1024,-37/256`.

The equality of the full spectral polynomial also follows structurally from
orthogonal conjugation. More strongly, because Q fixes the all-one root,
`1^T A^r 1=1^T B^r 1` for every nonnegative integer r, and the corresponding
all-one rooted resolvents coincide. Thus that single marked scalar response
also cannot restore H. The missing information is vertex-exclusion/disjointness,
not merely omission of an ordinary root moment.

At order five the exact repair is

    H=1+2 tau+2 C5=1+2 tau+2 det(A)+2 M23.

The first equality follows by taking expectation of the binary five-vertex OCF
identity; all monomials use distinct unordered edge variables. The determinant
identity is the permutation split into cycle types (5) and (2,3). Consequently
one M23 sidecar suffices together with the spectrum at n=5. At n<=4, H=1+2tau,
so the loss first occurs at n=5 in the weighted domain.

### Historical scope

`01-canon/theorems/THM-499-the-spectral-reframe-boundary-H-is-finer-than-the-spectrum.md`
already proves the *binary* spectral boundary at n=6, with disjoint triangle
pairs as the missing coordinate; binary H is spectral through n=5. This does
not conflict with the weighted n=5 boundary. The packing hierarchy continues
in `THM-505-the-ocf-non-spectral-defect-H-equals-spectral-skeleton-plus-witt-defect.md`.
No claim of global novelty follows from this bounded historical inspection.

## 2. Odd-fourset cocycle and sharp nonzero support floor

Let c be the F2 indicator on triples of a cyclic triangle and define epsilon
on foursets by epsilon(X)=sum_{Y subset X,|Y|=3} c(Y) modulo two. Local four-vertex
classification gives epsilon=(Pf(S_X)^2-1)/8. Thus epsilon=delta c and
`delta epsilon=delta^2 c=0`: in each five-set every triple is counted twice.
This proof is all-order and does not depend on finite extrapolation.

If epsilon(X)=1 on one fourset X, then for every outside vertex v the five-set
X union {v} must have another odd face Y_v. The Y_v are pairwise distinct,
because Y_v contains v and no other vertex outside X. Together with X these
produce at least `1+(n-4)=n-3` support faces. The bound is sharp: take one
cyclic triangle as an ordinal block and all other vertices singleton blocks.
Exactly the foursets containing that triangle have odd cyclic-triangle count.
This argument is the generic nonzero cocycle support bound specialized to
four-element faces; it is not a claim that all cocycles are tournament-realized.

The verifier independently determines cyclic triangles from outdegrees and
Pfaffians from signs, checking their identity on all labelled tournaments of
orders4,5,6, all five-face cocycle equations, and the support floor. Distributions:

* n4: `{0:48,1:16}`;
* n5: `{0:384,2:640}`;
* n6: `{0:3840,3:7680,4:11520,5:4608,6:5120}`.

Ordinal C3 witnesses are additionally checked for every n=4..20. Direct Paley7
construction gives exactly14 of35 odd foursets, not35. A singleton four-face
support on five vertices explicitly fails the cocycle condition.

### Historical scope and correction target

`THM-1455-skew-pfaffian-parity-and-a-mod-16-law.md` already records the five-set
evenness census and the all-odd-order mod16 propagation. The separate file
`THM-1455-why-eleven-spectra.md` proves epsilon's cyclic-triangle interpretation,
but still calls k=2 at n7 unexplained and erroneously says Paley7's k=14 means
all foursets are odd. The cocycle support proof resolves that stated local
question (n7 requires k>=4), and the direct Paley check refutes that sentence.
This is a recovery/correction; the bounded inspection is not a repository-wide
priority certificate. No Riemann-hypothesis implication has been inferred.

## 3. Structural audit of PROOF_N748.md

The proof at `/tmp/session139_research/stability/PROOF_N748.md` was audited
independently at the level of implications, quantifiers, and equality. Subject
to its explicitly inherited Session138 inequalities, its argument is sound:

1. For a candidate Delta<=kappa V with V>0, the retained full greedy matching
   gives B<=kappa D, while the nonlinear lower bound keeps B and yields
   `kappa V>D^2/(2n)+c`. Completing the square gives
   `V-D>c/kappa-n kappa/2`; no correlation or sign assumption is omitted.
2. Removing L>=1 light matching edges leaves positive heavy cardinality
   r_h>n/15. Thus n>2r_h and the stronger matching coefficient is applicable.
   The replacement r_h<=R-L is valid because all factors are positive.
3. The epsilon estimate and V<=n/7 imply `L<CD`; the replacement of the
   additive error by error times L uses exactly L>=1, not an unstated limit.
4. Distinct assigned edge groups contribute disjoint wedge pairs. Therefore
   heavy mass and light mass can be charged simultaneously:
   `W>=eta D_h+D_l^2/(4L)>=eta D_h+D_l^2/(4CD)`.
   This remains valid when D_l=0. With x=D_h/D and D>0, minimizing on [0,1]
   is legitimate; the checked condition C eta<1/2 places the critical point
   x=1-2C eta in [0,1]. The minimum is h=eta-C eta^2.
5. h>kappa contradicts B<=kappa D and B>=W, ruling out light selected edges.
   Then all selected weights exceed kappa; D>0 would again contradict
   B<=kappa D. Consequently D=B=0, and inherited matching equality applies.

The proof therefore establishes all n>=748, not n10..747. It retains the
same fair-pair equality classification as the inherited matching theorem.
A separate audit is checking endpoint arithmetic; this section is specifically
the structural audit, not a claim to independently reprove every Session138 input.

## Reproduction and hostile controls

    python3 audit_bridges.py
    python3 -O audit_bridges.py

The two output files agree byte-for-byte. Four controls (normal and -O for
each change) deliberately substitute the false H coefficient1/55 or the false
Paley count35. Every one raises an explicit RuntimeError, even with assertions
disabled; see `negative_controls.out`. No check relies on Python `assert`.
