# Independent audit: categorical, honeycomb and MRB bridge

Date: 2026-10-22. Verdict: **PASS after the four scope repairs below**.
The author applied the repairs to the source note and checker during this
audit; their repaired files were then reread and replayed successfully.
No source-worker files or Git state were modified by this auditor.

Audited source: `session31-categories/CATEGORICAL-HONEYCOMB-MRB-AUDIT.md`
and `check_categorical_parameters.py`. Source hashes are recorded in
`audit_bridges.json`. The original checker replay is saved separately as
`category_original_replay.json`; `audit_bridges.py` imports none of its code.

## 1. Repairs and their exact failure boundaries

1. **Parity is essential for the supertrace claim.** One fixed finite-order
   operator preserving every finite-dimensional grade has ordinary traces
   in a fixed cyclotomic field. For supertraces it must also preserve parity.
   Otherwise, on a grade of dimension `1|1`,
   
   `G(a) = [[a,1],[1-a^2,-a]]`
   
   has `G(a)^2=I`, ordinary trace zero, and supertrace `2a`. Taking
   `a=(r^(1/r)-1)/2` separately on each grade defeats the unrestricted
   supertrace assertion with one global involution. The ordinary-trace
   obstruction, and the usual parity-preserving supercharacter obstruction,
   survive unchanged. This is a genuine missing-hypothesis repair, not a
   numerical discrepancy.
2. **Dense/dilute branch labels require `m>=4`.** The expressions
   `g_dilute=2-4/m` and `g_dense=4/m` occupy the standard ranges `[1,2)`
   and `(0,1]` only there. At `m=3` their values are `2/3` and `4/3`, so the
   labels exchange. They meet at `m=4,g=1`. At `m=8` the stated charges
   `c_dilute=0,c_dense=-2` are correct. Negative loop fugacity for `m>8`
   is an analytic/statistical-mechanics parameter, not a positive loop
   probability measure. This does not affect the scalar algebra.
3. **Unitary-target positivity means the canonical positive spherical
   pivotal structure.** A pivotal tensor functor into that specified target
   preserves dimensions; hence a zero-dimensional TL strand cannot have a
   nonzero image. Positivity must not be attributed to an arbitrary pivotal
   choice merely from the existence of a unitary structure.
4. **State the affine KL parameter exclusion explicitly.** In the rank-one
   convention of the note, it is `k+2 not in Q_{>=0}`, with the appropriate
   KL/ordinary affine category and Lusztig's big quantum group. Positive
   integer WZW level six is instead handled by the rational integrable
   category and the semisimplified tilting construction.

All four points are now present in the repaired source. The repaired source
also makes the fixed finite order independent of grade explicit and strips
off the initial character monomial.

## 2. Algebraic obstructions: accepted with independent proofs

For integer `m>=2`, `2m^(-1/m)sin(pi/m)` is nonzero algebraic. If its quotient
by `pi` were algebraic, `pi` would be algebraic. Thus `A(m)` is transcendental.
For an object in a complex pivotal fusion category, its dimension is an
eigenvalue of an integral fusion matrix, hence an algebraic integer. The
claimed exclusion follows. It excludes literal dimensions, not analytic
normalizations or character special values.

For `beta=pi A(8)`, eliminate `y=sqrt(2)` from

`y^2-2=0`, `y(beta^4+2)-3=0`.

The resultant is `-(2 beta^8+8 beta^4-1)`. Its reciprocal, up to sign,
is `x^8-8x^4-2`, Eisenstein at two. Hence beta has degree eight and field
norm `-1/2`; an algebraic integer cannot have that rational nonintegral
norm. This verifies the obstruction without calling the source's
minimal-polynomial routine. The independent checker performs the resultant,
coefficient divisibility and norm checks exactly.

For every prime `p`, the shift `g_p+1=p^(1/p)` has Eisenstein polynomial
`x^p-p`, so `[Q(g_p):Q]=p`. If all coefficients belonged to a fixed number
field `K`, then `Q(g_p)` would be a subfield of `K` and its degree could not
exceed `[K:Q]`. This proves the fixed-field obstruction. Finite interleavings,
fixed power substitutions and grade shifts do not enlarge the coefficient
field beyond one finite compositum. The claim should not be read as covering
arbitrary analytic changes of variable with unrestricted coefficients.

The Gelfond–Schneider argument for `exp(-pi sqrt(2))` also has the correct
branch: it is the value of `(-1)^(i sqrt(2))` obtained using `log(-1)=i pi`.
Both arguments of this algebraic power satisfy the theorem's hypotheses.
Transcendence separates it from the root-of-unity parameters; it does not
exclude defining a different generic quantum group at that real parameter.

## 3. Scalar and categorical transport: accepted scope

The normalized `S00` identity and Chebyshev identity are substitutions, not
tensor functors. The positive square-root simplification of the Nienhuis
candidate uses `cos(2pi/m)>=0`, exactly explaining its `m>=4` restriction.
The source's 90-digit finite checks cover `m=4,...,40`; they are controls,
not proofs of criticality or of the general identities.

The TL strand-dimension obstruction is valid for complex-linear pivotal
strong tensor functors. At loop value zero its identity is negligible;
tensor-ideal closure kills every positive tensor power in the standard
semisimplification, leaving the unit sector. This statement does not replace
the dilute lattice algebra, which must retain vacancies and marked sectors.
The vacuum partition value at `n=0` loses all nonempty closed-loop
configurations, whereas the single-loop derivative or a marked open strand
retains information. Those observables are additional data.

## 4. Primary-source checks

- Neshveyev–Tuset, Theorem 3.3 and the paragraph immediately following it,
  give the generic braided equivalence, preserve highest weights and the
  underlying-vector-space functor, and state the extension to
  `hbar not in Q*`. This is compatible with `hbar=i sqrt(2)` and does not
  identify the resulting infinite category with `SU(2)_6`.
  [Primary text](https://arxiv.org/pdf/0711.4302).
- Creutzig–Ridout–Rupert §1.2 gives the excluded nonnegative rational shifted
  levels and the lacing-number convention; for `sl2` the lacing number is
  one. Gaitsgory's §0.1 distinguishes the established `G(O)`-integrable KL
  equivalence from the conjectural Iwahori-integrable extension.
  [Parameter statement](https://arxiv.org/html/2112.13167v1#S1.SS2),
  [extension distinction](https://arxiv.org/html/1810.09054).
- Huang's theorem requires simplicity, no negative weights, one-dimensional
  weight zero, self-contragredience, complete reducibility of every
  `N`-gradable weak module, and `C2`-cofiniteness. This supports the source's
  warning against unconditional VOA modularity.
  [Primary text](https://arxiv.org/pdf/math/0502533).
- Creutzig–Jiang–Orosz Hunziker–Ridout–Yang, Theorem 4.2.6, gives a braided
  tensor structure on `O_c^fin`. Its identification with lower-bounded
  `C1`-cofinite generalized modules is made in §3. Neither result alone
  asserts finite semisimple modularity at `c=0`.
  [Primary text](https://arxiv.org/html/2002.03180v2).
- Zhou–Batchelor §2.1 supplies the branch ranges. Read–Saleur §5 retains
  vacancies and marks the proposed dilute-lattice Morita equivalence as a
  conjecture. Jacobsen–Ribault–Saleur §6 explicitly replaces the dense state
  space by one including the trivial representation and the planar
  algebra by a Motzkin/dilute analogue. None is used here as a rigorous
  proof of the honeycomb SAW scaling limit.
  [Branch ranges](https://arxiv.org/pdf/cond-mat/9611156),
  [dilute polymers](https://arxiv.org/pdf/hep-th/0701117),
  [vacancies and diagram algebra](https://arxiv.org/pdf/2208.14298).
- Gukov–Vafa §5 states the two-torus rationality condition using both
  complex structure and complexified Kähler parameter in the same imaginary
  quadratic field. The note correctly refuses to infer it from one CM
  complex structure alone. [Primary text](https://arxiv.org/pdf/hep-th/0203213).

## 5. What this audit does not establish

No coefficientwise MRB character, MRB-preserving tensor functor, analytic
family for the CM defect, or dilute-polymer categorical equivalence is
constructed. No arithmetic conclusion about the MRB constant or either sign
of its difference from `A(8)` follows. These remain explicit open construction
problems rather than consequences of matching numerical parameters.
