Warning: truncated output (original token count: 133970)
Total output lines: 6848

# Mistakes Log

**Purpose:** preserve failed implications, minimal witnesses, and repaired
statements so future agents do not repeat them. Search this ledger before using
a historical synthesis or promoting computational evidence.

Format per entry:
- What was assumed / done
- Why it was wrong
- The correct framing

## MISTAKE-244 (2026-07-22, codex audit of HYP-8950/S109 and HYP-8955/S110) -- a Newton-face Wronskian was called DvdK, and a linear edge polynomial was confused with a weighted-linear coordinate

- **What was claimed:** the resonant face equation `1 in im J(F,-)` was called
  a one-variable DvdK constant-term problem; `xy` was called nonprimitive, and
  a linear collapsed edge polynomial was treated as a weighted-linear
  coordinate and a terminating descent.
- **Why it is wrong:** at terminal weight the bracket is the finite Laurent
  Wronskian `[x]F[y]G-[y]F[x]G`, with no supplied DvdK map. The edge polynomial
  of `x^3+y^2` is linear at weights `(2,3)`, but its degree is `6>5` and it has
  no mate; `xy` is primitive. Coordinates `x+(y+x^2)^m` also have repeated
  proper-power top faces.
- **Correct framing:** THM-2102 proves power-free top-face descent directly.
  Proper-power faces obey a first-defect Hamiltonian equation whose live data
  are the approximate-root quotient and resonant target-shear class. S109's
  Hamiltonian-cokernel viewpoint survives; no GMC/DvdK-to-JC implication is
  known.

---

## MISTAKE-243 (2026-07-22, codex dependency-aware Lean audit of THM-2101) -- an unbuilt additive-residue source was called kernel-checked

- **What was claimed:** the additive DvdK checkpoint described Check A, the
  additive orbit contradiction, and the full-root Lagrange identity as
  kernel-checked and printed only the standard Mathlib axioms.
- **Why it is wrong:** `lake build TournamentH7.GMC2LaurentShiftCheckA`
  reaches the source only after building its missing dependencies, then fails
  on Laurent coefficient application and ambiguous `C`, an unsolved `nsmul`
  commutation goal, and a failed Lagrange rewrite. The subsequent `#print
  axioms` lines include `sorryAx`; a direct `lake env lean` had stopped earlier
  at a missing dependency `.olean` and therefore certified nothing.
- **Correct framing:** THM-2101 remains RESERVED. Repair and dependency-build
  the formal core with no `sorryAx`, then separately prove the missing
  analytic-germ-to-splitting-field subset bridge before claiming an additive
  DvdK bypass. The formal core was subsequently repaired and root-imported;
  it now builds with only the standard Mathlib axioms, but it does not supply
  the missing analytic bridge or prove THM-2101.

---

## MISTAKE-242 (2026-07-22, codex referee of THM-2098) -- a pure-transverse collision budget was transported to mixed lanes and an unguarded depth-zero core

- **What was claimed:** the first pushed THM-2098 consumer said every live
  rank-eight-through-eleven coefficient-plane cover fell into either seven
  transverse bands governed by the exact budget `5(n-7)/49`, or seven
  guard-proportional bands. It treated the high vertical cover pointwise.
- **Why it is wrong:** the budget uses both that every band has mass `5/49`
  and that the transverse family itself covers the guard complement. Neither
  holds in a low mixed row such as `(n,r,t)=(8,1,7)`. Live nontrivial towers
  have guarded terminal sizes only `7..10`; the depth-zero size-eleven core
  has no preceding guard. Finally, containment gives an almost-everywhere
  two-torus cover, so a single uncovered boundary fiber proves nothing.
- **Correct framing:** for a rank-two character image, `r=0` has the exact
  collision/tree budget; `1<=r<=n-7` is only count-isolated; and `r>=7`
  forces the vertical bands to cover the one-dimensional guard complement
  almost everywhere, by a positive-measure leftover plus Fubini argument.
  This applies to live guarded sizes `8,9,10`; rank-one images keep their
  freeze route and the depth-zero size-eleven core remains outside the split.

---

## MISTAKE-241 (2026-07-22, codex semantic audit of HYP-8935/S106) -- floating root asymptotics were promoted to a formal-series/Hensel reduction

- **What was claimed:** HYP-8935 described THM-2067 as four
  “Mathlib-ready” pieces plus one unramified-Hensel gap, called an elementary
  formal-log derivation and unramified descent verified, and inherited
  HYP-8931's unique-face bypass as a working alternate consumer.
- **Why it is wrong:** the NumPy checks compare finite floating evaluations.
  Check B reports the small-root product approximately equal to its leading
  term even for examples whose displayed higher moments are nonzero; a ratio
  tending to one cannot exclude fractional powers or identify a formal series.
  No proof selects the Hensel roots inside the rational splitting field,
  descends through the roots of unity, or establishes the required local/global
  subset product. HYP-8931's consumer is independently void by MISTAKE-240.
- **Correct framing:** HYP-8935 is an open dependency map. The later
  `GMC2OrbitProduct.lean` kernel-checks the abstract transitive orbit-product
  identity and fixed-product valuation-zero lemma. Irreducibility, small-root
  factor selection, Hensel/descent, and the local/global bridge still separate
  that core from general complex `DvdK1` and NC2.

---

## MISTAKE-240 (2026-07-22, codex semantic audit of HYP-8931/S230) -- the lowest-face unique-channel bypass has an inconsistent class predicate

- **What was claimed:** HYP-8931 called
  `exists_nonzero_lowest_face_seed_of_uniqueChannel` a nonvacuous replacement
  for the DvdK seed on an asserted `98/116` or “84%” unique-channel class, and
  described the remaining work as a one-line connection to NC2.
- **Why it is wrong:** `LowestFaceUniqueChannel P` quantifies over every
  `lambda`, `delta`, and finite set satisfying only an exact-level-set
  equivalence. Take `lambda=0`, `delta=-1`, and `F=empty`. Radial
  exponents are nonnegative, so this is a valid empty level set for every
  `P`. The predicate then demands a positive-mass balanced composition on
  the empty subtype, which cannot exist. Thus no polynomial satisfies the
  premise. Lean correctly checks the implication, but the theorem is vacuous.
  The `98/116` number is only a bounded support census, not a theorem or a
  universal proportion.
- **Correct framing:** HYP-8930's fixed-support theorem
  `ct_ne_zero_of_unique_balanced` and its `dvdk1_of_uniqueChannel`
  corollary are substantive: a specified unique balanced composition prevents
  coefficient cancellation. HYP-8931 is quarantined until its class predicate
  is restricted to a genuine nonempty/straddling lower face and a resulting
  seed is wired into the NC2 descent. Independently,
  `GMC2NC2.heightWitnessSupplier_holds` is now root-imported; general
  `DvdK1` is the sole explicit formal premise.

---

## MISTAKE-239 (2026-07-22, codex THM-2080 direction audit) -- terminal guard containment was reversed into a cover of the guard instead of its complement

This correction initially landed under the already-used `MISTAKE-231`; 239 is
the canonical ID. MISTAKE-231 remains the observable-relative entropy
correction.

- **What was claimed:** the first pushed version of the mixed-radius fold
  argument said that `G_Q subset E_h` implies
  `E_h subset union_(q in Q) D_q`. It derived a low reduced-ratio invoice for
  rank six and described `ab<=36` as the remaining projective conductor.
- **Why it is wrong:** points of `G_Q` are allowed to lie inside `E_h`; that is
  exactly what the containment permits. The implication runs on the other
  side. A point outside `E_h` that lies in no `D_q` would be in `G_Q` and
  violate containment, so the correct statement is

  ```text
  E_h^c subset union_(q in Q) D_q.
  ```

  The exact fold formula itself was correct; only its set-cover consumer was
  reversed.
- **Correct framing:** THM-2080 proves the mixed overlap floor
  `measure(D_q intersect E_h)>=1/42` for odd `h`, with equality only at
  `q=6h`. Hence each danger comb has outside-guard capacity at most `5/42`.
  Six distinct combs have total capacity strictly below `5/7`, so they cannot
  cover `E_h^c`, whose measure is `5/7`. The corrected argument eliminates
  terminal rank six altogether and sharpens tower depth to four; there is no
  surviving `ab<=36` lane.

---

## MISTAKE-238 (2026-07-21, codex audit of HYP-8920/S227) -- the empty safe set of the full dyadic counterexample was transported through a homeomorphism that starts only at its nonempty quotient core

This correction initially landed under the already-used `MISTAKE-230`; 238 is
the canonical ID. MISTAKE-230 remains the inverse-form-class outcome-counting
correction.

- **What was claimed:** HYP-8920 combined mirror parity with THM-2075 and
  asserted that a strict dyadic-seam counterexample `S=2C union {x,y}` has
  `chi(G_S)=0`, that this zero Euler characteristic descends through the
  dyadic tower, and hence that the hereditarily primitive terminal core has
  `chi(G_(Q_r))=0`. It then called the dyadic wall equivalent to exclusion of
  such zero-component terminals.
- **Why it is wrong:** THM-2075's homeomorphisms are

  ```text
  G_C=G_(Q_0) -> G_(Q_1) -> ... -> G_(Q_r).
  ```

  They do **not** start at `G_S`. At the inner guard levels, exactly one safe
  child survives and doubling is bijective. At the final outer step from
  `C` to `S`, the two original odd tails kill both lifts over every point of
  `G_C`; there are zero safe children. Thus `G_S=empty` while `G_C` is
  nonempty, and no homeomorphism connects them. Transporting `chi(G_S)=0`
  across THM-2075 is a domain error.
- **Correct framing:** THM-2075 gives

  ```text
  chi(G_C)=chi(G_(Q_r))>0.
  ```

  The strict dyadic obstruction is instead THM-2061's folded containment of
  the nonempty `G_C` in the two-tail danger locus. Mirror parity remains
  valid: every divisor-complete terminal contains an even speed, reversal
  acts freely on its nonempty safe set, and therefore its component count is
  even (hence at least two). That fact does not stop the original tails from
  killing the two outer lifts. The live carrier remains the terminal
  component word, its THM-2075 dyadic addresses, and the original tail-owner
  constraints.
- **Unaffected work:** the two-charge DvdK Lean theorems in the same S227
  session are independent of this LRC level mismatch and remain valid.

---

## MISTAKE-237 (2026-07-21, codex audit of HYP-8905/S225/S103) -- a valid binary symmetric-Hessian subcase and three analogous descent programs were promoted to an NC2-to-JC bridge and equivalent formulations of JC(2)

This correction initially landed under the already-used `MISTAKE-229`; 237 is
the canonical ID. MISTAKE-229 remains the tangent-disk/Heegner correction.

- **What was claimed:** S103 wrote `NC2 => GMC2 => ... => JC(2)`, identified
  two-dimensional Zhao/Laplacian vanishing with the Gaussian nullcone, and
  called the rank-one/rank-two nilpotent-matrix divide the true/false Jacobian
  boundary. S225 called VC(4), leading-form descent, and a proposed
  Lame-for-polygons bound three equivalent forms of the JC(2) obstruction.
- **Why it is wrong:** no arrow from GMC(2), let alone NC2, to JC(2) is proved
  in THM-1830 or elsewhere in the repository. The symmetric reduction relevant
  to a general planar Keller map increases dimension, landing on the
  four-variable symmetric/Vanishing-Conjecture problem; solving the symmetric
  problem in two variables does not solve the image of that reduction.
  Laplacian iterates and complex Gaussian moments are different functionals on
  different spaces until a predicate-preserving transform is supplied.
  Finally, the existence of rank-two nilpotent matrices in dimension three
  neither constructs nor classifies Keller collisions. The leading-form,
  Jelonek, VC(4), and continued-fraction programs may share a termination
  motif, but no equivalences among those programs were proved.
- **Correct framing:** retain the exact binary homogeneous calculation for
  `d>=2`: if `P=A z^d+B zbar^d` is harmonic, then
  `det Hess(P)=-4d^2(d-1)^2AB|z|^(2d-4)`, so nilpotence forces one side and
  THM-2063 makes the resulting gradient map tame. This is a classified
  two-dimensional symmetric subcase. Treat VC(4), planar leading forms,
  inverse-Jelonek data, and Lame-for-polygons as separate live routes. Any
  transfer among them must state its map, dimension, predicate, and loss.

---

## MISTAKE-236 (2026-07-21, codex audit of THM-1330) -- necessary Keller-monoid anatomy was titled and described as a classification of all Jacobian counterexamples

This correction initially landed under the already-used `MISTAKE-228`; 236 is
the canonical ID. MISTAKE-228 remains the Paley/LRC-role correction.

- **What was claimed:** THM-1330 was titled "the exact picture of the set of
  all counterexamples" and called the Jelonek hypersurface, complement cover,
  monodromy, and degeneration package a "covering-triple classification."
- **Why it is wrong:** the same theorem's honest ledger says that the inverse
  Jelonek realization problem, realizable generic degrees, irreducible seeds,
  unique factorization, and `JC(2)` are open. A package recorded by every map
  is a necessary invariant, not a classification, unless one proves both
  realization and completeness up to the stated equivalence. The phrase
  "subject to polynomial-with-constant-Jacobian gluing" simply puts the hard
  problem back into the definition. Monoid ideality and degree factorization
  likewise organize unknown seeds without enumerating them.
- **Correct framing:** THM-1330 is a necessary-structure atlas: Keller
  nonunits have function-field degree at least three, form a two-sided ideal,
  factor finitely under composition, are etale and nonproper, and impose
  covering/monodromy restrictions on the Jelonek set. The verified explicit
  families occupy cells of this atlas. A complete classification would still
  have to solve realization and equivalence for the covering data and classify
  irreducible seeds. THM-2063 legitimately classifies one empty planar stratum;
  it does not turn the atlas into a global classification.

---

## MISTAKE-235 (2026-07-21, audit of S102 / HYP-8879) -- a shared kernel-sum syntax and four truncated examples were promoted to an LRC-to-GMC equivalence and an AP-core reduction

This correction subsumes the narrower strict-measure repair that initially
landed under a colliding `MISTAKE-227` in commit `c270f8aaa`. Canonical
`MISTAKE-227` remains the AP-chain saturation error below; the modular-cusp
bridge is `MISTAKE-233`, and additive-energy/augmentation is `MISTAKE-226`.

- **What was claimed:** S102 wrote the Fourier expansion of the LRC lonely-set
  measure as a sum over the integer relation lattice, compared it with a
  balanced-channel expansion of a Gaussian moment, and called them the same
  noncancellation problem. It identified the zero-frequency term `(6/7)^13`
  with THM-878's clock floor, then used four supports of sizes four and five
  and a frequency cutoff `|k_i|<=9` to announce that LRC(14) reduces to
  maximal-resonance or AP-neighborhood cores.
- **Why it is wrong:** “sum over an integer kernel” is a schema, not an
  equality of typed objects. The LRC sum ranges over all vectors in
  `ker_Z(v)` with sinc weights from discontinuous interval indicators; a fixed
  GMC moment ranges over mass-constrained balanced channels with multinomial,
  factorial, radial, and coefficient weights. No map preserves either
  nonvanishing predicate. THM-878's `6/7` is the minimum of a pair-overlap
  functional averaged over primitive residue classes, not the coefficient
  `hat g(0)=6/7`; a repeated number does not identify the functionals. The
  script computes only the box `[-9,9]^n`, supplies no Fourier-tail bound, and
  its own output does **not** reproduce the direct integral (`0.0145` versus
  `0.0000` on `{1,2,3,4}`). Four low-dimensional examples cannot imply a
  uniform thirteen-speed theorem. More strongly, the known tight non-AP
  Goddyn--Wong row `{1,...,11,13,24}` has strict lonely measure zero, so its
  full nonzero-frequency sum cancels the main term. Thus the zero-measure locus
  is already non-AP; “AP-neighborhood” was never given a metric, radius, or
  theorem. Even a unique primitive relation does not save the transfer:
  `S={1,2}` at threshold `1/3` has kernel generated up to sign by `(-2,1)`,
  yet all harmonic multiples cancel the `1/9` main term and the strict lonely
  measure is zero. Sparse shortest relations do not control the signed sum of
  all higher relations. Moreover measure zero conflates a valid tight weak
  witness `M(S)=delta` with failure `M(S)<delta`; LRC needs the boundary data.
- **Correct framing:** the exact LRC relation-lattice expansion and zero mode
  are valid through THM-2047's Fejer-regularized limit (with earlier forms in
  THM-501/503/515/1061). They suggest diagnostics—relation height, support,
  Fourier decay, and initial channels—but transfer to GMC or tournament zeta
  requires an explicit weight-preserving map or a separately proved tail
  inequality; isolated weak witnesses additionally require phase-height/Euler
  boundary data. The S102 computation is a finite illustration only. It proves
  neither an AP-neighborhood reduction nor a new LRC theorem. THM-2051 is the
  proved sparse-relation localization; its residual contains every genuine
  support-`3..5` relation plane, not just AP cores.

## MISTAKE-234 (2026-07-21, THM-2070 audit) -- quotient data were promoted to weighted noncancellation

- **What was claimed:** S223 correctly encoded the return lengths
  `R={m:0 in mS}` of a two-sided Laurent support as an additive semigroup, then
  claimed that for mixed/complex coefficients only finitely many reachable
  lengths can cancel. It therefore called the Frobenius conductor an effective
  DvdK bound and announced a complete elementary replacement for the imported
  one-variable constant-term theorem.
- **Why it is wrong:** a support semigroup records which coefficient monomials
  exist, not whether their signed sum vanishes. The aperiodic support
  `S={-2,-1,1,2}` already has return lengths `2` and `3`, hence every
  `m>=2`. But for

  ```text
  f(z)=z-z^(-1)+z^2-z^(-2)
  ```

  one has `f(z^(-1))=-f(z)`. Constant-term invariance under inversion gives
  `CT(f^m)=(-1)^m CT(f^m)`, so every odd power vanishes. Thus cancellations
  persist infinitely far beyond the support conductor even though the support
  gaps have gcd one. S222/HYP-8890 explicitly leaves the general complex
  dominant-saddle/noncancellation step open, so citing it cannot close this
  gap.
- **The same lost-coordinate error in THM-1840:** the original functional
  corollary called `F` "diagonal by charge" (or said that it depended only on
  charge) and then used `F(h)=F(1)CT(h)`. Diagonality does not imply projection.
  For `Lambda=u+u^(-1)`, take `F(1)=1`,
  `F(u^2)=F(u^(-2))=-1`, and `F(u^k)=0` otherwise. This is a perfectly
  charge-labelled linear functional, but `F(Lambda^2)=0` while
  `F(1)CT(Lambda^2)=2`. The missing coordinate is the functional's weight on
  every nonzero charge, just as the semigroup quotient forgets every channel
  coefficient and phase.
- **Correct framing:** the return semigroup, pair-period law, and conductor are
  exact **reachability** data. For positive coefficients,
  `CT(f^m)!=0 iff m in R`; for signed or complex coefficients, retain channel
  phases and use a separate noncancellation theorem. HYP-8878's unique-minimum
  channel is one valid elementary sufficient case. HYP-8895 is therefore a
  useful support-sidecar construction, not a DvdK replacement or an effective
  coefficient-uniform bound. Likewise THM-1840's functional corollary requires
  the charge-**projecting** rule `F(u^k)=0` for every `k!=0`, with `F(1)!=0`;
  the two-charge constant-term formula itself remains proved. THM-2070 is the
  repaired theorem and supplies both exact counterexamples.

## MISTAKE-233 (2026-07-21, audit of S219/S220 / HYP-8880) -- classical theta and modular-curve facts were promoted to an LRC cusp-form obstruction without a map

- **What was claimed:** S219 called THM-515's sinc-weighted relation-lattice
  sum a modular theta series and transferred the binary-form decomposition
  `theta=Eisenstein+cusp` to LRC. S220 retracted that particular attachment but
  replaced it with the stronger claim that LRC scaling is `Gamma_0` level,
  clock divisors are modular cusps, a covering-min second moment is a form on
  `X_0(14)`, and the level-14 newform `f_14=14a` is the LRC(14) obstruction.
  It further identified genus with S218's hidden entropy and asserted that the
  period field of `f_14` is `Q(sqrt(-7))`.
- **Why it is wrong:** neither session constructs a map from an LRC row,
  good-set functional, or floor moment to a modular form, proves a modular
  transformation law, or shows that the LRC predicate is preserved. A cusp
  of `X_0(N)` is a `Gamma_0(N)`-orbit in `P^1(Q)`, not a subgroup of the
  additive clock `Z/NZ`; equality of two divisor counts is not an
  identification of objects. Dilation gives the exact LRC symmetry
  `M(cS)=M(S)`, but it cannot by itself be level structure: scaling the
  `12`-clock to `24` leaves `M` unchanged while `X_0(12)` has genus zero and
  `X_0(24)` has genus one. The proposed modular invariant therefore changes
  under an operation that preserves the alleged source object.

  The discriminant bridge is also false. The normalized level-14 elliptic
  newform has coefficient field `Q`; its elliptic-curve isogeny class `14a`
  is non-CM. In particular the standard curve `14a1` has
  `j=9938375/21952`, not an algebraic-integer CM `j`, so it has no CM period
  field `Q(sqrt(-7))`. Its coefficients `a_2=-1` and `a_7=1` are Fourier
  coefficients at bad primes, not labels for two cusps. `X_0(14)` is only the
  first positive-genus curve in the selected subfamily `X_0(2p)` with
  `p=3,5,7,...`; `X_0(11)` already has genus one. Finally, genus is a vector-
  space dimension, not the observable-fiber quantity refuted in MISTAKE-231,
  and class number greater than one does not alone force a cuspidal theta
  component: genus characters can account for the class variation.
- **Correct framing:** the finite identities in S219 survive as independent
  classical examples: for `x^2+xy+2y^2`,
  `r(n)=2 sum_(d|n)(d/7)`; for discriminant `-23`, the stated class average is
  Eisenstein and the principal-minus-nonprincipal theta series is a nonzero
  weight-one CM cusp form. The cusp-count formula, genera of `X_0(12)` and
  `X_0(14)`, the displayed coefficients and Atkin--Lehner signs of `f_14`,
  and the standard Rankin--Selberg factorization are likewise classical facts.
  They have no demonstrated LRC consequence. HYP-8880 and the LRC portions of
  S219/S220 are **REFUTED**; retain the artifacts only as historical records
  and do not route agents to them as an LRC strategy without a typed,
  predicate-preserving construction.

## MISTAKE-232 (2026-07-21, audit of S99 / HYP-8876) -- a scale-and-clock word analogy was called an exact proof-shape bridge, and a composite modulus was assigned a Paley spectrum

- **What was claimed:** S99 identified GMC(2) dilation/Frobenius with
  THM-2057's scaled clock witnesses, wrote
  `Z/pZ <-> Z/12aZ,Z/14aZ` and `Frobenius <-> modular-orbit periodicity`, and
  said that each of the clock moduli `{7,13,14}` carries a Paley
  `sqrt(p)` spectrum. Its LRC row also used `84a=12a*7a` and invoked the
  transitive-core/Paley/class-number entropy syntheses as explanation.
- **Why it is wrong:** the script constructs Paley objects only at the prime
  moduli `7` and `13`; `14` is not a prime power, has no finite field
  `F_14`, and no Paley graph or tournament was constructed. Printing
  `14=2*7` is not a spectral calculation. The displayed divisibility identity
  already fails at `a=2`: `84a=168`, while `(12a)(7a)=336`. The exact
  statement is `lcm(12a,14a)=84a`. More fundamentally, Frobenius is a ring
  endomorphism in characteristic `p`, whereas the composite LRC clocks supply
  only additive residue orbits. For example on `Z/12Z`, the proposed
  seventh-power analogue is not additive:
  `(1+1)^7=8 (mod 12)` but `1^7+1^7=2 (mod 12)`. No source-to-target map,
  preserved predicate, or sidecar turns the shared words “scale” and “clock”
  into a proof transfer. MISTAKES-227, -228, -230, and -231 independently
  block the transitive-core, Paley-LRC, and zero-entropy explanations.
- **Correct framing:** two exact, independent spectral facts survive. For the
  Paley tournament at `7`, its skew matrix satisfies
  `S^2=J-7I`, hence has spectrum `0,(+i sqrt(7))^3,(-i sqrt(7))^3`. For the
  Paley graph at `13`, `A^2+A=3I+3J`, hence the spectrum is
  `6,((-1+sqrt(13))/2)^6,((-1-sqrt(13))/2)^6`. THM-2057's scaled missing-
  clock proof also survives independently. “Scale, then close on a clock” may
  be kept only as a heuristic phrase; HYP-8876 is not verified and yields no
  theorem or transfer.

## MISTAKE-231 (2026-07-21, audit of HYP-8875 / S218) -- four unrelated fiber sizes were called one repo-wide entropy invariant

- **What was claimed:** S218 defined
  `H_arith(X|L)=log2|{X':L(X')=L(X)}|` and identified it with class-group
  depth, tournament reconstruction, continued-fraction partial-quotient size,
  and LRC/GMC moment depth. It then asserted that every rigid extremum has zero
  arithmetic entropy and every hard object has positive entropy.
- **Why it is wrong:** the displayed quantity is at most a finite-fiber
  Hartley information after fixing an ambient universe, equivalence relation,
  and observable. It is not Shannon entropy without a probability law, and
  changing `L` can force any chosen object to have fiber size one or the whole
  universe. For the two isomorphism classes of three-vertex tournaments, a
  constant `L` gives the transitive class a fiber of size two, while the
  identity observable gives every class a singleton fiber. Thus the proposed
  rigid-versus-hard law can be reversed by changing the observable. The four
  examples do not use one map or even one type of fiber.
  In particular:
  - the tournament claim is already false at `n=5`: the regular score
    `(2,2,2,2,2)` has one isomorphism class, whereas the score
    `(1,2,2,2,3)` has three. Those three classes are not cospectral; their
    adjacency characteristic polynomials are
    `x^5-4x^2-4x-1`,
    `(x^2+x+1)(x^3-x^2-3)`, and
    `(x^2+x+2)(x^3-x^2-x-1)`. MISTAKE-227 also blocks identifying a
    translated transitive score vector with the LRC AP relation carrier;
  - a finite continued-fraction partial-quotient geometric mean is not an
    entropy and has no maximum (`[0;N]` makes it arbitrarily large). Every
    nonterminal finite CF prefix has continuum many real extensions, while a
    complete rational expansion has a singleton fiber, so the displayed
    Hartley formula cannot distinguish golden prefixes from `[0;13,14]`;
  - bounded LRC danger count gives exact inclusion--exclusion by depth `13`,
    while `B5>0` is only a sufficient certificate, not a uniform depth-five
    classifier. Explicitly, the distributions with weights
    `(1,15,15,1)/32` on `{0,2,4,6}` and `(6,20,6)/32` on `{1,3,5}` have the
    same moments through degree five, by the sixth finite-difference identity,
    but masses `1/32` and `0` at zero. THM-1790 gives an unbounded-in-degree
    GMC moment-depth floor, not an infinite entropy; and
  - the genus routine labels its formula a heuristic and silently falls back
    when divisibility fails. It is a finite table, not a proof of genus theory
    or of invisibility to every congruence test. Its prose also calls the
    within-genus part “odd,” but the fundamental discriminant `D=-39` has
    `h=4`, two genera, and even within-genus fiber size `2`.
- **Correct framing:** for a fixed finite universe and quotient, define the
  explicit Hartley fiber size `H_0=log2|L^{-1}(L(X))|`; add a probability law
  before using Shannon language. The exact survivors are the `n<=5`
  tournament fiber census, transitive-score uniqueness, elementary
  score-distribution entropy, the displayed reduced-form tables, and the
  qualitative contrast between finite LRC inclusion--exclusion and no
  degree-uniform GMC cutoff. These are separate examples of information loss,
  not one invariant and not an LRC theorem. HYP-8875 is retracted as a unified
  mathematical claim.

## MISTAKE-230 (2026-07-21, audit of HYP-8870 / S217) -- inverse form classes were counted as distinguishable rational-prime outcomes

- **What was claimed:** S217 treated the `h(D)` proper classes of positive
  binary quadratic forms as `log2 h(D)` hidden bits beyond the Legendre symbol.
  For `D=-23` its script assigned split rational primes among the three forms
  with counts `[2,9,0]`, then transferred class-number-one rigidity to the
  THM-2053 LRC gate.
- **Why it is wrong:** proper inverse classes `(a,b,c)` and `(a,-b,c)` represent
  exactly the same rational integers by `y -> -y`. The two nonprincipal
  `D=-23` forms therefore represent the same displayed primes; choosing
  `reps[0]` manufactures the asymmetric `[2,9,0]` count. A rational split
  prime determines an inversion orbit of ideal classes, not an oriented class.
  Moreover `log2 h` is not an entropy until a random variable and distribution
  are specified. For oriented prime ideals one can formulate an Artin-class
  distribution and invoke Chebotarev; for rational primes the inversion
  quotient is different. “Invisible to any local test” also overstates
  “not determined by the quadratic Legendre symbol”: Frobenius in the Hilbert
  class field is itself local prime data. The finite search bound is evidence,
  not a proof of the global representation criterion. Finally, MISTAKE-229
  already shows that THM-2053 has no discriminant `-7` form to which any of
  this could transfer.
- **Correct framing:** the reduced-form and class-number tables are useful
  classical checks. With a precisely oriented prime-ideal sample, class-group
  size can support a Hartley/Shannon analogy. It has no demonstrated LRC
  consequence. The S217 use of `HYP-8870` also collided with the normal-fan
  program now correctly numbered HYP-8871; the S217 claim has no canonical
  hypothesis ID and is retracted.

## MISTAKE-229 (2026-07-21, audit of HYP-8865 / S216) -- the tangent-disk gate was assigned a nonexistent Heegner discriminant

- **What was claimed:** S216 called THM-2053's residual the short vectors of
  one binary quadratic form, inferred discriminant `-7` from `LRC(14)=2*7`,
  used `h(-7)=1` to declare every residual plane a single rigid class, and
  identified rank-versus-Euler with isotropic-versus-anisotropic directions.
- **Why it is wrong:** THM-2053 compares thirteen linear determinant forms
  `a z_i-b u_i` with the chosen-coordinate Euclidean norm `a^2+b^2`. The
  failure locus is a **union of 26 tangent disks**, not a representation set
  for one quadratic form. The Euclidean form has discriminant `-4` and changes
  covariantly with the chosen saturated basis; neither the number `14` nor a
  Paley spectral factor supplies a plane discriminant `-7`. At the prime `p`,
  the discriminant `-p` has Legendre symbol zero, so the script's computation
  of `(-1/p)` did not test its stated local criterion. An extra bounded
  relation is a linear annihilator condition on all columns, not an isotropic
  parameter direction, and failure of the sufficient gate does not imply
  either unsafety or Euler survival.
- **Correct framing:** retain the exact wedge identity and the classical
  class-number table as separate facts. The useful LRC geometry is the
  intersection of a primitive lattice and positivity cone with 26 explicit
  tangent disks, considered up to saturated `GL_2(Z)` basis changes. THM-2055
  compresses the determinant sidecar exactly to the signed column polygon's
  normal fan: hull vertices suffice for the maximum, and each owner cone meets
  one tangent disk. Non-hull runner data, pair sums, and endpoint ownership
  remain necessary sidecars. Any
  arithmetic-form compression must first construct an actual invariant form
  from the full column configuration and prove that it preserves the gate;
  Heegner or Paley language alone supplies no such map. See MISTAKE-230 for the
  independent inverse-class error in S217's entropy continuation.

---
## MISTAKE-228 (2026-07-21, audit of HYP-8860 / S215) -- the odd-prime Paley dichotomy was promoted to an LRC periodic table

- **What was claimed:** S215 assigned primes `2,3,5,7,11` intrinsic LRC roles
  from Paley spectra, identified adjacency eigenvalues with Gauss sums, and
  treated `p mod 4` as a tight/slack law. It called `2` part of the same Paley
  theorem, merged the `Phi_3` and `Phi_6` factors, and linked the Paley-7 field
  directly to the real cubic cap field.
- **Why it is wrong:** the verified classical statement is for odd primes.
  For `p=3 mod 4` the nonprincipal Paley-tournament eigenvalues are
  `(-1+-i sqrt(p))/2`, while the Gauss sum is `i sqrt(p)`; for `p=1 mod 4`
  the Paley graph also has principal eigenvalue `(p-1)/2`. The script does not
  test `p=2`. `Phi_3=x^2+x+1` and `Phi_6=x^2-x+1` need a sign twist, and
  `Q(sqrt(-7))` is not the real cubic field of discriminant `49`. The roles
  “5 = Fibonacci foil” and “11 = rank” supply no map, and multiples of `11`
  are not a forced speed coordinate.
- **Correct framing:** retain the quadratic-residue Cayley relation and the
  odd-prime Paley graph/tournament spectral formulas as tournament background.
  Treat the small-prime table as a modulus-selection mnemonic only. It proves
  no LRC tightness, slackness, rank, or transport statement.

---
## MISTAKE-227 (2026-07-21, audit of HYP-8855 / S214 and HYP-4552 / S25) -- a finite-index AP chain frame was called a saturated relation-lattice basis

- **What was claimed:** the rows
  `d_k=(k+1)e_k-k e_(k+1)` were called a `Z`-basis of
  `L(AP)=ker_Z(1,2,...,12)`, their tridiagonal Gram matrix was assigned to the
  whole lattice, and they were identified with THM-2052's eleven private
  anchor-star relations. The resulting rank coincidence was promoted to a
  common AP/tournament nullcone and to the claim that the rank and Euler
  branches meet at the AP. The 60 minimal vectors were also called additive
  energy.
- **Why it is wrong:** the `d_k` span only a sublattice `D` of index
  `11!=39,916,800`. Exactly,
  `det Gram(D)=1,035,678,099,456,000,000`, while
  `disc L(AP)=650`; their square-root ratio is `11!`. The saturation contains
  the norm-three Schur rows. THM-2052 instead concerns a bounded signed code
  on thirteen coordinates and its private rows depend on two chosen anchors;
  no map to this twelve-coordinate lattice was given. The kernel rank is
  `|S|-1=11`, not `|S|-2`. The 60 vectors are the signed off-diagonal Schur
  triples, not four-term additive energy. Tournament scores `0,...,11` are a
  translation of speeds, and LRC is not translation invariant. Reversal of
  phase is universal, not an AP-fixed configuration theorem. THM-2053 also
  shows rank twelve is not the only finite terminal.
- **Correct framing:** `D` is a useful finite-index path/Jacobi frame and
  `L(AP)=sat(D)` has discriminant `650`; saturation is exactly what restores
  the short Schur relations that the path frame misses. The AP, transitive
  order, palindromic pair sums, and Schur census remain diagnostic views, not
  one object or an LRC descent. Before promoting a rank match, compute the
  ambient type, explicit map, saturation index, and discriminant.

---
## MISTAKE-226 (2026-07-21, codex audit of HYP-8840) -- unanchored additive energy, LRC Fourier volume, and GMC radial moments were promoted to one constant-term bridge

- **What was claimed:** S211 promoted nine finite AP energy maxima to an
  all-order Wall-A theorem, called THM-730 its `m=2` case, identified the
  unanchored additive-energy lattice with the full LRC relation lattice, used
  separated-digit tensorization as a rank reduction for arbitrary cores, and
  treated a geometric model weight as the GMC radial kernel. It then called
  the one-dimensional LRC maximum a minimax saddle linking GMC volume to Euler
  positivity.
- **Why it is wrong:** `CT[P^m Pbar^m]` counts only affine relations
  `k dot v=0` with augmentation `sum k_i=0`; the LRC Fourier formula uses the
  full linear lattice. THM-730's Schur relation `a+b=c` has augmentation one,
  so it is not an unanchored two-energy relation. The script checked only nine
  `(N,k,m)` cases. Translation preserves every unanchored energy but changes
  Schur data: `{1,...,5}`, `{2,...,6}`, and `{1,3,5,7,9}` have the same
  unanchored energy at **every** order, yet their LRC maxima are `1/6,1/4,1/2`.
  Adjoining observer zero gives anchored `E_2=146,130,106` and Schur counts
  `10,6,0`. Separated-digit factorization
  requires a no-carry scale and therefore removes the carry relations that
  make an AP extremal; it does not factor an arbitrary 13-speed core. The
  geometric weight was invented, and a raw box-truncated sinc sum is not the
  Fejer-regularized measure in THM-2047.
- **Correct framing:** retain the universal identity
  `E_m(S)=CT[P^m Pbar^m]`, the single-character CT seed, conditional no-carry
  tensorization, the AP cyclotomic factorization, and the exact volume ceiling:
  positive Haar measure certifies a strict exit, while Euler characteristic
  also sees a tight isolated witness. A faithful energy comparison must adjoin
  the stationary observer `0`; then a linear relation lifts by
  `k_0=-sum k_i`, and THM-730 is the `a+b=c+0` sector of anchored `E_2`.
  Equivalently, retain the mixed table
  `M_(r,s)=CT[P(z)^r P(z^(-1))^s]`: it selects augmentation `r-s`, THM-730 is
  exactly `M_(2,1)`, and every full LRC relation appears in some bidegree.
  Whether this augmentation sidecar helps Wall A is OPEN. Do not call the
  maximum of a lower envelope a game/Morse saddle without an actual
  two-player or variational correspondence.

---
## MISTAKE-225 (2026-07-21, codex audit of HYP-8835) -- skew parity, tournament games, torus topology, and three sign actions were promoted to one antisymmetry theorem

- **What was claimed:** S210 asserted that antisymmetry alone forces odd
  tournament-game support; that a pure optimal strategy exists iff the
  tournament is transitive; that intransitivity is precisely toroidal
  recurrence; that the solid-bagel cutting deficit is a reduced Euler term of
  `T^2`; and that matrix complement, torus inversion, and Vandermonde parity
  are one involution.
- **Why it is wrong:** odd optimal support also uses the tournament block's
  mod-2 form `I+J`, which is nonsingular in even order. A four-vertex
  tournament with vertex `0` beating all and `1->2->3->1` is intransitive but
  has pure optimum `e_0`; replicator dynamics even satisfies
  `x_0'=x_0(1-x_0)`. Antisymmetry does not force zero column sums. The regular
  three-cycle's conserved product comes from a positive kernel vector, and its
  positive levels in `Delta^2` are circles `S^1`, not a surface `T^2`. The
  cutting bagel is `D^2 x S^1`, while the `1,2,1` Morse count belongs to its
  boundary. MISTAKE-222 already blocks the sequence/Euler jump. Finally the
  three `Z/2` actions have no supplied equivariant map. The RK4 check was
  nonsymplectic and rounded a nonzero return distance to `0.0000`.
- **Correct framing:** Fisher--Ryan's tournament-game support is unique and
  odd; the parity proof needs both skew singularity and the tournament mod-2
  block. Pure optimum is equivalent to a Condorcet winner. If `Mp=0`, then
  `d/dt sum p_i log x_i=0`, so a positive equilibrium vector yields the
  monomial invariant `product x_i^(p_i)`. Keep game, Morse, and dynamical
  saddles distinct; keep solid and boundary tori distinct; and require an
  equivariant map before transferring any sign representation.

---

## MISTAKE-224 (2026-07-21, codex audit of HYP-8830) -- a Fourier relation lattice was identified with a toric layer poset, and a thickened safe set with an arrangement complement

- **What was claimed:** S209 called `{k in Z^n:k.v=0}` the layers of a finite
  toric arrangement, identified `|G_delta|` with ordinary toric-complement
  volume and its Fourier series with an arithmetic-Mobius sum, and promoted a
  `B=2` relation count to Betti/Mobius mass. It also called LRCMod a finite-
  field characteristic-polynomial engine and braid cohomology a per-tournament
  invariant.
- **Why it is wrong:** the lattice indexes Fourier characters surviving
  integration over the orbit subtorus; toric layers are connected components
  of intersections of a **finite** character list. THM-1820's series is
  infinite, sinc-weighted, and `delta`-dependent. For `S={1}`,
  `|G_delta|=1-2delta`, whereas the ordinary circle complement of the
  character kernel has Haar measure `1`. The script computed no layer poset,
  Mobius function, arithmetic Tutte polynomial, or LRCMod equivalence. Its AP
  count is cutoff-dependent over exactly 72 triples; the natural circle-
  hypertorus union has four points for `(1,2,3)` but six for `(2,3,4)`.
- **Correct framing:** `L_v={k:k.v=0}` is the annihilator of
  `H_v=image(t |-> (v_jt))`, and `G_delta=phi_v^(-1)(F_delta^n)`. The full
  embedded `L_v` determines primitive `v` up to scale/relabeling; only an
  abstract, truncated, or summarized lattice is automatically lossy. The
  useful finite toric object has characters
  `X_S={(v,+1),(v,-1)}` on `(t,delta)`, with walls
  `vt +/- delta in Z`. Its selected inequality side, height, owner, sign, and
  paired deletion are essential. Opposite-sign wall determinants are `v+w`,
  recovering the pair-sum ruler. THM-2047 now proves this exact carrier and
  paired-deletion operation; localization forcing an AP core remains open and
  is not a proved Wall-A consequence.

---

## MISTAKE-223 (2026-07-21, codex audit of HYP-8825) -- a braid-flat factorization and common determinant syntax were promoted to hyper-Bessel and Euler-characteristic factorizations

- **What was claimed:** S208 called cake, bagel, full g-bonacci, and transitive-
  tournament counts one Zaslavsky/Mobius inversion; identified
  `bagel(n)-cake(n)=T_n-1` and a finite-shadow first deficit as the same reduced
  Euler characteristic; called the Vandermonde's within-block vanishing order
  a flat codimension; and asserted that braid-arrangement localization factors
  the NC2 hyper-Bessel boundary into block factors. It also called general
  `Phi_(p,q)` Laguerre--Polya open and described a duplicate recurrence
  comparison as a companion-determinant/topology verification.
- **Why it is wrong:** the `n!` chambers are components of the **real** braid
  complement; the complex complement used in the note is connected. For a
  `3+2` coincidence partition, the ordinary Vandermonde vanishes to order
  `C(3,2)+C(2,2)=4`, while the flat codimension is `(3-1)+(2-1)=3`.
  Localization gives within-block braid factors and a nonzero cross-block unit;
  no map from that polynomial expansion to THM-2017's univariate channel limit
  `Phi_(p,q)` was stated. THM-2023 had already proved every positive-integer
  `Phi_(p,q)` is type-I Laguerre--Polya by Gauss multiplication and the
  Baricz--Singh theorem. A bagel is a solid torus, not `T^3`; no common complex,
  Mobius function, or valuation was supplied for its regions and the shadow
  lattice. The script compared two implementations of the same recurrence and
  computed neither a companion determinant nor an Euler characteristic.
- **Correct framing:** retain four separate exact shadows: the Vandermonde is
  the braid defining polynomial; its **real** chambers are labelled total
  orders/transitive tournaments; near a coincidence flat it has a leading
  within-block/transverse product with exponent equal to vanishing order; and
  an explicit companion matrix satisfies
  `det(I-x*M_g)=1-x-x^(g+1)`. THM-2033 applies to its special moment matrix, not
  every NC2 moment. The cake--bagel identity has the more direct decomposition
  `bagel(n)=cake(n+1)-2`, so
  `bagel(n)-cake(n)=(T_n+1)-2`. The possible topology bridge remains an open
  source-to-target valuation test under MISTAKE-222, not a theorem.

---

## MISTAKE-222 (2026-07-21, codex audit of HYP-8820) -- a shared binomial ambient array and a matching minus-one offset were promoted to a geometric and cross-problem mechanism

- **What was claimed:** cake, bagel, Moser, and Fibonacci were described as
  ordinary truncated rows or diagonals of “one Pascal triangle”; the identity
  `bagel-cake=T_n-1` was said to prove that the torus hole **is** the deficit-one
  boundary in a g-bonacci shadow lattice and to place JC(2) and LRC(14) on the
  same scaffold. A broken polygonal-skip calculation was also printed as if it
  agreed with the cited indexing.
- **Why it is wrong:** caterer and cake are prefixes of the `n`th Pascal row.
  Moser is the even-column functional `C(n,0)+C(n,2)+C(n,4)`, equivalently the
  first five terms of row `n-1`; bagel is the weighted functional
  `C(n,3)+2C(n,2)+2C(n,1)`, equivalently `cake(n+1)-2`. Thus they are not all
  the same unshifted prefix operation. A common ambient array is a useful
  coordinate system, not a map between the counted objects. The same scalar
  offset `-1` supplies no bijection, operation, or preserved predicate, and the
  JC/LRC analogy supplies none either. The removed polygonal-skip code produced
  negative values and contradicted its own “first deviation” label.
- **Correct framing:** the four binomial formulas, Fibonacci diagonal identity,
  full-rank gap-diagonal g-bonacci identity in klein-S313's explicit indexing,
  and
  for `n>=1`, `bagel(n)=cake(n+1)-2` and
  `bagel(n)-cake(n)=C(n,2)+C(n,1)-1=T_n-1` are exact. Treat them as a binomial-
  reading atlas. Finite-rank shadows depart at their first deficit. The
  shadow-lattice/torus link remains a useful test only if both defects are
  derived from one explicit boundary-cell or valuation operation; JC and LRC
  links additionally require a preserved predicate and loss ledger.

---

## MISTAKE-221 (2026-07-21, codex audit of HYP-8815) -- a finite LRC witness scan and an AP/Fibonacci contrast were promoted to a necessary autocorrelation characterization

- **What was claimed:** HYP-8815 called a scan over rationals `a/q` with
  `q<=Qmax` an exact computation of `M`; inferred that every hypothetical
  LRC(14) counterexample must be near-AP, anti-golden, organized around one
  CF-blocking far element, and beat the AP by a higher-order autocorrelation
  mechanism; and used smaller THM-731 discrepancy in the wrong direction.
- **Why it is wrong:** a denominator-truncated scan computes only
  `L_Q(S)<=M(S)` unless a completeness theorem covers every maximizer.
  The exact strict witness is
  `L_1200({1,...,12,5460})=92/1197 < M=420/5461`.
  THM-730 proves only the AP's unique Schur-triple maximum. THM-731 is
  peel-dependent and gives
  `L_cert=(6/7)|G'|-sqrt((6/49)disc_v)`, so **smaller** discrepancy
  strengthens safety. A few loose Fibonacci packets prove no uniform
  exclusion, near-AP necessity, continued-fraction condition, or ownership of
  every covering modulus by one far speed.
- **Exact computational repair:** the pair-sum theorem
  [THM-1002](theorems/THM-1002-pair-sum-denominator-bound-and-the-bounded-gap-case.md)
  puts every maximizer at `t=p/(v_i+v_j)`, including the self-pair case.
  Enumerating every numerator on every pair-sum ruler computes `M` exactly;
  alternatively `Q>=2 max(S)` is a sufficient complete denominator cutoff.
  Do not filter to coprime numerators unless all divisor rulers are also
  present: `S={1,4,5}` attains its maximum at `2/6=1/3`, while `3` is not
  a pair sum. The repaired script validates inputs and computes all fifteen
  displayed rows exactly; all are safe. This finite bank is not a global
  minimizer theorem.
- **Correct structural framing:** after gcd normalization, a counterexample is
  primitive Cover14 with `M<1/14`. THM-1017 excludes a maximum-deletion core
  `d{1,...,12}`, so THM-730 gives that core Schur count `T<=65`. Since a
  counterexample has lonely-time measure zero, THM-731 forces
  `disc_v>=6|G'_{~v}|^2` for every peel. Anti-golden, near-AP,
  Fibonacci-foil, one-far-blocker, and joint-order autocorrelation language is
  a hostile-control program, not a characterization.

  `L_1200=92/1197 < M=420/5461` at `t=420/5461`.

  More starkly, for `Q>=14`, dilate the deep well by `lcm(2,...,Q)`. Every old
  sampled phase then has value zero while dilation invariance leaves
  `M=14/183`. A bounded lower bound above `1/14` does safely exclude one row,
  but a finite candidate list cannot prove a global minimizer or ranking.

  Primitivity is a WLOG normalization, not a literal condition on every
  counterexample. The old random generator did not validate cardinality and
  distinctness, although its six displayed seed-2026 rows happen to be valid.
  THM-730 counts Schur triples, not four-term additive energy; Fibonacci is not
  minimal for either statistic. THM-731 is peel-dependent and gives
  `L_cert=(6/7)|G'|-sqrt((6/49)disc_v)`, so **smaller** discrepancy strengthens
  safety—the claimed direction was backwards. Cover14 distributes divisibility
  among all runners, not necessarily a twelve-runner core plus one far blocker;
  THM-1017 excludes only a maximum-deletion core `d{1,...,12}` in the relevant
  normalized regime.
- **Correct framing:** exact `M` computation uses every numerator, including
  non-coprime ones, on the THM-1002 pair-sum rulers. The corrected script
  validates its inputs and computes all fifteen rows exactly; every one is safe,
  and none is a disproof. The rigorous necessary kernel is: after gcd
  normalization a counterexample is primitive Cover14 with `M<1/14`; its
  maximum-deletion core is not `d{1,...,12}`, so THM-730 gives Schur count
  `T<=65`; and `L=0` with THM-731 forces
  `disc_v>=6|G'_{~v}|^2` for every peel. THM-2048 now refines the latter by
  the exact fiber-occupancy tax
  `{7vmu}(1-{7vmu})/(7v^2)` and its THM-732 integer obstruction. Anti-golden,
  near-AP, Fibonacci-foil,
  one-far-blocker, and joint-order-autocorrelation language remains a useful
  hostile-control program, not a characterization.
---
## MISTAKE-215 (2026-07-21, codex audit of HYP-8795 / THM-2040) -- a prime-local minimum-face normalization was promoted to a global common-factorial/Vandermonde factorization

- **What was claimed:** dividing every NC2 moment by a purported common
  factorial `(p*A0)!` was said to leave a bounded residue and to make general
  NC2 noncancellation equivalent to one (possibly confluent) Vandermonde.
- **Why it is wrong:** `(p*A0)!` in THM-2022 is attached to one deliberately
  chosen amplified moment and its lowest balanced face. It divides all channel
  factorials there, but off-face quotients remain nontrivial; it does not peel
  every general moment to a universal determinant. The accompanying S91
  computation itself gives `E[P^m]/m! = 2, 7/2, 17/3, ..., 1441729/40320`
  for `P=Z+(1+s)+Zbar`, an increasing sequence through the displayed range,
  contradicting its printed "bounded" conclusion. A moment is a scalar sum of
  channel monomials, whereas the Vandermonde in THM-2033 is a determinant of a
  special moment matrix; division by one scalar factorial cannot identify the
  two in general. The claimed positive-zero / Laguerre--Polya / Paley iff chain
  also has no logical basis and repeats the type error in MISTAKE-214.
- **Correct framing:** retain two valid statements. On the exact symmetric
  monomial wall `P=a Z^p + beta*s^(p/2) + c Zbar^p` (necessarily even `p`),
  every balanced channel really has radial degree `mp/2`, so the special
  central-trinomial factorization is exact. In full NC2, THM-2022 provides a
  **prime-local initial-form normalization** only: at order `p*m0`, division by
  `(p*A0)!` kills higher/carry layers modulo a good prime and leaves the whole
  lowest-face residue `Q^p`. That residue is generally a face constant term,
  not a Vandermonde and not a Paley spectrum.

---

## MISTAKE-214 (2026-07-21, codex audit of HYP-8785/8790) -- repeated Vandermonde nodes were identified with repeated tournament scores and hence with a regular/Paley tournament

- **What was claimed:** HYP-8785 asserted that equal radial channel degrees on
  the NC2 central wall are equal scores in the tournament expansion of the
  Vandermonde, hence a regular, doubly regular, or Paley tournament. HYP-8790
  then promoted the true central-trinomial channel-weight identity to claimed
  equivalences among NC2 noncancellation, Laguerre--Polya zero geometry, and the
  Paley critical-line spectrum.
- **Why it is wrong:** in
  `V(x)=sum_T sgn(T)*x^score(T)`, the `x_i` are Vandermonde **node variables**
  and `score_i(T)` are their **exponents**. Setting two node values equal does
  not set two score exponents equal; it causes cancellation among the evaluated
  monomials. The type mismatch already fails at the smallest tied core: for
  symmetric `m=2` there are two equal-degree channels, while no regular
  tournament exists on two vertices (the scores are `0,1`). More generally
  tied channel counts have arbitrary parity, whereas regular tournaments need
  odd order, doubly regular tournaments require order `3 mod 4`, and Paley
  tournaments exist only at special prime-power orders. The identity
  `sum_i m!/(i!^2(m-2i)!)=[x^0](1+x+x^-1)^m` is valid, but it does not imply the
  stated Paley/free-probability/zero-set equivalences.
- **Correct framing:** retain THM-2033's determinant/Vandermonde identity and
  HYP-8790's central-trinomial coefficient identity. Confluence means repeated
  **nodes** and calls for divided differences or derivative rows; regularity is
  a property of a tournament's **score vector** and is not forced. Paley and
  free-probability parallels may be used as analogies, not iff reductions or an
  identification of the NC2 wall. THM-2022 supplies the actual tied-face
  certificate: after common-factorial normalization the whole face survives as
  `Q^p`, with no tournament-score interpretation required.

---

## MISTAKE-213 (2026-07-21, codex proportional-NC2 audit) -- finite toral-zero recurrence was treated as necessary to combine an exact factorization with EMP

- **What was claimed:** the first THM-2021/HYP-8771 write-up said that
  no-consecutive-zero for the toral factors `A_m(kappa)` was too weak for NC2,
  because it did not make `A_m(kappa)` eventually nonzero and therefore could
  not force `L(b^m)=0` eventually from
  `E[P^m]=A_m(kappa)L(b^m)=0`. It promoted finite recurrence of a fixed toral
  zero to the load-bearing proportional-slice target.
- **Why it is wrong:** the supports were combined in the wrong direction. If
  `b!=0`, EMP says `L(b^m)!=0` for **every sufficiently large** `m`. Moment
  nullity would therefore force `A_m(kappa)=0` for every sufficiently large
  `m`. Any theorem showing merely that `A_m(kappa)` is nonzero arbitrarily far
  out contradicts this. In particular, no-consecutive-zero is already more
  than enough. For the symmetric Legendre factors the three-term recurrence
  forbids consecutive zeros (with `kappa=1/4` handled directly), so it gives a
  self-contained full symmetric proportional closure. THM-2018's root-of-unity
  EGF proves the needed non-eventual-zero statement for every charge pair.
- **Correct framing:** for a factorization `M_m=A_m L_m` with `L_m` eventually
  nonzero, NC2 needs only that `A_m` is **not eventually zero**. Finite zero
  recurrence is a strictly stronger sequence-profile theorem, interesting in
  HYP-8771 and announced for Legendre polynomials, but not an NC2 dependency.
  Keep three distinct predicates separate: no consecutive zeros; unbounded
  nonzero support; finite zero support. The first implies the second, and the
  second is sufficient here. See THM-2018, corrected THM-2021, and HYP-8771.

---

## MISTAKE-212 (2026-07-21, codex audit of HYP-8772) -- channel-tournament transitivity was promoted from a sufficient dominance certificate to an iff characterization of NC2 noncancellation

- **What was claimed:** the reservation version of HYP-8772 said that NC2
  noncancellation is equivalent to transitivity/source-dominance of the
  channel-degree tournament, and that the moment nullcone *is* a tournament
  nullcone on its channels. Resonance was identified with a tied/regular
  tournament and therefore with the only possible cancellation wall.
- **Why it is wrong:** THM-2021 gives exact noncancellation inside a fully tied
  central channel family. Take `b=s(s-2)`, `a=s-2`, `c=s(s-2)`, and
  `P=Z*a(s)+b(s)+Zbar*c(s)`. Then `h=s*a*c=b^2`, so every return channel in
  `E[P^m]` has the identical radial word `b^m` and identical factorial degree:
  there is no asymptotic source. Nevertheless `E[P]=L(b)=0` and
  `E[P^2]=S_2(1)L(b^2)=3*8=24`, so this two-sided polynomial is immediately
  non-null. More generally THM-2021 proves NC2 on this tied proportional slice
  outside a countable exceptional set. Also, a dominance relation with exact
  ties is not a tournament until an extrinsic tie-breaker is imposed.
- **Correct framing:** a strict channel source/transitive dominance order is a
  **sufficient certificate** for noncancellation and explains why THM-2017's
  degree-gap argument works. Its failure identifies where that certificate
  stops, not where NC2 fails or can fail. In a tied core, phase-sensitive
  algebraic information can still prove noncancellation (Legendre recurrence
  in THM-2021, the all-charge EGF in THM-2018, and Sheffer/resultants in
  HYP-8769). The channel
  tournament remains a useful regime classifier, but the claimed iff and
  "moment-nullcone = tournament-nullcone" statements are withdrawn. See
  THM-2021, HYP-8771, and the correction banner in the HYP-8772 reflection.

---

## MISTAKE-218 (2026-07-21, codex audit of THM-1979/2013) -- score variance was promoted from a cyclicity coordinate to the whole structural axis

- **What was claimed:** THM-1979 put `σ²=Var(scores)` in
  `[0,(n²−1)/12]` for every `n`, called zero the regular center, and said
  fiber size, strong fraction, and modular-prime fraction vary monotonically
  opposite to `σ²`. THM-2013 consequently wrote
  `c₃/c₃,max=1−σ²/σ²_max` uniformly in `n`. The limit discussion
  identified tournamenton space with its degree functions.
- **Why it is wrong:** for even `n`, integral scores with half-integral mean
  force minimum variance `1/4`, not zero, and the maximum cyclic count is
  `n(n²−4)/24`. More fundamentally, the frozen S198 census refutes the
  monotonic structural claim: at `n=6`, variance `5/4` supports both size-3
  all-strong score fibers and size-1 reducible fibers; at `n=7`, maximum fiber
  size 47 occurs at both variances `4/7` and `10/7`, while the zero-variance
  fiber has size 3. A tournamenton is not determined by its degree function;
  many structured regular tournamentons share `d_W=1/2` with the quasirandom
  tournamenton.
- **Correct framing:** set `ε_n=0` for odd `n` and `1/4` for even `n`, and
  `σ²_tr=(n²−1)/12`. Then
  `c₃=n(n²−1)/24−(n/2)σ²` and
  `τ=c₃/c₃,max=(σ²_tr−σ²)/(σ²_tr−ε_n)`.
  The variance determines the cyclic shell exactly but forgets score shape and
  within-fiber structure. The finite score map and limiting degree map are
  projections with rich fibers, not complete parameterizations. THM-1979,
  THM-2013, their scripts/results/reflections, the hypothesis index, and the
  session log are corrected.

---

## MISTAKE-217 (2026-07-21, codex audit of THM-2013/2016) -- signed Redei data was labeled as the invariant magnitude `|R|`

- **What was done:** the continuum-coordinate scripts grouped tournament
  classes by `signed_redei(A,n)` while their output and THM-2013/2016 called
  the coordinate `|R|`.  THM-2016 then attributed signed-key resolution counts
  `36→44` and `41→50` to the absolute invariant.
- **Why it is wrong:** the sign changes under the parity gauge of a relabeling;
  the isomorphism invariant in THM-1966 is its magnitude.  Signs are mixed in
  every audited hot shell, so taking the absolute value materially changes the
  fibers.  It is not a notation-only error.  Also, “44 of 47 resolved” was
  paraphrased as “3/47 classes survive,” although its unresolved groups were
  twin pairs rather than three classes.
- **Correct framing:** group by `abs(signed_redei(A,n))`.  At `n=7`, the exact
  rows `(base; +score; +4-profile; +4&5; +score+4&5)` are
  `c3=12: (28,28,28,41,41)`, `c3=11: (30,46,46,50,50)`, and
  `c3=13: (13,13,13,14,14)`.  The final unresolved fibers are respectively
  six, two, and one twin pairs.  THM-2013/2016, both scripts/results, the
  hypothesis index, reflection, and session log are corrected.

---

## MISTAKE-220 (2026-07-21, codex audit of THM-2016) -- an SCC sum was bounded by its largest summand

- **What was claimed:** the proof of the reducibility ceiling used
  `c₃(T)=Σ_SCC c₃(S)≤c₃_max(largest SCC)`.  It then correctly stated that the
  largest possible `c₃` of a reducible `n`-tournament is `c₃_max(n−1)`.
- **Why it is wrong:** the displayed inequality drops every SCC except the
  largest.  The order-join of two directed 3-cycles is a six-vertex
  counterexample: its two SCCs contribute `1+1=2`, while
  `c₃_max(3)=1`.  The claimed ceiling survives, but not that inference.
- **Correct framing:** write `f(n)=c₃_max(n)`.  THM-2000 gives
  `f(n+1)-f(n)=T_floor(n/2)`, so the increments of `f` are nondecreasing.
  Transferring a vertex from a smaller SCC-size part to a larger one cannot
  decrease the sum of `f` over the partition.  Concentrating any partition
  with at least two parts at `(n-r+1,1,...,1)` proves
  `Σ_i f(n_i)≤f(n-r+1)≤f(n-1)`, with equality at `(n-1,1)`.  Also distinguish
  the attained reducible temperature `f(n-1)/f(n)` from the first all-strong
  discrete shell `(f(n-1)+1)/f(n)`.  THM-2016 and its script/reflection are
  corrected.

---

## MISTAKE-211 (2026-07-21, codex GMC(2) audit) -- distinct return atoms were treated as separate equations inside one scalar Gaussian moment

- **What was claimed:** THM-1770(B) said that because distinct minimal balanced
  charge multisets give distinct coefficient monomials, the equation
  `E[P^m]=0` forces every atom monomial to vanish separately.  Part (D) then
  promoted this to closure of every pair-only/star support, and later work
  treated twice the primitive return length as a universal detection cutoff.
- **Why it is wrong:** a polynomial identity in *free variables* separates
  monomials; a scalar equation at one chosen coefficient point does not.  For
  `P=aZ^6+bW^2+cW^18`, the first return is `m=4` and
  `E[P^4]=4*6!*a*b^3+4*18!*a^3*c`.  With all coefficients nonzero, choosing
  `c=-6!*b^3/(18!*a^2)` cancels the two distinct primitive atoms exactly.
  Hence first-return minimality removes composite balanced words but does not
  remove cross-atom cancellation.  Exact radial-channel examples also have
  moments through `2R` zero and first die at `3R`, so `2R` is not a uniform
  cutoff.
- **Correct framing:** THM-1770(A) survives: every balanced word at the least
  return level is primitive.  Parts (B)--(D) are retracted.  To deduce
  one-sidedness one must prove that every positive-negative coefficient product
  lies in the **radical of the full multilevel moment ideal**, using an actual
  resultant/Hankel/cumulant tower.  The faithful state is at least bivariate
  `(charge, radial height)`; charge balance alone forgets the factorial weight
  and coefficient phase that permit cancellation.  See the resolved court case
  `CASE-gmc2-first-return-cross-atom-cancellation`, THM-2014, and HYP-8765.

---

## MISTAKE-209 (2026-07-21, codex audit of THM-1985/THM-1990) -- a harmonic SUBSET was computed as a term MULTISET, and the claimed linear-growth iff misses the whole Bertrand boundary

- **What was done:** THM-1985 and THM-1990 correctly found the simplex identity
  `sum 1/C(n,k)=k/(k-1)` and usefully organized many repo sequences by reciprocal
  mass.  But both files say that a sequence is a **subset** of the harmonic
  numbers while their scripts sum one reciprocal per indexed term.  Repeated
  values are therefore counted repeatedly.  THM-1990 also states that the
  reciprocal series diverges iff the sequence grows at most linearly;
  THM-1985 calls rapidly convergent unnamed census constants transcendental.
- **Why it is wrong:** a subset forgets multiplicity.  The labeled-tournament
  values `2^C(n,2)` are `{1,2,8,64,...}`, while the switching values
  `2^C(n-1,2)` are indexed as `1,1,2,8,64,...`; their **supports are identical**,
  so their support-harmonic masses are equal, not separated by `+1`.
  Similarly the subset masses of `0!,1!,2!,...`, Fibonacci, and Catalan lose
  one duplicated `1` relative to the termwise values `e`, the reciprocal-
  Fibonacci constant, and the displayed Catalan constant.  On convergence,
  for the strictly increasing support enumeration, `b_n ~ n log n` is already
  superlinear but `sum 1/b_n` diverges; more
  generally Bertrand series have iterated-log boundary exponents.  Finally,
  fast convergence alone proves neither irrationality nor transcendence.
  THM-1985 has two further independent errors: the repo contains conflicting
  duplicate `THM-1370` IDs, and the tournament file proves the omissions
  `7,21` plus finite coverage, while completeness of the Hamiltonian-path
  value spectrum (all other odds) remains conjectural; and
  `C(n,3)` counts all triple slots, not the maximum cyclic triples.  The latter
  is `n(n^2-1)/24` for odd `n` and `n(n^2-4)/24` for even `n`.
- **Correct framing:** distinguish
  `sigma_set=sum_m [m in image(a)]/m` from
  `sigma_multi=sum_m multiplicity(m)/m`; their difference is the nonnegative
  **collision tax** `sum_m (multiplicity(m)-1)_+/m`.  Decide support convergence
  from its counting function by Abel--Stieltjes summation, equivalently from
  summable relative occupancies of multiplicative blocks.  Use Bertrand's
  iterated-log test at near-linear growth; label unnamed convergent constants
  numerical unless an arithmetic proof is supplied.  THM-2000 is the repair
  packet.  Its exact tournament replacement is
  `sum_(n>=3)1/M_3(n)=75/4-24 log 2`.  Any use of the proposed complete
  Hamiltonian-path spectrum must be labeled conditional.  The exact simplex
  telescoping theorem in THM-1985/1990 is unaffected.

- **Two further audit corrections:** the later, explicitly named file
  `THM-1370-h-spectrum-omits-7-21-all-n.md` proves that `7,21` are omitted
  for all `n` and that every other odd value through `609` occurs by `n=8`; it
  explicitly labels global completeness of `odds minus {7,21}` a conjecture.
  Thus THM-1985/1990's claimed positive-density H-spectrum and reciprocal
  divergence are open, not corollaries of THM-1370.  Conversely the tournament
  triangular-number series was undersold: Gauss's identity gives
  `sum_(r>=0)q^(r(r+1)/2)=(q^2;q^2)_infty^2/(q;q)_infty`, equivalently a
  `theta_2` value.  At `q=1/2` this is the exact form of `1.641632560655...`;
  identifying it is not an open problem.

- **Legacy data corrections:** the old S447 A000568 list ended with
  `29305744576145`, which is not A000568(13); the correct value is
  `48542114686912`.  THM-2005 now reconstructs A000568(0)--A000568(20) by
  Burnside's all-odd-cycle formula and certifies the deduplicated reciprocal
  tail.  The old `H=1+2^(n-2)` neighbor sum also began at the unrealizable
  value `H=2`; its realized support begins at `n=3`, with mass
  `0.764499...`, not `1.264499...`.
- **Namespace warning:** that H-spectrum file and
  `THM-1370-elliptic-jacobian-conjecture-all-dimensions.md` both declare the
  identifier `THM-1370`.  Bare references are therefore ambiguous; use the
  full filename until the duplicate theorem namespace is repaired.

---

## MISTAKE-210 (2026-07-21, codex sequence-corpus audit) -- tournament counting growth was reversed, making two everywhere-divergent series look entire

- **What was claimed:** `03-artifacts/drafts/tournament-theory-comprehensive.tex`,
  `04-computation/complex_tournament_s339.py`, and
  `04-computation/beyond_finite_s339.py` state that
  `sum A000568(n)/n^s` converges (even entirely), and that the labeled EGF
  `sum 2^C(n,2) z^n/n!` is entire or has positive convergence radius.
- **Why it is wrong:** the orbit count satisfies
  `A000568(n)>=2^C(n,2)/n!`; this grows faster than every fixed power of `n`,
  so `A000568(n)/n^s` does not even tend to zero for any fixed complex `s`.
  For the labeled EGF, consecutive absolute terms have ratio
  `2^n|z|/(n+1)`, which tends to infinity for every `z!=0`; its radius is
  exactly zero.  The old comments reversed the comparison: `n!` beats a
  fixed exponential, not the quadratic-exponent growth `2^(n^2/2)`.
- **Correct framing:** reciprocalization creates the legitimate analytic
  objects.  `Z_V(s)=sum A000568(n)^(-s)` converges exactly for `Re(s)>0`
  (abscissa zero), and `sum z^n/A000568(n)` is entire.  The normalized Burnside
  correction `A000568(n)n!/2^C(n,2)-1` also decays and has a legitimate
  index-Dirichlet series (`tournament_dirichlet_s291b.py`).  THM-2000 records
  the proof; the old scripts remain historical negative evidence and must not
  be cited for analytic continuation.

---

## MISTAKE-219 (2026-07-21, codex reciprocal-atlas audit) -- finite prefixes, offsets, and proxy supports were promoted to sequence constants or all-n laws

- **What was done:**
  `04-computation/reciprocal_sums_of_repo_sequences_opus_S447.py` labels short
  prefix sums “tail-converged,” contains the wrong entry
  `A000568(13)=29305744576145` (the correct value is `48542114686912`), and
  computes the odd numbers minus `7,21` as though that were the proved global
  H-spectrum.  `04-computation/reciprocal_sums_of_our_sequences_kps_S128c144.py`
  mixes offsets and term multiplicities; its frozen output itself reports
  `match=False` for finite simplex approximations at `k=2,3`.  Separate
  sequence syntheses also propagated `2 selfK(n)=SC(n)` as all-`n` after the
  coincidence at `n=5,6,7`.  THM-1985 similarly promoted the prime-core
  identity `m({1,...,k})=H_k/C(k+2,2)` from the flagship `k=12` case to all
  `k`.
- **Why it is wrong:** a finite prefix supplies neither a tail bound nor an
  arithmetic type.  Offset changes can add duplicate `1`s, so they change the
  term-multiset mass while leaving the support mass fixed.  The H-spectrum file
  proves only two omissions plus finite coverage.  Finally THM-849 computes at
  `n=8` that `2 selfK(8)=404`, whereas `SC(8)=176`; the all-size self-line law
  is false, not a blue-line counting artifact.  For the deep-well measure,
  THM-819 gives the primitive sum
  `2/[(k+1)(k+2)] sum_(u<=k,(u,k+1)=1)1/u`; it reduces to
  `H_k/C(k+2,2)` when `k+1` is prime, including `k=12`, but not uniformly.
- **Correct framing:** every reported number must carry (i) its exact index
  range, (ii) support versus multiplicity semantics, and (iii) either a proved
  tail estimate or the label “finite prefix.”  Use
  `THM-1370-h-spectrum-omits-7-21-all-n.md` by filename and keep global
  H-spectrum divergence open.  Use A002785 for self-converse tournament
  classes and A051337 for strongly connected tournament classes.  THM-2000's
  optimization-safe referee replaces the two reciprocal scripts as the
  canonical analytic certificate; THM-849 remains the decisive self-line
  counterexample.  THM-819, not the prime-special display in THM-805, is the
  uniform deep-well law.

---

## MISTAKE-194 (2026-07-19, klein-S322, against my own THM-1290 runs AND flagging mac-mini-S54's census template) — THE UNGUARDED PAIR-COUNT MASK PRUNE IS UNSOUND: "missing unit-pairs > slots ⟹ prune" ignores that ONE future element that is a MULTIPLE of q satisfies pinning at q outright; the prune is valid only when maxnext < q (no future multiple possible)

- **What was done:** the in-branch bitmask prune (inherited from `lrc_gap_census40_S54.c`, mac-mini-S54, and generalized in my S319/S320 harness) pruned a DFS subtree at modulus q whenever the mult-bit was unset and (unhit unit-pairs) > (slots left). All THM-1290 runs (B=55 v1, B=64 S320, the LRC-mode run, the bottom-spectrum censuses) used it.
- **Why it was wrong:** a family can acquire a multiple of q AFTER the prune point — in descending-element DFS the modulus q ITSELF (q = 19, 23, 25 as a speed) arrives late — and then leaf pinning at q is satisfied via the multiple regardless of pairs. Detected by redundant spectroscopy: the v1 run counted 50 pinning survivors at B=55, the 14-mask rerun 49, and the PATCHED binary finds SIX survivors at w_max=54 alone where v1 had 4 and the rerun 3 — both prior counts were undercounts. The missed families are exactly tiny-anchor + mid-band clusters CONTAINING mask moduli as speeds.
- **The correct framing:** prune on pair-counts only when `maxnext < q`. Patched in `lrc14_subgap_census_klein_S319.c` (S322); all gates byte-identical; FULL re-certification suite launched (B=55 both parts, B=64, bottom census). **Until it completes, THM-1290's exhaustiveness claims and the S321/S322 "complete table" claims are UNDER-REVERIFICATION** (banners posted). ⚠ **mac-mini-S54's original harness has the same unguarded prune at q = 23, 25** — the n=12 gap census "EMPTY to height 48" (HYP-4117) and any downstream consumer should be re-run with the guard; flagged to mac-mini by broadcast, not unilaterally overridden.
- **Lesson:** an in-branch prune must be justified against EVERY leaf-escape clause, not just the one it models (here: the pair clause but not the multiple clause). Redundant re-runs with different pruning configurations are not waste — the 49-vs-50 diff was the only visible symptom.
- **RESOLUTION (klein-S323 harvest):** the patched suite re-ran every affected computation (~275.8B nodes): gap B=55 and [56,64] both 0-survivor/0-hard, LRC-mode B=55 0/0, bottom spectrum h ≤ 45 the identical 8 families. **Every THM-1290/S321/S322 conclusion re-certified sound**; the unsound prune had skipped real subtrees (222.0B vs 112.7B nodes on [56,64]) but hid no gap member, no counterexample, and no off-ladder family. mac-mini's n=12 re-run flag remains open on their side.
- Impact on existing results
- Source (who found it, when)

---

## MISTAKE-193 (2026-07-19, boxeph-S136; renumbered from 192 — kind-pasteur first-pushed 192 at 19:09) — pushed conflict-marker debris to main because the marker-grep GATE had inverted logic: `grep -c markers && continue` proceeds exactly when markers ARE found

**What happened.** Resolving a rebase conflict on `agents/.session-state.json` (nested conflict hunks), my single-hunk extraction left two `<<<<<<<` markers behind. The safety check I ran — `grep -c "<<<<<<<" file && git add ... && git rebase --continue && git push` — printed the count (2) and then PROCEEDED, because grep exits 0 precisely when it finds matches. The mangled file reached origin/main for ~1 minute before repair (valid JSON restored, verified marker-free, repushed).

**Lesson (MISTAKE-153 genus, new root cause).** A detector wired into a `&&` chain must FAIL on detection: gate with `if grep -q '<<<' f; then abort; fi` or `! grep -q '<<<' f && push` — never `grep -c f && push`. And after ANY hand-resolution, verify the file parses (here: `json.load`) before committing; a syntax check would have caught what the marker check was mis-wired to miss.

**Fix.** Repair commit on main; this entry. The same session-state file has now conflicted three times today (S133, S136 twice) — it should move to per-agent paths or get a merge driver; flagged to the fleet.
## MISTAKE-191 (2026-07-19, boxeph-S132, self-report on boxeph-S126; renumbered 189->190->191: opus-S400 first-pushed 189, opus-S402 first-pushed 190) — the frozen S126 output "990/271" is NOT reproducible from the committed script+seed (which gives 770/196); the "271" coincides with S124's DIFFERENT 271

**What happened.** The owner asked for a rerun of "the 271 mod-19 evaders" (a figure my own S131 letter carried from S126). Exact reproduction of the committed `lrc14_mod19_spread_kernel_boxeph_S126.py` part (B) — same lcg seed 999, same 6000 trials, same covering/spread predicates — yields 770 band+covering families and 196 mod-19-rung passers. The frozen `lrc14_mod19_spread_kernel_boxeph_S126.out` states "band+covering:990; +spread:271" — and that line's print FORMAT does not match the committed script's print statement, so the committed script is not the one that produced the frozen output (edited post-run, or the .out came from a variant). Suspiciously, "271" equals S124's different count (families CONTAINING 19 among its 1066-family bank) — likely transcription contamination between the two banks.

**What survives.** All qualitative S126-(B) claims stand under reproduction: the q=38 intra-modulus conditions are FEASIBLE and no sampled family reaches 3/38. Only the constants were wrong/unreproducible.

**Lesson (extends MISTAKE-183 and THM-1251's carried-constant rule to one's OWN artifacts).** Commit the script and its output in the same state that ran; when citing a count — even your own from three sessions ago — re-derive it from the committed artifact, not from the .out header or memory. A print-format mismatch between a script and its supposed output is a one-grep reproducibility check; run it before propagating any constant.

**Fix.** S132 reruns and all downstream statements use the committed-script figures (770/196); canon references to "271/990 evaders" should be read as "the S126 evader bank (committed-script reproduction: 196/770)". Files: 04-computation/lrc14_evaders_23_25_rungs_boxeph_S132.py + .out.

## MISTAKE-183 (kind-pasteur-2026-07-19-S128c84) -- I spent a session proving a theorem that was already PROVED in canon, because I searched for the METHOD and never searched for the STATEMENT

- **What happened.** Across cont.83 and cont.84 I worked death-star's residue `#runs <= 2D/3` / `D >= 3P`, whose entire purpose is to establish `sojourn <= 2/21` with equality only at `(1,2,3)`. That statement is **THM-1203** (codex-2026-07-18-S77): `mu(BAD_d) <= 2/21`, equality exactly at four-term APs, PROVED, computer-assisted, with a Lean kernel certificate on the finite core. It was in canon before I started.
- **Why I missed it.** I was chasing a *mechanism* -- runs, mirror pairs, gap coordinates, lattice tubes -- and I searched the repo for those words. codex's proof reaches the same conclusion through an object with no overlapping vocabulary at all: delete the four phase points to a **non-arithmetic additive triangle** `(p,q,p+q)` and bound six torus triangles of leg `1/7`. Nothing in my search terms would ever have surfaced it. I also inherited the target second-hand, from death-star's letter, and never went back to ask whether the *goal* was open -- only whether their *route* was.
- **It was worse than a near-miss.** THM-1203 section 10 is titled "Corrections to the preceding maximiser sketches" and refutes my own six-box draft **by name**, with the explicit witness `d = (1,6,7)`, `u = 3/4`. That is MISTAKE-181, pre-registered against exactly the reasoning I was extending. The refutation of my line and the proof of my target were in the same file, and I read neither.
- **Two attribution errors fell out of the same gap.** `|B| = (1/7)^3 = 1/343`, which I presented as a closed-form correction in THM-1211, is already codex's THM-1183 ("ambient volume 1/343", "AP X-ray density 2/21"). And my THM-1153 `n/(12n+1)` law is **THM-633** (mac-mini-S33), proved *and* Lean kernel-pure, ladder included -- boxeph-S123 had already credited it back to THM-633 while I was still calling it mine.
- **The rule I am adopting.** Before working any inherited residue: **grep canon for the TARGET STATEMENT, in the form of its inequality and its constants, not for the method.** Here, `2/21` alone would have found it in one search. A residue handed over in a letter is a claim about the state of canon, and claims about canon get verified against canon, not against the letter.
- **What survives, and it is the smaller half.** The witness law `P(D) = 3*floor(D/14) + e(D mod 14)`, the three gap-coordinate lemmas, the lattice-tube identity, and the branch-2 detection-floor finding are net-new. The ceiling they were built to prove is not.

## MISTAKE-176 (kind-pasteur-2026-07-18-S128c83, auditing my own THM-1150/1151) -- |B| = 0.003367 was a COARSE-GRID ARTIFACT; the true value is (1/7)^3 = 1/343 = 0.0029155 exactly, so every "28.28x concentration" I published should read 98/3 = 32.67x

- **What was claimed (THM-1150 step I, carried into THM-1151, THM-1152, THM-1154):** the bad set B is "six isolated boxes" at the permutations of (1/4,1/2,3/4) with measure |B| = 0.003367, and the maximising direction (1,2,3) concentrates the sojourn by a factor 28.28.
- **Why it is wrong (two errors, one geometric and one numerical):**
  1. **B is not six boxes -- it is six 3-SIMPLICES.** Writing the four defining slacks as `f1 = 2/7 - g_min`, `f2 = 2/7 - (g_mid-g_min)`, `f3 = 2/7 - (g_max-g_mid)`, `f4 = g_max - 5/7`, one has the identity `f1+f2+f3+f4 = 6/7 - g_max + g_max - 5/7 = 1/7`. A set cut out by four nonnegative affine functions whose SUM IS CONSTANT is a simplex, never a box. I never formed the sum.
  2. **The measurement was a midpoint grid**, whose error on a polytope is `O(surface * h)` -- on this set that is tens of percent even at `180^3`. The successive grid estimates 0.00333 (60^3), 0.00283 (120^3), 0.00267 (180^3) do not even converge monotonically; I read the first one as the answer. A Kronecker (low-discrepancy) sequence converges properly: +0.33%, +0.05%, -0.07% at 2e5, 1e6, 4e6 points against 1/343.
- **The correct framing.** In gap coordinates `x = g_min`, `y = g_mid-g_min`, `z = g_max-g_mid` (unimodular, volume-preserving) the constraints are `x,y,z <= 2/7` and `x+y+z >= 5/7`. Substituting `u = 2/7-x` etc. gives `{u,v,w >= 0, u+v+w <= 1/7}` -- the corner of the cube `[0,2/7]^3`, a simplex with leg 1/7 and volume `(1/7)^3/6`. Six ordering regions give **|B| = (1/7)^3 = 1/343**, exactly. Both boundary checks are strict (`x >= 1/7 > 0`, `g_max <= 6/7 < 1`), so there is no clipping and no wraparound.
- **Impact.** The sojourn values are UNAFFECTED (0.095230 = 2/21 is exact and was computed by exact affine-cell arithmetic, not by grid). Only the RATIO to |B| was wrong: the concentration at (1,2,3) is `(2/21)/(1/343) = 686/21 = 98/3 = 32.67`, not 28.28. THM-1154 has been annotated with a correction banner. More usefully, the identity that exposes the error is also the engine of THM-1211: it yields the general single-run bound `L <= (1/7)/d_exit` with no centre-hit hypothesis, and the incentre of the simplex -- all four slacks `1/28` -- is exactly the permutation of (1/4,1/2,3/4), which is where death-star's 1/28 comes from.
- **The lesson, and it is one I keep relearning.** I measured a quantity I could have computed in closed form. The grid told me a number; the algebra told me the SHAPE, and the shape came with a proof, a correction, and a new lemma. When a set is defined by finitely many affine inequalities, integrate it, do not sample it -- and always add up the constraints before assuming the shape.
- **Source.** kind-pasteur-2026-07-18-S128c83, self-caught while building on death-star-S58b's corrected maximiser bound. See THM-1211, `04-computation/simplex_slack_kakeya_kps_S128c83b.py`.

---

## MISTAKE-174 (opus-S388/S389, caught by codex-S78) -- the zero-cell tent inequality points in the opposite direction

- **What was claimed:** THM-1195 said that on every cell between consecutive
  zeros of `g(t)=min_v ||vt||`, the graph lies *under* the triangle through a
  cell maximum.  It inferred `cell area <= H*L/2`, summed this over the zero
  partition, and advertised `integral g >= 1/28` as an LRC(14) certificate.
  THM-1200 then compared that proposed threshold with an independence
  heuristic and interpreted their algebraic equality as an explanation of
  LRC tightness.
- **Why it is wrong:** every individual triangular wave is concave on a cell
  containing no one of its zeros.  The pointwise minimum is also concave
  (intersect the hypographs).  A concave function lies **above** the chords
  joining its endpoints to a maximum, so the universal inequality is
  `cell area >= H*L/2`.  The old proof reversed the geometry.
- **Exact in-scope counterexample:** for the actual 13-speed tight family
  `{1,...,11,13,24}`, the consecutive-zero cell `[1/24,1/13]` has active
  affine pieces `24t-1`, `t`, and `1-13t`, with slopes `24,1,-13`.  Its area
  is `185/100464`, whereas `H*L/2=11/8736`; the excess is `3/5152` and the
  ratio is `370/253>1`.
- **Correction and impact:** THM-1195 is withdrawn as a certificate; its
  arrangement routine still computes exact integrals, but the labels
  “CERTIFIES” do not follow.  THM-1200's parity result is independent and
  survives, while its tent/independence interpretation is withdrawn.  The
  correct summed law is the lower bound
  `integral g >= (1/2) sum_C |C| H_C`; it supplies no `1/14` pointwise
  witness without additional control of cell shapes or height distribution.
  See THM-1201 and the exact ordinary/`-O` replay
  `lrc14_tent_cell_direction_audit_codex_S78.py`.

---

## MISTAKE-175 (death-star-2026-07-18-S58b) [RENUMBERED from MISTAKE-173 by opus-S390: collision with MISTAKE-173 (opus-2026-07-17-S389), which pushed first at 23:12:22 vs 23:17:40. Content untouched.] (death-star-2026-07-18-S58b, auditing kind-pasteur THM-1150-six-boxes) -- the maximiser "centre-hitting criterion" dropped the mod-4 freedom AND the permutation symmetry; the standoff reduction is refuted (117 non-proportional hitters), but the bad<=2/21 CEILING survives with the maximiser being a SIX-RAY ORBIT
## MISTAKE-182 (death-star-2026-07-18-S58b) -- the centre-hitting criterion dropped mod-4 freedom and permutation symmetry

> This entry was originally also numbered `MISTAKE-173`. That identifier
> belongs to the earlier opus reduced-denominator/representation error below;
> the death-star audit is canonically `MISTAKE-182`. Its centre-incidence
> correction is sharpened by MISTAKE-180 and exact replacements THM-1206/1210.
> The historical language below predates THM-1181's correction from six boxes
> to the cyclic-gap polytope and must not be used as a proof dependency.

- **What was claimed (THM-1150-six-boxes, kind-pasteur S128c78, steps III & V):** the bad set B is six boxes around the six permutations of (1/4,1/2,3/4); a geodesic of direction (d2,d3,d4) hits the centre `(1/4,1/2,3/4)` iff `d proportional to (1,2,3)` ("elementary and complete"); hence the whole maximiser claim reduces to ONE Diophantine standoff: every NON-proportional integer direction keeps sup-distance `> rho ~ 0.0412` from all six centres (verified on 7 sampled directions, thin margin 0.0458 vs 0.0412).
- **Why it is wrong (two independent misses, both the project's classic "sampling misses structured families" trap):**
  1. **Dropped mod-4 freedom.** For the SPECIFIC centre `(1/4,1/2,3/4)`, hitting needs `d2 u = -1/4, d3 u = -1/2, d4 u = -3/4`. At `u=3/4` this is `d2=1, 2*(3/4)=3/2=1/2, d4*(3/4)=1/4`, and `d4*(3/4)=1/4 mod 1  <=>  3 d4 = 1 mod 4  <=>  d4 = 3 mod 4`. So `d4 in {3,7,11,...}` ALL hit, not just `d4=3`. Verified: `(1,2,7),(1,2,11),(2,4,14)` hit the centre exactly (min sup-dist 0), with positive bad measure `0.0246 / 0.0147 / 0.0246` (stable to N=50400). The step "substituting forces d4/d2=3" ignored the modular solutions.
  2. **Dropped the permutation symmetry.** B is symmetric under coordinate permutation (six boxes = six permuted centres), so the six PERMUTATION-RAYS of (1,2,3) -- `(1,3,2),(2,1,3),(2,3,1),(3,1,2),(3,2,1)` and multiples like `(3,9,6)=3*(1,3,2)` -- each hit their own box. All six give sojourn EXACTLY `2/21 = 0.095238` (N=1e6). (III) only examined the single centre `(1/4,1/2,3/4)` and so found only one of the six maximiser rays.
  A full scan found **117 non-proportional directions in [1,10]^3 that hit a centre** (the claim was 0), so the standoff reduction (V) is dead: standoff is NOT bounded below on non-proportional directions.
- **The correct framing.** The maximiser is the **coordinate-permutation ORBIT of (1,2,3)** (six rays), all at `2/21`, NOT the unique ray (1,2,3). The bad-measure ceiling `bad(d) <= 2/21` is UNBROKEN (measured max over all non-orbit directions in [1,10]^3 is `< 2/21`; the orbit attains it). But the proof cannot be a standoff/miss argument -- it must be a **SOJOURN-MAXIMISATION**: `bad(d) = meas{u : (frac(-d2 u),frac(-d3 u),frac(-d4 u)) in B} <= 2/21`, equality exactly on the six-ray orbit. That is an equidistribution/geometry-of-numbers statement (the resonant AP-ray maximises sojourn), the SAME "equal-spread AP is extremal" motif as the inverse theorem -- not a Diophantine standoff.
- **Impact.** THM-1150-six-boxes steps (III) and (V) are REFUTED; step (I) six-box geometry and (II) flow picture stand; the `bad <= 2/21 < 0.164 <= |S(P)|` CONCLUSION (uniform-r=5 analytic tail) is likely still true but its stated proof route is invalid. THM-1149-maximiser-proof-sketch step (III) (Moebius "fixed ratio forces (1,2,3)") has the same gap -- it forgot that a fixed ratio can be held on positive measure by any orbit member, and the wrapping case-analysis it defers is exactly where the mod-4 solutions live. Uniform r=5 remains OPEN, and the honest remaining piece is the sojourn ceiling, not a standoff bound.
- **Source:** death-star-2026-07-18-S58b. Scripts `lrc14_maximiser_permutation_hitters` / `lrc14_maximiser_sixray_sojourn` (+out). -> kind-pasteur.

## MISTAKE-169 (kind-pasteur-S128c68, caught by codex-2026-07-18-S67) -- a one-period window was said to contain a full safe gap

**What happened:** the four-comb gap-recursion scout asserted that an arbitrary
interval of length `1/k` contains a complete safe gap of the danger comb
`D_k`, of length `6/(7k)`.  It then used the adjacent-ratio threshold `7/6`.

**Why it is false:** after scaling to `k=1`, the interval `[1/2,3/2]` cuts the
closed safe arc `[1/14,13/14] mod 1` into two pieces, each of length `3/7`.
There is no safe subinterval of length `6/7`.  The scout's own measured minimum
`7*k4*L=4.9493` already contradicted its printed prediction `>=6`.

**Correction:** a one-period window contains at most two pieces of total safe
length `6/(7k)`, so it always contains one of length `3/(7k)`, and this is
sharp.  The sound recursive adjacent-ratio threshold is consequently `7/3`.
THM-1137 uses that exact lemma to prove a genuine multiplicatively spread
`r=6` cone.

**Impact:** only the explanatory claim in the S128c68 reconnaissance is
invalid; its sampled exact residual lengths remain data.  No pre-existing
canon theorem depended on the `7/6` recursion.

## MISTAKE-168 (death-star-S58, caught by codex-2026-07-18 tail audit) -- midpoint cells, a minimum branch, and an AP slice were promoted to exact uniform horn facts

**What happened:** the original `r6_Gsigma_exact_bands_deathstar_S58.py`
formed a list of centre/arc collision points, evaluated `G(sigma)` only at the
midpoint of each intervening cell, and labelled the whole cell good or bad.
THM-1132 then called the resulting endpoints exact.  The same theorem scaled
`T_old=min(N/(6mu),1/(3L))` by `3/7` as though its component branch were always
active, and described the equal-step phase AP as a reduction of the five-killer
problem.  It also said a component equal in length to one danger tooth could be
covered by that open tooth.

**Why it is false:** inside a fixed cyclic order, `G` is the maximum of affine
gap lengths and can cross `1/7` without any arc endpoints colliding.  For
example, `G(1/7)=2/7` and `G(17/50)=31/175>1/7`, contrary to the reported band
ends.  The printed “breakpoint-midpoint minimum” `1/14` even contradicted the
same output's exact value `G(1/5)=2/35`.  Algebraically, the exact identity is
`R_sharp=(3/7)R_comp` with `R_comp=(1/(3L))/k5`; when the counting branch wins,
`R_old` is a different quantity.  Geometrically, only offsets
`{b,b+d,...,b+4d}` form the asserted one-variable AP.  Topologically, the
residual is closed and a closed arc of equal length cannot fit inside an open
danger tooth.

**Correction:** with `u=min(sigma,1-sigma)`, the largest centre gap is

```text
max(u,1-4u) on [0,1/4],  u on [1/4,1/3],
max(3u-1,1-2u) on [1/3,1/2],
```

and `G=H-1/7`.  Therefore `G>1/7` exactly on
`[0,5/28) U (2/7,5/14) U (3/7,4/7) U (9/14,5/7) U (23/28,1]`,
of measure `9/14`; its minimum `2/35` occurs at all four nonzero
fifths.  The horn condition is `L>=1/(7k)`.  The sampled exact row and the
local horn survive, but global maximality, the AP landing/drift bridge, and
arbitrary-offset uniformity remain open.

**Affects:** THM-1132, HYP-7605, the sharp-horn search wording, the exact-bands
script/output, the auxiliary horn endpoint tests, and the S58 reflection.
Source: codex tail audit, 2026-07-18.

---

## MISTAKE-167 (boxeph-S111, caught by codex-2026-07-18) -- five named rows were reported as a universal covering-family floor

**What happened:** `lrc14_gap_theorem_test_boxeph_S111.py` evaluated five
hand-picked covering families—two AP examples and three non-AP examples.  One
non-AP row attained `M=1/12`.  The script, output, reflection, and hypothesis
ledger then said that *all* covering non-AP families have floor exactly `1/12`
and that `(1/13,1/12)` is empty.  There was no enumeration or proof supporting
either universal quantifier.

**Correction:** the computation is a five-row probe.  It supplies candidate
sharpness evidence only.  The proposed universal gap theorem remains open and,
on the fully covering stratum, is stronger than the corrected `INVcov` target.
The sieve still proves the stated non-covering fragment independently.

**Rule:** a minimum over a displayed witness table is not a family-wide floor.
Every claimed exhaustive minimum must state the enumerated domain and row count,
or be labelled finite evidence.

**Affects:** S111 script/output/reflection and HYP-7645.  The logical correction
“gap theorem implies the inverse target” survives; the alleged verification and
sharpness do not.

## MISTAKE-166 (boxeph-S108/S109, caught by codex-2026-07-18) -- the inverse target omitted Covering(2..14)

**What happened:** `LRCMSplit.lean` stated the open dominance premise as

```text
no Lonely13 time  =>  some speed dominates every other speed 13-fold.
```

Its comment said that `no Lonely13` already entails the `Covering` premise.  It
entails divisibility coverage only through **13**, not through 14.  The
distinction is real: the exact family

```text
{1,2,3,5,7,8,9,10,11,12,17,19,104}
```

has `M=8/105<1/13`, covers `2..13`, misses 14, and has
`104/19<13`.  It therefore falsifies the stated open premise.  S109 faithfully
transported that too-strong premise to the Finset target.  Lean had checked a
valid implication *from* the false/overstrong premise; the kernel was not at
fault, but the named mathematical target was not the intended crux.

**Correction:** the open statement is now explicitly

```text
INVcov: Covering(2..14) and no Lonely13 time
        => some speed dominates every other speed 13-fold.
```

In the `no Lonely14` branch, `counterexample_needs_all_divisors 14` supplies
`Covering(2..14)`; band monotonicity supplies `no Lonely13`; then `INVcov` and
`ap_core_bridge` close the branch.  The repaired kernel chain is
`LRC14_of_INVcov` / `LRC14_finset_of_INVcov`.

**Rule:** always put the threshold on a sieve-derived covering assertion.
`M<1/13` sees moduli only through 13; a hypothetical LRC(14) counterexample
(`M<1/14`) sees moduli through 14.  Never silently exchange those strata.

**Affects:** corrected `LRCMSplit.lean`, `LRCFinsetBridge.lean`, their root
comments, S108/S109 reflections and ledger entries.  This is the formal-target
instance of MISTAKE-161's covering-overload warning.
## MISTAKE-166 (boxeph-S110, caught by codex-2026-07-18-S67) -- five examples and a parameter identity were promoted to a global PG/tournament gap

**What happened:** S110 observed `183=13^2+13+1` and `14=13+1`, evaluated
five hand-picked integer families, and promoted their ordering to a spectral
gap `[14/183,1/12)` and a Singer/doubly-regular-tournament transport.  The
script contained no Singer set, no Covering filter, and no claimed 100-row
random census.  THM-724 was also cited beyond its proved single-killer scope.

**Exact refutation:** primitive Covering rows have exact values

```text
M(2*{1,...,12} union {13}) = 1/13,
M({2,3,5,8,9,11,12,13,14,15,17,20,23}) = 3/37 in (1/13,1/12),
M({1,...,12,364}) = 28/365 < 1/13.
```

Moreover `{1,...,12,182}` and `{2,...,12,182,184}` have the same residue
subset modulo 183 but global maxima `14/183` and `13/93`.  Thus the ambient PG
residue support forgets the integer lift and all other rational charts.

**Tournament correction:** a Singer `(183,14,1)` set has size 14, whereas a
tournament connection set on 183 vertices has size 91 (doubly regular
parameters `(183,91,45)`); a regular tournament on 14 vertices is impossible.
The augmented deep-well residue set and a Singer set are merely opposite-looking
additive examples, not vertices connected by a proved LRC-preserving map.

**Rule:** a shared parameter pair is not a functor.  Before transporting an
extremum or tournament class, name the vertex set, quotient, preserved LRC
predicate, and forgotten fibre; then test two distinct lifts of one quotient
point.  A short example table supports only those rows, never a population gap.
Source: THM-1131, HYP-7604/7635, corrected S110 script/output/reflection.

## MISTAKE-165 (opus-S375, caught by codex-2026-07-18-S67) -- a unit-residue kill count was applied to nonunit speeds

**What happened:** after correctly refuting the extended `q>14` sieve lemma,
THM-1110 proposed the repair `s*k_q<phi(q)`, where `k_q` counts nonzero
**unit** residues in the forbidden window.  Its hypothesis required only that
no speed be divisible by `q`, not that every speed be coprime to `q`.  A
nonunit speed can hit a larger gcd stratum and kill many more unit numerators.

At `q=90`, unit speeds kill `2` numerators, but speeds of gcd `5` kill `8`.
The three residues `5,25,35` partition all `24` unit numerators into their
kill sets; adding the redundant speed `1` makes `{1,5,25,35}` primitive while
all reduced `p/90` remain blocked.  This refutes the advertised arbitrary-row
`s<=11` reach already at four speeds.

**Correction:** for `g=gcd(v,q)<q`, the exact kill count is

```text
K_(q,g) = [phi(q)/phi(q/g)] * #{r in W_q : gcd(r,q)=g}.
```

The sound union bound is `sum_i K_(q,gcd(v_i,q))<phi(q)`.  The old
`s*k_q` form survives only when all speeds are units modulo `q`.  The exact
referee checks the stratified formula over `28,680` `(q,v)` rows.

**Rule:** “nonzero modulo `q`” and “a unit modulo `q`” are interchangeable
only when `q` is prime.  Every numerator count over `(Z/qZ)^*` must stratify
the speed action by `gcd(v,q)` before using a uniform orbit size.

## MISTAKE-164 (codex-2026-07-18-S67) -- Route B crossed the `1/14` threshold in the wrong direction and treated an affine carry chart as a quotient

**What happened:** the S101 Route-B gap reduction wrote a maximizing margin as
`M=s/q<1/13` and then worked in `13s<q<14s`, eventually specializing to
`q=13s+1`.  But `q<14s` is equivalent to `M>1/14`: it is the already-safe
band, not the hypothetical LRC-counterexample band.  Writing `q=13s+e`
makes the reversal explicit: safe means `0<e<s`, the LRC boundary is `e=s`,
and a counterexample has `e>s`.  The associated fourteen-gap defect is
`sum(g_i-s)=e-s`, so the small-gap pigeonhole changes sign at the boundary.

Two further steps were silently strengthened.  Even at `e=1`, twelve
`s`-separated points can have one gap `s+1`, so they need not be exactly
`s*{1,...,12}`.  And because `gcd(s,13s+1)=1`, the subset
`s*(Z/qZ)` is the whole cyclic group: offsets modulo `s` form an affine chart
with wrap carry `-e`, not nontrivial cosets.  A residue AP still determines
speeds only modulo `q`; integer lift coordinates carry the cross-modulus
cover obligations.

**Scope corrections downstream:** S102's `q=183` two-gap picture is a
suboptimal local maximum for its example (the global maximum is at `q=24`),
so it refutes an arbitrary-local-picture implication, not a premise explicitly
conditioned on the global maximizer.  S103's continued-fraction argument forces
the `13` carrier at `14/183`, but the `14` carrier comes from Cover14/AP-core,
n…83970 tokens truncated…ct measure in the last Farey cell). Dilation/scale invariances are free consistency tests -- run them on every new invariant. Secondary workflow lesson: `timeout N python | tee out` loses ALL buffered output on timeout with exit 0 (tee masks the kill); use `python -u` + per-section flush for long sweeps. -> THM-593, HYP-3840.

## MISTAKE-100 (2026-07-03, kind-pasteur-S37, caught by mac-mini-S31 THM-612 next session) — claimed "the tight locus is a SINGLE family (the AP), NO GW" from an under-powered search

**What happened:** searching for tight families (M=1/14) over APs {a,a+d,..}, dilates c·{1..13}, and ~2500 RANDOM 13-subsets of {1..30}, I found only dilated APs and concluded (HYP-4062) "primitive tight = {1..13} UNIQUE, no GW family." mac-mini S31 (THM-612) then verified **GW = {1..11,13,24} = AP[12→24]** is tight (M=1/14 EXACT, primitive, NOT a dilated AP, 6 tight points at units/14, non-covering) — a genuine second tight family.

**Why it was wrong:** the search was STRUCTURALLY BLIND to the GW shape. Random 13-subsets of C(30,13)≈1.2·10⁸ at 2500 samples has ~0 chance of hitting the specific one-residue-moved family {1..11,13,24}; and my structured candidates (APs + dilates) exclude it by construction (it is an AP with a single entry moved, neither an AP nor a dilate). Same class as mac-mini's own weak-adversary traps (MISTAKE-097/098): concluding "does not exist" from a search whose generator cannot produce the object.

**The correct framing:** the tight locus is at least {AP, GW}, both non-covering and small-speed — so the S37 REDUCTION and the 14-grid REPULSION stand (they only need tight ⇒ non-covering, which both satisfy); only the uniqueness/no-GW classification was wrong. LESSON: before claiming "no X exists" from a search, ask what shapes the generator CANNOT produce, and construct the adversarial shape explicitly (here: perturb the AP by moving one residue). A random/structured search is evidence of abundance, never of absence.

## MISTAKE-099 (2026-07-03, mac-mini-S30, caught in-session ~15 min in) — spent the session's opening re-deriving the Φ6 covering-min construction that I MYSELF refuted in a court case I filed 3 sessions earlier (S47)

**What happened:** tasked to "work the covering-min core", I conjectured `covering-min(n) = n/Φ6(n)` with tight family `{1..n-2,(n-1)n}`, verified `M(T_n)=n/Φ6` exactly + `apex≡-1 mod Φ6` + `g=2`, and started building a general theory around it — before recalling that `CASE-convergent-not-covering-min` (which **I filed**, mac-mini-S47, opus CONFIRMED S32) already REFUTED exactly this: the construction is NOT the covering-min for `n≥7` (beaten by 2/13, 2/15, 4/33 at n=7,8,9). The true covering-min is `1/n` over ALL covering families (opus-06-30 even block `2·{1..13}`, imprimitive) and empirically `~7–11%` above `1/n` for primitive ones.

**The resolution:** reading `02-court/active/` + `definitions.md` surfaced the correction within ~15 min; the session then produced real results (THM-610 deep-hiding lemmas; n=11 counterexample 3/31; margin map). Not a wasted session, but the opening detour was avoidable.

**Lesson:** the CLAUDE.md Step-5b "scour for leads" scan must include **`02-court/active/` for the exact quantity you're about to work** — especially cases you filed yourself. A named target ("the covering-min", "14/183") should trigger a `grep -rl` of court + hypothesis INDEX for that number BEFORE deriving. My own S47 filing sat one `ls 02-court/active/` away. Related recurring pattern: the "weak/wrong adversary" family (MISTAKE-090/098) — here the dual error, re-deriving a known-refuted *upper* construction. -> CASE-convergent-not-covering-min, THM-610, MISTAKE-090.

## MISTAKE-112 (kind-pasteur-2026-07-05-S15/S17): the consecutive multi-fold law M(D_l) = 14/(14(13-l)+1) is FALSE at l = 4, 5

The S15 instance verified the law on a sub-range and claimed it for l = 1..6 (HYP-4177). Exact merge-grid enumeration (S17): M(D_4) = 17/155 (claimed 14/127) and M(D_5) = 19/155 (claimed 14/113) -- the actual values are BELOW the law (the binding pair migrates to the 154+1 pair grid mid-ladder; l = 1, 2, 3, 6 match). The 2/25 FLOOR survives at every rung and is what HYP-4212's domination assembly consumes (LRCMultiFoldRows.lean certifies all six).

LESSON (the recursion of MISTAKE-102): a closed-form law verified at the ladder's ends and one midpoint is NOT verified -- the binding-pair structure can migrate at interior rungs. Enumerate the FULL parameter range exactly before stating a law; the merge grid makes this cheap (THM-592/HYP-4108: denominators are pair sums).


## MISTAKE-113 (mac-mini-S7, recurring from sibling-S6): adversarial free-fraction minimization on a FIXED grid masks sub-grid free sets

**What:** S5 (HYP-4282) concluded "consec[1..11] combs TILE the circle at radius 2/25 (phi_worst -> 0)" from a coordinate-descent that MINIMIZED the free fraction measured on a grid of 1600 points.  S7 exact arithmetic shows the free set has measure **0.000529** -- NONEMPTY -- but smaller than one grid cell (1/1600 = 0.000625), so the grid read it as 0.  Distinct-freq combs do NOT tile; kps's CircleClearFloor is correct.

**Why it's instructive:** minimizing a grid-estimated measure REWARDS configurations whose true free set hides between grid points -- the optimizer drives toward sub-grid-resolution free sets and reports 0.  The sibling hit the SAME trap in S6 (v1 adversarial "gridmax 0" artifacts).  Two instances, same lesson.

**Fix:** for covering/tiling questions use EXACT rational arithmetic (test midpoints of the elementary intervals cut by all arc endpoints -- finitely many, exact) OR verify any claimed phi=0 at >=100x the search resolution before believing it.  Never trust a minimized grid-measure at its floor.


## MISTAKE-114 (kind-pasteur-2026-07-06-S36, guard-railed by opus-S118 HYP-4506 and self at S38): the "window is too narrow / Dx<D" width narrative is a per-family SYMPTOM, not the obstruction -- the root is the ARITHMETIC of 3N+2

**What:** kps-S36 (HYP-4517) framed the n=12 first-gap emptiness via a metric mechanism -- the resonance-ladder crossing width Dx is smaller than the resonance spacing D, so the grid "skips" the gap. Read as *the* reason, this is misleading. opus-S118 (HYP-4506) proved first-gap emptiness is NON-MONOTONIC in N: N=13 is NONEMPTY (mediant 3/41 attained by {1..11,13,36}) while the WIDER-window N=12 is EMPTY. A width/Selberg story is monotonic (narrower => harder), so width cannot be the deciding quantity.

**Why it's instructive:** the deciding quantity is arithmetic -- the mediant 3/(3N+2) is achievable iff 3N+2 is PRIME (N=7->23, N=13->41 prime, nonempty; N=12->38=2*19 composite, empty). The metric Dx<D is a true computed fact about the single-outlier ladder subfamily, but it is downstream of the arithmetic: whether a base's resonance grid ALIGNS to the mediant is decided by the factorization of 3N+2, not by a width budget. kps-S38 verified the reconciliation (opus's nonempty witnesses ARE ladder families at N=7,13) and refuted a too-clean spin-off hypothesis ("F(N)={1..N-2,N}+3(N-1) gives the mediant iff 3N+2 prime" -- FALSE: F gives the mediant only at N=7,13, not at primes N=5,9,15).

**Fix:** keep Dx<D as the constructive/symptom side (the ladder is the WITNESS builder) and treat the ARITHMETIC of 3N+2 (mac-mini HYP-4562/4572 mod-19 clearance; opus O-arith / Fan-Sun gcd template) as the obstruction side. Do not present width/Selberg as the root cause of (G); it is non-monotonic and thus cannot be.

-> HYP-4282 (S5, the artifact), HYP-4312 (S7, the resolution), HYP-4292 (sibling S6, same trap), kps CircleClearFloor.


---## MISTAKE-114 (2026-07-06, mac-mini-S20) -- the FAST-EXACT-M helper (S16-S19) skipped non-coprime witness numerators, UNDERESTIMATING M.

**What was wrong.** My fast exact-M `Mfast` (offered to the fleet S16 as an O(n^2*max) exact-M via the witness-denominator lemma q|(v_i+-v_j)) had `for a in range(1,q): if gcd(a,q)!=1: continue`. The lemma is CORRECT (M's reduced denominator DIVIDES a pairwise sum/diff), but a witness at a SUB-denominator q' dividing (v_i+-v_j) appears over q=(v_i+-v_j) as a NON-coprime numerator (q/q')*a', which the skip DISCARDED -- so Mfast could MISS the true (larger) M.

**Caught by.** n=6 set {1,3,4,5,18}: Mfast gave 4/23=0.174 (looked like a gap member in the second gap (1/6,2/11)); the independent fine grid gave 2/11=0.1818 (the BOUNDARY, loose). The q=11 witness appears as 2a/22 with 22=4+18, and was skipped -- a FALSE POSITIVE.

**Fix.** Remove the gcd skip; check ALL a in [1,q). Re-verified vs grid: AP=1/13, doubled-apex=2/25, block=2/25, single-lift {1..11,23}=1/12, n=7 {1,5,6,11,16,17}=5/33 all correct.

**Impact (assessed, mostly benign).** (a) This session's 'n=6 gap member' = FALSE (it is 2/11, loose). (b) S16 targeted near-AP search RE-RUN with fixed M: still 0 in gap (15,976 families) -- conclusion HOLDS. (c) S17 n-specificity: the n=7 gap member 5/33 is GRID-confirmed independently -- STANDS. (d) S18 equioscillation AP=phi(n): counts UNCHANGED with fixed M (AP witnesses at a/13 coprime, 13 prime) -- STANDS. (e) S19 Fekete: direct energy, unaffected. The bug underestimates M so it risked FALSE NEGATIVES in searches; the n=13 emptiness rests on the fleet's correct-M exhaustive work (concurrent lift census), not my buggy searches. Files fixed: lrc_fastM_highscale_probe / lrc_leaveoneout_alignment / lrc_witness_denominator_dichotomies / lrc_equioscillation_count _macmini_S1x.py.

## MISTAKE-115 (opus-2026-07-06-S122, self-correcting opus-S120/S121): 'gap member = (N-2)-AP + exactly 2 defects' is FALSE -- 3-defect gap members exist; the defect count does NOT govern (the ORDER does)

**Claim (S120, wrong):** every LRC first-gap member is an (N-2)-term dilated AP + exactly 2 defects (longest-AP = N-2); the crux residual is a '>=3-defect Freiman-stability exclusion'.

**Refutation (S122, exhaustive at N=7):** {1,3,4,5,7,13,18} is a gap member (M=3/23 in (1/8,2/15)) with longest-AP {1,3,5,7} = 4 = N-3, i.e. 3 DEFECTS.  It coexists with the 2-defect member {1,2,3,4,5,7,18} at the SAME value M=3/23 (order 2).  So (a) 3-defect gap members exist, refuting the 2-defect signature; (b) the defect count is NOT the governing parameter -- two families of the same order (2) have different defect counts (2 vs 3).  The S120 signature was over-fit to 3 examples (N=6,7,13 first members) that happened to be 2-defect.

**Correct frame:** the governing parameter is the ORDER k of the value s/(Ns+k) (opus S116/S117), and the crux is kps's ACHIEVABILITY GAUNTLET (HYP-4557): in-gap values exist at every order for every N, and (G) at N=12 is that EVERY order's value is unattained -- a uniform-over-orders exclusion, NOT a bounded-defect one.  {1,3,4,5,7,13,18} is kps's 'no-isolated-runner species'.

**Lesson:** a structural signature read off 3 examples is not a theorem; exhaustive enumeration at the smallest nonempty cases (N=6 gave 1 member, N=7 gave 2 of different defect count) is cheap and would have caught it.  The proof map (00-navigation/LRC14-PROOF-MAP.md) crux line is corrected.

## MISTAKE-127 (klein-2026-07-08-S193, self-caught): the S192 "large-spread half via arc-count pigeonhole #arcs < rho*Vmax, prove #arcs <= c*spread with c<1" is the WRONG TOOL -- VACUOUS on the extremal near-dilated-AP family (where c ~ 1.17 > 1 AND rho* ~ 0.6), even though good periods abundantly exist

**Claim (klein-S192, wrong route):** the large-spread half of THM-527-A closes via the pigeonhole #{good ruler periods} >= rho*Vmax - #arcs(G*), so it suffices to prove [#arcs <= c*spread with explicit c<1] + [rho* >= rho0 > c]. Verified "zero failures" over 25 RANDOM primitive clusters/spread (c~0.2, rho*~0.99), and flagged c<1 as the residual.

**Refutation (S193, exact):** on the EXTREMAL near-dilated-AP family E_d=d*{0..9}u{p} (the low-rho* shape, longest-AP=10), #arcs ~ 1.17*spread (block-like, (k+1)/(k-1) > 1, so c<1 is FALSE) and rho* ~ 0.60, so the pigeonhole bound rho*Vmax - #arcs ~ -1545 is VACUOUS (deeply negative). Yet a good period abundantly exists: #good ~ 1612 ~ rho*Vmax (d=300), and #good/Vmax -> rho* = 0.594 with |#good - rho*Vmax| <= 7 -- the TRUE discrepancy is O(1), NOT #arcs (=3170). The random-cluster test missed this because random e is generic (small c, large rho*); the near-AP extremal was never tested.

**Root cause:** the crude bound |#good - rho*Vmax| <= #arcs is Koksma-Hlawka with the grid discrepancy 1/(2Vmax) x the total variation 2*#arcs of 1_{G*} treated as an ARBITRARY union of #arcs intervals. That is blind to the fact that the arcs of G* and the ruler grid {(j+1/2)/Vmax} share the SAME Vmax-arithmetic (both come from the phases frac(e_i x)), so they are not adversarially aligned and the real discrepancy is tiny. Only the Erdos-Turan analysis of the STRUCTURED indicator 1_{G*}=F(frac(v.x)) -- discrepancy driven by resonances a.v = 0 (mod Vmax) -- sees the cancellation.

**What survives / correct frame:** the CONCLUSION (a good period exists for large spread) stands, verified. The correct route is Erdos-Turan: |#good - Vmax*rho*| <= Vmax*D*, D* <= C_m(1/(H+1) + sum_{Vmax|a.v, ||a||<=H} 1/r(a)). PROVED sub-result: for near-AP the low-height resonances are d-INDEPENDENT (AP-supported a give a.v = d*sum(i*a_i), and Vmax=9d+14 forces sum(i*a_i)=0 for bounded a -- identical resonance set d=5..100), so the discrepancy stays spread-uniform => #good/Vmax -> rho*. THM-527 part H + THM-663 corrected.

**Lesson:** (1) an arc-count / total-variation (Koksma-Hlawka) bound is the WRONG instrument for "does an arithmetic grid hit a set built from the SAME arithmetic" -- it ignores the grid-set correlation and is vacuous exactly at the structured (extremal) cases; use Erdos-Turan on the resonances instead. (2) NEVER validate an extremal-sensitive claim on random samples only -- test the known extremal family (here near-dilated-AP, cf. MISTAKE-126). The extremal is where c is largest AND rho* is smallest simultaneously. Files: lrc14_nearAP_gridhit_klein_S193, lrc14_resonance_reduction_klein_S193.

## MISTAKE-126 (opus-2026-07-08-S155, court case filed): "block+outlier is the k=11 tail D3-minimizer / D3 >= D3_10 = 0.4646 / fixed-window cluster-monotonicity" (LEM-009, klein-S186/S187, kps-S86) is FALSE -- a dilated AP + interior point goes lower

**Claim (klein-S186/S187, kps-S86, wrong):** over the k=11 prim-diam>=25 tail, D3 is minimized by the block+outlier {0..9,D} (value 0.4587 / limit D3_10=0.4646), and every tail shape satisfies D3(E) >= D3_{c(E)} >= D3_10 = 0.4646 where c(E) = max points in a length-9 window ("cluster size"), which is decreasing in c.

**Refutation (S155, exact, by klein's OWN D3 code):** A = (0,3,6,8,9,12,15,18,21,24,27) = the AP 3*{0..9} (common difference 3) + interior point 8, primitive, prim-diam 27 (in the tail), has D3(A) = 88747403972619401646021583/195916463945506515076905312 = 0.452986. This is < D3_10 = 0.4646 (the claimed global bound) AND < 0.4587 (the claimed minimizer). Verified identically by klein-S184 exact Farey moments and opus-S148 moments_exact. A thorough search (56840 shapes) gives true tail min ~0.4530 at A (and its reflection). A has R2 = 590 = SAME as {0..9,25} but different D3 -- so D3 is not a function of R2 and the max-R2 shape is not the min-D3 shape.

**Root cause:** D3 is DILATION-INVARIANT (W_{cE}(x)=W_E(cx) => equal moments) and so is prim-diam; but the fixed-window "cluster size" is NOT. A has window-cluster 5 (predicting D3 >= D3_5 ~ 0.6) but contains a length-10 AP (its dilation-invariant "cluster" is 10). A is the tail analog of the EXHAUSTIVE minimizer 2*{0..9}u{9} = (0,2,4,6,8,9,10,12,14,16,18) (D3=0.4356, prim-diam 18) -- both "AP_10 (energy 570) + 1 point (+20) = R2 590", AP at a different scale.

**What survives:** the k=11 CLOSURE is not threatened -- true tail min ~0.4530 >= bar 0.3312 (margin +0.12); klein's block-decorrelation LIMIT values (D3_10=0.4646 etc.) are correct FOR THEIR FAMILIES; THM-662's R2 BOUND (<=590) stands (A satisfies it, though the uniqueness-of-maximizer sub-claim over-extends past the exhaustive range).

**Correct frame:** the dilation-invariant axis is the LONGEST AP in E, not the fixed-window count; min D3 is monotone in longest-AP (0.76/0.67/../0.467/0.453 at longest-AP 2..10); the extremal family is "AP_10 + 1 point" at any scale, tail min ~0.4530 at scale 3 + interior. The k=11 tail closes IF tail-inf >= bar (strongly evidenced) but via the AP-extremal picture, not window-cluster monotonicity.

**Lesson:** any extremal/monotonicity claim for a dilation-invariant functional (mu, D3, PZ) MUST be stated on a dilation-invariant axis. A "cluster/window" count is scale-dependent and silently misses dilated copies of the extremal structure -- exactly the trap that put a shape 0.012 below the claimed global bound. Verify candidate extremizers against their dilates before claiming global minimality. Files: lrc14_cluster_monotonicity_opus_S155, lrc14_tail_true_min_opus_S155; CASE-tail-D3-min-is-not-block-outlier-dilated-AP-counterexample.

## MISTAKE-134 (mac-mini-2026-07-09-S65 cont.14/15) — duplicated opus-S190's moment-floor discharge; grep-before-Lean skipped under context pressure

**What happened:** In cont.14 I wrote `LRCMomentFloorDischarge.lean` discharging `hsize`
(`clusterSize_shapeOf_le`) and re-exporting a five-parameter moment-floor assembly. In cont.15,
finally grepping before the NEXT leg, I found `LRCMomentFloorConcrete.lean` (opus-S190, hours
earlier) already contained an identical `clusterSize_shapeOf_le`, the discharged `hbonf`
(`bonferroni_concrete`), concrete `nuShapeConcrete`/`measGPConcrete`, and a STRICTLY BETTER
four-parameter terminal assembly `lrc14_from_momentfloor_concrete {hMoment, hB, hsmall,
hpartA}`. My file was removed (nothing imported it beyond the root; root import dropped).

**Root cause:** the MISTAKE-131 lesson (grep `TournamentH7/*.lean` for the statement shape
BEFORE writing Lean) was skipped in cont.14 — I grepped for *witnessG2 definitions* but not for
*existing discharges*, under end-of-context time pressure. The stale docstring in
LRCWitnessMomentFloor ("opaque … cannot be a theorem") reinforced the false impression that no
discharge existed; I flagged the stale doc in the same session yet did not draw the inference
that a discharge file likely already existed.

**Rule reinforced:** before ANY new Lean file: `grep -rl <key-identifier>` over the whole tree
— for the DEFINITIONS *and* for the DISCHARGES/consumers; a stale doc is itself evidence that
someone has been working the area. Context pressure is exactly when the rule matters most.

---

## MISTAKE-136 (klein-2026-07-09-S232, self-caught at S233): THM-684(I)'s orthogonality box object was MISIDENTIFIED -- the character layer sum equals the COMMON-MULTIPLIER (partial-live) count A_t(U) = #{c : c*u in B for all u in U}, NOT the product count M_t(U) = #{y in B^t : prod y = prod u}

**What happened:** S232 canonized "layer sum = M_t(U)/(q-1)" with M_t the product box count,
and ran the whole raw-vs-connected scale analysis on M_t. S233's convention check (Mobius peel
against the centered pair object) failed by O(q); re-deriving the orthogonality showed the
tuple sum over prod-chi = chi_0 forces ALL u_l*y_l^{-1} equal (not their product), giving
A_t(U). Direct character sums at q=61 confirm A_t exactly (1e-9) at t=2,3; the product form
M_2 = A_2 fails at 60/78 pairs at q=139.

**Why it slipped through:** (1) both objects share the main term b^t/q-ish, and at SMALL q the
integer deviations are small enough that the two counts frequently coincide (q=61 test
supports agreed exactly); (2) S232's numerical verification exercised the CS cascade (a true
statement about M_t) rather than the identity itself; (3) at t=2 the two objects are cousins
(ratio vs hyperbola parametrization of pair correlation) and both had appeared legitimately in
the program (THM-683 I vs S230's hyperbola), blurring the distinction.

**The cost:** S232's "raw M_3 devs ~ q" attribution and S233's first script run (wrong Mobius
normalization, garbage P_3 at scale 0.62*q). The qualitative conclusion of S232-III (raw
counts contain lower layers; connected form needed) survived by luck -- it is true for both
objects.

**The gift in the correction:** A_2 = THM-683's ratio object verbatim, and A_13 = LM(q)
itself -- the character program's box counts ARE the partial live counts. The corrected
cascade then yielded the relation-triple law with exact torus constants (THM-684 S233
addendum).

**Rule:** when an identity is "proved by orthogonality", VERIFY THE IDENTITY ITSELF
numerically (not just downstream bounds), at t >= 3 and at q large enough that deviations
resolve the candidate objects apart. Two counts sharing a main term are indistinguishable
exactly where verification is cheapest -- push one deviation-scale beyond agreement before
canonizing.

## MISTAKE-137
**Session:** mac-mini-2026-07-15-S109 (caught by S110's THM-868 GF referee, same day)
**What happened:** the S109 figurate note printed a "Fibonacci defect" line as F(n+2) − G(n)
= 0,1,1,2,3,5,8,13,22,38,... and narrated it as "itself Fibonacci-like until its own holes open".
The correct deficit is F(n+1) − G(n) = 0 (n < 8), then 1, 4, 13, 33, 76, 159, ... — opus-S317's
independently computed "deviations" — with exact GF x^8/((1-x)^5 (1+x)^2 (1-x-x^2)) (THM-868).
The fake pattern was pure index shift: F(n+2) − F(n+1) = F(n) makes any off-by-one look
"Fibonacci-like" while G still equals Fibonacci.
**Lesson:** an off-by-one against a linearly recurrent sequence MANUFACTURES a plausible pattern
(the recurrence reproduces itself under shifts). Before narrating "the deficit looks like X",
check the index against the exact GF — the GF referee catches in one line what eyeballing cannot.
**Status:** corrected in the S109 draft + THM-868; no downstream theorem used the wrong line.

## MISTAKE-152 — "Fejes Tóth equality iff regular" overclaim (klein-S313 cont.4b, corrected cont.5 same day)

**What happened:** LEM-020's cont.4b addendum claimed the pair-energy floor S₂ ≥ 6/7 has
equality "iff the 13 points are regular." FALSE: the tent kernel is convex but NOT strictly
convex, and the equality set is the full 12-dimensional polytope P = {g_i ≤ 1/7,
g_i + g_{i+1} ≥ 1/7} (adjacent-only overlap) — S₂ = 13/7 − Σg = 6/7 identically there.
Caught next continuation by sampling P directly (300 exact points, all at 6/7).

**Lesson:** for piecewise-linear/merely-convex kernels, energy minima have FLAT BOTTOMS;
"unique minimizer" claims require strict convexity. Always sample the suspected equality
set's neighborhood before claiming uniqueness. The flat bottom was the better theorem: the
covering adversary's playground is the polytope, and covering analysis = Kronecker line ∩ P.

**Affects:** LEM-020 addendum (corrected in place), reflection
the-coverage-spectrum-one-grammar-four-instruments.md (statement softened by this entry).

## MISTAKE-153 — K_{7,7}, K_{7,8}, and K_{8,8} were called open Zarankiewicz cases (death-star-S29/S30; corrected by codex-S20)

**What happened:** THM-922 labeled `K_{7,7}` open, and HYP-7106/S30 described the
class-coloring computations `108` and `144` as meeting open ordinary crossing-number
cases.  The numerical values and restricted class-coloring minima were correct, but the
literature status was not checked.

**Correction:** Woodall's cyclic-order computation proved `cr(K_{7,7})=81` in 1993.
Deleting one vertex from an eight-vertex part gives
`(8-2)cr(K_{8,n}) >= 8cr(K_{7,n})`; hence `cr(K_{7,8})>=108` and then
`cr(K_{8,8})>=144`.  Zarankiewicz drawings attain both bounds.  Thus all three values
are theorems.  S30 still proves a useful restricted result: its cyclic parallel-class
coloring minima equal the ordinary optima.

**Lesson:** distinguish an open general conjecture from a finite case already settled by
computation and deletion averaging; and distinguish ordinary crossing number from a
book/class-coloring restriction.  Check the primary case-status literature before
calling a small Zarankiewicz value open.  Source: D. R. Woodall, *J. Graph Theory* 17
(1993), 657--671, doi:`10.1002/jgt.3190170602`.

## MISTAKE-138
**Session:** mac-mini-2026-07-16-S127, caught S128 (same machine, next session)
**What happened:** S127 declared FragmentationCount.lean and TieSplitWalk.lean "kernel-verified,
zero errors" based on `lake env lean FILE 2>&1 | head -N; echo EXIT: $?` — but in a pipeline, `$?`
reports the LAST command (head), so the "verdict" was head's exit code, always 0. The actual
`lake build` in S128 surfaced five real errors in FragmentationCount (renamed `le_or_lt`, a fragile
nlinarith, scoped notation, and a statement that was FALSE without 0 < lam / 0 ≤ L in the empty
branch). The files were then repaired and now genuinely build.
**Lesson:** (1) never read `$?` after a pipe — use `${pipestatus[1]}` (zsh) or, better, use the
ARTIFACT as the verdict: the .olean's existence after `lake build` is the only build proof that
cannot lie. (2) A hypothesis-free inequality that "compiles" can still be FALSE as stated — the
empty-branch counterexample (negative lam) was caught only because linarith refused it: when a
prover balks at an "obvious" branch, check the statement before blaming the tactic.
**Status:** all three ladder files (FragmentationCount, TieSplitWalk, KillerBudget) now build with
oleans emitted; the S127 session log's "kernel-verified" claim corrected in the S128 entry.

## MISTAKE-158 -- a canon lemma (THM-523's q-witness) was independently re-derived and presented as new

**What was claimed:** kind-pasteur-S128 cont.50 presented a "sieve-margin lemma" as new
rigorous content in THM-995(IX): *if some q in {2..13} divides no speed, then t = 1/q gives
M >= 1/q > 1/14.*

**Why it is wrong:** this is verbatim **THM-523's q-witness lemma** (mac-mini,
2026-06-16), which states and proves exactly *if S contains no multiple of q in {2..14},
then tau = 1/q is lonely and M(S) >= 1/q >= 1/14*, together with the covering-set
necessary condition.  The re-derivation was correct mathematically but the novelty claim
was false, and the file did not cite the prior result.

**Correct framing:** THM-523 owns the lemma.  The only increment in THM-995(IX) is the
**strictness split** (q <= 13 gives the STRICT inequality 1/q > 1/14 with margin >= 1/182,
whereas q = 14 gives only M >= 1/14) and the consequent **pinning** of the tight locus to
the "covers 2..13, misses exactly 14" stratum.  Corrected in-file at cont.51.

**How to avoid:** before claiming a lemma as new, grep canon for its STATEMENT SHAPE, not
just its name -- here `grep -l "1/q" 01-canon/theorems/` or a search for "q-witness" /
"covering set" would have surfaced THM-523 immediately.  The covering/sieve reduction is
old and central; assume any elementary statement about small-divisor witnesses already
exists.

## MISTAKE-160 — the empirical covering floor M ≥ 1/9 (THM-995 X) undershoots; it contradicted a proved theorem (boxeph-2026-07-18-S85)

**What happened:** THM-995 (X) reported, from 3000 samples + local descent, an empirical
covering-family floor `M ≥ 1/9` (min at V = [3,4,11,12,13,15,18,20,24,42,55,64,67]). This is
WRONG: it undershoots and, more tellingly, it **contradicted the already-PROVED THM-724**
(primitive covering-min `= 14/183 = 0.0765 < 1/9`, attained uniquely at the deep well
`{1..12,182}`). Independent brute force (all pair-sum denominators, THM-999) confirms covering
primitive families with `M = 14/183` (deep well) and `M = 1/13` (the near-dilated-tight family
`2·{1..12}∪{13}`, covering + primitive + `ρ=1.08`), both `< 1/9`.

**Correction:** the true covering minimum is `14/183` (THM-724). The COMPACT sub-case
(`ρ = v_max/v_2nd < 13`) floor is conjecturally `1/13` (16k-family adversarial hunt, zero
counterexamples; extremal `2·{1..12}∪{13·odd}`; consistent with THM-726's proved `1/13` for
≥2 outliers). The dependent claim "12-subset floor `M(V∖{v_max}) ≥ (1/14)(1+1/ρ)`"
(the S84 reduction target) is likewise **FALSE** — the near-dilated-AP families violate it
(`M(V')` down to `1/13 < (1/14)(1+1/ρ)`) while `M(V) ≥ 1/13` still holds by their dilation
substructure. Elementary descent/sieve/measure all fail on these families.

**Lesson:** an EMPIRICAL floor from sampling + local descent can miss a measure-zero structured
stratum (here: dilations of the tight AP made primitive by one swap). ALWAYS cross-check an
empirical extremum against every PROVED bound in the repo — a sampled floor above a proved
minimum is a red flag that the sampler missed the extremal family. Structured adversarial
families (dilated APs, near-tight perturbations), not random/descent samples, probe the true
floor. Source: `lrc_12subset_floor` / `lrc_covering_infimum` / `lrc_compact_1over13_hunt`
_boxeph_S85 (.py + .out). Affects THM-995 (X) [corrected in place], the S84 12-subset-floor
reduction [refuted], HYP-7355.


## MISTAKE-161 (death-star-2026-07-18-S57) — conflated "covering 2..13" with "covering 2..14"; the 1/13 inverse theorem needs covering 2..14

**What happened:** In cont22 I "corrected" THM-1029/1038 to say the far-element candidate is the
smallest multiple of **13** (26, 39, …), "not 182, because covering 14 is never required (missing 14
only gives M≥1/14<1/13)." This is **WRONG**. It silently switched the covering hypothesis from **covers
2..14** (what THM-724, boxeph THM-1017, and the compact-floor-1/13 conjecture all use) to the weaker
sieve notion **covers 2..13**, and analyzed the cover-gap at level 1/13 over that wrong class.

**The witness (this is the valuable part):** `V = {1,2,3,5,7,8,9,10,11,12,17,19,104}` has
`M(V) = 8/105 = 0.076190` (verified exactly, stable to denominators ≤ 4000; witness `t=8/105`, min
attained by `v=1` and `v=104`). Then `1/14 < 8/105 < 14/183 < 1/13`: it is a **primitive, covers-2..13,
ρ=104/19=5.47<13** family with **M < 1/13** whose 12 non-max speeds are **NOT a dilated AP** — and it is
**below the deep-well covering-min 14/183**. So the statements "covering ⟹ (M<1/13 ⟹ AP core)"
(THM-1017), "compact ρ<13 floor = 1/13" (MISTAKE-160), and "M<1/13 ⟹ ρ≥13" are ALL FALSE if "covering"
means **2..13**. They are TRUE (and boxeph THM-1017 line 27's `14∣v_max` is valid) only with covering =
**2..14**: `V` **misses 14** (M ≥ 1/14 by the sieve, so no LRC(14) violation), and every covering-**2..14**
family with M<1/13 has an AP core (0 non-AP found in 138,129 perturbations, v_max∈{182,364,546}).

**Correction:** the LRC(14)-relevant class is **covering 2..14** (equivalently threshold 1/14: `M<1/14 ⟹`
covers 2..14 by the sieve). There the core misses 13,14 ⟹ `13∣v_max` AND `14∣v_max` ⟹ `182∣v_max`
(boxeph, correct). THM-1038's original candidate `182` was RIGHT; the cont22 "mults of 13" retraction is
withdrawn. My cover-gap enumeration and the "candidate correction" in cont22 analyzed the wrong (covers-
2..13, level-1/13) class and are void; the cover-gap *technique* (exact criterion, soft-Weyl bound) is
threshold-agnostic and survives, but must be applied at covering-2..14 / far-element-182.

**Lesson:** "covering" in this project is overloaded — the sieve-margin lemma uses **2..13**, but every
LRC-reduction theorem (THM-724/726/1017, the compact floor) needs **2..14** (= all residues 2..n, n=14).
The gap between them is a real, populated stratum (M ∈ (1/14,1/13), covers 2..13, misses 14), witnessed by
`V` above. Always state which covering is meant; at threshold 1/13 the two differ and non-AP "false
alarms" like `V` appear. Source: `lrc_covergap_uniform`, `lrc_covering214_test` _death-star_S57 (.py+.out).
Affects cont22 [synthesis §6/§7, cover-gap reflection — corrected], THM-1029/1038 [182 restored], sharpens
THM-1017 [covering must be 2..14].

## MISTAKE-155 (opus-2026-07-17-S367) — "it filled every slot so far" when only one slot existed

In THM-1065 I proposed extending the Bonferroni ledger to B7 by the same means that filled the S2 slot, arguing that "every slot filled so far was filled by containment and counting rather than by new analytic machinery." That sentence was literally true and substantively misleading: **exactly one slot had been filled**, at k=2 — the single value of k where the technique is provably sharp (THM-1012/1025). Generalising from one data point read as generalising from a track record.

S367 measured it: the containment floor loses a factor of ~5 per additional speed (exact/floor 3.5, 24.5, 114, 200, 2101 at k=2..6) and the fragmentation upper bound is ~1190x loose at k=7. Both are VALID; both are useless for a ledger needing O(1) relative accuracy.

**The lesson, which generalises past this project:** when proposing that a method extend, count how many times it has actually succeeded and check whether those successes share a special feature. Here they shared *all* of it — k=2 is where a single alignment assumption is the only assumption, so no compounding can occur. **A method that is sharp at the boundary case is evidence about the boundary case.** The mechanism of the success has to be checked for scale-dependence before the success is extrapolated. See THM-1070.


## MISTAKE-156 (opus-2026-07-17-S369) — scanning a SLICE of a dilation-invariant family, for the third time

While testing 13-term arithmetic progressions I scanned (a,d) = (1,d) for d = 1..89 plus a handful of others, found every d ≥ 2 gave uncovered ≈ 0.116–0.129, and wrote the conclusion **"every 13-term AP with d ≥ 2 stays uniformly away from 0."** That is FALSE. A wider scan — in the same session, by a script I had already written — returned min uncovered = 0 at (a,d) = (2,2), i.e. {2,4,…,26} = 2·{1,…,13}, a DILATE of the tight family. The whole diagonal a = d is tight.

The correct statement: **a 13-term AP is tight iff a = d**, i.e. iff it is a dilate of {1,…,13}; among primitive APs (a,d ≤ 16) only {1,…,13} itself.

**This is the third time dilation invariance has bitten in this program** — MISTAKE-154 (proposing a min-speed threshold), THM-1055 (the stratum run that looked like a threshold), and now this. Each time the mechanism was identical: sampling a SLICE of the parameter space of a quantity that is constant on dilation orbits, and reading the slice as the population. The fixed slice a = 1 meets the tight orbit in exactly one point (d = 1), which I had already set aside as the classical case.

**Rule going forward:** when a quantity is invariant under a group action, a parameter scan must range over ORBIT REPRESENTATIVES (here: primitive families, gcd = 1), never over a coordinate slice. I now have three instances of this and should treat any un-normalised scan in this project as suspect by default. See THM-1080, THM-1055, MISTAKE-154.


## MISTAKE-157 (opus-2026-07-17-S374) — proposing a route without spending five minutes attacking it

In THM-1100 I proposed the bounded-denominator conjecture — an absolute Q₀ with every primitive 13-family admitting a lonely p/q, q ≤ Q₀ — as the successor to the retired ledger route. I hedged the *evidence* carefully (explicitly refusing to read the sampled maximum 25 → 32 → 39 as a bound), and that hedging was right. What I did not do was spend five minutes trying to **construct** a counterexample from the definition.

The counterexample is one line. Blocking modulus q needs only one speed divisible by q, since that runner then sits at the origin for every p. So a single speed divisible by lcm(1..Q) blocks every q ≤ Q at once, and V = {lcm(1..Q)} ∪ {12 coprime speeds} refutes any absolute bound.

**The lesson is specifically about direction of effort.** I was careful in the right way about the wrong thing: I audited my *sampling* (having been burned three times by dilation) but never audited the *claim*. Escalating search maxima are evidence that a supremum is not being reached — which is exactly what an unbounded quantity looks like — and I read that as 'my search is weak' rather than 'the conjecture may be false'. **When a search keeps finding worse, try to prove it can always find worse before proposing the bound.** See THM-1105, and MISTAKE-152 for the earlier form of confusing a sampled maximum with a population bound.

## MISTAKE-170 -- a sufficient compact residual and an equality input were promoted to equivalences

**Sessions:** boxeph-2026-07-18-S113/S114; corrected by
codex-2026-07-18-S74 (THM-1099/1149).

**What happened:** S113 called the compact floor

```text
primitive + Cover14 + rho<13  =>  M>=1/13
```

equivalent to LRC(14), whose actual target is only `M>=1/14`.  S114 then
identified a maximum-deletion AP-core assertion with the twelve-speed
equality classification and with that compact/LRC residual.  The cited
reduction proves only a forward sufficient chain; no reverse implication or
reverse embedding was supplied.

The forward use also omitted a real hypothesis.  Equality rigidity can
classify `V\{v}` only after `M(V\{v})=1/13` has been established.  THM-1149
proves an exact alternative: every deletion can be loose, producing thirteen
pairwise-disjoint private essential regions.  A primitive compact row with
`M=8/105<1/13` realizes this all-loose crown while covering `2..13`; it
misses `14`, showing exactly where the missing cross-modulus input must act.

There were three accompanying source-level scope errors: the S113 script's
`range(2,14)` omitted modulus 14, its `QMAX=250` value was a bounded search,
and `V[:-1]` removed 24 although deleting 13 exposes the dilated AP.  The
claim that THM-1013 handles every dilated-AP-core compact row also omitted
THM-1013's condition on the extra speed.

**Correction:** retain four typed arrows:

```text
Cover14 compact strict row
  -> tight deletion                         [OPEN crown collapse]
  -> d[12] deletion                         [OPEN n=12 equality]
  -> 13d divides extra speed                [PROVED THM-1149]
  -> primitive 14-carrier ratio conflict    [PROVED THM-1149].
```

The compact `1/13` floor is a stronger sufficient route to the `1/14` LRC
target, not a proved equivalent.  S113/S114, HYP-7665/7675, THM-1013, and the
single-row script/output are corrected in place.

**Lesson:** before declaring two inverse statements equivalent, write every
arrow with its threshold, extraction hypothesis, and reverse realization.
Classification of an object does not extract that object from a larger one.

## MISTAKE-171 -- distinct four-comb moduli do not force a 4/3 actual-mean gap

**Sessions:** kind-pasteur-2026-07-18-S128c69/70; corrected by
codex-2026-07-18-S74 (THM-1148).

**What happened:** THM-1141 proposed that four distinct danger combs on a
core-safe interval should have longest survivor gap at least `4/3` times the
actual mean survivor gap.  The exact legal row

```text
P={1,...,8}, J=[1/14,13/112], killers=(108,109,110,111)
```

has five survivor components with

```text
L/mean=638/573=1.1134...<4/3.
```

Yet it easily satisfies the desired metric inequality:
`7*111*L=319/72>1`.  Nearby teeth coalesce; the component count falls, so a
nearly uniform gap word can have a large mean.  The two-comb linear endpoint
law does not force its full descent to be sampled by the surviving
multi-comb word.

**Correction:** retain overlap-cluster count together with the labelled
metric gap word.  THM-1148 replaces the false actual-mean lemma by a sharp
four-residue multiplier cone, the exact Q4 mass/component gate, and the
corrected THM-1137 `Phi` transfer.  Uniform `r=5` remains open, with
`m(3,4,5,6)` an explicit infinite proof-method residual.

**Lesson:** “maximum exceeds a baseline mean” and “maximum exceeds the
actual component mean” are different claims.  Overlap can improve the first
while refuting the second; always name the denominator being averaged.
## MISTAKE-170 (codex-S67, caught by concurrent S74 audit) -- `Covering` did not repair the inverse premise

**What happened:** MISTAKE-166 correctly observed that `no Lonely13` supplies
divisibility only through 13, so the proposed dominance premise was amended to
assume `Covering(2..14)`.  The amended proposition `INVcov` was then described
as a genuine noncircular open target.  Dilation was not audited.

**Exact refutation:**

```text
W=2*{1,...,13}={2,4,...,26}
```

is positive and covers every modulus `2..14`: use `2q` for `q<=13` and speed
`14` for modulus 14.  Dilation preserves the exact AP maximum, so
`M(W)=1/14<1/13`; hence there is no `Lonely13` time.  But the largest two
speeds are 26 and 24, so no speed 13-dominates the rest.  Therefore literal
`INVcov` is false.

**Correction:** retain `LRC14_of_INVcov` and its Finset bridge only as valid
conditional implications from a refuted premise.  `ResidualINV` remains the
exact counterexample interface and is equivalent to working LRC(14) under the
cited AP bridge, so it is diagnostic rather than a smaller theorem.  Any live
noncircular inverse supplier must include primitive normalization and a proof
that after gcd reduction Covering is rederived from the no-`Lonely14` branch.

**Rule:** after adding a divisibility hypothesis, test it under common dilation.
`Lonely` is dilation-invariant; `Covering(2..14)` is not.  A corrected domain
can still admit new extremal dilations that falsify the intended conclusion.

**Affects:** `LRCMSplit.lean`, `LRCFinsetBridge.lean`, the formalization
manifest/picture, THM-1131, HYP-7615/7625/7675, and the S108/S109/S111/S114
reflections.  See THM-1158 and `LRCINVcovCounterexample.lean`.

## MISTAKE-171 (kind-pasteur-S128c69/c70, caught by S74/codex audit) -- three different gap means were conflated

**What happened:** THM-1141 compared the observed longest gap with the
uniform-interleaving benchmark `3/(7 sum k)` and called that benchmark the
actual mean gap.  It proposed the universal lemma
`L_max >= (4/3) mu_actual`.  THM-1147 then averaged a complete linear
two-comb descent and treated that average as a source of the required
four-comb dispersion.

**Exact refutation:** for core `P={1,...,8}`, core component
`[1/14,13/112]`, and killers `(108,109,110,111)`, the final gaps have lengths

```text
319/55944, 305/55944, 291/55944, 277/55944, 13/3024.
```

Hence `L_max/mu_actual=638/573<4/3`, while
`7*111*L_max=319/72>1`.  The first four gaps are themselves the exact
THM-1147 pair-law values for `(108,111)` and `j=8,9,10,11`; their max/mean is
only `319/298`.  The full positive descent continues to `j=30`, but the core
window exposes only a short nearly constant slice.

**Correction:** distinguish (i) the mean over a complete pair-law branch,
(ii) the actual surviving-component mean, and (iii) the uniform-interleaving
benchmark.  With `D=L_max/mu_actual` and
`B=mu_actual/(3/(7 sum k))`, the exact target is
`D*B>(sum k)/(3k4)`.  Pair/end-point dispersion controls `D`; mass/overlap or
multiplier gains control `B`.  Neither may silently replace the other.

**Rule:** every use of “mean gap” must name its numerator, component count,
and sampled index set.  A closed-form full branch says nothing about the
accessible indices after core truncation and obstruction by other combs.

**Affects:** THM-1141, THM-1147, HYP-7676/HYP-7560, and the proposed r=5
nonuniformity route.  THM-1148 owns the exact guardrail row.

---


## MISTAKE-172 (opus-2026-07-17-S384) — excluding the regime where the answer lived

Testing whether beat frequencies certify loneliness, I restricted the search to q > 14 — reasoning that q ≤ 14 was "the classical sieve, already understood" and I wanted to probe the new regime. The run then reported the tight families {1,…,13} and {1,…,11,13,24} as having best margin **0**, i.e. failing to certify, which I briefly read as the certificate breaking exactly where it mattered.

It was my filter that broke. Those families certify at **q = 14 = 1 + 13** — a sum beat, sitting exactly on the boundary I had excluded. With q ≤ 14 restored they certify with margin 6.

**The lesson:** when excluding a regime as "already understood," check that the objects you most care about do not live in it. The extremal families are extremal *because* their witness sits at q = 14; excluding q ≤ 14 removed the witness from precisely the cases that test the claim hardest. A filter chosen to isolate the interesting regime silently redefined the experiment. See THM-1175.


## MISTAKE-173 (opus-2026-07-17-S389) — testing a claim on REDUCED denominators

Checking whether THM-401's setup ("M(S) is attained at a pair-sum time") is complete, I computed the maximising point t, took its denominator **in lowest terms**, and asked whether that denominator is a pair sum. Four of 25 families came out "difference-only", and I briefly believed I had a counterexample to a PROVED canon theorem and was preparing a court case.

The test was wrong. A point t has many representations, and reducing the fraction can turn a pair-sum denominator into something else. For V = {9,15,16,23,25,27,31,35,37,41,43,46,51} the optimum is 1/6, whose reduced denominator 6 is a difference — but **1/6 = 4/24 and 24 is a pair sum**, so pair-sums attain the identical value with deficit 0.00000000. All four cases were this artifact. THM-401 stands.

**The lesson:** when a theorem says "the optimum lies at a point OF THE FORM m/q with q in some set", the test is whether the POINT admits such a representation — not whether its reduced denominator belongs to the set. Reduction is not representation. I caught this only because I ran an explicit per-case verification before filing; the aggregate count alone said 4/25 and looked convincing. See THM-1200.
## MISTAKE-180 (codex-2026-07-18-S79) — discarding the integer branches in a torus congruence

The superseded six-box draft claimed that the closed geodesic
`u -> ({-d1 u},{-d2 u},{-d3 u})` hits the centre
`(1/4,1/2,3/4)` only when `d` is proportional to `(1,2,3)`.  In solving
`d_i u = n_i-r_i/4`, it effectively compared the three rational ratios after
discarding the independent integers `n_i`.  Those wrap integers are the
problem, not a removable nuisance.

The exact counterexample is

```text
d=(1,2,7),    u=3/4,
({-d_i u})=(1/4,1/2,3/4).
```

THM-1206 proves the complete replacement.  If `g=gcd(d)` and `e=d/g`, the
geodesic hits a labelled centre `r/4`, with `r` a permutation of `(1,2,3)`,
iff `e=+r` or `e=-r (mod 4)`.  Thus `(1,2,4m+3)` supplies infinitely many
nonproportional exact hits.  In particular its proposed uniform
positive standoff for nonproportional directions is false.

This correction does not refute the conjectured `2/21` *measure* ceiling:
incidence at one phase is not positive sojourn.  The missing coordinate is
contact order/sojourn length, which a centre-only classifier destroys.

**Rule:** never cancel or suppress the lift integers in a torus congruence.
First primitive-normalize, then solve in the torsion quotient; here the whole
criterion lives in `(Z/4Z)^3/{+/-1}`.

The earlier death-star audit that first exposed this mistake also found 117
nonproportional labelled-centre hits in its bounded scan and emphasized the
six coordinate-permutation AP rays.  Those observations are now subsumed by
THM-1206's exact congruence and THM-1203's exact equality classification;
they are telemetry, not a separate “six-box” theorem.

**Affects:** the superseded six-box draft, HYP-7600, and THM-1181.  See THM-1206.

## MISTAKE-181 (kind-pasteur-S128c77, corrected codex-S77) — BAD does not force exact balance

**What happened:** the superseded maximiser proof sketch observed the
balanced four-gap point and then treated membership in the full BAD region as
if it forced that equality point.  From there it argued that positive BAD
measure requires a persistent `1:2:3` edge ratio and concluded that every
nonproportional direction has zero BAD measure.

The implication is false.  In max-gap coordinates the exact condition is

```text
Delta_1+Delta_2+Delta_3+Delta_4=1,
1/7 <= Delta_i <= 2/7.
```

This is a three-dimensional inequality region, not the single balanced point
`Delta_i=1/4` (nor `1/8`, which refers to survivor lengths in a different
coordinate).  The explicit nonproportional direction

```text
d=(1,6,7)
```

has exact BAD measure `5/147>0`; at `u=3/4` it even hits the labelled balanced
centre `(1/4,1/2,3/4)`.  Thus both the zero-measure conclusion and the
proportionality premise used to reach it are refuted.

**Correction:** THM-1203 keeps the whole inequality region.  BAD forces every
pair difference into three bands, then deletes to one non-arithmetic additive
triangle `(p,q,p+q)`.  Six exact torus triangles, a sheared shifted-grid tail,
and a 99-pair exact core prove the desired uniform ceiling
`mu(BAD)<=2/21`.  The four triangle-deletion obligations subsequently force
all three adjacent gaps equal and classify the equality locus.

**Rule:** an extremal equality configuration does not parametrize its entire
sublevel set.  Before transporting an equality ratio along a flow, write the
full defining inequalities and test a non-equality interior point.

**Affects:** the superseded maximiser proof sketch, HYP-7595, and
the centre/standoff continuation.  See THM-1203 and MISTAKE-180.

## MISTAKE-184 (codex-2026-07-19-S82) — a multiplied beat relation was mistaken for extra equidistribution

THM-864 allowed an arbitrary presentation

```text
qB-pA=+/-y
```

and claimed a localization gain proportional to `1/y`.  Multiplying all
three relation entries was therefore allowed to manufacture an arbitrarily
stronger estimate without changing the underlying relation.  The claimed
`y` starting points are a `1/y`-net only when the associated step numerator
is primitive modulo `y`; the proof explicitly discussed the nonprimitive
case but incorrectly called a repeated proper subgrid equidistributed.

The boundary case

```text
delta=1/13,       A=3744,       B=3745,
E=[1/3,1/2],      kappa=1,
p=q=y=12
```

meets every displayed hypothesis: `gcd(A,B)=1`,
`qB-pA=12=y`, and `A=26qy`.  If both danger conditions hold, then

```text
||t||=||(B-A)t||<=||At||+||Bt||<=2/13,
```

so the restricted overlap on `E` is zero.  Since `13|A`, the exact global
pair mass is `rho=4/169`, and hence

```text
error=|E|rho=2/507.
```

THM-864's clean right side is instead

```text
13rho/[y(p+q-1)] + (8kappa+10y+8)/(13B)
 =13129/3359265,
```

smaller by the positive amount `531/14556815`.  Thus the theorem is false as
stated; a battery containing only primitive relation presentations could not
test this failure.

The invariant repair is to primitive-normalize the transverse relation
before attaching a clock.  For coprime coefficients `p,q`, choose `u,v` with
`qu-pv=1` and put

```text
h=qa-pb,       k=ub-va.
```

Then `(a,b)=(pk+uh,qk+vh)` and the unimodular change of coordinates preserves
`gcd(a,b)=gcd(k,h)`.  A valid positioned estimate may depend on this primitive
`(p,q;k,h)` datum (and on the actual starting-point multiplicity), but not on
a scaled presentation.  In the counterexample the primitive relation is
simply `B-A=1`; replacing it by twelve copies creates no new phase samples.

More exactly, with `c=B^(-1) (mod A)` and
`k_*=(yc-sigma q)/A`, one has the exact identity
`gcd(k_*,y)=gcd(p,q)=d`.  After dividing the relation by `d`, Bezout applied
also to `Bc=1 (mod A)` proves the reduced step primitive.  The unperturbed
starts are therefore exactly `d` superposed copies of the uniform grid of
size `y/d`.  Thus a genuine `1/y`-net occurs exactly when `gcd(p,q)=1`.
This repairs only the net claim.  The published quantitative
bound also treats an unwrapped path `J` as a circle arc when it may wrap and
compresses several error ledgers without proved constants, so the displayed
estimate still requires a fresh proof even after primitive normalization.

**Rule:** relation height, transverse clock, and orbit multiplicity must be
defined after primitive normalization.  Never infer discrepancy decay from
a coefficient that changes when the same Diophantine relation is multiplied.

**Affects:** THM-864's theorem statement and proof, its clean and exact error
bounds, HYP-6925's localization summary, and any downstream argument that
uses the asserted `1/y` gain.  The static height-seven classification in
THM-605 is unaffected because it requires coprime pattern coefficients.

## MISTAKE-185 (codex-2026-07-19-S82) — static height-seven positivity was promoted to finite-interval independence

THM-598 and THM-602 split a pair or cluster using only relations of
coefficient height at most seven.  They then asserted that if every such
phase completes one cycle across a target interval, all pair intersections
are forced close to their global mean.  This conflates two different facts:
THM-605 says a *globally phased primitive pattern* can have zero overlap only
at height at most seven; it does not say a short segment of a high-height
exact orbit samples that positive global overlap.

The exact counterfamily is

```text
(a,b)=(64K,75K),
I_K=[407/(896K), 407/(896K)+449/(4928K)].
```

For the base pair `(64,75)`, exact rational endpoint enumeration gives the
largest zero-overlap gap

```text
[407/896,489/896],       length=41/448=451/4928.
```

Thus `I_K` lies strictly inside a scaled copy of that gap and

```text
D_(64K) intersect D_(75K) intersect I_K = empty.
```

Nevertheless every nonzero integer vector `(r,s)` with
`|r|+|s|<=7` satisfies

```text
|64r+75s|>=11,
|64Kr+75Ks| |I_K| >=11*449/4928=449/448>1.
```

So THM-602 declares its truncated resonance lattice zero, and THM-598 calls
the pair resolved at every listed low pattern, while the local pair overlap
is exactly zero.  The missing relation is the high primitive exact relation

```text
75*(64K)-64*(75K)=0,
```

of height `139`.  Its global fixed-phase overlap is positive, consistently
with THM-605, but the chosen interval is shorter than the common `1/K`
period and can miss it completely.

The original THM-598 audit had two additional warning signs.  Its
`PQ<=16` list contains twenty-one primitive rows, not the claimed thirteen,
and the sum of the displayed individual envelopes over that list is
`2788339/2162160>1`, not a small tail.  Its advertised hard examples also
had minimum resonance `19` at window length `0.01`, hence `19L<1` and were
frozen rather than resolved by its own definition.  The computation did not
test the theorem hypothesis it was cited to verify.

There is a second algebraic guardrail: an HNF basis of a truncated relation
lattice cannot automatically be completed to a unimodular matrix unless the
lattice is saturated.  A torsion-sheet sidecar is otherwise lost.  This is
the cluster-level analogue of MISTAKE-184's repeated improper subgrid.

**Correction:** THM-605 parts (i)--(ii), including the exact nine static
channels, remain proved.  The dynamic carrier must retain both (a) the full
primitive exact-relation/gcd period and its torsion sheets and (b) low-height
detuned relations with their actual interval phase.  A height-seven cutoff
alone is not a finite-interval inverse theorem.

**Rule:** global positivity of every phase fiber does not imply that a finite
orbit segment samples the fiber.  Before truncating a resonance lattice,
retain the exact integer kernel, its saturation index, and the ratio of the
target interval to the common orbit period.

**Affects:** THM-598 Parts B--D, THM-602's fully-resolved branch and claimed
HNF renormalization, and THM-605(iii)'s assembly restatement.  THM-599's
global `c`-averaged torus-band identities are unaffected.

## MISTAKE-186 (boxeph-2026-07-19-S129) — a `∀ m` loneliness hypothesis made the mod-19 spread lemmas VACUOUS (caught before wiring)

**What happened.** In S127/S128 I formalized the mod-19 antipodal-spread lemmas (`LRCMod19Spread.lean`)
with the closeness hypothesis stated as `∀ b, ∃ i, ∀ m : ℤ, |c_i·(b/19) − m| < 2/19`. The inner `∀ m`
is WRONG: for a fixed real `x = c_i·(b/19)`, `∀ m, |x − m| < 2/19` is UNSATISFIABLE (`|x| < 2/19` and
`|x−1| < 2/19` give `1 < 4/19`, false). So the hypothesis can never hold and `antipodal_spread` /
`antipodal_cover` were VACUOUSLY TRUE — they built kernel-pure and sorry-free, but said nothing.

**Why it slipped.** I copied the shape from `LRCMod13Blocking.no_middle_band_witness_of_tight`, which has
the same `∀ m` (there it is a one-off contrapositive helper used only at `m = 0`, so its vacuity is
harmless; as a MAIN hypothesis it is fatal). "Kernel-pure, sorry-free, builds" does NOT imply "non-vacuous"
— a false/unsatisfiable hypothesis passes the kernel silently.

**The fix.** The intended condition is "some runner is within `2/19` of SOME integer", i.e. `∃ m` (equal to
`dist(c_i·(b/19), ℤ) < 2/19`, i.e. `margin < 2/19`), not `∀ m`. Changed `∀ m → ∃ m` in `hclose` (and in
`no_middle_band_of_close`), using the witnessed `m` in the `sieve19_single` contradiction. The hypothesis is
now SATISFIABLE (e.g. `{1,…,12}`, `M = 1/13 < 2/19`) and the lemmas are meaningful. Verified by wiring
`hclose` to the ledger's `margin` framework (`LRCMod19LedgerBridge.antipodal_cover_of_margin`): `margin
v (b/19) < 2/19 ⟹ hclose` via `le_margin_iff`, which only type-checks because the `∃ m` form matches.

**How to apply.** When a Lean lemma's hypothesis is a loneliness/closeness condition, sanity-check that it
is SATISFIABLE (exhibit a witnessing family) before trusting or wiring it — a `∀ m` where an `∃ m` was
meant is vacuous and the kernel will not complain. Prefer stating closeness as `dist(·,ℤ) < c` or
`margin < c` (which is `∃ m`) rather than an unguarded `∀ m`. See [[lrc14-crux-state]], HYP-7812,
`LRCMod19Spread.lean`, `LRCMod19LedgerBridge.lean`.

## MISTAKE-187 (death-star-2026-07-19-S59) — searched my own thread's vocabulary, not the target's: 30 minutes re-proving THM-1004/1005 because "Hamming" ≠ "single-far"

**What happened.** I planned "the single-defect single-far stratum of the N=12 gap, closed by an
absorption lemma + finite check" as a new theorem. I DID run the MISTAKE-183 statement-grep first —
but with MY thread's vocabulary: "single-far", "single-outlier", "absorption". Those greps found only
empirical sweeps (mac-mini-S26) and THM-633 (i=12), so I built and ran the closure... which exactly
replicates klein's THM-1004 (Hamming-1 rigidity, 2026-07-17) and THM-1005 (Hamming-2), down to the
identical interval table (ℓ(i=5) = 7/1000 = THM-1004's L_5). The rigidity thread calls the same object
"Hamming radius ≤ 2 of the AP"; one grep for "Hamming" in 01-canon/theorems would have surfaced both
files instantly. Worse: my OWN S58g/S58h session entries cite "THM-1004/5/6" by number — I read those
citations the same morning and never dereferenced them, mis-filing them as n=14-kernel-specific.

**Why it slipped.** Vocabulary mismatch defeats statement-grep: two threads can name one object
("replace k elements of {1..12}") by disjoint terms ("Hamming-k perturbation" vs "k-defect k-far").
And a theorem number cited in your own notes feels "already integrated" — it isn't, until dereferenced.

**The damage.** ~30 minutes of redundant compute; the near-claim was caught at writeup time (the grep
for codex's n=12 Hamming banks surfaced THM-1005's title). Silver lining: the replication is an
independent double-witness of THM-1004/1005 (different code, same rationals), and the cross-N
extension (THM-1284) survives as the genuinely new content.

**How to apply.** (1) When statement-grepping, ALSO grep the canonical synonym families for the object:
for near-AP work that means "Hamming", "defect", "outlier", "replacement", "perturbation", "far" — not
just your own term. (2) A THM-number citation appearing in your own writing is a POINTER, not absorbed
knowledge: `ls 01-canon/theorems/ | grep <number>` and read the title before planning anything in its
neighborhood. See MISTAKE-183, MISTAKE-131, THM-1284 §5.

## MISTAKE-188 (kind-pasteur-2026-07-19-S128c86, correcting opus-2026-07-19-S396 / THM-1235) — a rung-realization negative scoped to a region where the rung cannot live: "D=2 not found" while {1..12,26} = 2/27 sat inside the scanned shape

**What happened.** THM-1235 (opus-S396) reported for the slack-1 ladder D/(14D−1):
"Testing which rungs are realised, only D = 1 and D = 3 turned up; D = 2, 4, 5, 6, 7, 8
were not found." Exact computation (gate-verified) shows M({1,…,12, 26}) = 2/27 exactly —
primitive, pair (1,26), s = 27, D = 2, slack 1 — and {1..12, 26} lies inside the scan's
own shape family {1..12, x}, x ≤ 400. The D=2 rung IS realised.

**Why it happened (most probable reading).** The session's search pipeline was aimed at the
interval (1/14, 3/41); 2/27 = 0.0741 > 3/41 = 0.0732 lies OUTSIDE that interval, so an
in-interval filter applied before the rung-realization check silently excluded the only
place the D=2 rung can live. A negative about "which rungs are realised" was thereby
scoped to a region where the D=2 rung is impossible BY DEFINITION — a vacuous negative of
the same genus as detection floors (MISTAKE-162, HYP-7870 IV), but arising from a scope
filter rather than a weak searcher. The interval-emptiness claim itself is unaffected.

**The one-evaluation miss.** The canonical family for the D=2 rung is the direct THM-633
transfer K₂(13) = {1..12, 2·13}: at q = 13m+1, a = m, the far element 26 sits at distance
exactly 2/q (13m ≡ −1) and the base is in-band. One three-gap evaluation of the ladder's
own canonical family would have found it. The K-ladder {1..N−1} ∪ {cN} attains c/(cN+1)
at EVERY (N, c) tested (N ≤ 24, c ≤ 8, 0 violations — lrc14_ladder_realization_crossN).

**Rule.** A "value X is not realised" claim must (i) state the region actually searched
and check X lies inside it, and (ii) evaluate the canonical/constructive family for X
(here: the ladder shape that realises the neighbouring rungs) before the negative is
recorded. A realization survey inherits every scope filter of the pipeline it ran in.

**Affects:** THM-1235 (amendment banner added; slack-1 status now D=1,2,3 realised, D ≥ 4
open), downstream discussions of "isolated floor vs accumulation" (three consecutive
realised rungs), HYP-7840's framing. Scripts:
04-computation/lrc14_ladder_realization_crossN_kps_S128c86.py (+.out).

## MISTAKE-189 (opus-2026-07-19-S400, against my own S399) — headlined a "zero-uptake" lead whose statement-level half was already answered in canon (HYP-4096, fourteen days earlier)

**What happened.** The S399 history synthesis promoted boxeph-S114's proposal — mine the
Sungkawichai–Trakulthongchai paper for its equality case, "Wall A may be implicit in a paper we
already cite" — as the highest-leverage unclaimed lead, recommending "one session: fetch, read
the equality analysis." kind-pasteur-S1 (2026-07-05, HYP-4096, reflection
`the-tight-locus-rigidity-of-lrc13-kps-S1.md`) had already read the paper for exactly this and
recorded: it "does *not* characterize extremizers (it only cites Goddyn–Wong for the word
'tight'), so this is genuine open mathematics, not a literature citation." A one-grep check
("Sungkawichai" in 07-reflections/) would have surfaced it; the S399 session — whose own §3/§8
preached MISTAKE-183's grep-the-statement rule — ran six subagent sweeps and still missed it,
because none of the sweeps targeted the lead's own statement.

**Why it is instructive (beyond MISTAKE-183 again).** (1) Synthesis sessions are not immune:
aggregating 175 mistakes into a taxonomy does not apply the taxonomy to the synthesis's own
outputs. Rule: before RECOMMENDING a lead, grep its statement exactly as if about to work it.
(2) The recovery matters: the S400 session proceeded anyway because the PROOF-level half
(what the paper's internal structure pins, vs what it states) was genuinely unmined — and that
half yielded the S-T tightness cage + height-258,276 rigidity (HYP-7920). A lead can be
half-dead and still be the right thing to work; the mistake is mis-stating WHICH half is open.
(3) Corollary for backlog hygiene: a lead entry should cite the canon items that ALREADY bear
on it (HYP-4096 should have been in the S399 backlog entry).

**Cost:** none in session time (the deeper mining was the correct move regardless); the cost is
a mis-calibrated recommendation broadcast to the fleet for one day.

## MISTAKE-190 (opus-2026-07-19-S402, against my own S401/HYP-7930) — propagated a mis-extracted inclusion chain into a HYP headline: the "1/13 floor" route rode a garbled star-index; the proven chain ends in gridmax-land, not lines-one-dimension-down

**What happened.** S401's HYP-7930 headline asserted "acc(13-speed M-spectrum) ⊆ [1/13, 1/2]
via their acc(S(n)) ⊆ S*(n−1) + settled LRC(≤13)". The pinning session (S402) shows the
proven chain in the PUBLISHED v4 is acc(S_k(n)) ⊆ S_{k+1}(n) ⊆ S*_{k+1}(n) = S*₀(n−k−1):
for lines the terminal set is S*₀(n−2) — D-values of FINITE proper subgroups = grid-loneliness
values of the (n−2)-torus — which has NO 1/13 floor (gridmax((1,…,11); 14) = 1/14 exactly,
verified; interior window values realizable). The "= S*(n−1)" k=1 line in my S401 source
fetch was a mis-rendering. Consequences of the error: corollaries (C2)/(C3) (window
finiteness; 12-speed gap = finite list) are CONDITIONAL on their Conjecture 1.5, not
derived; only (C1) floor-isolation survives unconditionally — and by a DIFFERENT (stronger)
route: the proven "only upper accumulation points" phrase (now THM-1289). Also the
abstract-vs-body contradiction (the /abs fetch served the v1/v2 abstract whose equality
claim was WITHDRAWN in v3 and demoted to Conjecture 1.5 in v4) went unnoticed for one day.

**Instructive because:** (1) two extractions of the SAME source can disagree; the tie-breaker
is internal cross-validation (here: Theorem 1.3's mirror statement validated the 1.4 phrase,
and the version history + footnote-2 error mechanism explained the discrepancy) — one fetch
is never a pin. (2) A DERIVED-status headline still propagates: agents read headlines, not
status fields; keep conditional clauses IN the headline. (3) The garbled star-index was
precisely the paper's own historical error (their footnote 2: subtorus∩subspace is a
subGROUP — disconnected); when an extraction smooths over exactly the subtlety that broke
the source's earlier version, suspect the extraction. Cost: one day of a mis-scoped
headline; the S402 pinning was already scheduled, and the designed process caught it.

## MISTAKE-192 (kind-pasteur-2026-07-19-S128c91, against my own HYP-7955 H-G2) — hunted the "extinction prime" c(p) > 13 while my own gate, two sessions earlier, had already proved the sea never empties: the AP inhabits I(13,p,1) at EVERY p coprime to 14

**What happened.** HYP-7955 seam G proposed the extinction principle "c(p) > 13 ⟹
I(13,p,1) = ∅" and HYP-8000 (ex-7975) hunted the crossing with compiled code (gated 9/9, greedy
triage over 130 primes to 1200, parallel exhaustive sweeps). The triage found every
prime ALIVE — because the answer was a one-line theorem I already possessed: the
acceptance-test gate (HYP-7930, two sessions earlier) verified "(1..13) improper mod
every tested p", and the general proof is immediate from the pinch/located-maximizer
theory (THM-401/1002): every maximizer of the AP {1..13} sits at a sum-14 pair, hence
at denominator 14, so at any t = a/p with gcd(p,14) = 1 the min-distance is STRICTLY
below 1/14 and the AP is improper. **c(p) ≤ 13 identically; the extinction prime does
not exist; H-G2 was vacuous as posed.** Same for k=12: {1..12}'s maximizers live at
denominator 13, so I(12,p,1) ∋ AP12 always — the "S-T retrodiction control" had a
predetermined answer too.

**Why it happened.** The extinction idea was born inside the covering-number frame
(seam G's c(p) growth data) and never re-checked against the improper-tuple frame's
own gates — the two frames are one object, and the gate FACT ("AP improper at every
tested p") was recorded as an instrument check, not read back as a structural theorem
quantified over all p. A vocabulary-internal blind spot: MISTAKE-183's grep rule
would not have caught it (same session family, same files); what catches it is asking
"does a known object inhabit the set I want to empty, for structural reasons?"

**Rule (the eternal-inhabitant check).** Before hunting for the emptiness/extinction
of a set, list the known structured objects and check each for membership AT THE
GENERAL PARAMETER — a gate verified "at every tested p" whose proof is p-uniform is a
theorem, not a check. (Kinship: THM-1289's only-upper-accumulation logic; the
tight families are eternal inhabitants of every level-1 sieve.)

**What survives (real yields, not salvage-spin).** (1) c(p) ≤ 13 for all p coprime to
14 — now PROVED, with the AP cover as the uniform witness; c(p) data pins the band
{≤12 → 12 → ≤13} with the 12-onset at p = 181. (2) The ragged greedy-frontier (15
greedy-hard primes among 130) survives as a SIEVE-COST predictor: greedy-difficulty
ranks how thin I(13,p,1) is per prime. (3) The instrument lesson: exhaustive-MRV
order is terrible at FINDING covers (minutes) where randomized greedy finds them in
0.1s — finding and refuting need different search orders. (4) The REFRAMED question
(the one the hunt was groping toward): not "does I(13,p,1) empty" (never) but "does
it COLLAPSE TO THE TIGHT-FAMILY SHELL as p grows" — i.e. does c(p) reach and stick at
13 with few, classifiable minimal covers — the large-p ansatz-completeness question,
now sharply posed with instruments in hand.

**Affects:** HYP-7955 (H-G2 struck), HYP-8000 (ex-7975) (resolution recorded in place), the
S128c90 reflection (seam-G consequence paragraph amended), CONSTANTS-INDEX (no
extinction-prime entry to add — the non-existence noted at the c(p) line).

## MISTAKE-195 — Claim-check and claim-write executed in one compound command (death-star-S59m)
Ran `grep -c HYP-8070` and the stub-write in the SAME Bash call; the grep revealed the collision (kind-pasteur had pushed 8070 minutes earlier) but the write had already executed before I read the output. Cost: one renumber cycle (8070 -> 8075) and an INDEX edit race. Rule sharpened: the fetch-grep-claim sequence must be SEQUENTIAL TOOL CALLS — inspect the check's output BEFORE the write ever runs; never batch a namespace check with the write that depends on it.

## MISTAKE-196 — Emptiness claims from un-validated search instruments (death-star-S59n, caught in-session twice)
(a) The general-k ansatz carried y^(k+1) where the weight arithmetic demands y^k; the row-solver then "proved" absurd rigidity (dim 0 spaces EXCLUDING the known witness). Caught because the known solution must lie in every claimed solution space — an impossibility check that costs one substitution. (b) The first numeric Newton hunt returned 0/800 hits INCLUDING at k=2 where the witness exists: the search was drowning in the degenerate attractor (C = delta = 0 gives det = 0 = "constant"), and the filter silently rejected everything. The k=2 validation gate exposed it; the fix (Keller-forcing equation c0(0)*nu = 1 deleting the degenerate branch) is the general pattern. RULE: a search that cannot rediscover the known witness has ZERO evidential weight for emptiness; run the witness-rediscovery gate BEFORE believing any negative result, and prefer surgically deleting degenerate branches (reciprocal-variable equations) over post-hoc filtering.

## MISTAKE-197 — Presenting a classically-trivial identity as a novel theorem (death-star-S59o/p)
THM-1320's headline factorization det = −E₀(0)A(0)C(0) is det JF(0) in disguise: any Keller map's constant det equals the det of its linear part, and the six-function form's linear part is the antidiagonal of the unit constants. The identity was derived through the c₀-formula and felt like discovered structure; one normalization (G = L⁻¹F, next session) exposed it as linear-part invertibility. Caught in-session-chain and amended. RULE: before filing a "factorization/necessity" theorem, evaluate it at the simplest distinguished point (the origin, the identity, the trivial case) and check whether the claim reduces to a textbook fact there; the FRAME can still be canon-worthy — file it as a frame, not a theorem.

## MISTAKE-198 — Three-instance recurrence fit claimed as a verified law (S419 diagonal sums)
- **Who/when:** opus-2026-07-19-S419, caught by opus-2026-07-20-S420 (same agent, next day).
- **What:** S419 stated the owner-triangle diagonal sums 1,2,4,7,12,21,37 "satisfy a(n)=a(n-1)+a(n-2)+a(n-4) (verified on all three available instances)". The statement was literally true and still a mirage: the law breaks at the very next term (m=8: predicts 65, actual 68 — both on the pure Faulhaber grid and on the owner's corrected triangle). The true sequence (1,2,4,7,12,21,37,68,129,...) has NO constant-coefficient recurrence up to order 8 (kind-pasteur-S128c103 lists it OEIS-new).
- **Genus:** category E (small-case extrapolation; Moser's circle / width-of-G_n). A 3-term-deep fit of an order-3-support recurrence has zero spare verification instances.
- **Rule:** a recurrence of order r fit on a sequence needs >= 2r+2 terms with the fit verified on the surplus; for families with super-exponential drivers (Faulhaber columns), ALWAYS compute at least one continuation term beyond the data used, before claiming any law.
- **Repair:** S420 reflection section 4 (museum of impersonations); banner added to the S419 reflection and HYP-8155.

## MISTAKE-199 (Nth recurrence, death-star-S61b/c/d) — three ID collisions in one session cluster on fleet-wide owner prompts
The owner's arborescence and odd/even/Babai-Cameron prompts went fleet-wide, and I collided THREE times by filing before re-pulling: (1) THM-1445 (switching H-sum) vs opus/kp THM-1445, pushed 11 min earlier — renumbered to THM-1460 then again; (2) THM-1460 (arborescence det-shadow) vs mac-mini THM-1460, pushed 23 SECONDS earlier and carried further (two poles, ordinal-sum log) — renumbered mine to THM-1467; (3) THM-1465 (canonical member / Babai-Cameron 7.4 = 0 at every odd n) vs kp THM-1465 (5 min earlier) AND opus THM-1460 (10 min earlier), BOTH identical (all-even anchor n≡1, all-ODD anchor n≡3, via klein's score-parity law) — CEDED entirely, my file deleted. Net: of my three S61b-d "theorems," two were independent rediscoveries of same-day fleet work and one (THM-1467 switching-sum) is the only distinct survivor, plus the 3/8-mass confirmation of boxeph HYP-8295. HARDENED RULE (again): on any owner prompt, `git fetch && rebase` IMMEDIATELY BEFORE the checkpoint that claims an ID, not just at session start — the fleet moves in minutes, not sessions. And when a prompt is visibly fleet-wide (Babai-Cameron had klein+opus+kp already), default to CONFIRMATION/synthesis and do not file a competing theorem number at all. The distinct-contribution test must be applied BEFORE filing, not after the collision.

## MISTAKE-204 (boxeph-2026-07-20-S182) — THM-1635's ladder was verified by injecting its own ansatz: a numeric that models β_eff(m) = β(1 + c/m) tests the ladder's ARITHMETIC, not its {m^{-k}} scale PREMISE

**What happened:** THM-1635 §3 built a 1/m-Vandermonde ladder on the premise
that tied jump coefficients vary on the scale set {m^{-k}}, and the S182
machine check "confirmed" it — by SYNTHESIZING data of exactly that form
(β(1 + c/m)) and watching the tied sum come out (c₁−c₂)/m. Circular: the
check could not have failed on the premise, only on the algebra.

**The truth (referee, machine-verified in
`04-computation/tie_ladder_scale_referee_boxeph_S182r.py` + frozen `.out`):**
the scale premise is FALSE for generic tie arcs. Inverting a generic branch
r(t) = C²/t² + A/t + ... produces HALF-STEP Puiseux corrections
t = C r^{-1/2}(1 + u r^{-1/2} + ...), and Laplace over the germ turns those
into β_eff(m) = β · e^{-u√(2m)}(1 + O(m^{-1/2})): the honest scale set is
{e^{a√m} · m^{-k/2}} (measured: log β_eff/√(2m) → −u, −0.2957 at m = 512 for
u = 0.3). Even INTEGER-step corrections v/r, which do give a true 1/m ladder,
dress the base constant to e^{-2v}β (0.5488 measured) — the ladder's "β" was
never the raw fold coefficient. Re a ≠ 0 breaks "slowly varying" entirely;
imaginary a puts √m-drifting phases under the Vandermonde's nose.

**Genus:** the numerical cousin of MISTAKE-203. That one: representation-level
ARGUMENTS dodge (only function-level invariants survive). This one:
ansatz-level NUMERICS dodge — a check whose test data is generated FROM the
model's own form verifies closure of the model under its own assumptions,
nothing more. Also kin to MISTAKE-186 (vacuous hypothesis: the check was
true, but of a vacuously self-referential statement).

**Rule:** before building (or "verifying") an asymptotic ladder, DERIVE the
scale set from the geometry (here: invert the branch, Laplace the germ, read
off which powers of m appear in the exponent AND the prefactor). A
verification numeric must be fed data from the GEOMETRY (integrate the actual
germ), never from the model's ansatz. If the only available test data is
synthetic, say "arithmetic check" in the writeup, never "machine-verified."

**Repair:** THM-1635 §6 (verdict + the √m-graded repair route), THM-1630 §6
(dressed-constant amendment; domination survives since e^{O(√m)} never
crosses (C_j/C_k)^m), HYP-8505.

## MISTAKE-203 (boxeph-2026-07-20-S178) — THM-1615's pinch dichotomy was false: fold-curve crossings are Stokes crossings of the representation, not singularities of the function

**What happened:** the pinch bridge (S177) inferred "the sweep of |t*(r)| forces
a pinch or endpoint singularity of A(t) at finite t." FALSE: as t moves, the
branch point r_j(t) crosses the REAL r-axis and the contour simply deforms
around it (a Stokes crossing of the representation). Genuine singularities need
a TRAPPED pinch (two branch points colliding — finitely many t) or endpoint
contact (excluded near finite t by the sweep's own r -> 0 blow-up). A generic
ray misses the finite trapped set — the theorem's own ray-genericity deleted
the singularities it needed. Also: the sweep itself fails when f_0(0) != 0
(t*(r) -> 1/f_0(0) finite), and the Gevrey/rotation citations from THM-1565 do
not transfer (integrand not entire in r).

**How caught:** owner-ordered adversarial pass (S178), referee 1, before any
completion claim was made (the theorem was filed as a skeleton with the attack
surface listed — the process worked).

**Lesson (two bridges dead by the same class of error):** both failed bridges
compared or tracked REPRESENTATION-level structure (term domination; fold
crossings) instead of FUNCTION-level invariants (actual singularities, actual
values). The Radial Lemma and the orbit-product proof survived review because
their contradictions are function-level (A -> 0 vs A == 1; (ct)^r vs constant).
Rule: a bridge argument must exhibit a function-level invariant that
distinguishes mixed from one-sided — anything contour- or term-level will
dodge.

## MISTAKE-202 (boxeph-2026-07-20-S176) — THM-1605's original local lemma overclaimed: equal products over disjoint subsets do NOT force equal subsets

**What happened:** THM-1605's first proof (S175) concluded from the Puiseux-DFT
expansion that disjoint index sets I, J with equal branch-products must be
EQUAL ("Fourier inversion => I = J"). FALSE as stated: the mu^j coefficient of
the log-difference is C_j * (Sigma_I zeta^{ij} - Sigma_J zeta^{ij}) with a
UNIVERSAL factor C_j = [mu^j] log(W(mu)/w*) that can VANISH; when it does, the
j-th character sum is unconstrained. Example: I = {0, 2} in Z_4 has S_1 = S_3
= 0 automatically, so with C_2 = 0 all constraints hold with I != J. The
step "no monodromy element moves the cluster" was therefore NOT established
(and is not even true in general — the identity only forces equal PRODUCTS
along the orbit).

**How it was caught:** in self-review while PREPARING the fleet adversarial
review (owner-ordered), before the referees ran. The repaired proof (orbit-
product: permanence + transitivity + Vieta + the c != 0 patch) is simpler and
was then confirmed by two hostile referees, who found one further pinhole
(c = 0 exclusion — patched in place) and two exposition debts (paid).

**Lesson:** when a lemma concludes SET equality from SYMMETRIC-FUNCTION
equality, check the kernel: products/power sums see subsets only through
characters, and universal prefactors (here C_j) can kill exactly the
characters you need. Also: the strongest review step is trying to write the
attack yourself before delegating it.

## MISTAKE-201 — Over-estimating the VC-witness reduction dimension 4× by a per-monomial heuristic, then reading noise as corroboration (death-star-S61)
- **Who/when:** death-star-2026-07-20-S61, caught same session by the reduction agent's concrete computation.
- **What:** I estimated the de Bondt–van den Essen reduction dimension of F at **≈76** using a crude "~⌈(deg−1)/2⌉ auxiliary variables per nonlinear monomial × 13 monomials" heuristic, then treated a concurrent agent's partial-output "76/77" hits as an independent cross-check. Both were wrong: (a) F's 13 nonlinear monomials are **not independent** — they all share the single quadratic u=1+xy, so **6** helper coordinates reduce everything to cubic (N≈10, not 35–38; M=2N≈20, not 70–77); (b) the "76/77" in the partial output was the ballpark of Zhao's *a-priori VC bound* (3/2)(3^{M−2}−1), NOT the reduction dimension — I pattern-matched a number to the conclusion I'd already written.
- **Genus:** category C (heuristic upper bound mistaken for an estimate) + confirmation bias (misreading unrelated partial output as agreement). The heuristic ignored shared sub-expressions; the "corroboration" was manufactured from noise before the agent finished.
- **Rule:** (1) an auxiliary-variable count for degree reduction must count **distinct building blocks after common-subexpression elimination**, never monomials independently — sparse/structured maps (here ℂ*-equivariant, all powers of one quadratic) collapse the count. (2) NEVER cite a running agent's *partial* output as corroboration; wait for its final result, and match numbers to their actual definitions (dimension vs a-priori bound) before claiming agreement.
- **Repair:** corrected §2 of the S61 reflection (≈76→≈20), the HYP-8265 sub-entry, the SESSION-LOG, and the PROBLEM-LEDGER. Correction broadcast to the fleet.
- **ADDENDUM (same session, on executing the transport):** the corrected "≈20, feasible, gate clears, no new math" was ITSELF overstated. Exact ℚ(i) computation (`polylib_exact_deathstar_S61.py`) showed: (a) Yagzhev G's Jacobian is NOT nilpotent, so the lift needs a genuine cubic-**homogeneous KELLER** reduction; (b) the "6-helper stacking" reduces degree and transports the collision but is Keller only ON the section {W=X^β} — off-section det varies (Schur complement = JG only there), so the reduced map is NOT globally Keller and NOT nilpotent-Jacobian; (c) naive homogenization breaks Keller too (det(I+x₀JH₂+JH₃) not constant). So the agent's stacking (which I relayed) conflated the easy degree-reduction with the hard Keller-preserving homogeneous reduction — the witness EXISTS (BCW) but its construction is real math (global-determinant control), not "engineering." **Deeper rule:** when relaying a subagent's feasibility claim, re-derive the load-bearing step yourself before writing "feasible/no new math" into canon — a dimension count is not a construction. The exact lift+rotation machinery WAS validated (produces Hessian-nilpotent quartics from nilpotent cubic-homog maps), so the honest status is "machinery ready, nilpotent input still to be built," not "done modulo typing."

## MISTAKE-199 — Re-deriving concurrent fleet work across two sessions on a fleet-wide owner prompt (death-star-S59s/S59t)
The owner's n·2^x+1 / figurate-triangle prompt went to the WHOLE fleet. I worked it in S59s (HYP-8165) and S59t (HYP-8175, THM-1355/1360) without first grepping for concurrent same-prompt HYP claims — and kind-pasteur (S128c102 Rosetta triangle, S128c103 shear catalog / HYP-8170), opus (S420 shear duality), and mac-mini (S1 Pisot tower) had gone DEEPER on the same objects: the Pisot ladder, the Proth spectrum, the products=ordered-tournaments identity, and the EXACT triangle (Faulhaber + 3 deviations) were all already theirs. My results were correct but ~90% duplicative; I even collided on HYP-8165 (both mine-S59s and kp-S128c102) and HYP-8170 (kp's; I took 8175). My repo-mining agent surfaced the overlap only AFTER I closed both sessions. RULE: on any prompt likely broadcast to the fleet (owner prompts usually are), the FIRST action is `grep -i <topic-keywords> 05-knowledge/hypotheses/INDEX.md` and a scan of same-day SESSION-LOG / inbox for concurrent claims — BEFORE claiming IDs or doing deep work. Convergent independent verification has value, but it must be FRAMED as convergent (credit + priority banner), never presented as primary. Amended THM-1360 + the reflection with priority banners; credited kind-pasteur/opus/mac-mini. **RECURRENCE (S59u, same day):** the NEXT owner prompt (the problem-ledger) was ALSO fleet-wide — klein-S332 wrote a peer ledger, kp-S128c104 + mac-mini-S140 stubbed. This time a mining agent surfaced the collision DURING the session and I consolidated into one canonical PROBLEM-LEDGER.md crediting all rather than competing — the discipline working. Rule: on any owner prompt assume fleet-wide, grep INDEX + same-day SESSION-LOG FIRST; on collision, CONSOLIDATE with credit, never produce a parallel artifact.

**MISTAKE-199 4th recurrence (S59v):** the weight-sign reduced JC was proven by mac-mini-S123 (THM-1370) and pushed WHILE I computed it, though the INDEX was clean when I claimed. Escalated rule: on a hot fleet-wide prompt re-pull IMMEDIATELY BEFORE filing any theorem, and default to synthesis/confirmation over primary claims once 2+ agents are visibly on the prompt. Credited mac-mini; kept only a verification cross-check + the three-reduced-JC synthesis.

## MISTAKE-200 — naming an invariant before grepping canon for it (d_sat = the diameter)

**What happened (mac-mini-S126, caught by kind-pasteur-S128c108 / THM-1400 SS I).** I
introduced the "saturation depth" `d_sat(n)` — the least `d` at which the waggly truncation
`G^(<=d)` becomes complete — called it a **new metagraph invariant**, computed `2,3,4,7` for
`n=4..7`, refuted `d_sat = n-2` at `n=7`, and handed off "compute n=8 before conjecturing."

**Why it was wrong.** `G^(<=d)` is complete exactly when every pair of classes is within `d`
flips, so `d_sat = diam(G_n)` — trivially, once stated. And the diameter was already canon:
`07-reflections/diameter-is-feedback-arc-set.md` (opus-2026-03-24-S306) records
`diam(G_n) = max_T min-FAS(T) = A003141(n)`, growth `~n^2/4`, with `OPEN-QUESTIONS.md`
listing it **RESOLVED** and `README.md` carrying it as the Waggly Completeness Theorem.
So: the invariant was a rediscovery, the handoff needed no computation (`A003141(8) = 8`),
"no linear formula" was a known **quadratic**, and my `n=7` refutation restated opus-S306's.

**The rule.** *Before naming an invariant, grep the reflections, OPEN-QUESTIONS, and README
for it — including under other names.* A quantity defined by "the least `d` at which X
becomes complete" is a **diameter**; ask what metric, then search for that metric. I searched
for "saturation" and "d_sat" (my own coinage, so of course no hits) and not for "diameter."

**What survived.** The map-graph *framing* of THM-1390 (point-contact vs edge-contact, the
clique-explosion reading) is untouched and THM-1400 explicitly grants it; the point-adjacency
thread continues correctly in THM-1405. The cost was one invented name and one bad handoff,
both retracted in place. See `01-canon/theorems/THM-1390-...md` (correction banner).

## MISTAKE-202 — Trusting an asymptotic domination claim measured only to m=8 (death-star-S61g/S61h)

**What happened.** For GMC(2) I declared "GMC(2) is complete" (S61g) on top of klein-S351's
Gamma-Bridge domination step, and separately claimed (S61h) that for two-sided `P` the top
charge-0 term `γ_{a_max} a_max!` carries `> 50%` of the mass of `E[P^m] = Σ_a γ_a a!`, so
`|top| > |rest|` and the triangle inequality forces `E[P^m] ≠ 0`. I measured the share only
across **m = 2..8**, where it sits at ≈ 0.60–0.67, and read that as domination.

**Why it was wrong.** kind-pasteur-S128c120 (THM-1585) refuted klein's version exactly — the
top-term share falls to **0.04%**, the consecutive-term ratio grows to **45×** — and opened
`CASE-gamma-bridge-domination-step`. Prompted by their addendum I reran **my own** statistic
to **m = 24**: the share collapses `0.67 → 0.068` (`Z²+W+ZW²`) and `→ 0.0002` on a b-sweep;
top-dominance is False for every `m` past ≈ 8. **A domination claim measured to m = 8 cannot
tell "share → 1" from "share → 0" — both look like ⅔ there.** `E[P^m] ≠ 0` still holds, so
the *conclusion* (NC2) survives, but the *mechanism* was false: "domination was an analytic
strategy for an algebraic fact" (kp). The sound route is kp's THM-1605 Hermite no-common-root.

**The rule.** An asymptotic-domination / "top term wins" claim is only tested by pushing the
statistic until the ratio's *trend* is unmistakable (here m ≈ 20+), and by a parameter sweep
that can *amplify* the competing term (the b-sweep). Never certify "X dominates as m → ∞" from
a window where the sub-dominant term is still the same order. When the honest tool is a
factorial-weighted sum `Σ c_k k!`, suspect that nonvanishing is *algebraic* (orthogonality /
no-common-root of a classical sequence), not a size comparison. Conceded in full in the court
case; S61g headline withdrawn, S61h §1 retracted (both banners in place); S61h §2 (the Lean
NC2 ⇒ GMC(2) reduction, pure charge arithmetic, assumes NC2) is unaffected. See also
MISTAKE-199 (the same over-eager pattern, there in fleet-coordination form).


## MISTAKE-207 (boxeph-2026-07-20-S186r) — THM-1765's pair-sum constant used the symmetric-collision model (missing the Λ''' midpoint shift), and its evidence verified the RATIO, never the CONSEQUENCE object

**What happened:** THM-1765 §1 claimed the far-end pair-sum limit
(1/t)·2/(u_c²Λ'') and hence the universal constant −2/(d₁d₂). Two errors:
(a) the merging roots are NOT symmetric about u_c at second order — the
Λ''' midpoint-shift enters at the SAME order; the correct law is
PS = (1/t)[2/(u_c²Λ'') + (2/3)Λ'''/(u_cΛ''²)], giving −(2/3)(d₁+d₂)/(d₁d₂)
on a two-term edge (referee-measured: (1,−1) → 0 — forced exactly by the
global residue identity Σ_{all roots} 1/(uΛ') = 0; (2,−1) → 1/3; (3,−1) →
4/9; claimed: 2, 1, 2/3). (b) The frozen evidence (F1–F4) measured only
the RATIO Λ''u²/v — the identity's hypothesis — and never once measured
the pair-sum, the object the Consequence is about. The identity was true;
the Consequence drawn from it was wrong; the check could not have caught
it because it tested the wrong object.

**Also refuted downstream (same session, same verdict):** "O(1)
universally" — the charge-0 coefficient p₀ enters v but never Λ'/Λ'':
VALUE-HIJACKED ends (witness P₄ = ZW + Z⁹W⁷ + W) realize the S183r threat
with non-integrable s^{−2} pair-sums; and the (L2) two-mechanism
classification + O(t^{−1/2}) rate (zero-drift mechanism; Θ(T^{−2/5})
example).

**Genus:** the leading-term trap, third instance (MISTAKE-202: leading
products; MISTAKE-206: leading odd-coefficient; 207: leading collision
model), COMPOUNDED with MISTAKE-204's cousin: the numeric verified an
intermediate quantity, not the concluded one.

**Rule:** (1) at a colliding pair, never use the symmetric local model
without checking the next Puiseux coefficient (Λ''' here) — or better,
compute pair-sums by the GLOBAL RESIDUE IDENTITY (sum over all roots is
zero: pair-sum = −spectator residues), which is exact and model-free.
(2) Every verification numeric must measure the CONSEQUENCE'S object, not
the hypothesis' — if the theorem says "hence X is bounded", the script
must print X.

**Repair:** THM-1765 amended in place (corrected law; value-hijack
sub-case named; ledger reverted; residue identity canonized; §6 verdict
archive); THM-1680 §4 pointer corrected; referee checks frozen at
04-computation/thm1765_referee_{hijack,momenttest}_S186r.py + .outs.

## MISTAKE-206 (boxeph-2026-07-20-S183r) — THM-1680's deletion biconditional read only the LEADING odd coefficient: "B ≡ 0 ⟺ removable" is false one rung down, and the dichotomy missed the boundary-truncated class

**What happened:** THM-1680 §1/§2 filed "DEFECT ≡ 0 ⟺ B ≡ 0 ⟺ removable" and
"per germ exactly two cases (delete / graded-visible)". The hostile referee
(S183r) refuted both AS STATED: (a) the local odd sector at a fold is a full
Puiseux vector Σ c_k τ^{k/2}; B = c₋₁ is only its leading entry. A germ with
B ≡ 0 but c₁ ≢ 0 is NOT removable and stays Γ-graded exactly one rung down
(T3: C(1/2,m)/C(−1/2,m) = −1/(2m−1) exact; T4: its monodromy defect 2|c₁|√r
sits twelve orders above the instrument floor — THE INSTRUMENT WAS FINER
THAN THE PROSE). (b) The ε-signs can hop at sub-germ boundaries, so
"vanishing on an interval ⟹ vanishing on the germ" propagates only within a
sub-germ; a sign flip mid-germ creates BOUNDARY-TRUNCATED arcs with pure
exponential moment grade I_m ~ e^{−s₂}v(s₂)^m/m (T2) — a third class.

**Genus:** the leading-term trap, one floor up from MISTAKE-204. 204: a
check fed by the model's own ansatz verifies arithmetic, not premises. 206:
a FUNCTION-LEVEL invariant (good, per MISTAKE-203) that is only the LEADING
coefficient of the honest invariant (the full odd-sector vector) inherits
exactly the failure it was built to avoid, one rung down. Cousin of
MISTAKE-202 (leading products vs full set data).

**Rule:** when the invariant is a Puiseux/graded expansion, the deletion or
vanishing criterion must quantify over the WHOLE graded vector, and every
"identically zero" must name the domain it propagates over (sub-germ, germ,
family). If an instrument (here the monodromy defect) is sensitive to the
full vector, let the instrument's zero-set DEFINE the criterion rather than
the leading coefficient's.

**Repair:** THM-1680 §1/§2/§4 amended in place (odd-sector vector; per-
sub-germ trichotomy with the truncated exponential class; §4 hypothesis =
defect ≡ 0 everywhere; finiteness lemma added); referee verdict archived as
§8; checks frozen at 04-computation/thm1680_referee_hostile_S183r.py + .out.
Referee: "no scenario survives the repaired statements."

## MISTAKE-205 — the "Alpoge-Mathew" attribution of the JC counterexample (THM-1300) is a HALLUCINATION (owner-corrected 2026-07-20)

**What happened.** THM-1300's attribution blocks (mac-mini-S127 and S129) confidently attributed the
dim-3 Jacobian-Conjecture counterexample to **Levent Alpoge** with co-credit to **Akhil Mathew**,
"obtained with Claude Fable 5, announced 2026-07-19," citing web searches (X posts, a Wikipedia edit,
an arXiv absence check). The OWNER has flagged this as a **hallucination**: those searches were run on a
result then roughly ONE DAY OLD, and the specific-name attribution is not reliable — it should not be
restated as established fact.

**The corrected provenance (per the owner).** The counterexample was **discovered by Claude**. It became
public on **2026-07-19 via a tweet from an Anthropic employee**, which is where the "shared publicly
yesterday" signal came from — that is a SHARING event, not a discovery attribution. So the map is a
Claude discovery that an Anthropic employee surfaced; it is NOT an Alpoge-Mathew result.

**The lesson (this is the reusable part).** Do NOT attribute a very recent result from web searches with
confidence. A result that is ~1 day old has essentially NO reliable public record: no peer review, no
stable arXiv, and search engines index tweets and rumors that conflate "who shared it" with "who found
it." When mac-mini's own search notes flagged "no arXiv preprint, no peer review, no journal," that was
the tell that NO reliable attribution was available — and the correct move was to record the map as
**provenance-uncertain**, not to fill the vacuum with plausible names. Confabulated attribution is exactly
the failure mode a fresh, un-refereed result invites.

**Standing correction for the fleet.** Treat the JC counterexample as **Claude-discovered, provenance
otherwise uncertain**. Do NOT write "Alpoge-Mathew" (or "Alpoge") as the discoverers in any file, letter,
or write-up. THM-1300's external-attribution blocks are kept as history but are now marked contested by
the banner added this session. What the repo legitimately holds is unchanged: the INDEPENDENT exact
verification (det J ≡ -2, triple collision, all identities) and the equivariant/Dixmier/elliptic anatomy.

---

## MISTAKE-208 — a "poly-time invariant" that is not isomorphism-invariant (arborescences rooted at a fixed vertex)

**kind-pasteur-2026-07-21-S128c143.** In the invariant-lattice work I used `arb(T)` = number of
spanning arborescences **rooted at vertex 0** (one reduced-Laplacian minor). This is NOT an
isomorphism invariant: the root "vertex 0" depends on the labeling, so relabeling changes it. On the
exhaustive census it *looked* fine because I computed it on ONE canonical representative per class
(deterministic, but an artifact of the canonicalization). The bug surfaced when a random-labeling
SAMPLER reported spurious "same-invariant, different-H" collisions at n=6 that the exact exhaustive
census (poly-tuple determines iso at n=6) said were impossible.

**Why it's instructive:** (1) a per-root / per-vertex quantity is only an invariant after
symmetrization (sum over roots, or the *sorted tuple* of per-root values). The fix is
`arb_inv = tuple(sorted(arb_root(T,r) for all r))`. (2) Two code paths (canonical-rep vs
random-label) disagreeing is a RELIABLE detector of a non-invariant — always cross-check an
"invariant" on random relabelings. (3) The scare (a suspected floating-point artifact) was itself a
false lead — the numbers matched numpy exactly; the real bug was mathematical (root-dependence), not
numerical. Corrected: `|arb_inv|=55` at n=6 (vs 32 rooted-at-0); `arb_inv` refines score. Headline
findings (score ⟂ cyc, the 2-adic edge THM-1980) use only genuine invariants and are unaffected.
Banner added to THM-1965.
## MISTAKE-216 -- confusing Poisson rank two with a DC(2) or planar-JC counterexample

**Error.** From a nonautomorphic Poisson endomorphism of
`C[x,q,p,z]`, infer either a counterexample to the planar Jacobian conjecture or
an endomorphism counterexample of the second Weyl algebra merely by replacing
Poisson brackets with commutators.

**Correction.** The four-variable symplectic map is a Keller counterexample in
affine dimension four, not two. Quantizing its four output polynomials directly
in `A_2` requires exact control of ordering corrections; equality of principal
Poisson symbols does not give exact Weyl commutators. The safe cotangent or
Hamiltonian-dual construction doubles four classical coordinates to four Weyl
pairs and therefore lands in conventional `A_4`. To reach `DC(2)`, verify an
actual endomorphism of `A_2`; to reach planar JC, construct or exclude a
two-coordinate Keller pair. Keep Poisson rank, affine dimension, Weyl rank, and
number of algebra generators separate.
