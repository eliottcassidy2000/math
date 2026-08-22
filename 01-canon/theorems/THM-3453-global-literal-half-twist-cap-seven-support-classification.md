---
id: THM-3453
title: "Global transverse literal half-twist cap-seven support classification"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  The theorem classifies
  every modulus supporting a
  transverse strict literal half-twist cover by at most seven owners.  It is
  an exact Boolean-mask result whose labelled witnesses realize the fixed
  common centre c=1/(2Q) with zero mode cochain.  It has no arbitrary-centre
  classification, nonzero current, decrement, or LRC(14) consequence.
source: root-global-literal-half-rank-seven-2026-08-15
audit: >
  two independent immutable-package audits CLEAN: line-by-line theorem
  derivation checked period descent, transversality and deduplication, the
  mixed-prime and fixed-fibre dichotomies, predecessor directions, the
  17-adic invoice, exact universe, and scope; a separate clean-room strict-mask
  solver matched the divisor support for every 2<=Q<=200, independently
  recovered all atom ranks, and closed Q=366 and Q=578 by fibre, capacity,
  and invoice controls; normal/-O/stored replay, dependency pins,
  semantic/hash, AST/security, ID, routing, documentation, and diff gates
  passed
depends_on:
  - THM-3405-common-centre-gcd-gauge-and-boolean-half-twist
  - THM-3416-zero-mode-cochain-global-rank-six-support
  - THM-3434-seventeen-fibre-two-sided-mass-closure
  - THM-3435-dyadic-fibre-grid-decomposition-for-literal-half-twists
  - THM-3445-prime-even-half-twist-cap-seven-classification
  - THM-3451-prime-quarter-half-twist-cap-seven-classification
related:
  - THM-3429-prime-fibre-activity-descent-for-mixed-order-half-twist-seven-covers
script: 04-computation/lrc_global_literal_half_twist_cap7_support_thm3453.py
output: 05-knowledge/results/lrc_global_literal_half_twist_cap7_support_thm3453.out
script_sha256: fb0147ba4f5ee6e918ba8509667daa149e2a1ee835b51fe7a1ff7b35bae5f146
output_sha256: 4c22fe1f62507d40e3fb6b8bdb8e7bc14ad9c52cd3f228e8c823f803f06c84a1
semantic_sha256: 5f5d2c29d47a15abc8cf74d2b0f40769c220c3a646b4ba4c6bb21ca5d49cc1fa
hash_basis: LF-normalized bytes
---

# THM-3453 -- global transverse literal half-twist cap-seven support classification

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**  The proof and frozen
exact companion passed two independent immutable-package audits, including a
separate clean-room strict-mask reconstruction through `Q=200` and hostile
closures at `Q=366,578`.

## 1. Candidate statement

For `Q>=2`, a coefficient `r` modulo `2Q`, and a sheet `ell` modulo `Q`, put

```text
B_(Q,r)={ell: ||r(2ell+1)/(2Q)||<1/14}.                 (1)
```

Only **transverse** owners are admitted:

```text
Q does not divide r.                                    (2)
```

After sign normalization one may take `1<=r<Q`; the self-opposite class
`r=Q` is empty.  Owners are distinct, and repeated or sign-equivalent masks
may be deleted.  Let `rho_H(Q)` be the least number of transverse masks in a
literal cover of every sheet.

The candidate support antichain is

```text
D={8,9,10,11,12,13,14,15,23,25,29,38,51,68,148}.       (3)
```

The claimed all-modulus classification is

```text
rho_H(Q)<=7  iff  some d in D divides Q.                (4)
```

The exact atom ranks are

```text
rank 4:  8,9;
rank 5:  10,12;
rank 6:  11,15,23,25;
rank 7:  13,14,29,38,51,68,148.                         (5)
```

This is a Boolean strict-mask statement at a specified phase.  Directly from
the definitions,

```text
B_(Q,r)=D_(Q,r)(1/(2Q)).                                (5a)
```

After retaining the owner label, the containing THM-3398 mode has
`h=r mod 2Q`; its centre lift with `n=0` is `x=1/(2Q)`.  Thus every displayed
cover has a fixed common physical time and zero complete mode cochain.  The
theorem does not classify any other centre, supplies no nonzero
runner-current or decrement certificate, and does not prove LRC(14).

## 2. Inheritance and the positive direction

The proved inputs are:

- [THM-3416](THM-3416-zero-mode-cochain-global-rank-six-support.md), which
  classifies cap-six support as multiples of `8,9,10,11,12,15,23,25`;
- [THM-3434](THM-3434-seventeen-fibre-two-sided-mass-closure.md), which
  classifies odd cap-seven support and supplies the atoms `13,29,51`;
- [THM-3435](THM-3435-dyadic-fibre-grid-decomposition-for-literal-half-twists.md),
  which gives the exact dyadic fibre chart and direct witnesses at
  `14,38,68,148`;
- [THM-3445](THM-3445-prime-even-half-twist-cap-seven-classification.md),
  which classifies target-free `Q=2p`; and
- [THM-3451](THM-3451-prime-quarter-half-twist-cap-seven-classification.md),
  which classifies target-free `Q=4p`.

Every `d` in `(3)` has the direct witness frozen in the companion.  If `Q=hd`
and `(s_i)` covers `d`, then direct cancellation gives

```text
B_(Q,hs_i)=pi^(-1)B_(d,s_i),       pi:Z/QZ -> Z/dZ.     (6)
```

Transversality is preserved because `Q` divides `hs_i` exactly when `d`
divides `s_i`.  Hence every multiple of an atom is positive.  The exact-rank
claims in `(5)` follow from the displayed witnesses together with THM-3416
for the cap-six boundary and THM-3434 for the odd rank-seven atoms.  In
particular, none of `14,38,68,148` contains a cap-six atom, so their seven-owner
witnesses have exact rank seven.

The closest mechanism for the converse is THM-3435's labelled dyadic fibre
grid.  Its canonical descent hostile is the `Q=14` partition: bare projection
to the odd core `7` loses the doubled radius and the dyadic orientation.  The
repair below peels a whole factor `2^a p` on one exact fibre rather than
projecting that structure away.

## 3. Strong-induction period descent

For a selected owner define

```text
m_Q(r)=Q/gcd(Q,r),       L=lcm_r m_Q(r),       h=Q/L.   (7)
```

Every `m_Q(r)` divides `L`, so `h|r`; write `r=hs`.  Since `L|Q`, reduction
of sheets modulo `L` and cancellation in `(1)` give

```text
B_(Q,hs)=pi_L^(-1)B_(L,s).                              (8)
```

Moreover `Q` divides `hs` exactly when `L` divides `s`.  Thus descent
preserves transversality, literal OR coverage, and owner count; coincident
descended masks may be deduplicated.

Suppose, for contradiction, that `(4)` is false and choose the least positive
counterexample `Q`.  No member of `D` divides `Q`.  If the joint period in
`(7)` were proper, `(8)` would give a smaller counterexample; induction would
then put an atom of `D` in `L`, hence in `Q`.  Therefore

```text
L=Q.                                                     (9)
```

Odd `Q` is impossible by THM-3434.  Since `8` is forbidden, an even
counterexample has

```text
Q=2^a R,       R odd,       a in {1,2}.                 (10)
```

The direct strict masks on `Q=2,4` are empty, so take `R>1`.  THM-3416 and
the absence of its eight atoms force the family to have exactly seven owners.

## 4. One odd prime: the activity dichotomy

Fix an odd prime `p|R` and write `Q=pN`.  Call an owner `p`-active when
`p` does not divide `r`, and inactive otherwise.  Joint period `(9)` supplies
at least one active owner.

For an inactive coefficient `r=pu`, exact cancellation gives

```text
B_(Q,pu)=pi^(-1)B_(N,u),       pi:Z/QZ -> Z/NZ.         (11)
```

If at least one inactive owner exists, there are at most six of them.  Their
descended masks cannot cover `N`: otherwise THM-3416 would put a cap-six atom
in `N`, hence in `Q`.  Choose a base sheet `y` missed by every inactive mask.
On its `p`-point fibre `ell=y+Nt`, an active owner has phase

```text
r(2(y+Nt)+1)/(2pN)
 =r(2y+1)/(2pN)+rt/p.                                  (12)
```

The second term is a complete translated `p`-grid.  An open circular arc of
length `1/7` contains at most `ceil(p/7)` of its points.  At most six active
owners must cover the missed fibre, so

```text
p<=6 ceil(p/7).                                         (13)
```

The exact odd-prime solutions are

```text
p in {3,5,11,17,23,29}.                                (14)
```

This is the **mixed** branch.  Notice that `(12)--(14)` charge only active
complete `p`-grid sections.  They never substitute an ambient
`ceil(Q/7)` bound for an inactive lower-order pullback; that false move is
MISTAKE-392.

## 5. The exact odd-cofactor fixed fibre

It remains to treat the branch in which every owner is `p`-active.  Put

```text
D_p=2^a p,       M=R/p,       Q=D_p M.                 (15)
```

The cofactor `M` is odd, so `y_0=(M-1)/2` is an integer.  On the complete
`D_p`-point fibre

```text
ell=y_0+Mt,       t in Z/D_p Z,
```

one has `2y_0+1=M` and therefore

```text
r(2ell+1)/(2Q)=r(2t+1)/(2D_p).                         (16)
```

Thus restriction to this one fibre is **exactly** the literal mask
`B_(D_p,r)`, including its dyadic coset and orientation.  Since `p` does not
divide `r`, `D_p` cannot divide `r`; every nonempty local owner is transverse.
After harmless deduplication, the global cover therefore gives an at-most-
seven transverse cover on `D_p`.

Formula `(16)` repairs the apparent absence of a reflection-fixed `p`-fibre
when `N` is even.  The two nearest bare `p`-fibres have opposite residual
phases and do not give a prime literal half twist.  Enlarging the fibre by
the full dyadic factor makes the remaining cofactor odd and restores the
fixed sheet.  This is a componentwise fibre identification, not a claim that
a multi-component grid union maps bijectively; it respects the sign and
degree correction in MISTAKE-393.

## 6. Elimination when `v_2(Q)=1`

Let `Q=2R`.  If all owners were `p`-active, `(16)` would give a transverse
cover of `2p`.  The target-free hypothesis of THM-3445 is automatic: a factor
`p` among `5,11,13,23,29` would already put `10,11,13,23,29` in `Q`.
THM-3445 would then force `p=7` or `19`, putting `14` or `38` in `Q`.
Therefore every prime divisor of `R` is mixed.

Combining `(14)` with the forbidden atoms leaves

```text
prime_support(R) subset {3,17}.                         (17)
```

The atom `9` gives `v_3(R)<=1`, while `51` forbids simultaneous factors `3`
and `17`.  Hence either `R=3`, or

```text
R=17^b.                                                 (18)
```

At `Q=6`, the union of every transverse sign representative reaches only
the two sheets `{1,4}`.  The boundary `b=1`, namely `Q=34`, is negative by
THM-3445.  It remains to exclude `b>=2`.

## 7. The `2*17^b` two-sided invoice

Write

```text
Q=2*17^b=17N,       N=2*17^(b-1),       b>=2.          (19)
```

The `17`-coordinate is mixed.  On a missed fibre each active owner supplies
at most three points, by the exact complete-seventeen-grid law.  Therefore at
least six owners are active.  Since an inactive owner exists and there are
exactly seven owners, the profile is

```text
six 17-active owners and one 17-inactive owner.          (20)
```

Let `h` be the number of base sheets covered by the unique descended inactive
mask.  Its quotient order `m>1` divides `N`.  The possibilities are:

- `m=2`, whose strict block is empty;
- `m=17^c`, whose density is at most `3/17` by the quotient-order estimate in
  THM-3434; or
- `m=2*17^c`, `c>=1`.  After cancellation its coefficient is coprime to `m`,
  so its phases form a complete `m`-grid and its size is at most
  `ceil(m/7)`.  Since `m>=34`,

  ```text
  ceil(m/7)/m <= (m+6)/(7m) <=3/17.
  ```

Consequently

```text
h<=3N/17.                                               (21)
```

Every active seventeen-section has at least two points as well as at most
three.  On each of the `N-h` missed fibres, active incidence is at least `17`;
on each of the `h` covered fibres, the same six blocks spend at least `12`.
Thus

```text
17N-5h <= sum_(six active i)|B_i|.                      (22)
```

For an active owner, `gcd(Q,r)` divides two, so its quotient order is `Q` or
`Q/2=17^b`.  Charging the exact quotient and pullback multiplicity gives

```text
|B_i|<=2 ceil(17^b/7).                                  (23)
```

This order-sensitive charge is the correction-lineage safeguard from
MISTAKE-392.  Combining `(21)--(23)` and using `17^b=17N/2` yields

```text
274N/17
 <=12 ceil(17N/14)
 <=102N/7+72/7.                                        (24)
```

After multiplication by `119`, `(24)` implies

```text
184N<=1224.                                             (25)
```

But `N>=34`, a contradiction.  The failed boundary at `b=1` is genuine:
the scalar version of `(24)` is then inconclusive, which is why `Q=34` is
handled by THM-3445 rather than silently absorbed into the tail.

## 8. Elimination when `v_2(Q)=2`

Let `Q=4R`.  Every prime in the mixed list `(14)` already forces an atom:

```text
3 ->12,   5 ->10,   11 ->11,   17 ->68,   23 ->23,
29 ->29.                                                 (26)
```

Therefore every selected owner is `p`-active for every `p|R`.  Applying the
fixed fibre `(16)` gives a transverse cover of `4p`.  Its target-free premise
in THM-3451 is automatic: the excluded primes

```text
3,5,7,11,13,19,23,29
```

would respectively force one of `12,10,14,11,13,38,23,29` into `Q`.
THM-3451 therefore gives `p=17` or `37`, which forces `68|Q` or `148|Q`.
This final contradiction closes `(10)` and completes the proof of `(4)`.
**QED.**

## 9. Connection, loss, and hostile ledger

| field | exact content |
|---|---|
| source | a labelled transverse literal cover on `Q=2^aR`, `a=1,2` |
| target | for each `p|R`, either inactive masks on `Q/p` plus an active `p`-grid invoice, or one literal cover on `2^a p` |
| maps | `ell -> ell mod Q/p` in `(11)`; `t -> (M-1)/2+Mt` in `(16)` |
| preserved | strict endpoints, literal OR, owner count, `p`-activity, and on the fixed fibre the complete dyadic coset/orientation |
| destroyed | behaviour off the chosen fibre and the active affine lift character over the remaining base |
| required sidecars | the prime-activity vector and each inactive owner's exact quotient order/pullback multiplicity |
| cheapest strict hostile | `p=7`, where one active section has capacity one, not two |
| descent hostile | the `Q=14` partition versus the uncoverable bare odd-core layer `Q=7` |
| affine hostile | the `Q=51,p=17` lift character from THM-3429; local prime data do not reconstruct a global cover |
| tower boundaries | direct `Q=6`, prime-even `Q=34`, first tail row `Q=578` |

The faithful carrier is an activity-labelled fibre clutter.  There is no
intrinsic pairwise orientation and hence no reason to manufacture a
tournament.

## 10. Exact companion and audited status

Run from the repository root:

```bash
PYTHONHASHSEED=0 python3 -B 04-computation/lrc_global_literal_half_twist_cap7_support_thm3453.py
PYTHONHASHSEED=1 python3 -B -O 04-computation/lrc_global_literal_half_twist_cap7_support_thm3453.py
```

Both transcripts must equal
`05-knowledge/results/lrc_global_literal_half_twist_cap7_support_thm3453.out`
byte for byte.  The standard-library companion pins all five proved
dependencies, forbids `assert`, dynamic execution, and unexpected imports,
and freezes:

- a Fraction-based strict-mask reference through `Q=32`;
- period pullback through `Q=96`;
- the fixed odd-cofactor identity for both dyadic depths and every odd core
  through `149` (`103,488` owner rows and `11,582,552` sheet cells);
- all fifteen atom witnesses and independent no-cap exact ranks;
- the exact mixed-prime list and both dyadic residual ledgers;
- the order-sensitive `17`-adic density and invoice boundaries; and
- no-cap negative searches at `Q=6,34,578`, the last visiting `10,413`
  states and `10,412` branches.

These computations verify the stated interfaces and hostile boundaries; they
do not replace the strong-induction proof.  Two independent immutable-file
audits checked the proof, scripts, output, hashes, routing, correction
lineage, and a separate clean-room bounded/hostile path.  The theorem is
therefore **PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED**.  It does not
classify arbitrary centres, provide a nonzero current, or settle LRC(14).
