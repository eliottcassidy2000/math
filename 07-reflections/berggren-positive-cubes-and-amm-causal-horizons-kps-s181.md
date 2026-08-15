# Berggren positive cubes and AMM causal horizons

Status: **SYNTHESIS OF PROVED/VERIFIED RESULTS.**  Proof sources are THM-3370,
THM-3373, THM-3374, and THM-3375.  The cross-frontier comparisons below are
typed research heuristics, not theorem transfers.

## What the attached U-spine note contributed

The identities

```text
Q_r=(2r+1)^2+2=2C_r+1,
P_r=(2r+1,2r(r+1),2r(r+1)+1),
U P_r=P_(r+1)
```

were already present in audited THM-3334/3370.  The attachment was still a
useful prompt because it emphasized three coordinates that a scalar sequence
hides:

1. the parabolic branch parameter `r`;
2. the plane/height coordinate `C_r`; and
3. the positive-chamber coordinate of a second representation.

The third is the one that unlocked new mathematics.  Writing a cube pair as

```text
d=x+y,       e=x-y
```

turns positivity into the literal window `0<|e|<d`.  THM-3370's fixed-`d=11`
Pell orbit leaves that window forever.  THM-3375 repairs the object rather
than the technique: it lets `d` grow on a second Pell coordinate and keeps
`e/d` on an interior limiting ray.

## The moving-coordinate principle in two unrelated frontiers

The same research move appeared independently in AMM 12592.

- A width-five slack-kernel decomposition exists at `R=8` (THM-3373).
- The `R=512,D0=0` rule-A death has a 277-bit fatal constant.
- Five terminal rows have only 42 bits of universal Lucas-box coefficient
  capacity, and even 57 rows are insufficient (THM-3374).

Thus “width five” is a lawful operation scale but not a terminal repair
scale.  A successful flow must start moving slack by row `49`, long before
the death.  The analogy with the cube problem is methodological only:

| source | frozen coordinate that fails | restored moving coordinate |
|---|---|---|
| Berggren/cubes | fixed `d=x+y` | Pell-growing `d_j=h_j^2+2` |
| AMM | five-row patch at death | early accumulated slack/history |

No arithmetic map connects the two carriers.  What transfers is the hostile
question: *is the coordinate held fixed by the proposed proof exactly the
resource that must grow?*

## Typed connection ledger

### U-spine to positive cubes

```text
source:    parabolic depth r and scalar Q_r
map:       Q_r=x^3+y^3 -> (d,e,q) with d=x+y,e=x-y
kept:      scalar, depth, witnessed positive pair
lost:      order after quotienting by e -> |e|; Berggren ancestry from Q alone
sidecar:   0<|e|<d and the Pell divisibility h=29u
test:      THM-3375's exact orbit identities and ratio bound
```

### Gaussian U-spine fibre to Eisenstein cube fibre

```text
source:    C_r=N_Z[i]((r+1)+ri)
target:    q=N_Z[omega](x+y omega)
map:       only the affine scalar relation Q_r=2C_r+1
kept:      r and Q_r
destroyed: prime-factor labels, since gcd(C_r,Q_r)=1
sidecar:   (x,y,d,e,q), not a generic CM-field label
test:      the coprimality identity in THM-3370
```

The discriminant-`-8` form `[d,2a,q]` in THM-3375 is a third norm carrier.
Class number one gives a useful parametrization, but it does not identify the
Gaussian and Eisenstein toggles.

### Deletion-transform lesson

THM-3372 shows that diagonalizing a marked deletion deck loses factor order,
while a skew contraction restores one exterior response.  The cube analogue
is exact at the level of information loss: the scalar `Q_r` forgets the sign
and chamber of `e`; retaining `(d,e)` restores positivity.  There is no common
algebraic operation, so this is a representation-selection lesson rather than
an implication.

### LRC boundary

The Pell ray supplies a positive current in an arithmetic sense, but not an
LRC current.  It has no runner row, clock, selected arc, owner, phase, or
global exit.  THM-3354's typed-carrier no-go remains the correct guardrail.

## New exact targets

1. **Extend the slope compiler atlas.**  For
   `d=n^2u^2+2`, `a=mu*d+nu`, the remaining equation is

   ```text
   3W^2=n^2(4m^2-n^2)u^2+4(2m^2+2mn-n^2).
   ```

   THM-3376 now completes the primitive parity-correct range `n<=29`: exactly
   `(m,n)=(14,23),(26,29)` survive, and both give infinite positive rays.  The
   `(14,23)` seed is enormous, a hostile to bounded seed searches.  Continue
   beyond `29`, requiring four exact gates: an integral seed, a norm-one unit
   preserving odd parity/divisibility, an invariant positive ratio, and a
   modular obstruction certificate for every rejected slope.  Optimize the
   regulator only after those gates.

2. **Count all positive intersections.**  THM-3375 proves a logarithmic lower
   family, not an asymptotic.  The universal-torsor system

   ```text
   pt-rs=+/-1,
   3e^2=4(s^2+2t^2)-(p^2+2r^2)^2
   ```

   is the right starting point.  The ambient local singular series is not
   killed: THM-3375 proves the mod-`63` filter is already the full local
   obstruction for `X^3+Y^3=Q_r`.  The distinguished Pell ray is much thinner
   (`r mod 63` only `38,51,60` and `r=2 mod 5`), so orbit scars must not be
   mistaken for ambient obstructions.

3. **Test fixed-sum uniqueness.**  Fixed `d` supports finitely many positive
   candidates, but uniqueness is open.  Pell gap estimates should be tested
   against exact divisor data before a global conjecture is promoted.

4. **Move the AMM atom search earlier.**  Translate labelled width-five atom
   shapes only into matching degree-step words at rows `<=49`, retain every
   cap slack, and optimize accumulated entry-cone margin rather than immediate
   survival.  THM-3374 rules out spending more time on rows `102..106`.

5. **Preserve marked coordinates.**  In every proposed cross-frontier map,
   carry `(d,e)` for cubes, labelled slacks for AMM, and owner/phase for LRC.
   Scalar norms, total slack, and unlabelled periods are discovery views only.

## Bottom line

The U-spine attachment did not reveal a missing Berggren identity; it pointed
back to a missing coordinate.  Once the positive chamber was retained, the
fixed-sum obstruction became immediate and a moving-sum Pell ray proved
positive infinitude.  The simultaneous AMM calculation delivered the hostile
dual: a local operation may be real and exact while the target defect lies far
outside its causal capacity.  In both cases, the next proof starts by moving
the resource earlier, not by sharpening the final local step.
