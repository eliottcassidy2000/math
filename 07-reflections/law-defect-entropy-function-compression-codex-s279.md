# Law-Defect Entropy And Function Compression

Anchors: HYP-3201, HYP-3200, HYP-3162, HYP-3199, HYP-3161, HYP-3160, HYP-3154, HYP-3153,
HYP-3152, HYP-3151, HYP-3150, HYP-3147, HYP-3146, HYP-3142, HYP-3140,
HYP-3132, HYP-3122, HYP-3109, HYP-3092, THM-577, T1301, LTI-301,
LTT-201, OPEN-Q-108.

## One Sentence

An algebraic law is a zero-entropy compression theorem: it says the target
function survives a quotient, and every failure of the law is a named sidecar.

## The Upgrade

The previous function-compression lane had one clean example: commutativity.
Forget order in `(a,b)`.  Then `a+b` and `a*b` survive, while `a^b` and `b^a`
do not.  The new move is to make this universal:

```text
commutativity failure  -> ordered sidecar
associativity failure  -> bracketing sidecar
idempotence failure    -> multiplicity sidecar
distributivity failure -> context sidecar
group-action failure   -> representative/action sidecar
root-circle failure    -> off-circle/root-variance sidecar
moment-law failure     -> cumulant sidecar
```

This is more useful than saying "the quotient loses information."  It tells us
which kind of information was lost and whether the target proof function is
allowed to proceed.

## Information Theory Reading

For a compression `q:X->Y` and target `f:X->Z`, the theorem condition is:

```text
H(f | q) = 0.
```

When `H(f | q)>0`, the quotient has not merely become philosophically lossy;
it has a measurable residual.  In the scout, subtraction leaves `0.8` bits of
bracket debt after forgetting parentheses, exponentiation leaves `0.515625`
bits, and the K4 class-action quotient leaves `0.701205` bits of action debt.

So the LRC packet can now ask a precise question: is the next scalar being
used after a zero-entropy quotient, or is it silently spending sidecar entropy?

This does not conflict with the incoming HYP-3160 warning that ordinary row
entropy is not the k=8 extremal scalar, nor with HYP-3161/HYP-3200 correcting
the old `1/7` associativity-defect smell.  HYP-3201 uses conditional entropy
as a quotient-defect diagnostic, not as the quantity to maximize over rows.

HYP-3162 adds the sharper root-sidecar warning: the lawful cap is the
7th-cyclotomic/Joukowski ideal, while the dip is a real-rootedness and
rational-approximation defect.  So the root quotient must retain the
cyclotomic target as well as the off-circle spread.

## Lee-Yang / Pascal / Phi4

The Pascal cap is the good compression.  Pair-normalized Pascal mass is the
exchangeable, commutative, second-moment face; de Moivre-Laplace is what this
face looks like after the usual Gaussian compression.  That is why the
binomial row feels lawful.

The dip is the failed compression.  It is what remains when the pairwise
Pascal shadow cannot carry the target function.  In root language, a genuine
Lee-Yang circle quotient compresses six root radii to one `R` and gives
`q0=q6*R^6`.  Off-circle spread is then exactly the root-radius sidecar; in
HYP-3122 language it is the phi4 correction.  The quartic does not decorate
the pair mass; it pays the entropy bill that the pair mass cannot pay.

With HYP-3162 in the web, that bill has a number-field label too: the cap is
trying to approximate the 7th-cyclotomic de Moivre angles, and the dip is the
defect left when the rational/Pascal packet does not land on that cubic ideal.

## K8

The k=8 row is now less mysterious and more demanding.  The even part is
solvable:

```text
(t-1)(t-2)(t-4)(t-5) = u^4 - 5u^2 + 4.
```

That is the degree-4 / quadratic-in-`u^2` page.  But the odd Worpitzky side is
larger on the known consecutive row, so the proof cannot close by pointing at
the even page alone.  The even page is a legal compression only for the even
target.  For `L_y`, cap dip, or endpoint data, the odd sidecar has to be
shown zero, reconstructed, bounded, or routed to named debt.

## Monoid Correction

The user's correction about the flip action matters.  A true `V4` action would
have zero action entropy on the visible classes.  The fixed-path K4 quotient
does not.  The exact class relations are:

```text
flip a: T->{+}, +->{T}, (-)->{S}, S->{-,S}
flip b: T->{-}, (-)->{T}, +->{S}, S->{+,S}
flip c: T->{S}, +->{S}, (-)->{S}, S->{T,+,-,S}
```

So the quotient is a transformation/relational monoid packet.  The original
bit cube has a clean action; the class quotient has action entropy.  This is
exactly the difference between "looks like a group after compression" and
"has a canary/representative sidecar that makes the action legal."

HYP-3199 makes this guardrail concrete: the fixed-path table is a cover with
an `S` fiber, while the partial-score `x,y` Einheit chart is the legal section.
Law-defect entropy measures the cost of using the cover as if it were already
the section.

## Ear Decompositions

Ear language also becomes a compression law.  Strong connectivity compressed
to "has an ear decomposition" preserves reachability structure.  Factor-
critical compressed to "has an odd ear decomposition" preserves a parity
sidecar.  For LRC/Omega, the useful question is not whether ears are a nice
analogy, but which target survives the ear quotient:

```text
root-collision event?
odd-cycle/Omega membership?
zero-real component connectivity?
boundary crossing parity?
```

If the quotient forgets which root pair crossed the real segment, or which odd
cycle supplied the parity, then ear existence alone is another scalar shadow.

## Working Rule

Before using a compression in the proof packet, write:

```text
law_id
quotient_map
target_function
H(target | quotient)
sidecar_type
discharge_status
```

This makes the web of ideas less misty.  Pascal/binomial, Lee-Yang circle,
phi4 dip, k=8 quartic, Worpitzky odd side, K4 monoid action, and Newton/
Maclaurin moment inequalities are all asking the same thing: which quotient is
zero-defect for the function we are actually trying to prove?
