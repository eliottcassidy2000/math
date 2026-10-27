# First-hit guards and the two integer clocks

Session73 additive scope supplement. **PROVED elementary clarification**, not
an extension of the frozen core's claim to native ROOT positivity. The core
packet remains byte-for-byte unchanged. Its report is SHA256
`554e4ba747b77247242a78d2c914df575d694f0b00211694eac51cdb333fbd10`.

## Ordinary valuation cylinders are not killed first-hit ROOT cylinders

In core P2, U(n)=(3n+1)/2^v2(3n+1) is the ordinary odd-accelerated map on
all positive odd integers, including U(1)=1. A specified valuation word
k_1,...,k_h, with S=sum k_j, has the inherited exact cylinder

    n=c+2^(S+1)t,    U^h(n)=d+2*3^h t,    t>=0,

and its odd-index endpoint congruence modulo q has minimal source odd-index
period `2^S q/gcd(3^h,q)`, with the already stated image-divisibility guard.
This theorem preserves the valuation word and the *residue* of the endpoint.
It does not say that the word is legal under a protocol which stops at the
first visit to1, nor that the endpoint is literally1.

For a strict first-hit ROOT word from n>1, retain all extra conditions

    U^j(n) != 1 for 0<=j<h,   U^h(n)=1.              (FH1)

A killed protocol regards the operation after the first arrival as unavailable;
an absorbing convention may pad the trajectory, but those idle/padded states
are not further original ROOT moves. These conventions must not be exchanged.
The inequalities/equality in(FH1) are nonperiodic realization guards, excluded
from the periodic-congruence conclusion of core P1/P2.

**Exact fixed-word boundary.** The terminal equality U^h(n)=1 can occur in
this positive cylinder only for t=0 and d=1. Therefore a nonempty fixed word
has a strict first-hit realization iff d=1, c>1, and the finite orbit of c
has no earlier1; if so the unique realized source is c. A first-hit source
cannot be promoted to all congruent aliases by preserving the branch word.

**Proof.** The displayed endpoint is d+2*3^h t, with positive odd d>=1 and
t>=0. Equality to1 forces d=1 and t=0. Earlier first-hit checks on c are then
literal finite checks. More generally every prefix j satisfies

    U^j(c+2^(S+1)t)
      = U^j(c)+3^j*2^(S-S_j+1)t,  S_j=sum_(i<=j) k_i.

For t>=1 this is strictly greater than1. Thus within a full positive
valuation cylinder, an earlier visit to1 can only occur at its canonical
source t=0; removing that source, when necessary, is not a residue-class
condition on the infinite cylinder. This argument uses the entire specified
word to justify the common affine branch, not a guessed branch of a nearby
integer. QED.

**Every finite endpoint modulus has a same-word alias.** If the canonical
source c is a strict first-hit realization (so d=1), then for any q>=1 take
`t=q/gcd(3^h,q)>0`. The endpoint odd indices are0 and3^h t, hence equal
modulo q, while the latter literal endpoint is greater than1. Their source
odd indices differ by exactly `2^S q/gcd(3^h,q)`, the minimal P2 precision;
the valuation word and both integer clocks are unchanged. Thus a finite
endpoint modulus alone cannot certify exact first hit, even on a correctly
identified branch. An exact lift, a sufficient height bound, or the explicit
first-hit guard is needed. Any further finite source congruences can be
preserved by taking a suitable positive multiple of t. This is an endpoint
alias statement for a fixed word, not a claim that its higher sources never
reach1 at a later time.

### Literal hostile controls

- The one-letter word(4) has c=5,d=1. Sources5 and37 share its branch and
  endpoint odd-index residue0 modulo3, but their endpoints are1 and7.
  Only5 has first hit at the end of this one-letter word. No assertion is
  made that37 is globally unrooted.
- Keeping the native arithmetic source congruence n=3 modulo6 does not repair
  this: word(6) gives sources21 and405, both3 modulo6, with endpoints1 and19.
  They share endpoint odd-index residue0 modulo9 and source odd-index residue
  10 modulo192. The ordinary word clock is7, but only21 first-hits at that
  one-letter endpoint.
- Word(4,2) is an ordinary U-word at5 and133, with endpoints1 and19.
  At5 the first step already reaches1, so the second step is an illicit
  continuation in a killed first-hit protocol. At133 the intermediate value
  is25 and the two-step ordinary prefix remains pre-hit. The two sources
  have the same P2 word and endpoint odd-index residue0 modulo9.

If source weights encode a first-hit law, weighted cell sums must incorporate
these very guards before any Gram test or pushforward. Replacing them by the
larger ordinary valuation cylinder changes the measure/observable and can
invalidate a paid forcing argument. The core Gram theorem remains valid for
the actual guarded tests and original positive source masses; it supplies
neither the missing first-hit indicator nor its positive value.

## Shortcut time S is not literal operation time h+S

Define the shortcut parity-step map C(n)=n/2 for even n and C(n)=(3n+1)/2
for odd n. One odd U-step of valuation k takes exactly k C-steps: the odd
C-step performs the first halving, then k-1 ordinary even halvings remain.
The same U-step takes k+1 literal operations if `3n+1` and each individual
halving are counted separately. Hence a word has the distinct integer clocks

    shortcut parity-step clock = S,
    literal 3n+1/halving clock = h+S,
    odd-accelerated step count = h.                  (FH2)

All three add under an actual legal concatenation. Modulus refinement changes
none of them. A sign defined by one tick of the shortcut clock is (-1)^S;
using literal time instead multiplies that phase by (-1)^h. Thus importing a
signed ROOT current whose phase uses S requires retaining h as well; it is
incorrect to silently substitute (-1)^(h+S). This is a clock-conversion
identity, not a new theorem identifying every ROOT current with that phase.

Example5 reaches1 by shortcut path5,8,4,2,1 (four steps), while its literal
path is5,16,8,4,2,1 (five operations). The sign difference is real. At21 the
corresponding counts are6 and7. At1, a terminal/killed first-hit protocol has
already stopped; ordinary U(1)=1 is an artificial continuation relative to
that protocol, so(FH2) must not be used to charge additional ROOT moves there.

## Provenance, exact controls, and limits

This is a direct specialization of the source-pinned native-return N6 cylinder
and clock-even C2/C3 conventions already archived in the frozen core. It adds
no new assumption about the native signed-current formula. `check_scope.py`
checks the three explicit hostiles and all84 nonempty words of lengths1–3
with k_i∈{1,2,3,4}, their eight lifts t=0,...,7, every prefix, and both clocks
by separate literal and shortcut implementations. It tests5 strict first-hit
words in that finite universe and160 same-word finite-modulus aliases at
q=1,...,32. Normal and -O receipts are identical.

No first-hit existence theorem for every source, no ROOT weight positivity,
no stopping-time estimate, and no transfer from ordinary congruence legality
to first-hit realization is asserted.
