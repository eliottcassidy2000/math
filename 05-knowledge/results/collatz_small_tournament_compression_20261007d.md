# Small marked carry windows for the eight-bit receipt

**PROVED:** the marked four-/six-point decoder, ordered gluing law, and the
unexpanded two-parameter representation of the inherited eight-bit receipt.
**FINITE-EXACT:** the declared windows, motif collision, and literal controls.
**OPEN:** grounding the emitted child with exponent 924745897 modulo 2^32.
This package changes the proof interface, not the set of paid source guards.

## Inheritance and the information that must survive

The closest proved construction is the marked prefix-carry chart in
[frieze_collatz_reset2_20261007](frieze_collatz_reset2_20261007.md): positive
valuation words give increasing rational points, with the final clock retained
separately. [THM-4602, tournaments as chirotopes and frieze strips](../../01-canon/theorems/THM-4602-tournaments-as-chirotopes-pfaffian-criterion-and-legendre-frieze-strips.md)
explains why the signs of their ordered minors form a transitive tournament.
The present operation adds local four-/six-point ports and an exact ordered
composition interface to that inherited chart.

The hostile quotient is already visible in
[tournament_recursive_four_state_20261004](tournament_recursive_four_state_20261004.md):
the fixed-Hamiltonian-path masks `ab` and `c` represent the same strong
four-vertex class, but toggling the named arc `a` gives different classes.
The class label is not a sufficient state for future named flips. A faithful
four-symbol ordered-join codec does exist there; it is a storage convention,
not an identification of tournament operations with Collatz steps.

There is a second, six-vertex loss. The inherited canonical masks 344 and 345
in [HYP-3049, edge-perspective extension](../hypotheses/HYP-3049-a000568-edge-perspective-extension.md)
are distinct isomorphism classes. Both have score sequence `(2,2,2,3,3,3)`
and 43 Hamiltonian paths. The present exact check also gives the same complete
four-vertex induced-motif census, in the prior `(T,+,-,S)` convention:

\[
                         (1,2,2,10).
\]

Their ordered-pair sector-size and internal-sector decks agree, whereas the
cross-sector orientation decks differ. This reproduces only these two
inherited representatives; it is not a new enumeration of all six-vertex
classes. The missing datum is an oriented relation across the marked ports.

The corrected near miss is to treat an unweighted motif or tournament class
as the ordered arithmetic word. The least-used useful sidecar is the final
scale, together with the source-owned terminal parameter. The live concepts
are prefix carries, weighted minors, ordered gluing, repeat macros, exact
source guards, and terminal proof obligations.

## A small intrinsic chart and its decoder

For a positive valuation word `w=(a_1,...,a_r)` put

\[
 F_w(n)=\frac{P n+B}{Q},\qquad P=3^r,\quad Q=2^{a_1+\cdots+a_r}.
\]

For each prefix i, mark

\[
 z_i=B_i/3^i,\qquad z_0=0,\qquad
 \Delta_i=z_{i+1}-z_i=2^{a_1+\cdots+a_i}/3^{i+1}.
\]

Thus `z_0,...,z_r` are strictly increasing, `Delta_0=1/3`, and

\[
                    3\Delta_i/\Delta_{i-1}=2^{a_i}
                    \quad(1\le i<r).                 \tag{1}
\]

A window with r=3 or r=5 has four or six vertices. Its vertices are actual
chronological prefix carries, its pairwise observable is `z_j-z_i`, and its
orientation is positive precisely when i<j. There are no ties and no chosen
orientation gauge once chronology is marked. Equivalently, the columns
`(1,z_i)` have these differences as their ordered 2-minors.

**Decoder.** Retain the exact rational point tuple and the last valuation.
Check `z_0=0,z_1=1/3`, positivity of every gap, and that all ratios (1) are
positive powers of two. Their exponents recover the first r-1 letters; the
marked last letter completes the word. Re-encoding is an exact membership
check. This proves an iff, not merely a decoder on generated examples.

The final letter may instead be the typed expression `c+d`, d>=0, where c
is one shared deferred valuation. It is not silently given an arbitrary
value. The point tuple is independent of the last letter: for example
`(1,1,1)` and `(1,1,4)` give the same four points. Their final scales differ.
Deleting that sidecar loses an actual denominator and changes source guards.

All these sign tournaments are transitive. The useful arithmetic is in the
weights and marks, rather than a choice among four or six tournament types.
Four-/six-point windows are convenient small interfaces; they are not claimed
to be the unique or bit-optimal representation.

## Ordered gluing preserves the full affine map

Write a window summary as

\[
                  (z,s)=(B/P,Q/P),\qquad F_w(n)=(n+z)/s.
\]

If v follows w chronologically, direct substitution gives

\[
        (z,s)*(z',s')=(z+s z',s s').                   \tag{2}
\]

This operation is associative. It preserves the full carry and scale of
concatenation, and hence its affine endpoint and native congruence. The
ordered window list still retains the local word for exact export; the
summary is not used as permission to change that list. Reassociating (2)
does not change chronology. Permuting factors generally changes z even
when their lengths and total valuations agree.

For an ordinary concrete word the exact formal odd endpoint guard is
`Pn+B = Q mod 2Q`. A strict first-hit export must additionally exclude early
ROOT. The present application inherits that check from the authenticated
eight-bit receipt; an affine identity by itself is not a ROOT certificate.

## An unexpanded representation of the eight-bit family

Use the fixed heads from
[collatz_eight_bit_completion_20261007c](collatz_eight_bit_completion_20261007c.md):

```text
a = (2,2,2,1,2,9,16)
b = (2,2,4,1,3,3,1,3,3,1,1,1,2,2,3).
```

The source address is the exact pair `(e,t)` with e>=68 and positive odd t,
representing `n=2^(e+1)t-1`. Its proposed smaller child is
`h=(n+1)/256-1`. The actual path words are

\[
          1^e,a,c\quad\hbox{and}\quad 1^{e-8},b,c+2,    \tag{3}
\]

where c is read at the source's actual state after its fixed head. The
source tail `a,c` splits into letter lengths 3+5, hence vertices 4+6. The
child tail `b,c+2` splits into lengths 5+5+3+3, hence vertices 6+6+4+4.

Introduce typed formal parameters

\[
                         T=(2/3)^e,\qquad C=2^c.
\]

The repeat macro `1^e` has summary `(1-T,T)`, and the child macro `1^(e-8)`
has summary `(1-(3/2)^8 T,(3/2)^8 T)`. A last letter `c+d` makes only its
window scale depend on C, by a factor `2^d C`. Consequently every fixed
window and the full gluing identity can be stored as sparse exact rational
polynomials in T,C without expanding e repeated letters.

Applying (2) to the marked lists gives the exact formal identities

\[
             s_R=s_L/256,\qquad z_R=(z_L+255)/256.       \tag{4}
\]

Since `256h=n-255`, (4) proves `(n+z_L)/s_L=(h+z_R)/s_R`.
The terminal offset c+2 is essential to this equality.

Let `F_a(x)=(Px+B)/Q`, Q=2^34, and let
`r=(Q-B)P^(-1) mod 2Q`. The source-owned guard checked by the compressed
packet is exactly

\[
                   2\,3^e t-1\equiv r\pmod{2Q}.         \tag{5}
\]

This uses modular exponentiation, not construction of n or T's enormous
denominator. The inherited e>=68 proof makes every source prefix through
a grow above n. Its actual terminal c and the exact child integrality then
give strict first-hit paths (3), with an endpoint of 1 allowed only last.
The child is positive and strictly below the same supplied n.

For `(e,t)=(924745904,1)`, (5) holds and the symbolic child is exactly
`2^924745897-1`. The package audits (4), (5), and the ordered fixed heads;
it does not materialize this integer or assert that it reaches ROOT.
The packet's terminal status is explicitly **OPEN**. A forged ROOT status,
an altered exponent phase, a reordered head, or a request exceeding the
explicit materialization cap is rejected.

The mathematical compression is the repeated one-run plus a shared terminal
parameter: the recipe has O(log e + log t) source-address size, fixed head
data, and a deferred c. For fixed t, including the Mersenne choice t=1, the
address contribution is O(log e); binding c additionally stores that actual
integer. Exact
evaluation of T itself would again need a denominator with Theta(e) bits,
and literal export would restore the long word. There is no assertion that
arbitrary words have such a short expression, or that a finite tournament
alphabet removes the cost of grounding the open child.

## Exact API and controls

`encode_window` and `decode_window` expose the marked local iff.
`glue` composes exact polynomial summaries in chronological order.
`make_packet(e,t)` and `audit_packet` check the fixed D8 template, source
guard, and formal identities without expanding the source. `materialize`
binds c by the inherited actual receipt at the identical bounded source and
checks both decoded full words. Its default source bit cap is 4096; this
is an implementation bound, not a mathematical limit on (3).

Reproduce from the repository root:

```text
python -B -X utf8 04-computation/experiments/collatz_small_tournament_compression_20261007d.py
python -B -O -X utf8 04-computation/experiments/collatz_small_tournament_compression_20261007d.py
```

The 3,467 exact controls include all 1,088 words of lengths 3 and 5 over
letters 1..4; 64 ordered three-window compositions checked both ways and
against the independent full carrier; last-scale and reordered-word
hostiles; the two inherited six-vertex representatives and all their
four-vertex subgraphs; sixteen literal D8 receipts at runs 68,69,96,128;
the unexpanded astronomical packet; and nine malformed/forged requests.
Normal and optimized execution agree. The decoder and all-height symbolic
identity are proved above, rather than inferred from this finite universe.

The next arithmetic task is to supply a guarded continuation or a ROOT proof
for the marked child. The present interface preserves the information that
such a continuation needs and makes no additional coverage claim.
