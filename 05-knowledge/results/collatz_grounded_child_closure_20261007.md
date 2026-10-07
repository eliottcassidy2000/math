# A grounded family closed under guarded inverse child words

**PROVED:** the parameter isometry, unique guard refinement, finite-word
closure, and exact residue readers below. **FINITE-EXACT:** independent
literal controls and the two supplied large child certificates. **OPEN:**
covering arbitrary supplied integers by these grounded families. Closure of
a represented family is not closure of every obligation emitted by a selector.

## 1. Inheritance and the representation being extended

The closest proved mechanism is the completed terminal lift in
[the terminal-lift note, section 2](collatz_terminal_lifts_20261007.md).
It uses an actual supplied first-hit ROOT word, not an assumed child proof.
The new issue is whether further inverse child guards can be imposed without
losing that grounding. The canonical hostile is an illegal fixed source:
the child map `(2n-1)/3` cannot be applied to 3. Finding a different member
of a family where it is legal does not extend the supplied source 3.
The corrected near miss is to infer grounding from a reversible return to
an ungrounded larger core. The least-used sidecar is the **unchanged terminal
lift parameter**, retained alongside the whole inverse word and its carry.

Anchor: close a genuinely grounded family under the new child maps.
Niche: use ternary isometry to solve every finite guard exactly.
Wildcard: exploit stabilization to read long-certificate sources without
expanding them. The concept board is the seed certificate, parameter,
inverse word, ordered carry, exact phase, and first-hit terminal.

Let the supplied seed be an odd integer `a>1`, with first-hit ROOT word `w`
of length `L` and cost `A`. Write its carrier as

\[
F_w(n)=\frac{Pn+B}{Q},\qquad P=3^L,\quad Q=2^A,
\qquad Pa+B=Q.
\]

The inherited completed family is

\[
N_s=a+\frac{4Q(4^{Ps}-1)}{3P},\qquad s\in\mathbb Z_{\ge0}.       \tag{1}
\]

At `s=0`, the exact certificate is `w`. For `s>0`, its certificate is
`w` followed by the single valuation `2(1+Ps)`. Indeed,

\[
F_w(N_s)=\frac{4^{1+Ps}-1}{3}>1.
\]

The perturbation of the source is a multiple of `4Q`, so every native
valuation in `w` remains exact. Every intermediate state is its old state
plus a strictly positive even integer when `s>0`; there is no early ROOT.
The final displayed valuation reaches 1. This recalls the inherited proof
to make the grounding premise explicit, rather than rebranding it as new.

## 2. Exact ternary isometry and one phase for every word

**Theorem 1.** For distinct nonnegative integers `s,t`,

\[
v_3(N_s-N_t)=v_3(s-t).                            \tag{2}
\]

**Proof.** Assuming `s>t`, factor the difference in (1):

\[
N_s-N_t=\frac{4Q\,4^{Pt}(4^{P(s-t)}-1)}{3P}.
\]

The leading factors are ternary units. Lifting the exponent gives
`v3(4^(P(s-t))-1)=1+L+v3(s-t)`, and the denominator has valuation
`L+1`. This proves (2). In particular, for every `j>=0`, the map
`s mod3^j -> N_s mod3^j` is well-defined and injective, hence bijective.
Its unique compatible extension to the ternary integers is an isometry.
No assertion about the ordinary size of `N_s` follows from that metric. ∎

Let `v` be any finite tuple of positive valuation letters, with carrier

\[
F_v(n)=\frac{pn+b}{q},\qquad p=3^\ell,\quad q=2^{\sum v}.
\]

Allow the empty word; for nonempty `v`, restrict here to `q<p`, so its
inverse is a smaller-child operation. Its inverse family is

\[
H_v(s)=\frac{qN_s-b}{p}.                          \tag{3}
\]

**Theorem 2.** There is exactly one class

\[
s\equiv s_v\pmod {3^\ell}                       \tag{4}
\]

on which (3) is integral. Every nonnegative parameter in (4) gives a
positive odd integer, all actual valuations in `v` are legal, and

\[
F_v(H_v(s))=N_s,\qquad H_v(s)<N_s
\]

for nonempty `v`. Concatenating `v` with the certificate of (1) gives a
first-hit ROOT certificate for every such family member.

**Proof.** The necessary and sufficient congruence is
`N_s=b*q^(-1) mod p`; (2) gives its unique preimage. Final integrality
forces integrality at each backward step, since every power of 2 is a unit
modulo 3. Starting at the positive odd endpoint `N_s`, the formula
`(2^e y-1)/3` gives a positive odd integer whenever it is integral.
Induction therefore gives positivity, oddness, and exact forward valuations
throughout. An early state 1 would force all later actual odd states to stay
at 1, contrary to `N_s>=a>1`. Finally `q<p` and `b>0` give strict decrease
from `N_s` to the child. The known terminal suffix supplies grounding. ∎

This yields an exact ternary digit compiler: start with `s=0 mod1`;
at stage `j`, try its three lifts modulo `3^j` and keep the unique one
satisfying the target congruence at that precision. It uses no orbit search.

The implementation additionally retains the direct strict positivity test
`qN_s>b`. It is redundant on the legal native phase by the preceding
backward argument, but is useful as an independently derived boundary
check. It must not be advertised as an extra empirical gap in Theorem 2.

## 3. Closure by refinement, with the supplied point kept fixed

**Theorem 3.** If `u` is another contracting inverse word (or the empty
identity word), then

\[
\frac{q_uH_v(s)-b_u}{p_u}=H_{uv}(s),
\]

where `uv` means chronological **forward return** concatenation: first
`u`, then `v`. The parameter phase for `uv` reduces to the phase for `v`
modulo `3^len(v)`. Every such finite extension is a nonempty infinite
grounded family. Its source is strictly smaller than `H_v(s)` when `u`
is nonempty; the empty extension preserves the source exactly.

Here `(p_u,q_u,b_u)` is the carrier of the new forward return word `u`.
The proof is exact affine composition, Theorem 2, and the fact that the
new inverse operation contracts. In particular, a legal composite has a
legal old endpoint, so its phase refines the old phase.

The old elementary maps are included:

\[
G_1(y)=\frac{2y-1}{3},\quad
G_5(y)=\frac{8y-5}{9},\quad
G_{17}(y)=\frac{2048y-2363}{2187},
\]

with return words `(1)`, `(1,2)`, and `(1,1,1,2,1,1,4)` respectively.
See [the mixed normal-form package](collatz_child_normal_forms_20261007.md)
for their disjoint guards and full-carrier decoding. The theorem above
allows any finite contracting valuation word, not just these generators.

There are two distinct APIs:

- `refine_family` computes the exact subfamily on which an extension is
  legal. It is allowed to restrict the parameter set.
- `extend_point` takes a supplied family member and retains its **same**
  parameter. It returns no certificate if that parameter fails the new
  guard. It never selects another parameter to make the extension succeed.

For instance, seed 5 at parameter 0 takes `G1` to 3. Applying `G1` again
at that same point is rejected. The new family phase still exists, but
its members are different integers. This is the minimal useful hostile
against a mistaken universal-closure inference.

For every fixed word the legal parameter density is exactly `3^(-ell)`.
Adding `u` keeps exactly `3^(-len(u))` of that phase. These are parameter
frequencies, not the density of grounded source integers. For a fixed
seed and word, (1) and (3) grow exponentially in `s`, so their source set
has at most `O(log X)` members below `X`. Arbitrarily many finite positive
phase measures do not give a refinement-independent positive lower bound;
along unbounded word length they tend to zero. Nor may zero densities of
individual families be summed over a countably infinite family of words.

### Exact recursive reseeding, and its sharp nonzero-phase boundary

There is a stronger closure when the native phase is `s=0 mod p`, where
`p=3^ell`. Let `c=H_v(0)`. Its inherited ROOT word is `v+w`, with
`P'=Pp` and `Q'=Qq`. Direct substitution gives

\[
H_v(pt)=c+\frac{4Q'}{3P'}(4^{P't}-1).             \tag{3a}
\]

The **entire** child family is precisely the same canonical completed-lift
construction based at its newly grounded child. Thus it can be reseeded
recursively without changing any represented integer or its actual ROOT
word. The reversible parameter transport is `s=p*t`; the original point
and inverse word retain its provenance. `reseed_zero_phase` implements
this exact identity, rather than choosing a new fitting family.

For example, the seed-5 family after `G1` at phase `s=0 mod3` becomes
exactly the seed-3 family at parameter `t=s/3`. Both large children in
section 5 have zero phase, so both admit the same lossless reseeding.

The phase condition is sharp for **this canonical reseeding**. Suppose
instead its least phase representative is `r>0`, and put `c=H_v(r)`.
The inherited ROOT word at `c` includes the extra final valuation
`2(1+Pr)`, so `P'=3Pp`, `Q'=4Qq*4^(Pr)`. Writing

\[
k=Pp,\qquad K=\frac{4Qq4^{Pr}}{3Pp},
\]

the old child family and the canonical completed lift at `c` are
respectively

\[
c+K(4^{kt}-1),\qquad
c+\frac43 K(4^{3ku}-1),\qquad t,u\ge0.
\]

Equality would require `3*4^(kt)+1=4^(3ku+1)`. For `t>0` the left
side is 1 modulo 4 and the right side is 0; for `t=0`, equality forces
`u=0`. They therefore share **only their base point**, even though every
member of both families is grounded. This does not obstruct the general
phase-retaining closure of Theorem 3. It prohibits silently replacing its
nonzero phase by a different canonical lift.

The smallest seed used by the control is `a=3`, `w=(1,4)`, `v=(1)`.
Its phase is `1 mod3`, its first child is `828503`, and the actual ROOT
word is `(1,1,4,20)`. This is a concrete hostile to source-family recycling
and explains why the terminal boundary is part of the representation.

## 4. Long certificates need only the precision being queried

The parameter phase reader can be independent of the seed word length once
that length exceeds the requested ternary precision. Put

\[
W_L(s)=\frac{4^{3^Ls}-1}{3^{L+1}},\qquad w_L=W_L(1).
\]

Exact expansion gives

\[
w_{L+1}=w_L+3^{L+1}w_L^2+3^{2L+1}w_L^3.          \tag{5}
\]

Consequently, for `j>=1` and `L>=j-1`, `w_L` is already constant modulo `3^j` and
the binomial expansion gives

\[
W_L(s)\equiv s\,w_{j-1}\pmod {3^j}.              \tag{6}
\]

Thus the huge exponent `3^L*s` need not be materialized in a modular
exponentiation for such a query. For `L<j-1` the direct modular formula
retains the denominator precision `3^(L+1+j)`. A binary reader instead
inverts `3P` modulo the requested power of two and tests whether the
power of 4 has already vanished there.

For the inverse source (3), division by `3^ell` requires reading `N_s`
at `ell+j` ternary digits **before** division. Discarding those extra
digits loses precisely the required guard/quotient information.

An optional bounded materializer has an explicit bit cap. Symbolic source
residues and finite ROOT words do not require it. The word is finite but
may contain a final valuation of enormous magnitude; compact storage is
not a bound on the length of a literal halving trajectory.

For completeness, the defensive positivity cutoff can also be computed
without expansion. Let `C=4Qq` and `H=3P(b-qa)+C`. Then `qN_s>b` iff
`C*4^(Ps)>H`. If `floor(H/C)>=0`, this is equivalent to

\[
2Ps\ge\operatorname{bitlength}(\lfloor H/C\rfloor).
\]

The strict equality boundary is therefore exact. Intersect its least
parameter with (4); native legality already implies that this cut adds
no excluded legal positive member.

## 5. The current non-Mersenne children are inside a grounded instance

Use the already authenticated seed `a=2^1457-1` and its stored first-hit
word (`L=7347`, `A=13102`), from
[the preceding reset-two rules](collatz_reset2_rules_20261007.md).
At unchanged parameter 0, the inverse words `(1,2)` and
`(1,2),1^5,(1,2)` give respectively

\[
h=\frac{8\cdot2^{1457}-13}{9},\qquad
z=\frac{2048\cdot2^{1457}-29455}{19683}.
\]

The exact certificates have `(odd steps, total valuation)` equal to
`(7349,13105)` and `(7356,13113)`. The script directly replays both.
Their phases are `s=0 mod9` and `s=0 mod19683`, yielding infinite
grounded extensions of these two points. At positive parameters, the
one extra final valuation completes each inherited word.

This construction proves grounding directly from the 1457 seed. The
[paid-exit package](collatz_mixed_child_rules_20261007.md) gives a different,
stronger size-ranked interface: both points can be paid through the
strictly smaller seed `2^1451-1`, and the corresponding guarded receipts
hold on infinite exponent phases. Grounding selected family members and
paying arbitrary guarded obligations are deliberately kept separate.

## 6. Exact controls and the remaining target

[Program](../../04-computation/experiments/collatz_grounded_child_closure_20261007.py)
and [saved output](collatz_grounded_child_closure_20261007.out).

```text
python -B 04-computation/experiments/collatz_grounded_child_closure_20261007.py
python -B -O 04-computation/experiments/collatz_grounded_child_closure_20261007.py
```

The declared universe includes the supplied seeds 3, 5, 13; parameters
0 through 80; all contracting words of length at most three over letters
1,2,3 (plus the empty word); complete ternary residue universes through
precision five; all 39 nonempty generator sequences of length at most
three for each seed; and the two large supplied points. Independent
literal forward replay, the old completed-family recognizer, direct
modular exponentiation in 1512 cases, and exact difference valuations
test the symbolic compiler from different representations. The 306
symbolic refinement edges include huge parameters without expansion.
Malformed types, forged metadata, ROOT padding, source replacement,
expanding inverses outside scope, and allocation limits are hostile
controls. Checks remain active under `-O`; no unvalidated public object
is cached by numeric equality.

Reseeding controls compare eleven literal zero-phase members, all 25 pairs
of old/new parameters 0 through 4 in the nonzero-phase hostile, and a
symbolic large-child parameter `10^100` through binary and ternary residue
readers. Both empty identity extensions preserve their source exactly.

Normal and optimized runs agree on **322,863 exact checks**. The saved
normalized-LF output has SHA256
`ed4a05dc71769ee0650952905d04a1e72732fc5d1196f19eccdcfce11ad26e4e`.
Independent proof/code review and optimized replay passed.

The next global obligation is not another formal closure identity. It is
to put an **arbitrary supplied uncovered source** into a grounded family,
or give it a paid dependency in a well-founded union of such families.
The isometry solves every finite guard in parameter space; it does not
solve that source-membership problem or the unresolved ROOT obligations
in the arbitrary-exponent paid-exit family.
