# What 63 changes in the Chan Pak proof strategy

**PROVED scoped / FINITE-EXACT; Chan–Pak remains OPEN.** A numerical factor
of three need not be a removable piece of a tournament. The useful
replacement question is which same-count change creates the marked
capability needed for the next exact target.

The main new operation changes exactly three arcs while preserving H,
the named root, its sole predecessor and its entire star. At 63 it changes
the rooted state from (63,11,11,63) to (63,24,25,64), making the root
unit-ready. Every neutral move changing at most two arcs is trapped in a
small cube of module swaps; three is the sharp support threshold at this
order. The operation extends to an infinite family, and then to arbitrary
paired modules with an explicit square defect.

## The proved constructions and their boundaries

| Result | Exact statement | What remains outside it |
|---|---|---|
| [Retained-root reflection](history/APEX_REARM.md) | Reverse-complement an apex word; every word starting 10 and ending 00 changes defect 0→1 at the same H and same root star | Not every target has such a word |
| [Constant three-arc family](history/APEX_REARM.md) | For B=2^b, b≥1, H=3B²+4B−1 has an exact neutral rearm | The NR exit has feedback rank two and leaves the one-apex class |
| [Arbitrary paired-module mirror](MIRROR_MODULE.md) | Z=12p²+8hp−h², defect 0→h²; NR gives Y=(6p+2h)²−5h² | The module must physically realize (h,p); primitive Y≡23 modulo 24 |
| [Primitive two-slot bridge](primitive/PROOF.md) | H=13p+11h; NNR gives balanced (117p+29h,39p+9h), including 63→165→497 | Seven-adic residue coverage is not exact integer coverage |
| [Repeatable module growth](primitive/PROOF.md) | Append one vertex at a retained module cut: (b,t)→(2b−29h,2t−9h) | The module position is essential; this is not a universal rule on scalar b |
| [Incoming fork completion](INCOMING.md) | W(1001^k00) gains the exit 45·2^k+11, including 63→191 and 123→371=7·53 | This branch changes old arcs; no common renewed source is supplied |
| [Finite activation certificate](semigroup/PROOF.md) | For base 7 and multipliers (3,5,49), I7=(x²,y,z), with seeds 63,35,343 | Localization can hide infinitely many missing primitive targets |

Here h=H(M), and p=P_M(2) counts ordered two-slot path covers with empty
slots allowed. The exact primitive criterion for the norm output is
gcd(h,p)=1 and 3∤h. Thus nontransitive modules can produce primitive
balanced continuations even though their first repair defect is h²>1.

Exact growth of the 63 fibre is inherited: there are two strong cores,
both order eight, and every other H63 tournament is transitive padding
around one of them. The number of ordinary isomorphism classes at order n
is 2(n−7) for n≥8; see [primitive §1](primitive/PROOF.md) for the source.
This describes that fibre, not the full count spectrum.

## Why the outside connections were useful

Older graceful fixed-color fibres and Gold exchange matchings supplied a
specific question: which partition is preserved by every small move, and
what is the least larger move that crosses it? At 63, this produced the
sharp three-arc repair. Older planar finite-word conjugations suggested
how to realize converse by a local change retaining the required interface.
The exact operation transfers are in [HISTORY](history/HISTORY.md).

The arbitrary-module extension then produced a genuine quadratic norm.
Its congruence limits and finite inverse show why the physical response
pair (h,p) must remain attached to the arithmetic equation. The supplied
quantum and metagraph briefs make a related observation-versus-object
distinction; their use and the two square-sum size clocks are recorded in
[SOURCES](SOURCES.md).

The semigroup connection gives a positive finite-box compiler and a sharp
hostile to the hoped-for global inference. One can delete infinitely many
primitive allowed values while preserving multiplication, all odd squares
and cubes, every allowed multiple of three, and all principal
localizations. These are abstract countermodels, not tournament claims.
Thus even such algebraic information cannot supply missing prime targets
without an additional physical operation.

## The revised proof obligation

The inherited reduction [THM-4094](../../../01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md)
gives the efficient global target: realize every odd prime other than
seven and every 7p with prime p≥5, p≠7; the carry values 49,63,343 are
known. That reduction and the two 63 cores are inherited results.

The new work suggests several source families with physical transitions,
instead of another fixed-seed loop. Each family carries its root state,
module response and order bound. A target-directed inverse must select an
actual module or source, change family when needed, and decrease a stated
target height. The norm chart has a finite inverse description; the linear
chart has literal repeatable growth; the companion reflection adds a
complementary residue class. No theorem yet supplies the transitions and
source selection for every prime lane.

The decisive unresolved question is whether these local capabilities can
be supplied at every prescribed primitive target with a terminating
physical construction. The specific resources and their limits are in
[BOARD](BOARD.md) and [INCOMING](INCOMING.md).

## Reproduction and review

Run the commands in [AUDIT](AUDIT.md). They use exact integer subset DP,
direct permutation controls, two independent symbolic compilers and
literal adjacency rows. All-parameter results have combinatorial proofs;
finite observations are identified separately. Independent reviews covered
the reflection, feedback-rank law, two-slot wrapper, activation theorem,
countermodels and arbitrary-module norm proof.

No new theorem ID was reserved and no broad tournament census was run.
The output is a set of exact operations and a narrower remaining proof
obligation; no complete Chan–Pak proof is claimed.
