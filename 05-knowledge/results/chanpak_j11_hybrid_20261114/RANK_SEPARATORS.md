# Exact rank separators in the eleven-block Chan–Pak orbit

The source is wave 20 and `scripts/runterm.c` of
`05-knowledge/results/chanpak_waves_20261113`. The result below does not
extend its finite coverage bound: it gives an exact factorization test,
valid at arbitrarily large height, which distinguishes a prime branch of
the search from its composite branch.

Write the terminal state as (H,E,L). The eleven letters mean J followed
by e, O, OO, OZ, OOO, OOZ, OZO, OZZ, OOOZ, OZOZ or OZZZ; e means empty.
The first run uses the original source and is not part of the eleven-letter
recurrence. The auxiliary fourth coordinate is set by J to the old H.
Define each three-by-three terminal matrix by the following instructions:
J gives (H-E+L,H-E,H+L,h=H); O replaces (H,E,L) by (H,0,2L);
Z replaces it by (H+L-E,h,2L-h). This fixes both product order and coordinates.

## Induced rank detector and its arithmetic meaning

Take a terminal state after O, hence E=0, L even, H odd. In the valid orbit L>=1 and H>=1; E<H and E<L hold after any
leading O. The initial J of any next block strictly raises H on the
E=0 plane, and subsequent H never decreases. More generally all
terminal states lie in the nonnegative cone (H-E,E,L-E), as checked
by the eleven terminal matrices below. The first coordinate cannot
decrease inside a block: J and Z add L-E, O does not change H;
after J, L-E equals old L+E, and the leading O then clears E.
Each Z chain after that O has L>=h=old H, and preserves L>=h
by L→2L-h. This proves the claim for actual states and the formal
primitive seeds used below.

Follow any corridor of invertible letters (e and the six Z-ending
blocks). Their determinants are 1 or 2^length. Track the three 2-by-2
minors (a,b,c) of the linear map on the starting E=0 plane, in row
order (H,L),(H,E),(E,L), starting at (1,0,0). Each step is its second
compound matrix in that row-pair order. Ending with one of the four
O-ending blocks returns to E=0. The determinant of the resulting
two-by-two map to (H,L) is

* 2^k(b-c) for a pure O^k block (k=1,2,3);
* 6a+6b-8c for OZO.

Indeed the pure O^k rows are (1,-1,1) and (2^k,0,2^k); OZO rows
are (3,-1,3) and (6,0,8). Their compound rows give the formulas.
The pure empty-corridor case is the simplest separator.

Any such determinant zero is genuine rank one, not rank zero. In
the basis (H-E,E,L-E) all eleven matrices are nonnegative. Thus the
return matrix on (H,L), with E=0 at both ends, is nonnegative integral;
both columns are nonzero (track the coordinate vectors under the first
J and the subsequent updates, which leave L>0 and H>0 from that point).
Writing its columns as multiples of the same primitive integer vector
(h*,l*) gives literal output (uH+vL)(h*,0,l*) with u,v positive integers.
This is not an inference from a projective cone: it is exact in the
integer orbit. The factor uH+vL is at least 3 on a valid incoming
O-ending state. Any suffix preserves this factor by linearity. On the
primitive ray, (h*,0,l*) has h*>=1 and l*>=1; the subsequent J increases
h* strictly, and H never decreases on subsequent blocks in this cone.
Consequently:

* A primitive first coordinate h*>1 excludes prime targets at and
  after the separator.
* If h*=1 the separator may output a prime at that very block but
  excludes prime targets after a nonempty suffix.
* For a fixed separator and a prescribed target N, enumerate only
  divisors d of N with d=uH+vL. The residual target N/d must occur
  in the suffix from the primitive restart seed. These conditions,
  plus actual prefix reachability, are exactly sufficient, not
  just congruence tests.

This gives a multiplicatively weighted restart decomposition: cut
a word at successive rank-one returns; each new seed is the primitive
ray of the return and each edge weight is the integer scalar computed
by the previous seed and intervening corridor. There is no claim the
full set of rays is finite, nor that the residual prime language is
cofinite.

## New literal separator examples

Start with E=0. For the following returns, list the actual induced
(H,L) matrix as outer product (h*,l*) (u,v); read the words left to right.

| Corridor followed by its O-ending block | Primitive seed | Linear factor |
|---|---|---|
| O | (1,2) | H+L |
| e e OOZ e e OZO | (55,138) | 7H+9L |
| e OZZZ OOZ e OZO | (65,162) | 51H+52L |
| OOZ OZOZ e e OZO | (71,178) | 41H+45L |
| OZ e OZZ OOOZ O | (51,104) | 65H+81L |

The first example is the local reset; the others exhibit genuine
nonlocal resonances. All are certified with exact determinant
and literal path-state replays, not orbit coverage counts.

The structural takeaway is that the singular letters cannot simply
be removed (the wave-20 census already says they are load-bearing),
but they can be separated into two exact mechanisms. A singular
letter usually remembers a second coordinate; an actual rank loss
forces a divisor, potentially usable for composites but fatal to
a prime proof in the displayed situations.

## A finite suffix reduction for the prime subproblem

A useful exact normalization does not need the complete set of
rank-one returns. Let P be any actual initial run, then any sequence
of the eleven blocks avoiding the following certified patterns:

1. An O-ending block followed by a pure O^k block, k=1,2,3.
2. An O-ending block followed by any of the four nonlocal returns
   in the table.

All occurrences of pattern 2 are forbidden in a prime word,
including at the final block. An occurrence of pattern 1 is allowed
only as the final step of a prime word. Its pure exponent k does not
affect the target there: the final H is precisely the old H+L. Thus
the prime projection of the full eleven-block language is **exactly**
the primes in the union of:
* H(P);
* H(P)+L(P), for such a prefix P ending in O (including a legal
  initial run ending in O).

In the first bullet P obeys both avoidance rules. In the second
bullet append any single pure block to certify the target; if P
has length zero beyond the initial run, the same argument applies.
This gives a finite word-avoidance automaton with two exact absolute
readouts, not necessarily a finite numeric carry automaton. Adding
any further verified strict-interior rank separator only prunes
this language, without losing a prime. Readout union is essential:
throwing away the terminal pure reset entirely would wrongly remove
its H+L candidates.

The four nonlocal examples also exclude **all 7p with p prime**:
their scalar exceeds 7 (use H>=1, L>=2), and so does the primitive
first coordinate, even before further blocks. If qp is a product
of two integers greater than q, both would have to be divisible
by p, which is impossible if p>q (and if p<=q their product is
already greater than qp). The same argument applies to any fixed
q<min(h*,u+2v). For a local pure separator with H+L>7, any 7p
witness must place it within its final six blocks: the primitive
first coordinate starts at 1 and increases by at least one per
following block, so a suffix of length at least seven gives both
factors >7. The complementary case H+L<=7 has finitely many
incoming E-zero states. These are exact restrictions on the
prime and 7p targets, not claims that these targets are attained.

For this particular eleven-bank the seven-step bound sharpens to
**zero** for every 7p with p>7. There are two reasons.

1. The primitive seed (1,0,2^k) of a pure separator is itself
   physically reachable: use k initial O letters on the singleton.
   Any lawful suffix is a tournament construction from this seed.
   In particular it cannot yield the known forbidden count 7.
2. The incoming scalar H+L itself cannot be 7 on a valid E=0 state.
   In fact L>=H always in the underlying letter grammar. J gives
   L'-H'=E; O doubles L and leaves H unchanged. After any O, the
   first Z of a consecutive run gives
   L'-H'=L-H-h>=0 using L>=2H and h<=H; subsequent Z's preserve this
   difference because their incoming E=h. An initial run of Z's
   merely repeats the seed (1,1,1), if allowed. Consequently E=0
   implies L>=2H (the exception to an immediately preceding O
   is a J out of the singleton state, which also gives (1,0,2)).
   So H+L=7 would force H=1 and L=6. While H stays 1, the only
   possibilities are the singleton stutters and repeated O doublings
   from (1,0,2), reached by O or J. L=6 never occurs.

A proper-factor decomposition of 7p with p>7 must use the factors
7 and p. Thus a pure separator at any interior position excludes
7p exactly as it excludes primes. The two-readout reduction in the
previous section works verbatim for primes **and** 7p with p>7.
This is a hybrid use of the universal impossibility of 7 to prune
the orbit, rather than a claim that the pruned language attains
all the admissible prime and 7p targets.

## Exact combination with ordinal multiplication

The background follow-up strengthens the pure-reset rule beyond
the prime lanes. Let S be the strict eleven-block scalar language
and let R omit a pure block JO^k, k=1,2,3, whenever its entering
state has d=0 (immediately after an end-O block, or in an initial
source phase with d=0). Then the multiplicative closures of S
and R are exactly equal. Indeed, at a forbidden reset the return
factor H+W occurs already by placing J alone after that same
prefix. If the reset is terminal, use that replacement. Otherwise
its normalized output (1,2^k,0) is the actual singleton followed
by k letters O, an initial run within the allowed budget. The
suffix on that normalized seed yields a second physically admitted
odd factor. Both factors exceed one: d=0 and physical parity give
H+W>=3; a first block on (1,2^k,0) increases H to at least 3.
Induct on the positive count for a nonterminal split, and on the
number of forbidden occurrences for a terminal replacement.
Ordinal joins then realize the same scalar products without forcing
the original continuation states to agree.

For nonzero d a pure block outputs (H+W−d, 2^k(H+W),0), whose
coordinate content is gcd(H+W,d), since H+W−d is odd. It is a
valid prefix divisor filter, not a factorization into realized
tournament counts. As a separate state-preserving reduction at
d=0, the block word J,J,JO and the single block JOZO both give
(3(H+W),6H+8W,0) in (H,W,d) order. This substitution works under
all subsequent complete block suffixes. The multiplicative closure
statement does not imply individual composite membership in R.
