# The parabolic point is where a continuant becomes a congruence

**Session:** death-star-imw-pi-n-s2, 2026-08-02.
**Trigger:** an owner request to read Ivanov--Mikhailov--Wu, *On nontriviality
of homotopy groups of spheres* (arXiv:1506.00952; Homology Homotopy Appl. 18
(2016) 337--344), and connect it to the repository's frontiers.
**Outcome:** THM-3204, THM-3205, THM-3210,
`05-knowledge/reference/CORE-PAPERS-HOMOTOPY.md`, and one candidate
META-PATTERNS card.  Both `CORE-PAPERS.md` and `META-PATTERNS.md` are at the
shared startup byte budget and could not take a pointer or a card; a
compaction pass on the startup surface should add both.  This reflection records provenance and the
reasoning that produced the connection, not truth; the theorem files are the
truth surface.

## 1. What the paper actually runs on

The paper proves `pi_n(S^2) != 0` for every `n>=2` by filling Curtis' open
residue class `n = 1 mod 8` with odd-primary Toda elements.  Strip the homotopy
theory and the entire arithmetic input is **one determinant**: its Lemma 3
computes a `d^1` differential in the odd-primary lambda algebra, finds the
matrix `tridiag(1,-2,1)` of size `k+1`, and observes that its determinant
`(-1)^(k+1)(k+2)` is invertible mod `p` exactly when `p` does not divide `k+2`.
Everything else in the paper is classical homotopy input plumbed around that
gate.

That is a repository object.  THM-3183 and THM-3186 had just put a Jacobi
continuant at the centre of the factorial exterior transfer, with a
`cancellation wall` they could exhibit only by one numeric hostile.  So the
question was not "does topology transfer" -- it does not -- but "what kind of
arithmetic wall can a continuant produce at all?"

## 2. The answer, and why it is the useful form

Three facts, all elementary once asked (THM-3204):

1. **A unital Jacobi matrix has cyclic cokernel.**  Delete row 1 and column
   `n`; the remainder is triangular with the sub-diagonal on its diagonal, so
   `d_(n-1)=1` and `coker = R/(det)`.  A continuant therefore produces **one**
   divisibility, never a `(Z/p)^2`.
2. **The converse is an obstruction.**  Smith factors are equivalence
   invariants, so `d_(n-1) != 1` forbids continuant form outright.
3. **The vanishing period is `p` iff the step is parabolic.**  For
   `D_n = alpha D_(n-1) - c D_(n-2)` over `F_p` with `c` a unit, the zero set
   is an arithmetic progression whose difference is exactly `p` when
   `alpha^2=4c`, and otherwise `ord(x_+/x_-)` dividing `p^2-1`.  Since
   `gcd(p,p^2-1)=1` the two are exclusive.

Item 3 is the one I did not expect and the one that pays.  It says a gate of
the shape "`p` divides length plus a constant" is not a generic feature of
three-term recurrences: it is the **signature of the degenerate point** of the
Chebyshev family, `alpha=+-2`, where `U_n(+-1)=(+-1)^n(n+1)` and the
trigonometric solution collapses to a linear one.  IMW's gate is a congruence
because their matrix sits exactly there.

And it sits there *identically in `p`*: the diagonal entry is
`a(2,1)=binom(p-2,1)=p-2`, the off-diagonals are `+1`, and
`(p-2)^2-4 = p^2-4p = 0 mod p` for every odd prime.  The Adem-type structure
constants of the lambda algebra put the excess-one continuant on the parabolic
locus by construction.  That is the mechanism, and it was invisible while the
lemma was read as "a tridiagonal determinant happens to be `k+2`."

## 3. What transferred, and what refused to

**Transferred, positively.**  THM-3182 recorded two integral frames for the
same rational system with Smith types `(1,1,p)` and `(1,p,p)` and called the
discrepancy "one elementary modification of the output lattice."  Item 2 makes
that invariant: `d_2=p != 1` means the weighted Gauss--Manin reset is **not
equivalent to any unital Jacobi matrix**, over any coordinate change.  THM-3183
exhibited a Jacobi block for the scalar frame and none for the weighted one;
that was not a presentational accident.

**Transferred, negatively -- with an exact obstruction.**  The obvious next
move is to import item 3 into the factorial lane by choosing parameters that
make the transfer parabolic.  Each single step can be made parabolic, on the
conic `P_i=(i+1)(2i+1)^2v^2-4idv+i=0`.  Two consecutive steps cannot: `P_i` is
affine in `d`, and eliminating `d` gives

```text
(i+1)P_i - i P_(i+1) = -(8i^3+20i^2+12i-1) v^2,
```

whose cofactor is at least `39` for `i>=1`, forcing `v=0`, which the standing
hypothesis excludes.  So the single-congruence mechanism is unavailable there
by an exact integer, not by a failed search.  That also explains, without
apology, why THM-3176/3183's wall polynomials `H=24p+109` and the quartic `J`
are genuine arithmetic polynomials rather than congruence conditions, and why
looking for a `floor(s/2)` staircase law of the form "`p` divides depth plus a
constant" is looking for a signature the transfer does not have.

**Refused.**  There is no map to LRC(14).  I looked: no LRC canon object is a
Jacobi or continuant matrix, and manufacturing one would fail the connection
contract at the first line (no source object, no preserved predicate).  The
temptation was real -- IMW's endgame is a **covering by complements of two
singleton residue classes** (`k != 1 mod p` from route A, `k != 0 mod p` from
route B, and no `k` is both), which rhymes with LRC's residue bookkeeping.  It
is a rhyme.  Recorded as a tangent, not a bridge.

## 4. The unexpected dividend: building the engine paid twice

I built an exact odd-primary lambda-algebra engine to check Lemma 3
independently rather than trust a transcription.  It found two things.

**A convention that is not fixed by the displayed formulas.**  In each relation
the summand `j=k` reproduces the left-hand side with coefficient `-1`.  Read
literally, the paper's own `mu_0 lambda_2 = -mu_1 lambda_1` becomes
`-(1/2) mu_1 lambda_1`.  The self term must be omitted.

**A missing sign.**  With the unsigned Leibniz rule the paper's computations
all reproduce -- and `d^2 != 0`, already on `mu_2`, where
`d^2(mu_2)=4 lambda_1^2 mu_0`.  The Koszul sign `(-1)^|x|` is required;
3640 words fail without it and 0 fail with it.  The paper is not wrong: in
every computation it displays, `d` is applied to the right of an even-degree
`mu` prefix only, so the sign never fires.  Anyone extending those computations
by one step in a different direction would get silently wrong answers.  This is
the clearest case I have seen for the repository's rule that an imported
mechanism must be re-derived, not transcribed.

Once the engine ran, the right stratification became visible and the lemma
strengthened itself.  Grade `Lambda*lambda(2)` by index sum `S` and lambda
count `q`: `d` fixes `S` and raises `q` by one, so the `q=1` stratum receives
**no** differential and its `E^2` is a kernel, not a subquotient.  Excess
`e=S-l` then stratifies the `q=1` part, and the source/target dimensions are

```text
e=0 : 1  ->  0        (Toda elements are cycles)
e=1 : k+1 -> k+1      (SQUARE; the parabolic gate)
e>=2: strictly wide   (injective throughout the tested universe)
```

So Lemma 3 is not the first rung of a tower.  It is the unique square stratum,
and a determinant gate needs a square map.  Within it, `E^2` is exactly `F_p`
when `p | S`, generated by the explicit cycle `z_k = sum (i+1) u_i` -- the
discrete Dirichlet solution of the path Laplacian -- whose leading `mu_2`
coefficient is `k+1 = -1`.  The paper's Proposition 2 *assumes* an `E^2` term
of the form `C mu_2 mu_1^(k-3) lambda_1 + ...` with `C != 0`; the generator
supplies it canonically and pins `C`.

## 5. The dividend flowed back: THM-3210

Having asked "what is the *complete* wall, not one witness?" for the lambda
algebra, I asked it of THM-3186.  Its amplitude has `deg_v E_L = L-2`, so `E_3`
is affine in `v` and its zero set is a single rational graph,

```text
v = (n+3)(d-n-2) / ( 2d[ (4n+9)d - (n+3)(4n+7) ] ).
```

THM-3186's published hostile `n=1, d=5, v=4/105` is the `(1,5)` point of that
graph.  Pushing one step further: on the ray `d=n+4`,
`v=(n+3)/(3(n+4)(2n+5))`, **both** `E_3` and `E_4` vanish identically in `n`,
with `Delta`, `beta_n`, `c_(n+1)=2` and `E_5` all nonzero.  THM-3186's own
witness was already a double cancellation; only `V_3` had been tested.

The consequence is sharper than the original wall.  Visibility runs
visible / invisible / invisible / visible at lengths `2,3,4,5`.  It is **not
monotone in length**, so no "first visible depth" statistic and no search that
certifies visibility a step or two past the graph distance is sound.

## 6. Candidate META-PATTERNS card (not yet promoted)

`META-PATTERNS.md` is at its bounded startup byte budget, so this card is
parked here rather than silently compacting another session's cards.  It meets
the two-distinct-thread evidence bar.

> **Check the discriminant before building a period theory**
>
> **Trigger:** a lane's arithmetic gate is produced by a three-term
> recurrence, transfer matrix, continuant, companion form, or tridiagonal
> differential.
> **Action:** compute the step discriminant first.  On the parabolic locus
> the continuant is `(r+1)(alpha/2)^r` and the gate is one congruence of
> difference exactly `p`; off it the gate has period dividing `p^2-1` and can
> never be that congruence.  Separately, take the Smith form: a unit
> sub-diagonal forces a cyclic cokernel, so the obstruction is a single
> divisibility, and `d_(n-1) != 1` proves no continuant form exists at all.
> **Counterindication:** the closed form needs a *uniform* parabolic step.
> With index-dependent coefficients, check consecutive steps before assuming
> the degeneration is reachable -- eliminating the linear parameter between
> two neighbouring conditions is usually one line and usually fatal.
> **Evidence:** THM-3205 (odd-primary lambda algebra; the gate `p | S` is
> parabolic identically in `p`) and THM-3204 section 4 with THM-3182/3183
> (the factorial exterior transfer is nowhere consecutively parabolic, and
> its `(1,p,p)` reset admits no continuant form at all).

## 7. Concept board at close

| object | live question | cheapest next test |
|---|---|---|
| unital Jacobi cokernel | which repo gates are secretly single congruences? | Smith `d_(n-1)` of any transfer already in canon |
| parabolic locus | non-constant analogue of the period dichotomy | two-step transfer with slowly varying `alpha_i` |
| `q=1` lambda strata | is `e>=2` injective for all `e`, `delta`, `p`? | the dimension gap is suggestive, not a proof |
| `E_L` cancellation | is there a three-length invisible window? | resultant `Res_v(E_3,E_5)` beyond `n=3` |
| `p=2` lambda algebra | the paper's Conjecture 1 | no `mu` generators; the excess count does not apply |
| IMW two-gate covering | is it a method or a coincidence? | find a second lane where two singleton-complement gates cover |

## 8. One reservation race, resolved by moving

`THM-3209` was reserved and pushed by this session at 20:42:49 and by a
concurrent `[gmc3209-reset]` session at 20:42:52.  Three seconds of priority is
not a claim worth defending, and that session is named after the number, so
this session's theorem moved to **THM-3210** -- file, companion script, stored
output, declared hashes, results index, frontier line, and cross-links all
renumbered together, with the companion re-run so the hashes are the renamed
bytes.  Moving is the non-destructive direction: it cannot invalidate the other
session's work.  This is the hazard already recorded as MISTAKE-351/356, so no
new ledger entry is opened; what it adds is the concrete repair recipe --
renumber the *whole* evidence tuple at once and re-hash, never just the file.

## 9. Honest remaining frontier

Nothing here proves a new homotopy group, and nothing here touches LRC(14) or
GMC(2).  THM-3205 section 5 and THM-3210 section 4 are finite-exact in stated
universes, not theorems for all parameters.  The `floor(s/2)` Euclidean-depth
staircase remains open; what changed is that one candidate shape for it -- a
single congruence in the depth -- is now excluded by an exact obstruction
rather than by absence of evidence.
