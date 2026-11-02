# Independent adversarial audit: all-modulus width-three controller

**Verdict: ACCEPT.** The candidate `../PROOF.md` establishes the stated all-modulus realization theorem, with strong witnesses and arbitrarily large distinct Hamiltonian-path counts. I found no counterexample or missing mixed-prime compatibility condition. The conclusion is local/profinite universality, not integer spectrum completeness. No literature-priority claim is assessed here.

Reviewed final proof SHA-256: `b7da98667ef7547dd8814d522363cde9c7209f39283bfe4bce7a41cb434d4bf2`. The final version incorporates the quantitative growth justification below and explicitly declares the positive-width convention.

This review independently reconstructed the gate algebra and proof directions, rather than accepting the existing computation transcript. A separate scalar-controller implementation is included as a diagnostic; it does not pretend to replay actual graphs.

## 1. Carrier and gates

The polynomial lift has the right normalization for unordered path covers. In the OCF expansion, contracting a family of vertex-disjoint odd cycles produces `n-sum(|C|-1)` atoms. Partitioning these atoms into c unordered path blocks gives the Stirling number, and evaluating the falling-factorial sum produces the asserted power of x.

For transitive factors `x^r,x^s`, a cross-block cycle in `C3[L_r,L_s,1]` selects nonempty increasing subsets from both transitive blocks. Its length is odd exactly when the total selected size is even. The sum of the two finite-difference products selects precisely these terms, with coefficient two. The existing bilinear substitution kernel therefore extends this monomial computation to every actual factor profile. No assumption that arbitrary residue pairs are actual profiles is needed.

At x=1, the padded B context has first coordinate

```
(3u-v)h+(v-u)g.
```

At x=2, the multiplier after padding is `2^N=1` in each odd residue ring, and the remaining bracket is

```
(u-v)h+(3v-u-w)g+(w-v)f.
```

This independently confirms the important cancellation: on the n-part, the auxiliaries have `w=v`; on the 3-primary part, f is zero by reset. Uncontrolled higher coordinates really do disappear. Every relevant padded cyclic gate restores both reset sidecars. The V factor has enough padding because `N-3>=max(a,b)+1`.

The auxiliary evaluations and all four gate rows check exactly. In particular the n-part matrices are

```
D1=diag(2,1), A=[[5,-1],[1,0]],
A4=[[11,-3],[3,-2]], Vgate=diag(1,1/4).
```

Thus

```
A D1^2 A^(-1) Vgate=[[1,15/4],[0,1]],
D1^(-2) A=[[5/4,-1/4],[1,0]],
D1^(-3) A4=[[11/8,-3/8],[3,-2]].
```

The last determinant is `-13/8`, so treating E as invertible at 13 would be wrong. The candidate does not do so.

The same products have dyadic and 3-primary responses exactly as stated. Raising U to `r=2^a 3^b` annihilates all primary anchor actions: `16,-2` belong to `1+3Z`, and `-27` is odd. Therefore the shear step fixes both primary anchors, including the entire 3-primary pair. Moreover

```
gcd(tau,n)=gcd(5,n),   tau=(15/4)r mod n,
```

so its 5-adic valuation is exactly one whenever 5 divides n. No extra obstruction is hidden at powers of 5.

## 2. Unit-bank lemma and actual seeds

The unit-bank lemma is valid uniformly, including mixed dyadic/3-primary moduli. Its finite base can be made completely explicit. Under CRT modulo `8*9=72`, the following positive products generate the independent factors:

```
17          -> (1,8),
5^3 * 17    -> (5,1),
23^3 * 17   -> (7,1),
5^2         -> (1,7).
```

The middle two generate all units modulo 8; the first and last generate all units modulo 9. Hence the bank generates all 24 units modulo 72.

For the uniform lift, any proper subgroup of a finite abelian group of order `2^i 3^j` is annihilated by some nontrivial character to C2 or C3. Every such character of `(Z/(2^A 3^B))^*` factors through the corresponding reduction modulo 72 (or a divisor at small exponents). For the dyadic component this follows from the decomposition into -1 and powers of 5: the kernel of reduction modulo 8 consists of squares and therefore is killed by every C2 character. For the 3-primary component, the kernel of reduction modulo 9 consists of cubes and is killed by every C3 character. The cross-primary character vanishes on the component of coprime group order. The base calculation rules out the character. This proves generation at every height; it is not an extrapolation from finite enumeration.

The four seed pairs have the claimed values. Dividing the 3-primary target by `3^v` and independently multiplying by the inverse of `3^v` on the dyadic side gives a unit target modulo `2^a 3^(b-v)`. The unit lemma therefore supplies nonnegative exponents. Adding the finite order of 5 makes k positive without disturbing either anchor, even when the modulus is 1 (take order one).

Crucially, actual seed H has possible odd prime factors `{3,5,17,23}`, while actual seed g has possible odd prime factors `{3,7,37}`. Since n is coprime to 3, their joint gcd with n is one. Forced k>0 gives the deliberate axis preparation `h=0,g!=0 mod5`; it is not valid to replace this by an arbitrary nonzero pair.

## 3. CRT synchronization and singular-prime safety

The first shear independently chooses a residue exponent at each prime power of n0. Where p divides g, the old h is a unit and exponent zero suffices. Else tau and g are units, so h can be set to one. CRT combines these into one nonnegative exponent, and S preserves both primary anchors. Applying R then makes g a unit modulo n0. Modulo 5 the same operations turn the prepared seed into `(t,0)` with t nonzero.

Modulo 5, R is the swap and E is `[[2,4],[3,3]]`. I checked the full ten-state orbit and all 20 table endpoints directly. Each target first coordinate has a nonzero second coordinate. The fixed diagonal `(1,1)` is genuinely outside the usable axis argument; the proof appropriately avoids claiming global transitivity.

Before each E, the inserted shear sets h to zero modulo n0. E then sends g to `-2g`, a unit. In particular, at prime 13 this avoids the nonzero kernel vector `(2,3)`. There is no inverse of E anywhere in the construction. Before each R, the shear sets h=g modulo n0, so R preserves the unit diagonal there. Each inserted shear is invisible modulo 5, so it does not corrupt the chosen finite word. Consequently g remains a unit at every non-5 prime and becomes a unit at 5 at the endpoint.

The last equation has coefficient gcd exactly 5 in the c>0 case; the finite word has already enforced the sole necessary divisibility condition. It therefore lands at the desired full residue modulo n, not just modulo its radical. In the c=0 case the coefficient is a unit and ordinary inversion suffices. All preceding anchor congruences remain unchanged.

Edge conventions cause no gap. If n0=1, its preconditioning steps are vacuous and may simply be omitted. If n=1, there is no remaining odd-part constraint. If M=1, all anchor constraints are vacuous but the positive final block and strongness argument still apply.

## 4. Positive inverses are actual constructors

The primitive gates D1, A, and Vgate induce bijections on the full finite displayed controlled state because their relevant scalar factors are units and A's n-matrix has determinant one. The reset class is invariant and future gate responses depend only on this controlled state. Therefore a positive power of the same genuine gate realizes its inverse on every actual input in the reset class.

The bound `K!-1`, with `K=2^a(3^b n)^2`, is a valid common inverse exponent. It can be zero when K=1, which is simply the identity context, not a negative constructor. No recovery of higher polynomial coordinates or graph isomorphism is required. Every continuation uses only the controlled responses already proved independent of those lost coordinates.

All fixed factors, source padding, join gates, and repeated B contexts remain in the given K1/J/B grammar. The argument constructs a finite actual expression. None of its modular inverses signifies subtraction of graphs, a negative number of vertices, or an unrealized free state.

## 5. Strongness and a quantitative growth tightening

The candidate's final stripping argument is correct only at the terminal position, and it is used exactly there. Adding n to the last S exponent preserves the n-pair and the already-fixed primary anchors while ensuring at least one final S block. Ordinary composition order makes its final primitive gate A. Its output is exactly

```
P^(N-1) B(T,U2).
```

Deleting precisely these outer source vertices preserves exact integer H. The remaining graph `C3[T,U2,1]` is strong because all three blocks are nonempty and the quotient is a directed triangle. Internal padding is not removed.

There is an especially clean justification of distinct unbounded counts. Write `p=p_T(2)` and `q=p_U(2)`. Since `g=2(h+p)` and `v=2(u+q)`, the B formula becomes

```
H(B(T,U))=3hu+2up+2hq+4pq >=3hu >=3h.
```

Thus every genuine cyclic context multiplies its variable-input H by at least three. Joins and source padding cannot decrease H. Each U contains at least two D1 applications and its final A, even if its positive-power inverse is the identity. Therefore each U multiplies H by at least 27, and S by at least `27^r`. Increasing the final exponent by successive positive multiples of n consequently gives strictly increasing actual counts. Terminal stripping preserves those counts. This independently supplies the algebraic justification for the quantitative growth bound now present in the final proof.

## 6. Independent executable diagnostic

`crt_referee.py` was written independently of `../local_universality.py`. It checks only the scalar CRT controller and seed/unit arithmetic, not graph realization; the symbolic argument above supplies the positive-constructor and strongness audit.

Completed checks:

* Every parity-compatible target for every modulus 1 through 256.
* 64 targets each for three additional mixed moduli, including high 5-powers and simultaneous factors 7,13,17,23,37.
* **24,832 exact simultaneous scalar landings**, with all original dyadic and 3-primary seed anchors verified.
* **11,394 safely preconditioned R/E steps**, checking the maintained unit sidecar and invisibility of shears modulo 5.
* The full ten-state mod-5 orbit and all 20 table endpoints.
* The fixed-diagonal hostile modulo 5 and the singular-E hostile `(2,3)` modulo 13.
* Exact generation of the full unit group for each of the 36 unit moduli encountered.

All gates use explicit `require`/raise checks, not removable Python assertions. Normal and optimized outputs match byte-for-byte. The transcript prints the script SHA-256. These finite diagnostics supplement, and do not replace, the all-modulus argument.

No repository mutation or external state changes were made. This audit authorizes no publication action; the parent retains that responsibility.

Independent script SHA-256: `a64965acb5e60ae4e943d737d8f155471d47e9482837c8b5436518e4959ff86f`. Both independent output files SHA-256: `910b9ec30ff8fd527b98e4dd9ee874c364545725cfdcc392f373edec4b7f83e9`.
