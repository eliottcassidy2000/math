# Reverse-complement collisions: an explicit family of linear size that makes the row-collision bound sharp

Let p(0) = 0, p(2v + 1) = p(v), p(2v) = p(v) + v + (v mod 2), and consider a carrier row

    G(u) = c u + k p(u) + φ·(u mod 2),        c, k ≥ 1, φ ∈ Z.

[HYP-9752](../../hypotheses/HYP-9752-row-collision-lemma-a-carrier-row-cu-plus-kp-u-has-linear-energy-which-would-give-the-chan-pak-spectrum-positive-lower-density.md)
conjectures that the collisions with top split level l and prefix w number at most C·2^l/(1 + w). This note proves
the matching **lower bound** on a positive-proportion set of prefixes. It does so with one explicit algebraic family,
which also carries a large, computable share of the whole collision energy. The upper bound remains OPEN.

## 1. Two identities for the reverse-complement

For an a-bit window put rev_a(x) for the bit reversal and

    rc_a(x) = 2^a − 1 − rev_a(x)        (reverse, then complement).

**Lemma 1.** For every a ≥ 1 and 0 ≤ x < 2^a:
- (i) p(rc_a x) = p(x);
- (ii) x − rc_a x = rev_a x − rev_a(rc_a x).

*Proof.*
- (i) Unrolling the recursion gives p(x) = Σ_{i<j, b_i=0, b_j=1} 2^(j−i−1) + #{i : b_i = 0, b_(i+1) = 1}, where
  b_i are the bits of x in the window. Bits above the leading 1 of x are zero and contribute nothing, so the formula
  holds in any window containing x. Reverse-complement maps the bit string b to b′_i = 1 − b_(a−1−i).
  - A pair i < j with b′_i = 0 and b′_j = 1 is the pair i′ = a−1−j < j′ = a−1−i with b_(i′) = 0 and b_(j′) = 1.
    It has the same gap, so the first sum is invariant.
  - A pattern b′_i = 0, b′_(i+1) = 1 is the pattern b_(a−2−i) = 0, b_(a−1−i) = 1, so the second term is invariant
    too.
- (ii) Complement and reversal commute, and rev is an involution, so rev_a(rc_a x) = 2^a − 1 − x. Both sides equal
  x + rev_a x − 2^a + 1. ∎

Write δ(x) = x + rev_a(x) − 2^a + 1 for this common difference.

## 2. The collision family

**Theorem (RC collisions).** Let m = 1 + c/k. Take integers a ≥ 2 and j ≥ 1 with k | c·2^(j−1), an odd x with
2^(a−1) ≤ x < 2^a, and r with 0 ≤ r < 2^j. Put

    y = rc_a(x),   t = rev_j(r),   w = m·2^j − 1 − r,
    u = 2^j(2^a w + x) + t,   u′ = 2^j(2^a w + y) + t.

Then G(u) = G(u′), u ≠ u′, and distinct parameters (a, j, x, r) give distinct pairs. The pair has:
- lowest differing bit j;
- highest differing bit j + a − 1, so its top split level is l = a + j;
- prefix w;
- common low suffix t.

*Proof.* We apply the two-ended split of u, u′ at the common suffix t and the common prefix w (the lattice-line
equation (1) of [TWISTED_LATTICES.md](../chanpak_cone_collision_placement_20261115/TWISTED_LATTICES.md), which we
re-derived independently). The pair collides iff

    d w + δ r + [w odd] − ε [r even] = ((c/k) 2^j + 2^j − 1) δ + p(x) − p(y),

where δ = x − y, d = rev_a x − rev_a y and ε = (x mod 2) − (y mod 2). Here the inputs reduce as follows:
- x is odd with top bit 1, so y = rc_a x is even with top bit 0. The pair is a valid two-ended core with ε = 1.
- By Lemma 1, d = δ and p(x) = p(y).

So the condition becomes δ(w + r − m 2^j + 1) + [w odd] − [r even] = 0.

The choice w + r = m 2^j − 1 kills the first term. Since k | c 2^(j−1), the integer m 2^j is even, so w + r is odd
and w, r have opposite parities. Hence [w odd] = [r even], and the second term vanishes too. Distinctness holds
because (j, a, x, w, t) are read off the pair. ∎

- The case a = 1 is the single core (x, y) = (1, 0). Its ray w = 7·2^j − r for the 3-cycle row was found in
  TWISTED_LATTICES.md §3 (credited).
- The new points are:
  - every odd x with top bit 1 is a core, giving 2^(a−2) cores of span a;
  - the union over all spans has **linear** size.
- A single core gives only O(√N) collisions below N, as that note observes. The 2^(a−2) cores of each span make up
  the difference.

## 3. Consequences

**(a) The shape of HYP-9752 is sharp.** Suppose (m − 1)2^j ≤ w ≤ m 2^j − 1, with j as in the theorem and
a = l − j ≥ 2. Then every core of span a contributes exactly one collision at (l, w), with r = m 2^j − 1 − w.
Hence the number of ordered collision pairs at split level l and prefix w satisfies

    N(l, w) ≥ 2·2^(a−2) = 2^(l−j−1) ≥ (c/(2k))·2^l / w.

- In every dyadic block [W, 2W) of prefixes, these w occupy one interval of logarithmic length log₂(m/(m − 1)).
- So the conjectured bound N(l, w) ≤ C 2^l/(1 + w) cannot be improved in shape, and its constant satisfies
  C ≥ c/(2k).
- For the calibration rows this gives:

  | row | lower bound on C | largest count·(1 + w)/2^l measured at L = 20 |
  |---|---:|---:|
  | 3-cycle row 18U + 3p + 2[odd] | 3 | 13.99 |
  | 2U + p | 1 | 5.99 |

**(b) Linear energy from one explicit family.** Suppose k | c, and write m = 1 + c/k.
- *Exact count.* For L-bit inputs, the number of family members below N = 2^L is

      Σ_{j≥1} Σ_{a≥2} 2^(a−2)·clip(2^(L−j−a) − (m − 1)2^j, 0, 2^j).

- *Limit.* This tends to 2^(−1−⌈log₂ m⌉)·N unordered pairs.
- *Comparison with the measured energy.* The replay at L = 18 gives:

  | row | family members / 2^L | limit | share of all collisions |
  |---|---:|---:|---:|
  | 18U + 3p + 2[odd] | 0.0610 | 1/16 | 44.1% |
  | 2U + p | 0.1231 | 1/8 | 19.7% |
  | 12U + 3p | 0.0610 | 1/16 | 17.9% |

- *What the family explains.*
  - The linear growth of per-level totals in the prefix length L − l noted in HYP-9752: each admissible suffix
    length j adds about 2^(l−3).
  - Almost all of the "palindromic" lines found below.

**(c) The arithmetic dichotomy.** The family needs k | c 2^(j−1) for some j, that is odd(k) | c. When the odd
part of k does not divide c, the anti-diagonal mechanism survives only for cores with δ(x) = k/gcd(k, c 2^j)
exactly, at suffix lengths j in residue classes modulo the order of 2.
- *189U + 11p.* The family is empty through L = 18. The whole collision energy is 50–250 times smaller than in the
  two physical rows: 0.0027 against 0.138 and 0.62 per 2^L, unordered.
- *A general divisibility.* Equation (1) shows that every collision with a common suffix of length j ≥ 1 needs
  k / gcd(k, c 2^j) | δ.
- *The physical rows.* All physical carrier rows have k | c (for example 18 = 6·3), so they carry the full family.

## 4. Where the rest of the energy lives (NUMERICAL, L = 18)

Classify every collision by its two-ended core (x, y), as in §2.

| row | palindromic cores (d = δ) | single-point lines | other multi-point lines |
|---|---:|---:|---:|
| 18U + 3p + 2[odd] | 46.2% | 50.8% | 3.0% |
| 2U + p | 40.2% | 58.2% | 1.6% |
| 189U + 11p | 3.5% | 87.5% | 9.1% |

The palindromic cores are characterized exactly by the rigidity lemma: d = δ iff the polynomial Σ σ_i(x^i − x^(a−1−i))
vanishes at 2. This polynomial is antisymmetric under x ↦ 1/x up to the factor x^(a−1). So it also vanishes at
1/2, and it is zero; in other words the signed difference σ = bits(x) − bits(y) is a palindrome. Counting
mirror-pair choices gives exactly 6^(⌊a/2⌋−1) cross cores (times 4 for odd a). Among them, the cores with
p(x) = p(y) are:
- exactly the reverse-complement cores (2^(a−2) of them) for even a;
- the same plus two sporadic cores for odd a.

This is verified for a ≤ 20 (`hd_family.py`, `pal_struct.py`).

**What a proof of HYP-9752 must therefore do.** Subtract the explicit lattice-ray families: the RC family, and
the fixed-core rays of the concurrent note. Then bound the generic single-point lines. Their count is the number
of two-ended cores whose line meets the (w, r) box.
- Counting admissible cores alone gives only O(4^l/W² · log W) at level l and prefix scale W.
- What is missing is the probability, about W²/2^l, that a core's line actually carries a lattice point. That is
  an equidistribution statement for (p(x) − p(y))·d⁻¹ modulo δ, over the digital family of cores.

The row 189U + 11p, which has no family, shows how small the generic part can be.

## Interpretation

The reverse-complement is what the carrier construction does to a tape under the converse. In a carrier
tournament, the X-tape spells digits most-significant-first with mark = digit. The Y-tape spells them
least-significant-first with mark = complemented digit. Reversing every arc swaps these two encodings. Lemma 1(i)
is the digital trace of that converse symmetry (H(T) = H(T^op)). The theorem says that applying it to a *middle
block* of one tape, at the right prefix, preserves the count exactly.

This is offered as a heuristic reading. The proved statement is the algebra above.

## Replay

`python3 verify_rc_collisions.py` (about 10 seconds, numpy) checks:
- Lemma 1 for all a ≤ 16;
- for each of four rows at L = 18: that every family member is a literal collision, that the member count equals
  the closed-form count, and the share of all collisions;
- the sharpness constant c/(2k).

The classification tables come from `coll_classes.py`.
