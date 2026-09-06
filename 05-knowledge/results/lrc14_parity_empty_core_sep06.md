# Mixed parity breaks the universal odd triple ceiling

**Status: REFUTED extension + FINITE-EXACT bounded parity census.** The
`6/77` ceiling in
[THM-4434 — universal scale-three network projection bound](../../01-canon/theorems/THM-4434-lrc14-universal-scale-three-network-projection-bound.md)
cannot be extended to arbitrary primitive ternary-unit triples. Both the
selected network bound and the physical Haar bound fail. The two replacement
parity ceilings suggested below are **OPEN** at unbounded height.

## 1. Typed predicate and inheritance

Let `w=(a,b,c)` satisfy `1<=a<b<c`, `gcd(a,b,c)=1` and `3∤abc`, without
assuming oddness. Use exactly THM-4434's complete carrier set `Lambda(w)`,
raw projections `E_i(w)` and physical failure comb `F_w`. The residue,
strict-roof, Helly and interval-to-carrier identities are parity free; the
proposed extension is `min_i E_i<=6/77`, implying `mu(F_w)<=6/77`.

The inherited mechanism is
[THM-4032 — d3 affine defect boundary](../../01-canon/theorems/THM-4032-lrc14-d3-affine-defect-lattice-boundary.md),
Sections 3–4. Its distinct positive ternary-unit scope already includes
mixed parity, and its exact sidecars retain a common affine lift. It gives
an existence criterion for spoiled phases, not a uniform Haar ceiling.
The canonical hostile is the tail triple `(1,10,11)` in
[THM-4004 — three-detuned divisor combs](../../01-canon/theorems/THM-4004-lrc14-three-detuned-divisor-comb-profile.md),
Section 3. The corrected near miss is “phase-uniform failure implies full
failure”: that same typed row is safe at `x=4/33` although its cited divided
pack phase `y=1/11` is spoiled. The least-used relevant sidecar is the
small-coefficient exclusion in
[THM-4425 — all-height rank-one carrier closure](../../01-canon/theorems/THM-4425-lrc14-all-height-rank-one-carrier-closure.md),
Section 3, which explicitly removes `(1,1,1)` and `(1,2,2)` by oddness.

The live concept board for this probe is: owner residues; relation parity;
complete ray addresses; projection versus physical mass; body-safe phase.
All counterexamples retain the first four coordinates. None encodes the
last, so none is an LRC counterexample.

## 2. Exact first hostile and two parity peaks

The first failure ordered by `(c,a,b)` is

```text
w=(2,5,7),              Lambda={±(1,1,-1)},
E=(22/245,6/49,1/10),   mu(F_w)=min E=22/245,
22/245-6/77=32/2695>0.
```

These values are a short exact certificate: the two carriers each contribute
`(11/245,3/49,1/20)` to the projections. Thus physical overlap itself exceeds
the old ceiling; changing the projection selector cannot repair this row.

The complete height-79 universe has the following unique maxima for both
`min_i E_i` and `mu(F_w)`:

| Number of even coordinates | Number of triples | Maximum of both targets | Unique maximizing triple | Rows with `min E>6/77` |
| --- | ---: | --- | --- | ---: |
| 0 | 2910 | `6/77` | `(1,5,11)` | 0 |
| 1 | 9044 | `6/55` | `(1,10,11)` | 243 |
| 2 | 8694 | `11/140` | `(2,11,20)` | 1 |

For the mixed peaks the complete certificates are

```text
w=(1,10,11), v=(1,1,-1), Lambda={±v,±2v},
E=(6/55,12/77,23/154), mu(F_w)=6/55,
6/55-6/77=12/385;

w=(2,11,20), v=(1,-2,1), Lambda={±v,±2v},
E=(131/1540,11/140,3/35), mu(F_w)=11/140,
11/140-6/77=1/1540.
```

The first peak is precisely THM-4004's inherited selector hostile, now with
its exact complete network and physical masses attached. The second shows
that simply excluding the additive relation `a+b=c` does not restore the
odd ceiling.

## 3. First failed implication and strongest observed survivors

For odd `w`, `v.w=0` forces `sum_i v_i` and `||v||_1` even. Removing oddness
destroys this implication. The primitive ternary-unit direction `(1,1,-1)`
of norm three becomes possible; its speed relation is `a+b=c`. Every one
of the 243 one-even failures in the bounded universe has this single live
direction. These are complete ray supports, not selected carrier subsets.

The two-even failure has norm four, so its mechanism is different: the
odd-only finite head in the signed norm-four closure never included
`(2,11,20)`. The parity-free analytic tail of that proof does not validate
the omitted finite rows. Both the coefficient-pattern universe and the
finite speed universe must be changed before an extension is asserted.

The newly permitted norm-five pattern `(1,2,2)` is not a target obstruction
in this probe. Among rows with one live direction of norm five, the largest
selected projection is

```text
w=(2,19,20), Lambda={±(1,2,-2),±(2,4,-4)},
E=(48/665,11/140,46/665), min E=46/665<6/77,
mu(F_w)=173/2660.
```

The strongest bounded survivors are the sharp ceilings `6/55` for exactly
one even coordinate and `11/140` for exactly two. They are **FINITE-EXACT
through height 79**, and only **OPEN candidates** at all heights. No large
scan or all-height extension is part of this result.

## 4. Connection map and smallest next tests

**Proved continuation:** the [additive-family theorem](lrc14_additive_parity_empty_core_sep06.md)
now closes `a+b=c` at all heights with sharp bound `6/55`, uniquely at
`(1,10,11)`. The remaining mixed-parity circuit universe is separate.

* **Source → target:** THM-4004's three-detuned divisor branch → a
  body-specific scale-three Haar completion for mixed tails. The map divides
  the ten divisible speeds by three and retains the three exceptional
  speeds. It preserves the common sheet label and exact tail failure comb;
  discarding the body loses which phases belong to `G_H`. The needed
  sidecar remains `mu(G_H)` or a direct intersection argument with `G_H`.
  A mixed-tail ceiling, even if proved, would not supply that sidecar.
* **Smallest analytic test:** treat the full additive ray `a+b=c` using its
  exact address set `3∤k`, `|k|<min_i 3(w_j+w_k)/(14|v_i|)`, and prove or
  refute the candidate `6/55` bound without height sampling. This isolates
  every one-even failure observed here.
* **Independent correction test:** retain all parities in the existing
  signed norm-four finite head and test the candidate `11/140` boundary.
  The exact hostile `(2,11,20)` must be an equality control. This prevents
  an additive-ray repair from silently reusing the old odd-only base case.

## 5. Reproducible exact evidence

[Producer](../../04-computation/lrc14_parity_empty_core_sep06.py) and
[frozen output](lrc14_parity_empty_core_sep06.out). Reproduce from the root:

```sh
python3 -B 04-computation/lrc14_parity_empty_core_sep06.py
python3 -B -O 04-computation/lrc14_parity_empty_core_sep06.py
```

The universe is every sorted distinct primitive positive ternary-unit
triple with largest speed at most 79: 20,648 rows, including the 2910 odd
positive controls. The primary path literally constructs all shifted danger
intervals on a common integer denominator `42abc`, forms all six owner
assignments, and computes contact capacities and physical intersections for
each omitted owner. The independent path enumerates the full modular
integer-relation set under strict roofs and sums its rational formulas.
Every projection and every physical mass agree on every row; all three
native physical projections agree as well. The proof gates use explicit
exceptions and remain active under optimization. There are 388,026 gates;
normal and `python -O` outputs are byte-identical.

The semantic digest of all rows, all projections, physical masses and
complete carrier sets is
`2379c15b6e7371a35e6e86e57cd33c4c7fad427e96c1ceaefb93d2f7203c6663`.
The small witnesses prove the refutation independently of the census;
the census supplies only the stated bounded maxima and mechanism counts.
An independent referee separately enumerated full integer boxes for the
three hostile controls `(2,5,7)`, `(1,10,11)` and `(2,11,20)` and reproduced
their complete supports, all projections and physical masses. This bounded
referee audit passed; it did not rerun the full census or claim an
unbounded parity ceiling.

Raw SHA-256 hashes:

```text
source b99af309ff6f8643dedf923f5ee8d86d67b32ff2b0d6510209c565820894f399
output 5cd5d081aae910724948c5cd7341aeb09a6245b5fab44f265ac8a4020b63382e
```
