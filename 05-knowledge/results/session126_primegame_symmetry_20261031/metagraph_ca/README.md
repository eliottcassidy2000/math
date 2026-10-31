# The six-cell perimeter CA and the ten-node observer

## Results and scope

This is an exact finite-carrier experiment and proof, **not** an identification
of the metagraph with FRACTRAN, a Platonic solid, or a divisor lattice.

* Independently reconstructing all 64 fixed-path tournaments gives 12 directed
  isomorphism classes and 10 converse-merged nodes, exactly matching the bundled
  hash-pinned inherited atlas.
* On the geometrically ordered six-cell perimeter, rule 150 has four fixed
  points, all reached in at most two steps. The literal masks requested by the
  user satisfy `7 -> 26 -> 26`, `17 -> 63 -> 63`, and
  `55 -> 22 -> 37 -> 37`. They select the three nonzero fixed points in this
  stated bit gauge. Each fixed-point basin has 16 masks, so this is not a unique
  selection principle.
* Neither rule 90 nor rule 150 descends to a deterministic rule on the existing
  ten merged nodes. The pure-blue fibre `{11,46,53}` is already a counterexample.
* The minimal exact observer has 38 states for complement/reflection controls,
  29 for rule 90 alone, and 30 for rule 150 alone. With rule 150, complement,
  and reflection available, it has exactly 40 states: the reflection orbits.
* Allowing all six **coordinate-labelled** flips forces all 64 states. At
  `n=5`, even the current node plus its six one-flip observations is injective.
  The general proof is about this control/observer interface, not a general
  computational-complexity lower bound.

`verify_ring_ca.py` is self-contained, standard-library only, and does not read
an inherited atlas. It constructs tournaments, all CA transitions, and minimal
observer quotients. `independent_audit.py` instead uses pairwise arc-set
isomorphism, cyclic integer shifts, and pair-removal indistinguishability.
Normal and `python -O` runs agree byte for byte. The independent audit passes
**67,571** gates. All gates are explicit exceptions, not removable assertions.

## 1. Coordinate conventions — no hidden gauge

Vertices are `1,...,5`, with Hamiltonian path `1 -> 2 -> 3 -> 4 -> 5`.
The free cells, in explorer bit order, are

```
k:       0      1      2      3      4      5
cell:  (5,1)  (4,1)  (3,1)  (5,2)  (4,2)  (5,3)
```

Mask zero is the ascending transitive tournament. A set bit reverses the
corresponding nonconsecutive arc. Thus the integer mask is `sum x_k 2^k`;
`explorer_bits_lsb_first` in the TSV lists `x_0,...,x_5`, not ordinary printed
most-significant-bit-first binary.

The triangular boundary has perimeter order

```
[0, 1, 2, 4, 5, 3].
```

The original reflection `R` fixes bits 0 and 4 and exchanges 1 with 3 and 2 with
5. In perimeter coordinates `y_i = x_[0,1,2,4,5,3][i]`, it is precisely
`(Ry)_i = y_{-i}` with indices modulo 6. The complement `C` is `x xor 63`.

This chooses the **six-cycle as the CA neighbourhood graph**. It does not assert
that cyclic shifts of these cells induce automorphisms of the tournament
observer or the colored metagraph. Reflection has the inherited geometric
meaning; full ring-translation symmetry is an additional CA modelling choice.
In particular this is perimeter adjacency, not all adjacencies in the filled
triangular array.

The merged observer `pi(x)` numbers fibres by their least mask, starting at 0.
The individual-tiling blue flag is `b(x) = [Rx=x]`. A **pure-blue node** is a
whole `pi` fibre on which this flag is always true. They are different notions.
The three pure-blue fibres are

```
node 0: {0}
node 7: {11,46,53}
node 8: {16}.
```

The requested masks `7,17,55` lie respectively in nodes `1,9,6`. Only mask 17
is a blue tiling; its node is mixed, not pure-blue.

## 2. Why rule 150 is a structural candidate

Work over `F_2` and let `T` cyclically shift the six perimeter bits. Put

```
A = T + T^-1             (rule 90),
B = I + T + T^-1 = I+A  (rule 150).
```

These equations also define the rules without relying on their standard
numbers. The complete class of homogeneous reflection-symmetric radius-one
linear binary rules is `{0,I,A,B}`: the left and right coefficients must agree.
Of these, exactly `I` and `B` commute with bit complement, since a linear map
`F` satisfies `F(x+1)=F(x)+1` exactly when `F(1)=1`. Requiring genuine neighbour
interaction excludes `I`. Thus **B is the unique interacting rule in this
explicit class that preserves both reflection and complement**. This is a
reason to test B before inspecting particular masks, not a fit to 7/17/55.

Directly, `RA=AR`, `RB=BR`, `BC=CB`, but `A(Cx)=Ax`; rule 90 erases the global
complement difference rather than respecting it.

These are finite periodic-boundary rules. Infinite-line XOR/Fermat-product
identities are a separate model and are not being transferred to this ring.

## 3. Complete algebra and dynamics

Because `T^6=I` and the field has characteristic 2,

```
A^3 = T^3 + T + T^-1 + T^-3 = A,
B^3 = (I+A)^3 = I+A+A^2+A^3 = I+A^2 = B^2.
```

The cancellation of `T^3` with `T^-3` is specific to this six-cycle.

The equation `Ax=0` says `x_{i-1}=x_{i+1}` in perimeter coordinates, so its
kernel is the two-dimensional space of 2-periodic patterns. Hence `rank A=4`.
On `im A`, `A^2=I`; there is at most one transient step. The fixed-point
equation `Ax=x` is `Bx=0`, whose recurrence `x_{i+1}=x_i+x_{i-1}` has a
2-dimensional space of 3-periodic solutions. Consequently A has 4 fixed
points and 6 two-cycles on its 16-element image:

```
fixed:    0, 29, 46, 51
2-cycles: (5,24), (10,36), (15,60), (18,33), (23,57), (43,54).
```

Exactly 16 masks start on these cycles; the remaining 48 enter in one step.

For B the decisive formula is

```
(B^2 y)_i = y_i + y_{i-2} + y_{i+2}.
```

Thus B² is the projection obtained by XORing the three even positions and,
separately, the three odd positions. In explorer bits put

```
e = x_0 xor x_2 xor x_5,
o = x_1 xor x_3 xor x_4.
```

Then `B²x = 37e xor 26o`. The image is `{0,26,37,63}`, and each fibre has
16 masks. This image is exactly `ker A = Fix B`, so every trajectory becomes
fixed after at most two steps. There are 4 states at depth 0, 12 at depth 1,
and 48 at depth 2. `rank B=4`, while `rank B²=2`.

There is an explicit block decomposition, not just a finite-state count.
In the explorer-mask basis

```
[37, 26, 46, 10, 29, 5],
```

B sends these basis vectors to `[37,26,0,46,0,29]`. Therefore

```
B ~ I_2 direct_sum J_2(0) direct_sum J_2(0),
A ~ 0_2 direct_sum J_2(1) direct_sum J_2(1).
```

Equivalently `P=B²=I+A²` and `Q=A²` are complementary projections onto `ker A`
and `im A`. This explains the two persistent bits and the four transient bits.

### The requested triple, literally

| mask | B trajectory | even XOR e | odd XOR o | final merged node |
|---:|---|---:|---:|---:|
| 7 | 7 -> 26 -> 26 | 0 | 1 | 5 |
| 17 | 17 -> 63 -> 63 | 1 | 1 | 5 |
| 55 | 55 -> 22 -> 37 -> 37 | 1 | 0 | 6 |

This is an exact three-attractor statement in the stated gauge. It is **not**
a three-pure-blue-node statement: the three nonzero attractors occupy only two
merged nodes. All four attractors are blue *tilings*, but only fixed mask 0
lies in a pure-blue node. The choice of “nonzero” also uses the stated zero
mask; that three-point set is not complement-invariant because C exchanges
0 with 63 and 26 with 37.

## 4. Why neither CA is a rule on ten nodes

For a deterministic map `f` to descend through `pi`, equality of current
observations must imply equality of next observations:

```
pi(x)=pi(y)  =>  pi(fx)=pi(fy).
```

The same pure-blue node already violates this for each CA:

| source mask | pi | A image | pi(A image) | B image | pi(B image) |
|---:|---:|---:|---:|---:|---:|
| 11 | 7 | 46 | 7 | 37 | 6 |
| 46 | 7 | 46 | 7 | 0 | 0 |
| 53 | 7 | 36 | 5 | 17 | 9 |

Choosing one representative per node can produce a ten-state animation, but
it adds an arbitrary section and does not model the genuine carrier dynamics.
A faithful system keeps enough hidden state, quantified next.

## 5. Exact continuation quotients

For controls `F` and an observation `o`, define

```
x ~ y iff o(w(x)) = o(w(y)) for every finite control word w.
```

This is the coarsest observation-preserving deterministic quotient: equality
must survive every future control word; conversely this relation is itself
stable under every control. The producer computes it by partition refinement;
the independent audit computes the inequivalent pairs by iterative removal.
The full classes and transition functions are in `observer_quotients.json`.

| allowed controls | observe pi | observe (pi,blue flag) |
|---|---:|---:|
| C,R | 38 | 38 |
| A | 29 | 29 |
| B | 30 | 36 |
| A,B | 37 | 37 |
| A,C,R | 38 | 38 |
| B,C,R | 40 | 40 |
| A,B,C,R | 40 | 40 |
| each of the six coordinate-labelled flips | 64 | 64 |

These are exact `n=5` counts, not asserted formulas at other orders.

### Complement and reflection: an explicit 38-state memory

The identities `C²=R²=I`, `CR=RC`, and `pi R=pi` imply that the complete
continuation signature is just

```
q(x) = (pi(x), pi(Cx)).
```

C swaps its two entries and R leaves them unchanged. The first entry is the
current observation and the second is required to predict the next complement,
so this quotient is minimal, not merely sufficient. In this particular atlas
it also determines the blue flag; appending the flag does not refine it.

There are 16 singleton blocks, 20 doubletons, and two four-element blocks:

```
q=(3,3): {3,9,54,60},
q=(6,6): {19,25,38,44}.
```

Each four-element block consists of two reflection orbits. These are the only
additional identifications beyond the 40 reflection orbits. A preserves this
38-state quotient, as verified on every block.

### Adding rule 150: exactly the 40 reflection orbits

Reflection commutes with B and C and is invisible to pi, so reflection orbits
always give a sufficient 40-state carrier. The two possible coarser blocks
above are both separated by a single B observation:

```
3,9   -> under B -> node 5,    54,60 -> under B -> node 9,
19,25 -> under B -> node 2,    38,44 -> under B -> node 1.
```

Hence the minimal B,C,R quotient is precisely the reflection-orbit quotient.
There are 16 reflection-fixed masks and 24 nontrivial pairs, giving
`(64+16)/2 = 40` states. Appending A changes nothing.

### Each rule alone

Since `A³=A`, the complete A-only observation signature is
`(pi(x), pi(Ax), pi(A²x))`; subsequent observations repeat with period 2.
There are exactly 29 signatures.

Since `B³=B²`, the complete B-only signature is
`(pi(x), pi(Bx), pi(B²x))`; subsequent observations are constant. There are
exactly 30 signatures. If the observed output includes the blue flag, apply
that observation at all three times; the number becomes 36.

## 6. Why coordinate-labelled controls recover the full carrier

For any `n`, let the `m=(n-1)(n-2)/2` free bits specify a tournament containing
the fixed Hamiltonian path. Its converse-merged transitive fibre is the
singleton zero mask: a transitive tournament has a unique Hamiltonian path,
so every free arc must agree with that path order.

Suppose a quotient preserves pi and each individual labelled bit flip.
For distinct masks `x,y`, apply the word of flips whose XOR mask is x. It sends
x to 0 and y to `x xor y != 0`. The resulting observations differ because
only 0 lies in the transitive fibre. Thus x and y cannot be identified. Every
such quotient has all `2^m` states.

At `n=5` the exhaustive result is even stronger: already

```
(pi(x), pi(x xor 1), pi(x xor 2), pi(x xor 4),
 pi(x xor 8), pi(x xor 16), pi(x xor 32))
```

is injective. Every pair sharing its current node can be distinguished by one
labelled bit intervention. This stronger one-step fact is verified here only
for `n=5`.

The hypothesis matters. This is **not** a bound for complement-only control,
unlabelled/random flips, approximate observation, or arbitrary computation.
It is an exact observability statement for a specified control alphabet.

## 7. Portable reproduction and evidence

This folder is a self-contained publication packet. It uses only Python's
standard library; no network, installed packages, scratch siblings, or source
worktree are required. All paths are relative to the packet or a generated
temporary directory.

From this directory:

```sh
sha256sum -c SHA256SUMS
python verify_packet.py
```

`verify_packet.py` checks the complete file inventory and hashes, then copies
only this packet's required files into unrelated temporary directories. It
runs the producer and independent audit there in normal and optimized (`-O`)
modes, with Python isolated mode (`-I`) and no inherited `PYTHONPATH`. Every
producer artifact and both captured outputs must match the published bytes.
`reproducibility_report.json` records the verified result.

For a direct, in-place replay of the two mathematical checkers:

```sh
python verify_ring_ca.py > /tmp/ring-ca-producer-output.json
python independent_audit.py > /tmp/ring-ca-audit-output.json
cmp producer_output.json /tmp/ring-ca-producer-output.json
cmp audit_output.json /tmp/ring-ca-audit-output.json
```

Files:

* `ring_ca_results.json`: compact complete reconstructed census, cycles,
  basins, explicit linear decomposition, and summary checks.
* `transitions.tsv`: every one of the 64 masks, both bit orderings, observer,
  blue flag, C/R/A/B transitions, and A²/B².
* `observer_quotients.json`: compact complete minimal blocks and control
  transitions for both output conventions and every listed control set.
* `verify_ring_ca.py`, `independent_audit.py`: independent, optimized-mode-safe
  producer and auditor; the audit performs 67,571 explicit verification gates.
* `producer_output.json`, `audit_output.json`: deterministic replay outputs,
  including all-pair separation counts and minimal distinguishing depths.
* `reference_n5_atlas.json`, `source_pins.json`: one byte-exact inherited atlas
  and its repository snapshot/path/hash. The two research-lane atlas pins were
  identical; only one copy is needed. The independent census needs no atlas
  input, but the auditor verifies agreement and provenance.
* `verify_packet.py`, `reproducibility_report.json`, `SHA256SUMS`: portable
  relocated normal/optimized replay, receipt, and complete integrity manifest.
* `PEER_REVIEW.md`: the separate research-lane mathematical review. Publication
  changes only made outputs compact and the provenance dependency local; both
  checkers were rerun after packaging.

The proof references only files bundled here. The source repository path in
`source_pins.json` is provenance, not an executable dependency or relative link.
The packet's results are original finite derivations and exhaustive checks,
not claims inferred from outside citations.

No isomorphism to the FRACTRAN prime registers, F/U/S balance strata, divisor
complement, or Platonic actions is inferred merely from the small cardinalities.
The verified bridge is specifically **six free bits -> geometric perimeter
CA -> exact observer-memory requirements**.
These finite rules, with eventual periods at most two, cannot themselves be
an unbounded prime-generating FRACTRAN simulation. A separate unbounded-tape
CA construction is a different carrier and must be validated separately.
