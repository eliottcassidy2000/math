# Meta-Flashlight S95: Lifted Symmetric Sectors

**Session:** codex-2026-05-29-S95
**Script:** `04-computation/meta_flashlight_s95.py`
**Stored output:** `05-knowledge/results/meta_flashlight_s95.out`

## Constraint Discovered

The bare Cartan decomposition of a binary tournament adjacency matrix is inert.
For every 0/1 tournament adjacency `A`,

```text
S = (A + A^T)/2 = (J - I)/2
K = (A - A^T)/2
||S||_F^2 = ||K||_F^2 = C(n,2)/2
```

So the current Cartan/dark-mode analogy must not be applied directly to raw
tournament adjacency. If there is a nontrivial symmetric "dark sector", it must
live in a weighted object or a functorial lift:

- weighted attention/logit matrices,
- the OCF conflict graph `Omega(T)`,
- transfer matrices `M(T)`,
- heat kernels/resolvents,
- or another symmetric construction derived from the tournament.

This is a useful negative result: it prevents the Cartan thread from becoming a
pure metaphor.

## Probe Result

The OCF lift `T -> Omega(T)` is the smallest concrete variable symmetric sector
already in the repo. Exhaustive check through `n <= 6`:

```text
n=3: OCF mismatches=0
n=4: OCF mismatches=0
n=5: OCF mismatches=0, odd gaps=[7]
n=6: OCF mismatches=0, odd gaps=[7,21,35,39]
```

Conditioning on `t3` removes the first centered skew-chaos term, but it does not
determine `H` once `n >= 5`.

At `n=5`, the only ambiguous fiber is:

```text
t3=4: H range [11,15], |Omega|=[5,6,7], alpha2=[0]
```

At `n=6`, ambiguity becomes structural:

```text
t3=2: H range [5,9], alpha2=[0,1]
t3=4: H range [11,17], alpha2=[0,1]
t3=5: H range [15,27], alpha2=[0,1]
t3=6: H range [23,37], alpha2=[0,1,2]
t3=7: H range [33,37], alpha2=[0,1,2]
t3=8: H range [41,45], alpha2=[1,2,4]
```

The flashlight: the residual variation inside fixed `t3` fibers appears exactly
where the lifted symmetric sector `Omega(T)` gains independent sets.

## Three Cross-Thread Hypotheses

### 1. Lifted Dark Sector Hypothesis

The useful symmetric side of tournament theory is not `S=(A+A^T)/2`; it is the
first nonconstant symmetric functor applied to `A`.

In the pure tournament thread, that functor is `Omega(T)`. In the transformer
thread, it should be the symmetric part of a weighted attention/logit/hidden-state
matrix before binarization. In the adelic/ghost thread, it should be the local
heat/resolvent data after the conductor has split into prime places.

Same deeper pattern:

```text
antisymmetric competition -> symmetric lift -> scalar invariant
T                         -> Omega(T)      -> I(Omega,2)=H
attention logits           -> Cartan S      -> confidence/correctness
eigenvalue conductor       -> local factors -> global spectral invariant
```

This reframes Cartan not as "every tournament has dark modes", but as:
every useful tournament invariant is computed after a symmetry-producing lift.

### 2. Residual-Fiber Principle

The repo has many first-order laws: `t3`, score sequence, Paley spectrum,
quadratic Walsh term, conductor, or logit gap. Novel structure tends to appear
in the residual fibers after one of these laws is fixed.

Examples already present:

- `t3` fixed but `H` varies through `Omega` independent sets.
- Paley has optimal spectral gap but loses to Interval through higher Walsh terms.
- raw Cartan of binary adjacency is fixed, but lifted symmetric sectors vary.
- permanent gaps `{7,21}` survive after many local constructions fill the rest.

Search method:

1. Pick a first-order invariant with a beautiful theorem.
2. Condition on it.
3. Enumerate the fiber.
4. Ask which lifted symmetric object first separates the fiber.

This is a meta-hypothesis about where new hypotheses should be mined.

### 3. Conductor-Activation Principle

Several unrelated threads seem to be activation thresholds for new coordinates:

- k-periodicity: new correction layers activate at multiples of the tournament
  atom `3`.
- OCF: `alpha2` activates when two disjoint 3-cycles fit.
- Adelic space: new prime factors of `odd_part(C(n,2))` add local coordinates.
- Cartan bridge: symmetric information is invisible until the object is lifted
  beyond raw binary adjacency.

Possible common statement:

> New mathematics appears when the current coordinate system can no longer
> separate residual fibers; the system responds by adding either a new cycle
> layer, a new prime place, or a new symmetric lift.

This makes "where have we not looked?" concrete: look for variables whose
fibers have not been separated by the next natural lift.

## Immediate Tests

- For `n=7`, sample fixed-`t3` fibers and regress `H` against `alpha2`,
  `|Omega|`, Omega degree statistics, and transfer-matrix invariants.
- For circulant primes `p=7,11,13,17,19`, condition on quadratic Walsh energy
  or spectral IPR, then test which higher `Omega` independent-set statistic
  first predicts Interval beating Paley.
- For the Cartan probe, avoid binary adjacency. Use weighted top-k logit
  matrices and compare the symmetric sector against correctness labels or
  synthetic self-consistency labels.
- For adelic conductor tests, compare residual H-fibers before and after a new
  odd prime enters `odd_part(C(n,2))`.

## Working Slogan

Do not look for the theorem in the first invariant. Look in the fiber it fails
to separate.
