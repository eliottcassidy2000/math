# Mod-19 observers that retain the source and select certified inverse addresses

2026-10-04. **PROVED** elementary clock, affine-congruence, digit-lift and
certificate-selection statements below. **FINITE-EXACT** for the declared
script universes. No historical priority claim and no global Collatz
convergence or arbitrary-source coverage claim.

Artifacts: [script](../../04-computation/experiments/mod19_recursive_observers_20261004.py)
and [exact output](mod19_recursive_observers_20261004.out).

## 1. Inheritance and the useful question

The [sixth-clock note, section 3](sixth_clock_branches_20261004.md)
already proves `ord_(19^k)(2)=18*19^(k-1)` and its anchored synchronization
with a golden phase cycle at denominator `4*19^k`. The
[difference-family note](difference_families_20261004.md#beyond-76-there-is-an-exact-tower-with-19-branches)
proves that the golden two-coordinate system has 19 **child cycles** per
parent. An equality of periods is not a ring map or an integer itinerary.

The inherited [affine-axis note, sections 2–3](geometry_collatz_drift_carriers_20261004.md)
retains the ordered carry, exact dyadic source guard and real source/axis
comparison. The [Paley observer note, section 3](paley_geometry_observations_20261004.md)
proves that all quadratic orientation probes together see only the two
parities of word length and exponent sum. The
[cycle-gate note](collatz_mod6_20260921_pillai_convergents_cycle_gates.md)
retains the stronger ordered condition `Q-P | B` for an integer fixed
anchor; a small power gap alone does not imply that divisibility.

The [inverse-ray codec](inverse_ray_ternary_addresses_20261004.md)
already represents exactly the positive odd basin of 1 by unique guarded
first-hit codes `(parent, row, block)`, evaluates binary/ternary residues
without expanding the represented integer, and selects ternary addresses.
Its grammar's coverage of all positive odd integers remains open.

Closest mechanism: affine carry plus exact source guards. Canonical hostile:
formal modular iteration can continue after the integer word has become
illegal. Corrected near miss: the scalar inverse ray's 19 digits are not
the golden phase system's 19 child cycles. Least-used sidecar: the shared
block-index residue modulo 3 between ternary and 19-adic addresses.
The concept board is **source / word / carry / phase / address / certificate**.

## 2. Full affine observers preserve more than Paley signs, but still lose drift

For a nonempty positive valuation word `w=(a_1,...,a_j)`, retain

    P=3^j, Q=2^A, A=sum a_i,
    B=sum_(i=0)^(j-1) 3^(j-1-i)*2^(a_1+...+a_i),
    F_w(n)=(Pn+B)/Q.

Since 2 and 3 are units modulo `19^k`, this defines an affine permutation
with multiplier `P/Q` and translation `B/Q`. Composition retains both
coordinates: `(r_v,t_v) o (r_u,t_u)=(r_v*r_u,r_v*t_u+t_v)`.
Thus it can distinguish word order and carry that every Paley orientation
probe discards. For example words `(1,2)` and `(2,1)` have carries 5 and 7
and different affine permutations modulo 19, despite identical clocks.

The formal fixed-point equation is

    (Q-P)n = B mod19^k.

**PROVED exact gate.** Put `g=gcd(Q-P,19^k)`. There are no fixed residues
unless `g|B`; if it divides, there are exactly g, given by one residue
modulo `19^k/g`. Divide the congruence by g and invert its remaining unit
coefficient to prove both directions.

Words `(4,4)` and `(1,7)` have the same `P=9,Q=256`, hence the same
multiplier, real drift and projective discriminant `(P-Q)^2`. Their carries
are 19 and 5. The first affine reduction is the identity modulo 19; the
second is a nonidentity translation. At every precision `k>=1` the first
has exactly 19 fixed residues and the second has none. The exact rational
axes are respectively `1/13` and `5/247`: cancellation of 19 in the first
axis denominator is essential. This is a local refinement of the inherited
integer cycle gate, not a sufficient test for an integer cycle.

**PROVED all-precision hostile.** At each fixed precision k, words

    (1),       (1+18*19^(k-1))

give the same complete affine permutation modulo `19^k`. Their real
multipliers are respectively greater than and less than 1. This follows
immediately from the inherited order of 2; both carries are 1. Consequently
even the full finite affine observer cannot determine real drift in an
unbounded word universe. Extra precision can separate a specified pair;
no fixed precision solves the problem uniformly.

The two exact dyadic word cylinders are disjoint: one integer cannot have
two different first valuations. Each cylinder nevertheless projects onto
**every** residue modulo `19^k`, by coprime CRT. Thus any prescribed observed
input residue has distinct positive legal representatives for these two
words with identical affine observations. The first representative grows;
the second descends, since its axis `1/(2^(1+18*19^(k-1))-3)` is below 1.
This is a loss of information in the observation, not competing legal
choices at the same actual source.

## 3. The 112 chart has a complete formal shell decomposition

The marked word `(1,1,2)` gives

    F(n)=(27n+19)/16,         a=-19/11,
    y=11n+19,                y(F(n))=(27/16)*y(n).

The coordinate change is invertible modulo every power of 19. Write
`rho=27/16`. Exact arithmetic gives

    rho=10=2^(-1) mod19,
    rho^18=1+4*19 mod19^2.

The first identity gives order 18 modulo 19. Repeatedly taking nineteenth
powers of `1+19^s*c`, with `19` not dividing c, raises its 19-adic valuation
by exactly one; taking a power prime to 19 preserves it. Therefore

    ord_(19^k)(rho)=18*19^(k-1)                 for every k>=1.

**PROVED shell theorem.** Modulo `19^k`, the affine chart has its one fixed
point `n=-19/11`, and for each `0<=h<k` exactly one orbit on

    v_19(11n+19)=h,

of length `18*19^(k-h-1)`. Dividing y by `19^h` reduces to multiplication
by rho on all units modulo `19^(k-h)`. Its order equals their cardinality,
so its powers visit every one. These orbit sizes plus the fixed point sum
to `19^k`.

At the next precision the fixed point's 19 lifts consist of the new fixed
point and one 18-cycle. Each existing shell cycle has exactly **one** child
cycle, whose length is multiplied by 19. This is a one-coordinate phase
tree. The golden two-coordinate tower instead has 19 child cycles: its
parent phase has `19^2` lifts, not 19. Keeping dimension separates the two
recursions despite their matching primitive periods.

The source remains indispensable. The exact positive integer guard for
one 112 block is `n=7 mod32`. The source 7 legally follows
`7 ->11 ->17 ->13`; the next 112 block is illegal at 13. Its formal image
would be `185/8`, although every modular iterate still exists. The API
`chart112_observer(n,k)` returns coordinate, shell and formal period,
with this one-word source guard as a separate field. Repeated-block
legality must use the inherited higher-precision dyadic guard.

There is nevertheless an exact, useful source filter for any word. For
a requested endpoint residue e modulo `19^k`, intersect

    n=(Q-B)*P^(-1) mod2Q,
    n=(Qe-B)*P^(-1) mod19^k.                        (1)

The coprime Chinese remainder theorem gives one residue modulo
`2Q*19^k`. Its positive members have the exact supplied valuation word
and the requested endpoint address. Conversely every such legal source
lies in this cylinder. `endpoint_source_filter` implements this iff.
It can reject a proposed certificate endpoint's incompatible residue;
a matching residue still requires the endpoint's actual certificate.

## 4. Inverse channels admit recursive 19-adic address selection

Fix a positive odd certified hub u with `3` not dividing u. In source row
`r in {0,1,2}`, let `kappa` be the unique member of `{1,...,6}` satisfying
`2^kappa*u=1+3r mod9`. The inherited exact inverse channel is

    n_b=(2^(kappa+6b)*u-1)/3,        b>=0,
    U(n_b)=u,                      n_(b+1)=64n_b+21.

Let `v=v_19(u)`. If `k<=v`, every n_b has the same address `-1/3 mod19^k`.
If `k>v`, its exact block-index period is

    L_k=3*19^(k-v-1),

and every one of those `L_k` addresses occurs exactly once per block period.

**Proof.** The order of 64 modulo 19 is 3. For a nonzero block difference t,

    v_19(n_(b+t)-n_b) = v                         if 3 does not divide t,
                       v+1+v_19(t)               if 3 divides t.

Indeed factor the difference as
`2^kappa*64^b*u*(64^t-1)/3`, and use
`64^3=2^18=1+3*19 mod19^2` with the same elementary lifting argument.
This proves both the period and injectivity within one period. For negative
t use the unit `64^t` in the valuation identity, or swap the two indices.

**PROVED digit lift.** For `k>v` and `d in {0,...,18}`,

    n_(b+d*L_k) = n_b + d*19^k*((3n_b+1)/19^v) mod19^(k+1).   (2)

The normalized coefficient is a unit modulo 19. Expand
`64^(d*L_k)=1+3*d*19^(k-v) mod19^(k-v+1)` and substitute into the
difference formula. Thus the next desired source digit uniquely selects
d by one unit division. There are only three base addresses, tested at
precision `v+1`, followed by one carry computation per higher digit.

The API `select_ray19(hub_mod9,hub_mod19,row,wanted,k)` implements this
complete address test from finite registers, returning `(b0,L_k)` or
`None`. It does not require an expanded hub. When the hub register is zero
modulo the requested precision, it correctly uses the collapsed period-one
case. A false unit assumption would incorrectly divide by zero and invent
additional addresses.

For each fixed row, this channel visits one sixth of the nonzero shell
`v_19(3n+1)=v` when resolved beyond v; it does not traverse the full shell
visited by the 112 multiplier. That distinction is already visible in
orders `ord(64)=3` and `ord(rho)=18` modulo 19.

## 5. The ternary and 19-adic selectors share a block digit

The inherited ternary address selector gives

    b=b3 mod M,        M=3^(a-1)

for each desired source address modulo `3^a` in row r. Section 4 gives
`b=b19 mod N`, with `N=3*19^(k-v-1)` if `k>v`, or `N=1` in the
collapsed case. These are constraints on the **same block integer**.

**PROVED joint selector.** Both source addresses occur together iff

    b3=b19 mod gcd(M,N).                            (3)

When compatible, the complete block set is one progression modulo
`lcm(M,N)`. If `a>=2` and `k>v`, the gcd is 3: exactly one third of all
pairs of separately reachable addresses are jointly reachable. If `a=1`
or `k<=v`, that compatibility obstruction disappears. This is the ordinary
noncoprime Chinese remainder theorem, applied after the two exact address
bijections; the source moduli being coprime does not make the two clocks
independent.

**Small hostile.** In the root's row r=2, `n_0=5` and `n_1=341`. Address
`n=5 mod9` exists, and address `n=18 mod19` exists, but they cannot occur
together in that row. They require `b=0 mod3` and `b=1 mod3` respectively.

`select_certified_predecessors` composes (3) with the inherited root-pointer
codec. Every selected source has an actual edge to the supplied certified
parent. It preserves the unique first-hit certificate rather than merely
its address. There is exactly one root exception: parent 1, row 1, block 0.
If a progression contains that block, only its initial parameter is excluded;
the rest of that progression remains valid. For example the joint root
address `1 mod9`, `1 mod19` gives `b=0 mod3`; its first allowed member is
`b=3`, source `349525`, with its exact one-step odd return to 1.

To evaluate the certificate modulo `19^k`, use the node's exact exponent e
and update `(2^e*u-1)*3^(-1) mod19^k`. No 19-adic precision is consumed,
because division by 3 is invertible here. The kappa calculation still needs
the parent's mod9 register, and ternary source queries still consume the
inherited extra parent digit. These are different precision budgets.

The unexpanded demonstration constructs a certified parent of odd depth 6
with five block heights `10^100+i`. It selects a child with address
`5 mod27`, `84354 mod19^4`, block `63 mod61731`, and odd first-hit depth 7.
Its rigorous bit-length lower bound exceeds `10^100`; no represented integer
is expanded. The inherited structural proof establishes the route, while
the new observer chooses and verifies its requested addresses.

## 6. Transfer contract and finite audit

| Source -> target | Preserved predicate | Information lost / needed sidecar |
|---|---|---|
| Exact word -> affine map mod19^k | Endpoint congruences and composition | Real drift, full carry and integer domain; retain the word and source guard |
| Legal chart -> cylinder (1) | Exact word plus requested endpoint address, iff | Address does not identify a certified endpoint; retain endpoint witness |
| 112 source -> shell phase | Formal return period and shell | Repeated integer legality; retain dyadic precision and original source |
| Certified hub + joint address -> predecessor progression | Actual first-hit edge and both addresses, iff | Root self-return must be removed; retain parent pointer and block height |
| Base-2 clock -> anchored golden phase | Cyclic order and phase | Different dimension, orbit branching, ring operations and integer realization |

Reproduce:

```text
python -X utf8 -B 04-computation/experiments/mod19_recursive_observers_20261004.py
python -O -X utf8 -B 04-computation/experiments/mod19_recursive_observers_20261004.py
```

The script checks both primitive orders through `19^6`, excluding every
prime shortening; every formal residue modulo 19, 361 and 6859 in the
112 shell census; all 340 words of lengths 1–4 over `{1,2,3,4}`, with all
19 endpoint addresses at depth 1 and addresses `{0,1,18,19,360}` at depth 2,
three positive sources per filter (24,480 literal valuation replays).
The 340 words yield 181 distinct affine maps modulo 19.

The inverse address census uses hubs `{1,5,7,11,19,361,6859}`, all three
rows, precisions 1–3, and **every** requested residue at those precisions:
152,019 exact iff tests against independently enumerated block orbits,
plus 5,130 next-digit controls. All 3,591 pairs of addresses modulo 27 and
19 at those hubs are tested against literal sources; 189 are compatible,
and their guarded instantiated certificates are replayed to 1. The script
also checks the equal-clock fixed-point hostile at all residues through
`19^3`, the all-precision drift alias at four precisions, the isolated root
exclusion, malformed inputs, and the unexpanded depth-7 demonstration.
All checks survive optimization. No sampled density or universal source
coverage is inferred from these finite universes.

The useful outcome is a source-aware filter and an exact inverse-family
selector. Finding a legal, descending or completed certificate for an
arbitrary supplied source remains a separate obligation.
