# A monotone all-future envelope from finite stopped prefixes

Status: **PROVED** conditional envelope and monotonicity for every finite
source carrier; **FINITE-EXACT** native all-horizon certificates at the
declared three scales. Uniform 153N for every N>=4096 and arbitrary-source
Collatz entry remain **OPEN**. No unresolved trajectory is assumed to arrive.

## 1. Inheritance and the changed question

The native map P, 627-point stopped frontier V, Q-norm, source charges and
exact153N annulus budget come from
[session18_native_annulus](../session18_native_annulus/README.md). The
closest cancellation mechanism is the guarded one-edge synchronous
cancellation in [native phase lifts, section4](../session18_native_phase_lifts/README.md).
That result already erases a complete same-state zero packet before ROOT
entry. We generalize the **sufficient bound after a finite prefix**, not the
underlying cancellation identity or the specific64s+35/49 compiler.

The prior [fixed-horizon scale bound](../finite_horizon_budget_20261021a/README.md)
covers N>=2^(8d+12) through d; the present bound covers every future horizon
once a finite prefix passes a different check. These complementary
quantifiers must not be conflated into uniform all-N coverage.

Canonical hostile: N7 fails at native cut8. Corrected near miss: two charges
that sum to zero but occupy different states may arrive at different cuts.
Mixed old/new zero also retains an old-baseline correction. Least-used
sidecar: a **common unresolved current state at the same P-depth**, together
with separate signed old/new vectors and the original integer-clock paths.

Board: immutable source / synchronous current node / separate old and new
channels / unknown continuation / absolute pair envelope / original budget.
The anchor is long-horizon coverage; the niche is a finite-prefix stopping
certificate; the wildcard is a monotone quadratic envelope despite the
nonmonotonicity of the physical energy itself.

## 2. State and exact polarization

At native depth D, every original primitive source s<256N is either stopped
at a terminal in V or has been advanced exactly D primitive edges without
entering V. Sources s<4N carry the old charge A_N(s)=A_(64N)(s); the others
carry the new annular charge. Each has its actual sign (-1)^clock.

For each arrived terminal t, store c_t=(W_t,H_t). For each distinct pending
current node u outside V, store b_u=(w_u,h_u), the sums over its literal
source-disjoint fibre. A common pending node gives the same subsequent
trajectory, future first-arrival depth and extra clock to every fibre member.
At any later cut the entire fibre is absent or present, with one common
additional sign and one terminal. No terminal is assigned before arrival.

For four-coordinate blocks b=(w,h), c=(w',h'), define

\[
\delta(b)=Q(h)+\beta(w,h),\qquad
\kappa(b,c)=\beta(h,h')+\beta(w,h')+\beta(w',h).
\]

Here Q(x,y)=x^2-xy+y^2 and beta is its polarization. Then
delta(b+c)=delta(b)+delta(c)+kappa(b,c), and
delta(-b)=delta(b). The actual present increment is

\[
\Delta_D=\sum_t\delta(c_t)
=E(64N,D)-E(N,D).
\]

The map from source prefixes to blocks preserves exactly the observable
native Nyquist energy. It does not preserve all Fourier frequencies by the
four-vector alone, or all integer-clock cuts. The lossless source carrier
is `(N,D,s)` with the pinned deterministic rule and frontier: `decode`
recovers the complete edge-clock word, integer sum, mod15/mod35, stopped
node and depth. The receipt retains every pending current node, signed
channels, multiplicity, clock range, and a literal decoded witness. Source
rows, performed transitions and zero-group erasures also have reproducible
hashes; those hashes do not replace the checks.

If both w_u and h_u vanish, the block contributes zero at every future
native-depth cut and can be removed from this energy calculation. This
does not assert any of its sources reaches V. Only the numerical Nyquist
specialization is erased; its source address remains decodable.

## 3. All-future envelope

For each pending block put

\[
m_b=\max(\{0\}\cup\{|\kappa(c_t,b)|:t\in V\}),
\quad L_b=\max(0,\delta(b)+m_b),
\]
\[
\mathcal U_D=\Delta_D+\sum_bL_b+\sum_{b<c}|\kappa(b,c)|.       \tag{1}
\]

**Theorem.** For every finite native source carrier and every h>=D,

\[
\Delta_h\le\mathcal U_D.                                \tag{2}
\]

Proof: at h choose the actual subset of pending blocks that has arrived,
their actual terminal and continuation sign. An included b contributes
delta(b) plus a signed interaction with its one current arrived port;
this is at most L_b. An omitted b contributes zero, also at most L_b.
Two included blocks contribute a signed kappa only if they land at the
same terminal; that term is at most |kappa|. Every internal, old/new and
cross-block term occurs once in this polarization. Summing proves(2).
The bound even allows independently chosen signs, terminals and optional
nonarrival; actual dynamics is a subset of this larger finite choice space.

Thus `U_D<=153N`, together with the exact checks at cuts0,...,D, certifies
the original153N inequality at **every finite horizon** for that N. It is
not necessary to finish the unresolved trajectories. Equality in the
envelope need not be jointly realizable; separate absolute-value maxima
can conflict, so this is a sufficient upper bound, not an exact optimizer.

## 4. The envelope is monotone under actual prefix extension

Signed transport of a pending block preserves its delta and every absolute
kappa. Three operations remain.

* **Merge at one common current node.** Replacing a,b by a+b cannot raise
  the envelope. Its current-port maximum is at most m_a+m_b, so
  `L_(a+b)<=L_a+L_b+|kappa(a,b)|` by positive-part subadditivity. Every
  other pair cost decreases by the triangle inequality. The removed
  a,b pair cost pays the possible increase of their merged local term.
* **Arrival at a current terminal.** Move a into c_t with its actual sign.
  The present increment rises by at most L_a. Each remaining L_b can rise
  by at most |kappa(a,b)|, since only one port changes. Exactly those
  pair costs disappear when a leaves the pending list. Hence total U
  cannot increase.
* **Remove a both-channel zero.** Its delta and all kappa terms vanish.

Every synchronous primitive step is a sequence of these operations. Thus

\[
\mathcal U_{D+1}\le\mathcal U_D.                         \tag{3}
\]

This does not prove it ever crosses153N. A monotone sufficient bound can
remain above its target indefinitely. It also does not say Delta_D itself
is monotone. The potential retains adverse inter-block interactions
instead of deleting them because two local blocks are separately safe.

The implementation groups repeated vectors, up to an overall sign, only
to sum(1) faster: each type retains its multiplicity, including the
binomial number of internal pairs. This is scalar arithmetic compression,
not identification of distinct current nodes or their future terminals.

## 5. Finite native certificates with unresolved blocks

The verifier advances only to the listed D and reads no continuation beyond
that cut. It recomputes the source charges, exact source norm, every
performed edge and all earlier physical increments. At the cutoff it
computes(1) and independently decodes one source witness in every pending
block, including its full integer clock and both prime registers.

| N | Cut D | Pending nonzero blocks | Present Delta | Local terms | Pair terms | U_D | 153N |
|---:|---:|---:|---:|---:|---:|---:|---:|
|4096|105|195|421276|74675|61213|557164|626688|
|8192|110|348|793715|153217|191475|1138407|1253376|
|32768|130|760|2838302|651015|752165|4241482|5013504|

The first two scales were already covered by complete-arrival computations;
the new consequence is a certificate valid **without completing those
pending continuations**. N32768 is an additional supplied scale beyond the
older contiguous4096..16384 census. No assertion about all intervening or
larger N follows. The saved production replay verifies all three rows.
The largest earlier increments are respectively421512,798142,2858748,
each below its153N budget. The pending nonzero blocks contain respectively
373,689,1401 original sources; their largest performed integer clocks are
176,185,219. None of their continuations beyond the cutoff is read.

At the first two cutoffs,6991 and14045 complete synchronous zero groups
have been removed. At32768 the number is56767. These are groups formed at
different depths, not distinct new credit deposits. Their zero future
contribution is used once in the state representation; the source budget
stays153N throughout.

## 6. R-prime and ownership

Keep whichever authenticated ownership partition is selected, without
combining overlapping banks. Write

\[
\Delta_h=A_h+I_h,\quad A_h+M_h=153N,\quad I_h=R'_h+J'_h.
\]

The new estimate gives, for every h>=D,

\[
R'_h\le\mathcal U_D-A_h-J'_h.                           \tag{4}
\]

When U_D<=153N this proves the compensated bound R'<=M-J', exactly enough
for physical coverage. It does not claim the stronger sufficient R'<=M.
No previous negative interaction is deducted again. The proof estimates
the full energy directly; replacing a credit partition changes R' and J'
oppositely and leaves(4) valid.

In particular the incoming optimal permanent273/474 allocations in
[cross-scale credits](../cross_scale_frontier_credits_20261021/README.md)
can coexist with this bound as alternative certificates, with their
original addresses retained. Their mass is not added to the right side
of153N. None of the explicitly labelled recovered Session21 reconstruction
claims is needed as a proved dependency here.

## 7. Reproduction and remaining obligation

```text
python -B 05-knowledge/results/future_tail_blocks_20261021b/check.py --verbose
python -B -O 05-knowledge/results/future_tail_blocks_20261021b/check.py --verbose
python -B -O 05-knowledge/results/future_tail_blocks_20261021b/audit_receipt.py
```

`probe.py` is exploratory; `check.py` and its immutable receipt are the
production certificate. Algebra controls independently enumerate every
optional signed placement at two terminal ports for all pairs of blocks
in {-1,0,1}^4, and check merge/arrival monotonicity. Hostiles keep an
opposite pair at different nodes, a mixed-total-zero block, malformed
literal types, and the actual N7 cut8 failure. The production cases retain
actual old/new channel assignments, first-stop convention and source
filters. All checks remain active under -O.

The normal production replay completed4579989 body checks. Its canonical
LF receipt SHA256 is
`7731a9a9ddf8e5f0478b8d6455e53c6130870f2004353db2e0fefa345bb509c0`.
An independent full optimized replay passed4579990 checks including the
saved-receipt equality gate, with the same hash. Independent proof/API review
found no ROOT or eventual-arrival premise. A separate consumer recomputed
the four saved envelopes directly over
every literal pending pair, using only `Q(w+h)-Q(w)` and polarization;
it matched1153,557164,1138407,4241482 without the producer's type grouping.
That consumer is saved as `audit_receipt.py`; it checks the envelope arithmetic
independently and does not replace the production source-prefix replay.

The exact uniform target is now to bound U_N(D(N)) and all earlier cuts for
every N>=4096, perhaps using source bands and synchronous coalescence rates.
This packet supplies a valid stopping predicate and examples where it
fires before all sources are grounded. It supplies no universal bound on
D(N), no theorem that the predicate must eventually fire, and no proof
of arbitrary-source Collatz convergence.
