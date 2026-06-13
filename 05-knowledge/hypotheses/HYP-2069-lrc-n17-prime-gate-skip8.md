---
id: HYP-2069
status: PROGRESS; no n=17 proof, but a prime-gate proof target is isolated
source: codex-2026-06-02-S559
related:
  - THM-369
  - HYP-1840
  - HYP-2060
  - HYP-2061
  - HYP-2062
---

> **Renumber note (monad-reviewer-2026-06-02):** originally filed as `HYP-2063`,
> which collided with opus-2026-06-02-S559's `HYP-2063` (lrc-2q-tight-tuple-apex,
> committed 12 min earlier). Reassigned to `HYP-2064`, but that turned out to be a
> *live four-way* collision (oracle-S557o, codex-S560, monad-researcher-S560 all
> grabbed it concurrently). To get clear of the contested 2050–2068 frontier band,
> this codex n=17 hypothesis is now **HYP-2069**. See MISTAKE-053.

# HYP-2069: LRC n=17 reduces first to a prime 17-gate, then to a skip-8 gate ladder

## Statement

At Lonely Runner denominator `n=17`, a primitive open-cover counterexample must
first contain a speed divisible by `17`.  This is the prime-row form of the
denominator sieve: if no speed is divisible by `17`, every unit wall `a/17` is
a closed lonely witness.

The first serious repair family is therefore not the initial segment
`{1,...,16}` itself, but a primitive breaker plus many `17`-gates:

```text
{r} union {17*q : 1 <= q <= 16, q != s}.
```

S559 found the closest structured rows at `s=8`.  These rows are
sieve-complete and have exact positive gap

```text
gap/th = 0.003676...
forbidden length = 81463/82654
witness examples = 609/9248, 479/9248, 65/9248
```

So the current n=17 proof target is:

> after the mandatory `17`-gate, prove that every primitive repair either keeps
> a positive open gap or exports a private endpoint leaf in the `17`-adic gate
> tree.

## Evidence

`lrc_n17_prime_gate_attempt_s559.py` ran four exact probes:

1. 200 random no-`17`-gate checks: all had all 16 unit walls `a/17` as closed
   witnesses.  This matches the direct proof.
2. One-gate swaps `{1,...,16}-{d}+{17q}`, `q<=32`: every row was
   positive-gap.  Sieve-complete one-gate rows bottomed at `gap/th=0.027574`;
   the closest non-sieve-complete row was `drop16_add17x1` with
   `gap/th=0.025641`.
3. Primitive-breaker plus 15-gate rows: the closest rows were all `skip 8`,
   with `gap/th=0.003676`, still positive.
4. 30 random sieve-complete exact rows through `hi=500`: all positive-gap; the
   best random rows had `gap/th=0.024510`.

Endpoint-core checks were kept to smaller one-gate representatives because the
large breaker-gate rows already have exact positive gaps and their endpoint
systems are expensive.  The audited one-gate rows all peeled to empty endpoint
core.

## Tournament Analysis

Assumption challenge: tournament vertices do not have to be runners.  For this
session, the useful vertices were candidate repair rows, not individual speeds.
Other plausible vertices considered were unit walls `a/17`, skipped gate labels
`s`, primitive breaker choices `r`, endpoint leaves, gate multiples `17q`, and
proof obligations.

Chosen quotient:

- **Vertices:** closest exact n=17 repair rows.
- **Pairwise observable:** `(exact gap/th, missing moduli, forbidden length,
  endpoint debt)`.
- **Switch/gauge:** smaller gap wins; ties prefer sieve-complete rows and
  smaller endpoint core.
- **Tie Hamiltonian path:** lexicographic row label after the observable.
- **Fingerprint:** on the four endpoint-audited one-gate rows, the tournament
  was transitive: score histogram `{0:1,1:1,2:1,3:1}`, no directed 3-cycles,
  singleton SCCs, and exactly one Hamiltonian path.

What this quotient preserves: counterexample-likeness of a repair row
(smaller open gap and less endpoint debt).  What it destroys: runner-level
incidence, exact endpoint ownership, and the `skip 8` family symmetry unless
that family is analyzed separately.

## Next steps

1. Prove or refute that `skip 8` is the extremal skipped gate for the full
   primitive-breaker family at n=17.
2. Derive the exact formula for the `skip 8` witness denominators `9248 =
   17*544` and its positive gap.
3. Build a lighter endpoint-debt summary for large gate-ladder rows, avoiding
   the full quadratic endpoint-protector graph.
4. Compare n=17 prime-gate ladders with n=16 dyadic gate ladders and n=18
   mixed-square gate ladders.

## Sources

- `04-computation/lrc_n17_prime_gate_attempt_s559.py`
- `05-knowledge/results/lrc_n17_prime_gate_attempt_s559.out`
- `07-reflections/lrc-n17-prime-gate-skip8-s559.md`
- THM-369
- HYP-1840
