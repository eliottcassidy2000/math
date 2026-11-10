# Session243: prime gaps, Euclid numbers and physical tournament arithmetic

**COMPLETE, 10 November 2026. Independently audited partial results. Chan–Pak,
both Brocard questions, Polignac and prime-Euclid infinitude remain OPEN.**

This session merges the prior work without identifying arithmetic values
with physical Hamiltonian-path witnesses. The main new construction is a
compact literal one-vertex increase by any prescribed even amount. The
arithmetic bridge repairs an omitted offset-one test and makes the Brown
prime-power inverse exact at each specified factorial index.

## Proved results

- [THM-5110: primorial crossing gaps and Brown inverse](../../../01-canon/theorems/THM-5110-primorial-crossing-gap-repair-and-brown-prime-power-inverse.md).
  Two primality bits at P_n±1 restore the consecutive-prime gap lost by the
  Fortunate/lesser-Fortunate convention. In particular, 30030±17 are not
  consecutive primes: 30029 intervenes. Fixed integer offsets d>1 from cP_n
  eventually give composite endpoints. Brown solutions n!+1=m², n>=4,
  correspond exactly to adjacent coprime partitions of the full prime powers
  of n!/4. This is a finite inverse for each n, not an all-n completeness proof.
- [THM-5111: synchronized strong graft windows](../../../01-canon/theorems/THM-5111-synchronized-strong-graft-windows-and-prime-euclid-boundary.md).
  Four explicit source states give strong tournament counts
  6Q−H, 6Q+H, 6Q+3H and 6Q+5H for every nonempty actual gadget with counts
  (H,Q). Transitive gadgets give four consecutive odd counts around 6·2^m.
  The proof includes the exact constant-gap criterion under common padding,
  the small-prime twin-candidate wheel and its composite survivor, and an
  architecture-specific primorial-plus-one obstruction. A different graft
  realizing 31=30+1 prevents an invalid global obstruction.
- [THM-5112: compact strong vertex increments](../../../01-canon/theorems/THM-5112-compact-strong-vertex-increments-and-endpoint-gap-compiler.md).
  For every N>=1, there are strong tournaments T⊂T⁺, differing by one vertex,
  with H(T⁺)−H(T)=2N and H(T⁺)<=5N. For N>=2 their orders are
  ceil(log2 N)+3 and ceil(log2 N)+4. The bound is sharp within this family.
  A separate five-letter word comparison compiles any positive endpoint
  count E into gap 2E, but its lower output is always composite.

These are prescribed-gap results, not inverses for an arbitrary prescribed
path count. In particular they do not supply a 5-to-7 extension, divide a
49-path witness into a 7-path witness, or prove prime infinitude.

## What is inherited, rather than newly solved

Both Brocard meanings occur in the repository: primes between consecutive
prime squares, and the factorial Brown-number equation. Euclid numbers
P_n+1 and Euclid–Mullin factor iteration are distinct. The recovered source
map, precise meanings and primary literature are in [provenance](provenance/PROOF.md)
and [sources](SOURCES.md).

Sessions199/221 already proved all-even-gap existence and infinitude with
strong single-edge pairs. Session227 already gave fixed-root graft
product coordinates. THM-5000 already supplies arbitrarily long strong
odd-count windows. Session242 supplies the inherited prime/7p filters and
scoped physical descent. The present contribution is the compact induced
one-vertex pair, synchronized source family, exact restrictions and repaired
arithmetic predicates, not a new claim of global gap universality.

At the pinned audit baseline, kind-pasteur's THM-5109 was RESERVED and was
not used as a proved dependency. Concurrent checkpoint
`be5ce9a33c369e99327b29cd4b3417ddb087ca95` renamed it and labels its new
parts PROVED in a pre-audit publication. That incoming cone/chain-law work
and its ledger correction are preserved, but are not re-audited or used as
proof dependencies here.
The [synthesis](SYNTHESIS.md) records the six-concept board, transfer maps,
failed implications, surviving results and the remaining physical inverse.

## Proofs, independent audit and exact replay

- [Arithmetic proof](provenance/PROOF.md), [Brown inverse](brown/PROOF.md),
  and [independent arithmetic audit](arithmetic_audit/REPORT.md).
- [Physical proof](constructions/PROOF.md), [gap compiler](constructions/GAP_UNIVERSALITY.md),
  and [independent physical audit](construction_audit/AUDIT.md).
- [Aggregate receipt](verification.json) and [replay driver](verify_all.py).
  Run from this repository with Python 3 and Git available:

```sh
python3 -I -B 05-knowledge/results/session243_prime_gap_bridges_20261110/verify_all.py
```

The receipt records 16 successful normal/optimized replays and 26 deliberately
failing controls. The physical scripts run in isolated Python with only their
explicit trusted source directory inserted for sibling imports. Reviewed
producer/canon hashes and inherited Git blob pins are checked before PASS.
No third-party Python packages are needed.

Independent controls include 484,682 Brown partitions through n=60, all
4,396 source-bank grafts over labelled gadgets through order 5, 3,135 boundary
grafts, 10,000 compact arithmetic targets, 128 literal induced graph pairs,
and 36 endpoint-family examples. Finite tests supplement the all-order proofs;
they are not a new unrestricted spectrum-coverage cutoff.

## Publication and stopping point

The plan was preserved at checkpoint `4201c40ba3313f5326ccaf6e45681995030cfedc`;
the audited arithmetic bridge at `9332291ed45937bfd949b10111a59246cdc201f4`.
This final packet contains the audited physical extensions and replay receipts.
Concurrent unrelated work is preserved.

The remaining Chan–Pak obligation is target-dependent physical occupancy:
for each required prime or 7p target, exhibit an actual smaller state, a legal
constructor and a decreasing rank. Prescribed differences, local prime
support and scalar normalization do not establish this. No named conjecture
is closed. Session243 has ended.
