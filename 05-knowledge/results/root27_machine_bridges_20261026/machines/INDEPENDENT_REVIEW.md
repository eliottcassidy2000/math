# Independent review: native clocks and reduction types

2026-10-26 UTC. **ACCEPTED:** the elementary all-order identities and the
canonical marked-configuration encoding in `README.md`; the scoped
computability statements in `REDUCTIONS.md`. **FINITE-EXACT:** the independent
controls below. **CITED / not locally built:** the external Goldbach Lean
development and other primary literature. No claim about classical Collatz,
Goldbach, RH, P versus NP, or all-order black-core connectivity is promoted.

Reviewer: the separate `nonconsecutive_ports` agent. The producer files were
read without modification. This review is independent of the producer's
set-of-ones interpreter: its additional check below uses a dense tape and a
manually specified six-rule exponential machine. No producer module is
imported. This is a mathematical proof audit supplemented by bounded tests,
not an independent formal verification.

## All-order native proofs

The exact tables are load-bearing. In particular the exponential machine is
`1RB0LC_1RC1RA_1LA0RB`, not the cubic table
`1RB0LC_1LA1RB_0RA1LC`. Each displayed token is read before its one native
transition. Translation is a change of tape coordinates, not an extra step.

The bouncer's two sweeps move the alternating block left and then right;
their boundary tokens leave A at 2 and ones at
`-4n-2,-4n,...,0`. Counting the four transcript parts gives `8n+8`.
Summing over n=0,...,k-1 gives `4k(k+1)` and translates the head to 2k.

For the cubic machine the right pass reaches the first zero after k, writes
there, clears k on reversal, and returns through the k positions including
the temporary one at zero. The last bounce clears that temporary one.
Thus its stated intermediate invariant is preserved with k replaced by
k-1. The k=0 closing word fills position 1. The duration is
`sum_(k=1..n)(2k+6)+4=n^2+7n+4`; its cumulative formula is correct.

For the exponential machine, independently tracking the actual occupied
coordinates verifies the induction in (I). Write `h=2-2k` and
`H=2+2k`. The 2k right pairs shift the alternating occupied positions from
`h,h+2,...,H-2` to `h+1,h+3,...,H-1`. The first 2k left pairs shift them
to `h+2,h+4,...,H`. The final left pair clears h-1 and fills h; that last
one at h-1 exists precisely because k<n. The new head is h-2, giving
the stated invariant for k+1 with all other positions unchanged.

At k=n there is no solid interval left. The tail first creates the
alternating block at `h+2,...,2n+2`, then adds h and fills its gaps.
The resulting occupied interval is exactly `[-2n+1,2n+2]`, with A at
`2n+3`. Its length is `4n+2=2(2n+2)-2`, so the target parameter is
`2n+2`. The transcript length is
`2 + sum_(k=0..n-1)(8k+2) + 8n+1+4n+2 = 4n^2+10n+5`.
Iteration gives each of the four blank-tape formulas in the producer.
The positive durations make the direction of the native clock unambiguous.

The two rejected external claims fail independently. At n=1 the wrong
middle pair first requests C0 at the fifth read, whereas the actual read is
C1. The correct first 19 transitions create six ones, hence C4; the
recurrence with `+3` would require eight. This is a source correction for
the literal table and parameterization, not a correction to all bell
machines. The five displayed tables with no missing rule are syntactically
nonhalting on every tape; finite evidence about their patterns is not a
permanent dynamical classification.

## Marked configuration carrier and hostile

The proof of the encoding is sound on its stated space of **canonical**
marked transitive tournaments. Unique linear order, the unique origin,
bits, and the head-state marker recover every integer coordinate and the
complete finite-support configuration. One blank guard at both ends
suffices for one write/move operation, after which canonical trimming
restores the same convention. The special halted value must be retained
on both sides, as in the checker. These facts prove the actual commuting
one-step diagram, not just matching sizes.

The retained program sidecar is necessary: all unmarked tournaments here
are transitive, and each has exactly one Hamiltonian path. Even retaining
read labels, next states, and directions is insufficient without writes.
The producer's two-state hostile indeed returns to the literal blank
configuration in two steps for one table and halts on step three for the
other. The carrier is an exact simulation with unbounded size; it supplies
no bare-tournament halting predicate, black-detour rule, or A000568 formula.

## Quantifiers and external proof boundary

The reduction table consistently distinguishes a fixed sentence, a uniform
family, and a polynomial-time language. The fixed singleton language
argument is extensional. Integer-speed LRC has a finite exact event test;
no real-speed equivalence is inserted. Finite metagraph legality is required
to be decidable, and the proposed A000568 formula is required to be total.
The stated Collatz, twin-prime, and P!=NP forms are presentations, not
claims of hierarchy minimality or completeness for those fixed statements.

The two source theorems `halt_lemma'` and `halt_lemma_rev` in
[Leng's Content.lean](https://raw.githubusercontent.com/lengyijun/goldbach_tm/main/GoldbachTm/Tm25/Content.lean)
have exactly the claimed opposite implications between blank halting and
a failed even Goldbach instance. The source was inspected; it was not
compiled and its complete imported proof/axiom chain was not audited here.
The source's 45-step initialization and the producer's ten later recorded
checkpoints are different evidence types. The native table has 25 states,
49 defined instructions, and its missing instruction is state 24 on zero.
No native observation proves all Goldbach instances, and the source hash
in `README.md` pins the actual table used by the producer.

[Lagarias's primary paper](https://arxiv.org/abs/math/0008177) is the cited
mathematical equivalence used for RH. The further quantifier deduction is
correct: a strict violation of its computable-real upper bound admits a
finite rational interval certificate; dovetailing n and precision makes
failure recursively enumerable. This does not require deciding equality
of arbitrary computable reals, and does not prove RH.

The universal nonhalting compiler is properly restricted to a uniformly
effective total-predicate class; there is no claim to enumerate all total
deciders effectively. The Pi^0_2 totality obstruction applies to a uniform
compiler for every program, not to a possible special reformulation of
one Collatz sentence. Two reductions into nonhalting do not produce an
equivalence between their designated instances. Likewise a single Busy
Beaver value is distinct from a halting oracle, and an unbounded tableau
is distinct from a polynomial-size SAT reduction.

The inherited macro-to-native interface retains its needed hypotheses:
injective section, initialization, closed macro domain, positive finite
durations, absorbing halt, and halt reflection. A terminal predicate or
clock can descend through a quotient only when it is constant on fibres.
No missing precision coordinate is silently removed by this review.

### Generalized-Collatz source addendum

The later source paragraph was separately checked by `recurrence_actions`
on 2026-10-26. The [publisher's primary abstract](https://link.springer.com/chapter/10.1007/978-3-540-72504-6_49)
explicitly states Pi^0_2 completeness for Kurtz–Simon's generalized Collatz
problem. The recovered September implication atlas, row 31, already records
the family-level bridge. The added paragraph correctly keeps the compiled
map and its all-input predicate, labels the theorem CITED, and draws no
completeness or undecidability conclusion for the fixed classical 3n+1 map.
Different encodings into one complete family do not make their designated
instances equivalent or prove implications between named conjectures.
This is an accepted source/scope addendum, not a new audit of the paper's
full reduction proof: only the publisher abstract was independently opened
here, and the author PDF retrieval had failed. No producer code, numerical
claim or earlier executable receipt changed; no code rerun was needed.

## Reproduction and independent control

Producer `check.py` was run normally and with `-O`. Both outputs equal
each other and both stored `normal.json` / `optimized.json` byte for byte:
194 block tests, 6,640 configuration rewrite controls, and ten recorded
Goldbach checkpoints. Explicit exceptions keep those tests live under `-O`.

The following standalone code checks every intermediate exponential block
invariant for n=1,...,32 using dense tape storage. A second direct native
run, without prescribed read words, reaches the first seven blank
checkpoints k=0,...,6. It also checks the wrong fifth read and wrong first
block length. Run the fenced code from this file, for example:

    python -c "from pathlib import Path; import subprocess,sys; p=Path('05-knowledge/results/root27_machine_bridges_20261026/machines/INDEPENDENT_REVIEW.md'); s=p.read_text(encoding='utf-8').split('```python\n',1)[1].split('\n```',1)[0]; subprocess.run([sys.executable,'-c',s],check=True); subprocess.run([sys.executable,'-O','-c',s],check=True)"

```python
import json

RULE = {
    ('A', 0): (1, 1, 'B'), ('A', 1): (0, -1, 'C'),
    ('B', 0): (1, 1, 'C'), ('B', 1): (1, 1, 'A'),
    ('C', 0): (1, -1, 'A'), ('C', 1): (0, 1, 'B'),
}

def need(ok, name):
    if not ok:
        raise RuntimeError(name)

class Dense:
    def __init__(self, n, radius):
        self.radius = radius
        self.tape = [0] * (2 * radius + 1)
        for x in range(-2 * n + 2, 0):
            self.tape[radius + x] = 1
        self.state, self.head, self.time = 'A', 0, 0
        self.prefix = []

    def occupied(self):
        return {i - self.radius for i, v in enumerate(self.tape) if v}

    def tick(self, expected=None):
        need(-self.radius <= self.head <= self.radius, 'dense tape bound')
        pos = self.radius + self.head
        bit = self.tape[pos]
        read = self.state + str(bit)
        if expected is not None:
            need(read == expected, 'independent literal read')
        if len(self.prefix) < 5:
            self.prefix.append(read)
        value, move, state = RULE[self.state, bit]
        self.tape[pos] = value
        self.head += move
        self.state = state
        self.time += 1

    def pairs(self, first, second, count):
        for unused in range(count):
            self.tick(first)
            self.tick(second)

    def invariant(self, n, k):
        h = 2 - 2 * k
        expected = set(range(-2*n + 2, h)) | {h + 2*j for j in range(2*k)}
        need((self.state, self.head) == ('C', h), 'loop state and head')
        need(self.occupied() == expected, 'complete occupied set')

invariants = 0
for n in range(1, 33):
    run = Dense(n, 10*n + 20)
    run.tick('A0')
    run.tick('B0')
    for k in range(n):
        run.invariant(n, k)
        invariants += 1
        run.pairs('C1', 'B0', 2*k)
        run.pairs('C0', 'A1', 2*k + 1)
    run.invariant(n, n)
    invariants += 1
    run.pairs('C1', 'B0', 2*n)
    run.pairs('C0', 'A1', 2*n)
    run.tick('C0')
    run.pairs('A0', 'B1', 2*n + 1)
    need((run.state, run.head) == ('A', 2*n + 3), 'translated target')
    need(run.occupied() == set(range(-2*n + 1, 2*n + 3)), 'target block')
    need(run.time == 4*n*n + 10*n + 5, 'native block duration')
    if n == 1:
        need(run.prefix == ['A0', 'B0', 'C0', 'A1', 'C1'], 'fifth-token hostile')
        need(len(run.occupied()) == 6 and len(run.occupied()) != 8,
             'wrong +3 recurrence hostile')

run = Dense(1, 4096)
checkpoints = []
for k in range(7):
    time = 12*4**k - 18*2**k + k + 6
    while run.time < time:
        run.tick()
    parameter = 3*2**k - 2
    head = 6*(2**k - 1) - k
    need((run.state, run.head) == ('A', head), 'blank checkpoint state/head')
    need(run.occupied() == set(range(head - 2*parameter + 2, head)),
         'blank checkpoint exact interval')
    checkpoints.append([time, head, len(run.occupied())])
print(json.dumps({'status': 'PASS', 'blocks': 32, 'invariants': invariants,
                  'native_blank_checkpoints': checkpoints,
                  'hostiles': ['wrong fifth token', 'wrong +3 recurrence']},
                 sort_keys=True))
```

Independent normal and optimized executions passed with identical output:
32 complete blocks, 560 intermediate invariant checks, and checkpoints
`[time,head,ones]` equal to
`[0,0,0], [19,5,6], [128,16,18], [633,39,42], [2794,86,90],
[11723,181,186], [48012,372,378]`.

No new large census, external formal build, asymptotic dynamical claim,
or resolution of an open conjecture is contained in this audit.
