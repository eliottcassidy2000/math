# A program fits in a marked tournament; its semantics still need a decoder

2026-10-26. **PROVED elementary encoding; FINITE-EXACT controls.** This
is an explicit application of the existing marked-path bit carrier, not a
claim of a new universal-machine theorem or an equivalence between open
conjectures. Its purpose is to make the proposed tournament/machine bridge
exact and expose its operation boundary.

The closest mechanism is the repository's ordered-graph/marked-path
trijection, used in the preceding
[literal repair packet](../../lca_ports_reverse_20261025/PORT_FRONTIER.md).
The hostile is a single write-bit change that alters blank-tape halting.
The corrected near miss is mistaking finite control-graph equality for
configuration dynamics. The least-used sidecar is the marked Hamiltonian
order and the explicit interpretation of its nonpath arcs.

## Exact finite syntax

A binary one-tape machine with s>=1 nonhalting states has 2s transition
records. Each record contains a write bit, a direction bit, and a target
in0,...,s, where s denotes halt. Serialize a unary header1^s0 followed
by the records, with target width bit_length(s). The exact length is

    L=s+1+2s*(2+bit_length(s)).

Choose the least n with binomial(n-1,2)>=L. On vertices0,...,n-1 retain
the directed Hamiltonian spine P=(0,...,n-1). In lexicographic order on
all pairs(i,j) with j>=i+2, orient i->j for serialized bit1 and j->i
for bit0. Pad unused positions with zeros. There are exactly
binomial(n-1,2) available pair slots.

The marked pair(T,P) decodes the complete table uniquely, including
halt targets and writes. The map is injective, and relabelling vertices
while transporting P commutes with decoding. The orientation is a gauge
relative to P; forgetting P is not part of this theorem. This graph's
vertices are code positions, not physical tape cells or machine states.
No temporal meaning is assigned to its Hamiltonian paths.

## A one-arc hostile and a closure boundary

Let state0 be A, state1 be B, and state2 be halt. Use

| Read record | Write | Move | Next |
|---|---:|---|---|
| A0 | b | right | B |
| A1 | 0 | right | halt |
| B0 | 0 | left | A |
| B1 | 0 | left | A |

For b=0, the blank configuration repeats exactly after two steps. For
b=1, it halts after three steps. The two machines have the same control
states, next-state graph, read labels and head directions. Their tables
and encoded tournaments differ at one bit/pair. The write operation
changes the next observation; a preserved graph skeleton is insufficient.
The [native machine lane](../machines/README.md) independently checks
this same example; this is one shared hostile, not two discoveries.

For this serialization, the usual K_P operation complements every
program bit while retaining the spine. Its first header bit becomes0,
which the decoder reads as an invalid zero state count. Thus this exact
encoding's image is not even closed under that K operation. This is a
failure of this specified code language, not an impossibility theorem
for other encodings or guarded graph simulations.

The native lane's complementary configuration carrier uses the intrinsic
left-to-right order of tape positions as a transitive tournament, with
tape bits, head, state and origin retained as labels. Every such underlying
tournament has exactly one Hamiltonian path. Computation therefore cannot
be classified by that unlabelled path count: the information is in the
labels and permitted rewrites. Program syntax and running configuration
are separate carriers; together they specify an executable machine.

## State counts are not semantic invariants

For any table with s live states and any t>=s, keep all old live-state
indices, rename the old halt target s to t, and append t-s unreachable
states. Give each appended state arbitrary valid rules. No old rule
enters an appended state. Induction on native time proves that every
execution starting in an old configuration has exactly the same tape,
head, live state and duration, with only the terminal name changed.
This works for arbitrary initial finite tapes, not just blank input.

In particular both hostile programs above can have 5, 15 or 25 states
while retaining opposite blank outcomes. The specific 25-state Goldbach
construction is meaningful because of its verified transitions and
semantics; the numeral 25 alone is not a computational invariant.
Minimal state complexity would be a different, model-dependent question.
The checker compares all 64 one-state binary tables with their 25-state
padding on four tapes for 80 native steps, plus the two hostile programs
at all three sizes: 262 finite trace controls supporting the induction.

## A useful next theorem to ask for

For a proposed graph transformation F, demand an actual commuting square

    decode(F(encoded_configuration)) = native_step(decoded_configuration),

with an explicit finite clock conversion and preservation of the halt
predicate, ports and boundary extension. A simulation that merely agrees
on a final scalar count is weaker. Bidirectional simulation, injective
configuration conjugacy, and polynomial reductions of uniform input
families are different claims. None follows from using the same finite
tournament vocabulary for two programs.

The executable control checks every complete binary table at s=1,2:
20,800 tables, the K closure hostile for each, the two exact native
trajectories, and all40,320 joint frame relabellings of the order8 example.
It uses explicit exception gates and runs identically with Python-O.
