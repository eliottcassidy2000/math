# A fixed-label native-K obstruction to literal tournament tape writes

2026-10-26. **Verdict: PROVED, with a strictly fixed-label, fixed-order scope.**
This is one final gate on Session65's tape carrier, not a new general
Turing-machine impossibility claim. No repository files, commits, or pushes
were made by this subtask. Inputs are pinned in `SOURCES.json`.

## 1. Exact statement and proof

Let T be a tournament on the fixed labelled set [n], n>=1, and put

    U(T) = #{(i,j): i<j and i -> j in T},
    chi(T) = (-1)^U(T),
    C_n = binom(n,2) - (n-1) = (n-1)(n-2)/2.

For a directed Hamiltonian word Q=(q_0,...,q_(n-1)) of the current source,
the native operation is

    K_Q(T) = T^op xor E(Q):

retain its n-1 path arcs and reverse every other arc. The path pairs are
distinct because Q visits each vertex exactly once. Consequently exactly
C_n unordered pairs change orientation. Every changed pair changes U by
either +1 or -1, both odd. Therefore

    chi(K_Q(T)) = (-1)^C_n chi(T).                   (1)

For any finite sequence of L legal full-order native K moves, with each
word directed in its own successive source, induction gives

    chi(T_L) = (-1)^(L C_n) chi(T_0).                (2)

No hypothesis about blue/blackness was used; restricting the allowed moves
to black words cannot invalidate this necessary condition. Reducing n mod4
shows that C_n is even exactly when n=1,2 mod4. Hence at these orders chi is
constant on every fixed-label all-K component and every black-K component.
A single labelled arc reversal changes chi, so no such finite sequence
can implement that reversal. This excludes sequences of **any length**, not
merely sequences of even length or a bounded search horizon.

Equivalently the literal number of pairs on which the endpoints differ
must satisfy

    d_labelled(T_0,T_L) = L C_n mod2.                (3)

When C_n is odd, a single-arc target needs an odd number of moves; that
necessary parity condition is not a construction of legal words, still
less a black realization. A single full-order K move reverses more than
one arc at every n>=4, including odd-C orders.

## 2. Application to Session65's tape carrier

Use fixed consecutive three-vertex blocks. G0 has arcs 0->1, 1->2, 0->2;
G1 reverses only 0->2. For a binary word w of length m>=1 let

    E(w) = G_(w_0) join ... join G_(w_(m-1)),

with every inter-block arc pointing from an earlier block to a later block.
All pairs point upward except one distinguished pair in each 1 block, so

    U(E(w)) = binom(3m,2) - #ones(w).

A changed bit in cell j reverses precisely the labelled arc on
{3j,3j+2}; every other arc is unchanged. Such a write flips chi. Since

    3m = 1 or 2 mod4  iff  m = 3 or 2 mod4,

**a literal changed-bit write on E(w) cannot be implemented by any finite
sequence of full-order native K moves when m=2,3 mod4.** In particular it
cannot be implemented by a black-only sequence. The argument applies to
every tape word and every cell in these orders. An instruction that writes
the already-present symbol changes no arc and is not obstructed by this
argument. Head, state, program, origin, and window sidecars remain separate
from the tournament.

This strengthens the existing distinction between the controlled local
arc rewrite and native K without changing the tape encoding, its decoder,
or its exact isomorphism-edit/Hamming metric theorem. It does not require
the metric-product theorem: the present obstruction follows directly from
the literal labelled encoding and the definition of native K.

## 3. Relabelling hostile and exact limits of the conclusion

For a permutation pi acting by old label -> new label,

    chi(pi.T) = sgn(pi) chi(T).                      (4)

Proof: a pair {i,j}, i<j, retains its upward/downward convention when
pi(i)<pi(j) and reverses that convention when pi(i)>pi(j). Thus exactly
the inversion pairs of pi toggle the upper-arc indicators modulo2.
Their number has parity sgn(pi). An odd relabelling therefore gives an
isomorphic tournament with the opposite chi value. For n>=2, chi is not
an invariant of an unlabelled isomorphism class.

Accordingly, this packet does **not** establish an obstruction on the
unlabelled class graph, with arbitrary gauge/relabel steps, with a changed
vertex order, with padding, with a larger encoding gadget, or to all
Turing-machine simulations by tournaments. If relabellings pi_i are freely
interspersed with the L native moves, the endpoint multiplier is instead
(-1)^(L C_n) times the product of sgn(pi_i); the conserved sign need not
survive. This observation removes the stated obstruction; it does not
prove any desired unlabelled or gauge-enabled reachability result.

The operation is also **full-order K**, not K restricted to an individual
marked three-vertex factor. A local G0-to-G1 K on that factor would directly
perform its write, but allowing that different operation is outside this
obstruction and outside an asserted full-order native-K simulation.

### Padding hostile

At m=2 (n=6, C_n=10), append two common zero cells; at m=3 (n=9,
C_n=28), append one common zero cell. Both give m=4, n=12, C_n=55.
The source and target still differ in one designated arc, but native K
now flips chi rather than preserving it. Thus odd-length chains are
parity-compatible. This only **removes the conserved-parity obstruction**;
it does not produce a legal chain or a black chain. The check additionally
confirms that the identity Hamiltonian word at n=12 has the target sign
but is not the single-write target: it reverses 55 arcs, not one.

### Odd-C legal control, not a black control

At n=3, Q=(0,1,2) is directed in G0 and K_Q(G0)=G1 literally. This is a
positive control solely for legal unrestricted K at odd C_3=1. The word
is blue: its positional reversal is an anti-automorphism of G0. It is
not a black witness, and no claim of black realization follows.

## 4. Independently authored finite verification

`check.py` imports no producer or repository code. It represents tournaments
as adjacency matrices and performs the arc reversal explicitly. Checks
raise explicit exceptions and remain active under `python -O`.

For n=1..5 it exhausts **every labelled tournament and every vertex word**,
filters directed Hamiltonian words in the actual source, and checks exact
retained/reversed pair sets, the sign multiplier, continued path legality,
and K involutivity. The independent count is n! * 2^C_n legal pairs.

| n | Labelled tournaments | Legal K pairs | Retained arcs | Reversed arcs |
|---|---:|---:|---:|---:|
| 1 | 1 | 1 | 0 | 0 |
| 2 | 2 | 2 | 1 | 0 |
| 3 | 8 | 12 | 2 | 1 |
| 4 | 64 | 192 | 3 | 3 |
| 5 | 1,024 | 7,680 | 4 | 6 |

Total: **7,887 legal K transitions**. On the same universe the script
checks (4) for all **124,469 tournament/permutation pairs**, not only a
sample of odd permutations.

For all four length-two tapes and eight length-three tapes, it checks all
8+24=**32 directed one-bit writes**, including each cell and both symbol
directions. Adjacency-guided DFS independently exhausts their **16+64=80
legal Hamiltonian words** and confirms the source sign is preserved by
every native K. The DFS counts also agree with 3^(number of ones).
Each of the twelve tapes gets an explicit odd-relabel hostile. All 32
writes receive the common-zero-padding control described above.

The arbitrary-length chain obstruction is a proof consequence of (2),
not a claim that finite chains were exhaustively searched at n=6 or9.

Reproduce from this directory:

    python3 -B check.py > normal.json
    python3 -B -O check.py > optimized.json
    cmp normal.json optimized.json
    sha256sum -c SHA256SUMS

Both executions passed and their outputs agree byte-for-byte. The seal
records the report, source manifest, authored verifier, and both outputs.
The scoped gate is resolved; no further research lane is opened.
