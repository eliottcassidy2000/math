> **Curated publication copy.** Executable proofs and certificate data are in
> `code/`; independent audits are in `audit/`. Run `verify_all.py`.
> Historical exploratory scripts mentioned below are intentionally omitted;
> their finite searches are not dependencies of the certified theorems.

# Radix27: individual repeated-template obstruction

Status: exact certificate theorem, independently reconstructed by a
standard-library verifier. Scope: TWO FIXED EARS over a transitive spine.
Not a restriction on arbitrary tournament counts.

## Theorem

Let w be a length-three word over{0,1,2,3}, with at least one x-inneighbor
letter and at least one y-inneighbor letter. There are49 such words.
For every such w EXCEPT000, no actual seed state v and constant d can satisfy

    H(M_w^(k+1)v)=27 H(M_w^k v)+d   for every k≥0.

There is no bound on d in this theorem, and no lower bound on the seed count.
More strongly, no nonnegative formal36-vector with both all-prefix sums equal
to1 can satisfy the required scalar recurrence.

### Proof by exact separating identities

Put M=M_w, c=output, and R=c(M−I)(M−27I). A repeated constant-digit recurrence
implies R M^j v=0 for all j≥0. Define ell_xy and ell_yx to be the sums of the
three coordinates with Q=R_run=empty in their respective order blocks.
Every actual state satisfies ell_xy(v)=ell_yx(v)=1.

For each of the48 words, the supplied rational certificate specifies numbers
lambda_j, alpha, beta for which the coefficient row

    G = sum_j lambda_j R M^j + alpha ell_xy + beta ell_yx

is componentwise NONNEGATIVE, but alpha+beta is STRICTLY NEGATIVE.
For a nonnegative v satisfying the recurrence and conserved constants this
would give0≤Gv=alpha+beta<0, a contradiction.

`radix27_farkas_certificates.json` contains all48 exact rational identities.
The final three entries of each stored multiplier vector are legacy slots
for inequalities used during exploration; ALL ARE ZERO. The verifier checks
this explicitly. Thus neither a digit bound nor H≥3 is used in the theorem.

`verify_radix27_farkas.py` independently rebuilds each word matrix from
partition-state transitions, reconstructs the Krylov rows R M^j directly,
and verifies1,728 coordinate inequalities plus48 negative constants using
only integer arithmetic and Python `Fraction`. It does not import the
producer, SymPy, SciPy, or the saved candidate relation rows. Two hostile
mutations are rejected. Normal and optimized receipts match.

## Corollary: no stateless complete length-three radix27 compiler

Suppose a fixed physical template is assigned to every digit in
{-26,-24,...,24,26}, and the assigned action H→27H+d must hold at the seed
H>1 and at every descendant. Distinct digits require27 distinct templates.

A length-three word all of whose letters have x-bit1, or all of whose letters
have y-bit1, has at most O(2^m) growth when repeated. It cannot implement a
radix27 action on H>1 with |d|≤26, whose repeated counts grow at least
27^k(H−1)+1. Removing these15 slow words leaves exactly49 candidates.
The theorem excludes48; only000 remains, fewer than the27 needed. ∎

This is a fixed assigned-digit result. Do not silently generalize it to
phase-dependent changing physical banks or other carriers. The preceding
radix3/radix9 results and the spectral filtration have their own explicitly
stated scopes.

## Exploration and proof status separation

* `radix27.py` generated exact individual Krylov relations and tested87,381
  actual seed words of lengths0..8. No H>1 seed supported even ONE repeatable
  template with |d|≤26. Up to34 templates can nevertheless pass the first-step
  scalar test, illustrating why first-step matches are misleading.
* `radix27_relaxation.py` and `radix27_simple_lp.py` are FLOAT exploratory LPs,
  not proofs. One simple LP returned numerical uncertainty for203; the exact
  certificate nevertheless settles that word.
* `radix27_farkas.py` uses a floating solver only to propose multipliers, then
  accepts them only after exact rational verification. The standalone verifier
  is the proof gate and requires no numerical solver.

The 000 exception is real as a potential pure-ternary family mechanism:
H(2330^k)=14·3^k−17 for k≥1. This does not provide a complete digit bank.
At this stage length≥4 was unresolved. The following extension settles
length4; length≥5 remains outside this package. No arbitrary-length claim
follows from the uniformity of these finite word-family certificates.

## Separate length-four extension: radix81 fixed banks excluded

The same symbolic operator test (not a seed census) was applied to the225
eligible length-four words. Exact separating identities were found for220.
The only words not settled by this RELAXATION are

    0000,0303,1212,2121,3030.

This is NOT a claim that those five are realizable affine templates; it is
an explicit boundary of the current certificate. In particular the stronger
arbitrary-length assertion “only the all-zero word survives the relaxation”
is not supported by these data.

`periodic_length4_certificates.json` and the independent standard-library
`verify_periodic_word_certificates.py 4` establish7,920 coordinate inequalities,
220 strictly negative separating constants, and two mutation controls.
Normal/optimized receipts match. No digit or H-bound inequalities are used.

A complete fixed length-four radix81 bank would require81 distinct templates.
Words with one ear-bit constantly1 are excluded by slow repeated growth;
among the225 remaining words the certificates leave at most five. Hence no
such bank exists. This is another exact finite operator-family theorem, not
an arbitrary-length classification.

The generic producer `periodic_word_certificates.py L` uses numerical LP only
for proposals, followed by exact rational verification. The generic standalone
checker does not use numerical libraries or trust producer relation rows.

## Combined small-radix classification

On this two-fixed-ear/transitive-spine carrier, no STATELESS complete centered
scalar accumulator exists for ANY integer radix2≤B<243. Here a template is
assigned to each digit −(B−1),−(B−3),...,B−1 and must realize H→BH+d on every
descendant from a seed H>1. Even initially allowing unequal template lengths
does not help.

Indeed the common spectral filtration forces B=2^L or3^L for the lengthL of
each repeatable template, so all lengths agree. The prior binary theorem
excludes every2-power; the separate ternary certificates exclude3,9,27,81.
Every other B in the stated range is excluded spectrally. This does not
exclude partial digit systems, changing ears, nontransitive append blocks,
or unbounded/count-dependent physical template selection.
