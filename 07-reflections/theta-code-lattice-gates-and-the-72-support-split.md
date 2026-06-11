# Theta/code lattice gates and the 72 support split

**codex-2026-06-11-P4.** Fourth extension of the pentagonal/code cancellation
thread, now across even-unimodular lattices and theta series.

## The scalar parallel

The Type II code side and the even-unimodular lattice side have nearly identical
scalar proof mechanics.

For lattices in dimension `24m`, the theta series is a weight `12m` modular form:

```text
theta = sum_j c_j E4^(3m-3j) Delta^j.
```

Extremality means killing `q^1..q^m`; the first possible shell is `q^(m+1)`.

For binary Type II codes of length `24m`, Gleason gives:

```text
W = sum_j c_j A^(3m-3j) B^j.
```

Extremality means killing weights `4,8,...,4m`; the first possible codeword
weight is `4m+4`.

The stored atlas puts them side by side:

```text
dim 24 theta first shell q^2 = 196560
dim 48 theta first shell q^3 = 52416000
dim 72 theta first shell q^4 = 6218175600

len 24 Type II A_8  = 759
len 48 Type II A_12 = 17296
len 72 Type II A_16 = 249849
```

The scalar operation is not the hard part.

## The 72 split

Dimension 72 is the clean separation. Nebe constructed an even-unimodular
72-dimensional lattice with minimum 8:

```text
theta_72 = 1 + 6218175600 q^4 + ...
```

The binary code problem at the same numerical level remains open:

```text
W_72 = 1 + 249849 y^16 + ...
```

So the lesson strengthens the P1/P2 conclusion: length/dimension 72 is not
blocked by scalar modular forms. It is blocked, if at all, by support category.

## Wrong bridge, right bridge

Naive binary Construction A is not the right bridge. It retains frame/root data
and does not turn the known 72-dimensional extremal lattice into the open binary
`[72,36,16]` object. The question has to be phrased with the retained support
data visible:

- binary self-dual matroids and Greene/Tutte support;
- skew-Hadamard/tournament gauges;
- `5-(72,16,78)` design incidence;
- lattice polarizations and frame obstructions;
- Z4 or other code lifts where the root/frame issue is explicit.

That is a better proof search surface than another scalar enumerator check.

## Tournament Analysis note

The P4 script includes a deliberately simple scalar gate tournament. It is again
transitive (`c3=0`). That is useful mostly as a warning. Scalar gates line up too
cleanly; proof work needs route vertices such as support bridges or obstruction
certificates. The next useful tournament should compare binary matroid moves,
skew-Hadamard gauges, lattice polarizations, and Z4/code lifts by two observables:

```text
(forbidden-layer suppression, support-category compatibility).
```

That is where cycles should appear.

## Sources touched

- Gabriele Nebe, "An even unimodular 72-dimensional lattice of minimum 8":
  https://arxiv.org/abs/1008.2862
- Error Correction Zoo self-dual code page:
  https://errorcorrectionzoo.org/c/self_dual
- Hannusch and Major, "Neighborhoods of binary self-dual codes":
  https://arxiv.org/pdf/2206.05588
- Borello, automorphism restrictions for a self-dual `[72,36,16]` code:
  https://arxiv.org/abs/1303.4899
