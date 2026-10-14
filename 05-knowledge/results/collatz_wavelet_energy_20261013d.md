# Source preserving ternary refinement for Collatz energy

2026-10-13. **PROVED:** exact refinement, its Gram identity, a uniform
rounding budget, and the obstruction to a depth-blind sibling contraction.
**FINITE-EXACT:** 1,308,622,770 refinement checks and a depth-24 compressed
packet. **CONDITIONAL:** a packet decay criterion for power cancellation.
**OPEN:** the infinite dyadic energy bound, an asymptotic epsilon rate,
nonconditional amplitude, and universal ROOT entry. Self-audited, with an
independent arithmetic implementation for the stated finite controls.

The target inherited from [the amplitude certificate](collatz_amplitude_boundary_20261013c.md)
is `sum_(X<=n<2X) H(n)^2 <= 2 X^2` at all sufficiently large dyadic X.
The new rules retain the arithmetic source and the covariance that a scalar
energy or a fixed noise model discards. They do not close that bound.

## Inheritance and concept board

Closest proved mechanisms: the ROOT-supported signed transfer in the
amplitude note; the source-labelled carry polynomial of the
[five-mirrors note, Proposition 4](collatz_five_mirrors_20260929.md);
and [THM-4630, Q covariance model](../../01-canon/theorems/THM-4630-hofstadter-q-jitter-is-the-silver-ratio-arma-fixed-point-rho1-equals-1-minus-sqrt2.md),
with its [exact Pell conjugacy](golden_nested_trees_20261013c.md#6-incoming-q-covariance-work-a-genuine-two-step-pell-mechanism).
The exact class hierarchy in [THM-4629](../../01-canon/theorems/THM-4629-twisted-transfer-law-step-count-characters-of-collatz-orbits-hear-the-near-returns-by-congruence-class.md) already consumes binary digits and produces ternary digits. Formula (2) below is its centered, finite-prefix realization with exact cutoffs, not a claim to have first discovered that hierarchy. The Q theorem's sampling assumptions are not proved for Q, and are not
assumptions about Collatz. The incoming audit `cef0e4e3a96f` also rejects the
iid interpretation of the measured Collatz class discrepancies.

The hostile is a covariance changing sign with scale; the corrected near
miss is interpreting a fitted epsilon slope as an asymptotic rate. The
least-used sidecar is the *ordered* ternary address, together with its
endpoint carry. The six live concepts are ROOT phase, ternary wavelets,
Gram matrices, depth-dependent renewal, the Conway Pascal bridge, and
golden carry memory. The [Conway/golden companion](conway_golden_memory_20261013d.md)
closes the actual Conway maximum limit from classical inheritance, while
keeping the continuum approximation separate.

## 1. Exact wavelet refinement

Use `T(n)=n/2` for even n and `(3n+1)/2` for odd n. Define
`t(n)=(-1)^tau(n)` on the basin of 1 and zero elsewhere. Since the terminal
cycle `1 -> 2 -> 1` has even length, the homogeneous identity
`t(n)=-t(T(n))` holds even at 1. No universal convergence is presumed.

For q>=1 and 0<=r<3^q put

    psi_(q,r)(n) = 3 [n=r mod3^q] - [n=r mod3^(q-1)],
    H_(q,r)(N)   = sum_(1<=n<=N) psi_(q,r)(n)t(n).

The old discrepancy H is H_(1,2). Siblings satisfy, identically in N,

    H_(q,a) + H_(q,a+3^(q-1)) + H_(q,a+2*3^(q-1)) = 0.       (1)

Define addresses

    e_q(r) = 2^(-1)r mod3^q,
    o_q(r) = 2^(-1)(3r+1) mod3^(q+1).

For **every integer N>=0**,

    H_(q,r)(N)
      = -H_(q,e_q(r))(floor(N/2))
        -H_(q+1,o_q(r))(3 ceil(N/2)).                         (2)

Proof: the even sources n=2m contribute the first term, since
`t(2m)=-t(m)`. For odd sources n=2m-1, put j=3m-1. Then
`t(2m-1)=-t(j)` and the transported weight is psi_(q+1,o_q(r))(j).
That child weight vanishes outside j=2 mod3; summing through
`3 ceil(N/2)` includes precisely those j. The argument also covers q=1,
whose parent modulus is one. The empty prefix N=0 is the equality boundary.

This is an identity on the original sequence, not an independent choice of
sign at each refined terminal. Changing a source, root phase, or cutoff
changes the obligation.

## 2. The energy rule retains a full Gram matrix

Fix a source window X<=m<2X and write

    x_m = H_(q,e_q(r))(m),
    y_m = H_(q+1,o_q(r))(3m),
    v_m = H_(q,r)(2m) = -x_m-y_m,
    z_m = psi_(q,r)(2m+1)t(2m+1).

With sums over that *same window*, set
`xx=sum x_m^2`, `yy=sum y_m^2`, `xy=sum x_m y_m`,
`vz=sum v_m z_m`, and `zz=sum z_m^2`. Then

    sum_(2X<=n<4X) H_(q,r)(n)^2
      = 2(xx+yy+2xy) + 2vz + zz.                              (3)

This follows by adding `v_m^2+(v_m+z_m)^2`. The first three entries are a
positive-semidefinite Gram matrix, so `xy^2<=xx*yy`. Its off-diagonal entry
can have either sign. The odd-prefix carry supplies two more indispensable
terms. For q=1,r=2,X=1 the entries are `(4,1,2,-3,1)`, giving energy13;
dropping xy gives5, and dropping the carry gives18.

The [census](collatz_wavelet_energy_20261013d.json) checks every
`1<=m<2^24`, every residue r, and q=1,2,3 (39 channels), using the previously
authenticated parity bits. Selected H_(1,2) entries, normalized by X^2:

| log2 X | xx/X^2 | yy/X^2 | xy/X^2 | child energy/(2X)^2 |
|---:|---:|---:|---:|---:|
| 16 | 3.489802 | 2.239625 | -2.476625 | 0.388100 |
| 18 | 1.517509 | 0.233723 | -0.501365 | 0.374253 |
| 19 | 0.917934 | 0.104143 | +0.221983 | 0.733021 |
| 23 | 1.140095 | 0.711798 | +0.842815 | 1.768761 |

At X=2^23 the exact cross-sum is **59,307,814,126,550**. The correlation is
about +0.936; at X=2^17 it is about -0.949. Thus neither a universal negative
cross-term nor Q's constant `1-sqrt(2)` is available for these channels.
The old constant 2 cannot simply be required for every auxiliary residue:
H_(1,1) has normalized dyadic energy above 4.60 at X=2^25 in the inherited
all-residue census. A successful invariant must have channel-dependent
budgets and retain cross-information.

An appropriate proof object is `(q,r,source window, endpoint maps, Gram
matrix, odd carry)`. A verified matrix inequality can propagate an energy
budget. A measured correlation cannot supply that inequality at later
scales. Even abstract vectors x=1,y=1 and x=1,y=-1 have identical individual
energies and combined energies4 and0.

## 3. A compressed refinement packet

For an aligned source N=2^L M, iterate (2). A word with k odd branches
ends in `H_(1+k,r)(3^k M)`, with overall sign (-1)^L. Combine equal pairs
(k,r). Within each sibling triple subtract the median of its three integer
coefficients; (1) makes this exact. The median minimizes
`|c_0-a|+|c_1-a|+|c_2-a|` over integer a.

Start with `c_(0,2)=1`, refine, and apply this reduction after each level.
Let R_L be the resulting coefficient l1 norm. Both child maps send sibling
relations to sibling relations, so reducing between refinements is sound.
It follows that

    H(2^L M) = (-1)^L sum_(k,r) c_(L,k,r) H_(1+k,r)(3^k M),
    R_(L+1) <= 2 R_L.                                         (4)

The reduction preserves this observable identity. To recover the original
word history as well, retain the expansion/reduction DAG; the reduced
coefficients alone intentionally forget it. These are identities of sums,
not certificates that every contributing source reaches ROOT.

| L | surviving states | R_L | R_L/2^L |
|---:|---:|---:|---:|
| 3 | 7 | 7 | 0.875000 |
| 6 | 32 | 32 | 0.500000 |
| 12 | 1,360 | 1,383 | 0.337646 |
| 18 | 65,691 | 70,227 | 0.267895 |
| 24 | 3,193,187 | 3,585,706 | 0.213725 |

These exact reductions are real cancellation, without a probabilistic
premise. The ratio is nonincreasing by (4); a uniform decay rate is OPEN.
The saved odd-count grading gives a further finite diagnostic: at L=24,
the l1-weighted mean k/L is 0.577868, and only 0.211812 of that norm lies
where `3^k>=2^L`. Thus the remaining packet is not confined to the
noncontracting affine slopes. These are coefficient statistics after
quotienting, not probabilities or densities of integer trajectories. A
renewal argument that discards every smaller-slope term would lose most
of this particular residual obligation.

### Arbitrary source and a uniform terminal carry

Write N=2^L M+s, 0<=s<2^L. Along a word let k_j be its odd count after j
letters. Its endpoint after j letters has the exact form

    2^(L-j) 3^k_j M + s_j,
    s_0=s,
    s_(j+1)=floor(s_j/2)  or  3 ceil(s_j/2).                  (5)

For the normalized residual `u_j=s_j/(2^(L-j)3^k_j)`, an even step never
increases u. An odd step increases it by at most
`1/(2^(L-j)3^k_j)`, with j indexed before that step. Consequently
`0<=u_L<1+sum_(j=0)^(L-1)2^(-(L-j))<2`, and

    0 <= s_L < 2*3^k_L.                                     (6)

The wavelet's parent progression has spacing 3^k_L. The interval in (6)
therefore adds at most two supported terms, each of absolute weight at
most2. The error per *raw word* is at most4. Subtract the aligned packet
before reducing it; do not reduce terms with different endpoints as if
their carries agreed. There are 2^L raw words, giving the universal rule

    |H(N)| <= 2 floor(N/2^L) R_L + 4*2^L,  for all N,L>=0.   (7)

Here `|H_(1+k,r)(3^k M)|<=2M`, since exactly M sources lie in its parent
progression. This also covers M=0, including sources below the chosen
refinement scale. The estimate uses |t|<=1 and its exact recurrence only.

**Conditional power criterion.** If an independently proved uniform bound
`R_L<=C rho^L` holds with 1<=rho<2, choosing
`L=floor(log2(N)/(2-log2(rho)))` in (7) gives

    H(N)=O(N^beta),  beta=1/(2-log2(rho))<1.                  (8)

Also, merely `R_L/2^L -> 0` implies H(N)=o(N), by first fixing L in (7)
and then taking N to infinity. The inherited transfer has
`epsilon_M=-H(3M)/(6M)`, so (8) gives its corresponding power rate.
Qualitative H(N)=o(N) already follows from the
[preceding parity-mean theorem](collatz_signed_resonance_pair_mass_20261013b.md)
and this transfer identity, on its stated [THM-4581, Haar coalescence](../../01-canon/theorems/THM-4581-haar-coalescence-affinely-related-collatz-orbits-merge-almost-surely.md) input. The
new target is a quantitative packet estimate independent of that input,
not a stronger claim based merely on rederiving qualitative balance.
This is a new sufficient target on a completely explicit arithmetic packet;
it is not yet a proved estimate on that packet. For Mellin continuation to
the Mercator pole one needs beta below its real part0.9982092207..., plus
quantitative constants if using the finite noncancellation certificate.
Continuation alone does not certify a nonzero residue. Nor does this
triangle-inequality route supply the desired square-root energy bound.

## 4. Why a fixed contraction block cannot be recycled blindly

Here the old Collatz carry polynomial is exactly the missing memory.
For a word w of length L with k odd letters, starting from (q,r), its final
address satisfies

    r_w = 2^(-L)(3^k r+B_w) mod3^(q+k),
    B_empty=0,
    E: B stays B;  O at position j: B becomes 3B+2^j.         (9)

Positions j begin at zero. For fixed k,

    0 <= B_w <= 2^(L-k)(3^k-2^k) < 2^L 3^k.               (10)

Indeed the ith odd position is at most L-k+i, and summing its contribution
gives (10). Distinct words of fixed k have distinct B modulo2^L: the
congruence `3^k a+B=0 mod2^L` reconstructs the parity word of a, one bit at
a time. This is the usual finite parity-cylinder bijection.

If `3^(q-1)>2^L`, distinct words of the same k consequently have distinct
parent residues modulo3^(q+k-1). Their B differences are smaller than that
modulus by (10), so no pair can occupy a common sibling triple. The same
argument applies to every intermediate depth. Thus **for each prescribed
block length L there are arbitrarily deep starting wavelets for which the
sibling-only l1 grammar has exactly 2^L surviving unit coefficients**.
There is no strict uniform block contraction over all address depths in
that grammar. This does not exclude stronger rules using actual terminal
correlations or depth-dependent renewal.

Already after two letters from q=1,r=2,
`E after O` has address4 mod9, whereas `O after E` has address2 mod9.
Counts of the two letter types do not retain the source address.
The first triple-block reduction gives norm 7 at q=1 but norm 8 at every
q=2 residue. This is the cheap hostile against iterating a7/8 factor.

The five-mirrors work similarly excludes a uniform one-step Fourier gap
and retains carry, conductor, and renewal scale. The explicit bridge here
is (9), not an identification of its random Syracuse law with this ROOT
character. Fixed word length here and fixed odd count there condition on
different ensembles. A renewal estimate must retain that conditioning.

## 5. Reusable rules and remaining targets

An ordered block with values v_j has a useful exact signature:
`(length L, total S, U=sum prefix, E=sum prefix^2)`. Concatenation obeys

    S_AB=S_A+S_B,
    U_AB=U_A+U_B+L_B S_A,
    E_AB=E_A+E_B+2 S_A U_B+L_B S_A^2.                       (11)

Thus energy needs the incoming offset and first prefix moment, even when
the terminal total is zero. For vector-valued prefixes, retain the matrix
of cross-products. This is the deterministic alternative to fitting one
covariance number. The Conway max-plus control uses precisely such
context retention, with the denominator n carried through concatenation.

More formally, on a fixed source window let V collect terminal observables
and d their retained carries. If refinement gives `U=S V+d`, put
`W=(V,d)` and `A=[S I]`. Then `Gamma_U=A Gamma_W A^T`. An authenticated
upper certificate `B-Gamma_W` positive semidefinite transports to
`A B A^T-Gamma_U` positive semidefinite. The full matrix includes the
terminal/carry cross-block. This is a sound matrix proof rule for changing
controllers on unchanged sources; finding uniformly bounded certificates
B compatible with all refinements remains OPEN.

The BSD [Tunnell continuation](bsd_tunnell_terminals_20261013c.md) is a
different kind of source-preserving switch: two grounded tests cover 3,610
of the same 4,056 curves, leaving 446. It adds 724 formerly uncovered curves
under the note's stated arithmetic dependencies. That union operation
suggests multiple terminal *proofs* for a fixed source; it does not make
two additive Collatz channels probabilistically independent. Their overlap
still belongs in the Gram matrix. Golden digit normalization likewise
preserves exact value, while losing the original spelling unless the
rewrite history is retained.

The most specific next targets are:

1. Bound packet mass conditional on both odd count and conductor depth,
   using (9), instead of repeating a shallow7/8 contraction.
2. Construct a cone of source-compatible Gram matrices with a proved
   refinement inequality and an explicit odd-carry budget. Test the
   strongly positive X=2^23 row before claiming any sign restriction.
3. Compare such a cone with the exact carry-conditioned Syracuse renewal
   estimates; any transfer must preserve the length/cost conditioning.
4. Use the proved Conway limit and golden canonicalization as positive
   controls for a proof compiler, while rejecting free terminal signs and
   arbitrary source reencoding.

## Reproduction and controls

The C++ producer uses signed128-bit sums and explicit product-overflow
guards, checks 4,608 independent strict first-hit paths, and checks every
first-step and odd-carry identity in its declared universe. Its bit input
is the two-run authenticated output of the preceding amplitude census.
The new Python audit uses no bit input: it computes every strict orbit
through 6,144, checks 159,783 all-N identities and all sibling identities,
and independently reconstructs 429 Gram rows (all 39 channels, k=0..10).
It compares raw and reduced packets at L<=10 on four aligned sources,
and checks every word and residual through L=8 (87,380 carry checks).
The depth 24 packet is integer-exact; its future behaviour is not inferred.
Normal and optimized Python runs are required to agree.

```sh
clang++ -O3 -std=c++17 04-computation/experiments/collatz_wavelet_energy_20261013d.cpp -o /tmp/collatz-wavelet
/tmp/collatz-wavelet /tmp/math-amplitude-20261013c/parity26.bits 24 3 05-knowledge/results/collatz_wavelet_energy_20261013d.json
python3 04-computation/experiments/collatz_wavelet_audit_20261013d.py --census 05-knowledge/results/collatz_wavelet_energy_20261013d.json --depth 24 --output 05-knowledge/results/collatz_wavelet_audit_20261013d.json
```

Saved [audit and packet table](collatz_wavelet_audit_20261013d.json).
The bitfile is a reproducible scratch artifact, not a tracked dependency;
its producer and the independent jump-depth comparison are documented in
the preceding amplitude packet. No hypothetical non-ROOT source is counted
as tested by these finite strict-orbit controls.
