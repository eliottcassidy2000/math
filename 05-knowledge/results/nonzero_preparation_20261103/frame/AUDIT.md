# Independent audit of the fixed-core terminal response obstruction

2026-11-03. **FULL ANALYTIC / SOURCE / NORMAL-AND-OPTIMIZED REPLAY PASS.**
The candidate in [PROOF.md](PROOF.md) is accepted without correction, in its
stated fixed-core, common-prefix, retained-root frame. No producer file was
modified by this audit. Root owns status promotion.

## Independent derivation and types

For an arbitrary terminal tournament U, arbitrary nonempty actual port set W,
and the fixed marked six-core D, the ordered product gives
`h_total=17h`, `f_total=28h+224f`,
`j_total=17j+2s+44f`. Hence the actual attachment identities are

```text
F_r=17h+z(28h+224f), H=F_r-44f-17j-2s, z=2^k.
```

The source and target use the same physical root and fixed clock. W meeting
the terminal SCC is separately required for actual strongness. None of the
coordinate manipulations supplies a terminal graph or an implementation.
The proof allows weak U and does not silently assume it is a strong atom.

I independently checked the source sequence at every depth. C2 transport
preserves `h=f=3^m` at the retained u, and singleton marking forces j=0.
For the start coordinate there is an explicit deletion argument: for m>=2,
`K_m-u` is the original singleton y dominating an actual copy of K_(m-1),
whose root is the original x. Every subsequent C2 birth acts inside that
copy and is dominated by y. The m1 deletion is transitive. Thus
`H(K_m-u)=3^(m-1)`. The zero root gap gives s equal to this odd deletion
count. This verifies `s=N/3` without extrapolating the five finite controls.
Substitution gives the displayed H_m and F_m for every m>=1.

Cutting a Hamilton path after each specified port independently gives
`f>=|W|h`. The counted covers have different prescribed first endpoints, so
summing over ports introduces no overlap. Empty second slots remain included.
This is valid for every nonempty tournament. With one actual port, the two
nonempty prescribed slots in j would share that vertex, so j=0 exactly.

Comparing the slope of an all-clock H return gives
`h'+8f'=27N`. Setting the integer `b=f'-3N` therefore gives
`h'=3N-8b`, and constant comparison gives
`180b+17j'+2s'=2N-d`. The cut bound implies b>=0; nonnegative j',s' imply
b<N/90. Two or more ports would instead require `b>=3N/17`, impossible.
Thus the singleton conclusion has no order cutoff and does not rely on a
finite atom bank. The balance then gives
`s'=N-90b-d/2` and the exact root-response debt `F_r'-3F_m=-136b`.
Subtracting H yields the displayed physical defect change `-136b-d`.

The inherited gap decomposition is used with its actual directions:
`f'-h'=q-s'+delta`, where q is the Hamilton count after deleting the port,
q>=s', delta>=0 is even, and q is odd. The empty deletion has count1, so
there is no singleton exception. Since `f'-h'=9b`, substitution forces b
odd for d2 and even for d4. For m<=4 the balance leaves only b0, excluding
d2. At m5 the nonnegative possibilities b0,1,2 leave only b1 under parity,
giving exactly `(h',f',s')=(721,730,152)`. The displayed k1 source, child,
response debt136, and defect32424 are arithmetically consistent. The proof
correctly labels realization of that tuple OPEN.

For the stronger statement, require both H and F to triple as specified.
The constant and slope of F first force `h'=f'=3N`. The cut bound forces a
single target port. At each zero-gap singleton source or target, the two
nonnegative terms q-s and delta must separately vanish, hence s=q is odd.
The H equation demands `s'=3s-d/2`, which is even for d2. This excludes a
response-preserving +2 for every actual singleton tail h=f, including the
new transitive-clone output, independently of the special s=N/3 sequence.
The conclusion is not an obstruction to changing D, its marks, or the
prefix/root attachment; it is not an unrestricted Hamilton-spectrum result.

## Source and exact universe

I read `frame.py` in full. It imports no earlier mathematical implementation.
Its small path tables enumerate literal vertex permutations on every induced
subset, then form the seven coordinates by actual ordered subset partitions.
All 1,099 labelled nonempty tournaments through order5 are included, together
with all 32,767 nonempty port sets and all 5,405 single roots. Weak tournaments
and empty complementary path slots remain in the universe. The gap checks
retain actual indegree, deletion count, parity and singleton j.

The growing controls construct the literal C2 cuts and migrate the auxiliary
while retaining u. Endpoint deletion tables separately check the terminal
signatures through m5 and full attached graphs for m1/2, k1/2. The arithmetic
loop through m8 is clearly labelled a necessary-parameter universe, not a
realization census. Its integer bounds and parity filters match the proof.
Every check is an always-active `need`, so Python optimization removes none.

## Independent replays and frozen pins

Both commands were run independently with fresh output files:

```bash
python3 -B 04-computation/nonzero_preparation_20261103/frame.py \
  --output /tmp/nonzero_frame_audit.json
python3 -B -O 04-computation/nonzero_preparation_20261103/frame.py \
  --output /tmp/nonzero_frame_audit_optimized.json
```

Both pass **51,259 gates**. Both 17,972-byte receipts are byte-identical to
the frozen producer receipt; both 166-byte stdout streams agree, SHA256
`3be6e395ad2f6324b3d24f567f1fd8b7ef6e07e162d82ac2bd232550118f4e9a`.

| Audited artifact | Bytes | SHA256 |
|---|---:|---|
|Primary proof, pre-promotion|7,178|`22c30062fb3b40df5da20746fabe44e52274d4ba7571cc23f9ae664f38d034bd`|
|`frame.py`|7,260|`7ae438d343f0ae84e66c986003f6a6abe6756c225b00ddd5091b9934135a7a21`|
|`frame.json`|17,972|`b317cdafa85e07160b8f4a558819ac1a9d6cd86b714688e3b4207dc67465088a`|

**Final verdict: PASS.** The fixed-frame obstruction, forced singleton port,
root-defect debt, parity threshold and stronger response-preserving theorem
are proved by the argument audited above. The first surviving scalar tuple
remains a necessary target, with actual realization OPEN. No correction or
producer change is requested.
