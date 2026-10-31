# Eight bounded core speeds and five arbitrary tails

**Status:** PROVED by an analytic conductor and a FINITE-EXACT complete certificate, independently replayed in normal and optimized Python.

## F1. General local bad-set envelope

For 0<a=2r<1, an interval of length W and speed v>0, put x=vW.
The bad fraction is at most

    E_a(x) = (a floor(x) + min(frac(x),a))/x.

Indeed each full period contributes exactly a/v and an arbitrary residual
arc of length frac(x)/v contributes at most min(frac(x),a)/v.

Fix k>=2 with ka<1, and let

    N = ceil(a(k-1)/(1-ka)),       c = kaN.

Then E_a(x)<=1/k for all x>=c. To see this, its peaks at n+a have height
 a(n+1)/(n+a), which is <=1/k exactly when n>=N. The defining inequalities
for N imply N-1+a<c<=N, so on the remaining initial descending piece
E_a(x)=aN/x<=1/k. All later pieces lie below their endpoint peaks.

Thus k arbitrary (even shifted) speeds with min(v_i)W>=c have a simultaneous
safe point in the closed interval. Their open bad sets each have measure at
most W/k, so total bad measure is at most W. If they covered the entire
connected interval, at least two nonempty relatively open bad sets would
have a nonempty intersection, hence positive measure overlap. Otherwise
they would disconnect the interval; and a single bad set cannot cover it
because W/k<W. Positive overlap contradicts total measure<=W. This argument
keeps boundary equality, but claims existence rather than positive measure.

At r=1/14 the constants for k=1,...,6 are

    1/7, 2/7, 3/7, 4/7, 10/7, 30/7.

For k=1 use the separate maximum-open-gap argument. This file's finite
certificate uses only k<=5. The k=6 conductor does not by itself certify any
seven-core/six-tail family.

The companion conductor packet and independent audit give the general
analytic proof and controls. No literature-priority claim is made.

## F2. Padding an already present bounded tail

Assume every (m+1)-element subset of a pool P succeeds with k-1 arbitrary
tails. Let C be an m-element subset and consider k tails, one of which is
v in P.

If v is absent from C, use the inherited theorem for C union {v} and the
other k-1 tails. If v is already in C, choose any d in P minus C and use
the inherited theorem for C union {d} and the other k-1 tails. In the second
case the new condition is stronger than the original because v was redundant.
Thus **every bounded tail**, not just a missing bounded speed, discharges
an inherited branch. This uses constraint monotonicity, not a count of copies
or a statement that marginal safe times coincide.

## F3. The complete next layer

Apply F2 to the completed nine-core/four-tail theorem in `LRC.md` and use
F1 for the all-large regimes. For every C subset of {1,...,13} of size eight,
sort five arbitrary positive integer attachments. Maintain their exact
common cover and apply the same finite stopping tree as `LRC.md` P2, now
with constants (1/7,2/7,3/7,4/7,10/7).

Each finite candidate v<=13 is explicitly marked `inherit_nine_core` and
resolved by F2. Every other candidate retains full common closed intervals,
positive internal widths, and original-clock rational witnesses at the
finite leaves. No branch is omitted because it looks unlikely or because a
single marginal attachment succeeded.

The compact producer ledger enumerates all 1287 eight-element cores. Complete successful generation and independent replay prove that every
thirteen-speed row containing at least eight **distinct** speeds in
{1,...,13} is safe at radius 1/14. Attachments may repeat and may already be
in the core. It strictly enlarges the nine-core/four-tail class; for instance
{1,...,8} with any five speeds above 13 is now covered.

## Evidence and replay

    python3 -B five_tail_certificate.py --output regenerated.json
    python3 -B -O five_tail_certificate.py --output regenerated.opt.json
    cmp regenerated.json five_tail_pruned_receipt.json
    cmp regenerated.opt.json five_tail_pruned_receipt.json

The optional `--full-directory DIR` writes every committed tree. Each compact
record preserves exact original safe intervals and sufficient metadata for
an independent full-comb reconstruction of the entire tree and all leaf
witnesses. A SHA256 commitment alone is not mathematical evidence of success;
the supplied full independent reconstruction and normal/-O agreement are
required verification stages.

The earlier unpruned feasibility scripts and their outputs are research
scratch, not proof artifacts or additional claims. This certificate explicitly
uses the stronger padding rule and has its own separately named receipt.

**Remaining frontier:** the six-tail analytic conductor is available but a
complete seven-core/six-tail certificate is not claimed. Arbitrary core shapes,
optimal radius, and unrestricted LRC(14) remain open.

## Completed verification

All 1,287 eight-cores pass. Normal and optimized producer receipts are
byte-identical. Both independent reconstruction modes agree and bind all
1,287 canonical tree commitments, totaling 41,646,941 reconstructed bytes.

- 1,287 roots; 37,610 four-tail and 91,872 three-tail nodes;
- 315,277 two-tail and 373,687 one-tail nodes;
- 24,899 literal finite leaves, each with an original-speed witness;
- 16,731 inherited bounded-tail branches, including padding controls;
- maximum cutoffs for one through five tails: 84,80,81,79,100.

The independent controls pass 4,000 bad-measure-envelope comparisons, 3,012
boundary cases, seven hostile/endpoint cases, and reject eight mutations.
See `../audit_lrc5/AUDIT.md` and its normal/optimized receipts.

The 849,448-byte compact receipt has SHA256
`2dfbc2cb3e7e2969ffa8dedb5283637e4cd1db7bbe939d158dfed570770130cd`.

Among these bounded eight-cores, the smallest longest safe-interval width
is exactly 1/70, attained only by (1,2,6,7,8,9,10,11), as directly recorded
in the independently reconstructed root geometry. This bounded observation
is not a uniform lower bound for arbitrary eight-speed cores.

## Lossless receipt transport

The full compact JSON is distributed as a gzip/base64 text archive to keep
repository API transfers small. From this directory, first run
`python3 -B ../materialize_receipts.py` to recover the exact named JSON files.
The decoder checks byte counts and SHA256. The proof still requires full
regeneration and independent replay; compression and hashes are not proofs.
The packet-level `verify.py` decodes into memory automatically and compares
the regenerated receipts byte for byte.
