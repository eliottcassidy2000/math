# Independent audit: complete one-variable-run inverse

**Independent analytic/source/replay audit: PASS.** No mathematical correction
is required. This audit accepts [FAMILY_INVERSE.md](FAMILY_INVERSE.md) only in
its stated 34-entry, one-variable-run language. It does not turn the absence
of these templates into absence of adaptive preparations or spectrum values.

## Proof and completeness

I read the complete proof, `family.py`, and its literal dependency
`surgery_ready.py`, and recovered the inherited two-entry funnel formula in
`two_entry_funnels_20261103/TRANSPORT.md`. The actual bank consists of 34 named
core/port alternatives, not 34 isomorphism classes of every ready source.

The five-state recurrence retains precisely the observations used by the
ready-count formula. In particular its initial eta=1 is necessary for the
first1, and D for a nonempty word is determined by its first bit, not its final
bit. Directly, for a word w of length n,

    T(w)=2^n−sum_(j<n) w_j 2^(n−1−j),
    L(w)=1+sum_(j<n) w_j 2^j,
    B(w)=1+sum_(i:w_i=1) (T(w_<i)−eta(w_<i)).

These formulas independently reproduce the recurrence. Iterating either bit
k>=1 times gives exactly the affine states displayed in the primary; there
is no fit from a handful of sample lengths.

For p a^k q, the coefficient of 2^k is

    v(q)(c2^|p|+hT_p)                  if a=0,
    (1+v(q))(c2^|p|+h(T_p−1))          if a=1.

I independently obtained this by summing the suffix1 contributions to B and
L. The bank has c=G+p2>0 and h>0, while T_p>=1. Thus the stated lower bound
c2^|p|<=A is valid for every candidate and supplies a finite **proved** prefix
limit. The exact quotient determines v(q); the unique little-endian expansion
ending1 then determines q. Empty q occurs only in the permitted one-run case.
In particular, trailing or leading zeros are not lost by an unjustified word
quotient: trailing zeros would violate the required final1, and leading zeros
are retained by the little-endian expansion. All p are literally enumerated.

The intercept includes the actual initial-bit D term and all c,h,e,K data.
Hence the algorithm is complete for every finite p,q and a fixed run bit,
with the same unbounded variable k. Several presentations of one family are
allowed; none is silently identified with a physical source.

## Independent controls and replay

In addition to full producer replay, I ran a separate verifier which
read only the literal BANK assignment by AST. It did not import either producer
module. It enumerated all simple directed paths of each core by permutations
and reconstructed all seven core counts, including the ordered disjoint joint
term. It verified each declared gap b. It then used the closed word sums above
and an independently written slope-factor inverse; the constant term was
computed as H(p a q)−2A and checked at run lengths 2,4,7.

This independent path passed **6,687 always-active controls**. For the requested
slopes 756,2268,6804,20412, the per-query `(prefixes,divisor candidates)` are
respectively `(1014,65)`, `(3230,107)`, `(10910,143)`, `(30686,167)`. All fifteen
specified target intercepts have no presentation. The positive query (90,5)
has `(94,20)` and nine actual presentations, containing both gaps1 and2. These
agree exactly with the frozen receipt, including the individual positive
prefix/run/suffix records. The verifier is preserved as
[`family_audit.py`](../../../../04-computation/nonzero_preparation_20261103/family_audit.py),
with output [FAMILY_AUDIT.out](FAMILY_AUDIT.out). Its only adjustment on saving
was replacing the temporary worktree's absolute base by its own source
directory; the output bytes were checked unchanged.

```text
python3 -B 04-computation/nonzero_preparation_20261103/family_audit.py
```

Independently run normal and optimized producer commands were

```text
python3 -B 04-computation/nonzero_preparation_20261103/family.py --output /tmp/nonzero_family_audit_normal.json
python3 -B -O 04-computation/nonzero_preparation_20261103/family.py --output /tmp/nonzero_family_audit_optimized.json
```

Both pass **604,826 always-active gates**. Both resulting files are byte-for-byte
equal to the 5,502-byte frozen `family.json`.

| Frozen input | Bytes | SHA256 |
|---|---:|---|
|family.py|5511|`4e364088a89698dfd6e51e6b86682f9b295997762a28c425877f5c65bffb5d3a`|
|family.json|5502|`7bcff67f4d7c675c8b7f94b009892c45ec41bc7b6de0a524e31158b2c44798a4`|
|surgery_ready.py|8410|`024ac0b483d832865c6f12eabf0fe0f8676ac39188631e44910a0a4c9e66ca2f`|

The pre-promotion primary read for this audit had SHA256
`85fab314c6102c9f9323b9596452b6054b2b56f88935a820e0f1c778d73afb40`.
A status/link-only promotion does not change the audited mathematical scope.

## Boundary of acceptance

The returned targets are actual ready-funnel templates, but the inverse is a
scalar membership decision in this precise language. It does not supply an
old-name-preserving switch from the current graph to an equal-count template,
nor verify a joint retained-root response after such a switch. The primary
keeps those tasks separate. Multiple changing runs, other finite or growing
cores, and a word selected separately at each clock remain outside this
negative result. Those limitations are load-bearing and accepted.
