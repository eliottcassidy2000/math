# Independent audit: primitive branching and representation rebasing

Verdict: PASS on final frozen bytes.

PROOF.md SHA256: 458eb8c9ced79f763eeb9f9f52c4bb50b1a686d565a1eac5ada3de735d7d9001
verify.py SHA256: 0ba80ccc1304eb718865c005c7a80e2f45cf093db721a89c4a51add220362b09

## Density and the finite-menu theorem

All intermediate profiles in the explicit O^i J O^j Z J O^k family check.
The output has odd D and U=q(D+1), with q a power of two, so gcd(D,U)=1.
The displayed rational inverse recovers z,y,q on a nonempty open subset of
E=0, establishing injectivity of polynomial substitution. Independent infinite
parameter sets therefore give the required polynomial density on that plane.

Each fixed homogeneous source-extension rule minus its proposed affine child
has nonzero constant-c after restriction to the two E=0 planes. Its product
with finitely many other such nonzero polynomials is nonzero. Density supplies
an actual primitive source pair escaping every branch. Each pullback has degree
at mostone separately in each of six parameters, so the m+1 interpolation grid
in each parameter is sufficient and yields source order at most3m+7. The
argument does not claim arbitrary integer plane points are physical, and it
requires the stated fixed homogeneous rules and nonzero offsets.

## The71 all-word branches and physical repair

The positive-matrix inequalities correctly exclude every omitted word: N has
at mosttwo total letters; M has no extension; A/B/C/T already exceed2n+1;
J0/J1/J2/J3 can use only J, with at mosttwo per side. For n>=13 all numerical
strict inequalities hold. Source swaps are absorbed by transpose closure,
now explicitly stated. Singleton boundaries and freely rebased sources are
not extensions of the supplied nontrivial pair.

A separately authored checker evaluates the71 branches in original(D,E,U,d)
coordinates using the inherited scalar laws, reproducing exactly45 values
and the absence of81,83 from the H41 pair. It independently counts allfive
literal parent/rebase/child graphs by full permutation enumeration, confirming
three H41,C2 parents, H81,C2, and H83,C3. See
independent_branch_check.py and INDEPENDENT_BRANCH_AUDIT.json.

The last-letter proof of the twin equality fiber is valid. In particular,
e>0 sources satisfy f<=a because their last letter is Z or J; this is enough
to exclude the additional J branch. The primitive zero cases are precisely
t and the triangle. No general rebasing algorithm is inferred from them.

## Replay and scope

Final normal and optimized isolated replays agree exactly, and every release
hash verifies afterward. Controls include125 primitive parameter choices,
three literal primitive sources,71 branches/65 matrices,five literal rebasing
graphs, and the twin equality check on32,548 physical source profiles.

The theorem is representation-local: it does not exclude choosing a different
presentation for each target. The41→81/83 repair demonstrates that distinction.
No scalar nonrepresentability or all-odd coverage conclusion is warranted.

Final manifest: 170bba2ef38f4549fab7356ffcd0c966258f83899fcb34a47bdf75acd4f12ff0.
The final uncontaminated replay is in ../primitive-final-replay/; reviewer-added
scripts are kept outside that copy because the upstream verifier hashes every
Python/Markdown file beneath its own packet directory.
