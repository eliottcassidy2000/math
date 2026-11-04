# Independent review of the two-branch all-clock prototype

## Verdict

ACCEPT the all-clock affine lifting theorem and the literal paired two-birth witnesses. The construction is an infinite family, not a recursively closed grammar, and it does not return the actual 170807 suffix. The stated caveats are mathematically necessary and correctly retained.

The finite bank's 180 marked witnesses were all recomputed independently by permutation enumeration and exact cut positions. The 6 affine +2 and 7 affine +4 bank matches were independently reconstructed. These checks validate the witnesses, not by themselves the completeness of the bank census; complete labelled enumeration was also rerun in normal and optimized Python. Both census stdout and the complete JSON bank byte-match. The corrected verifier normal and optimized stdout also byte-match.

## Fixed-name control not supplied by response equality

The original verifier's two calls to `attach` used different label intervals for prefix/root when suffix size changed. That is a harmless relabelling for counting, but is not itself a test that actual births preserve old physical names. The independent checker instead starts with the labelled full source, fixes root index s+4, and appends x and y at indices s+5 and s+6. For s=1,2,3,4 and both branches it verifies:

- every old edge and name survives;
- the physical root still has the sole in-neighbor 0;
- the source, intermediate and target are strong;
- the source and target root satisfy E_root(2)=H;
- the counts satisfy the claimed 3H+d identity;
- direct permutation enumeration independently agrees for both full target tournaments at s=1,2.

The all-s proof follows by splitting covers over the uniform prefix cut, not by fitting these four tests.

## Exact recurrence boundary

The independent checker reapplies the proposed recipe at every in-neighbor of port0 in the actual displayed C4=[56,1,19,22,2,30]. Every candidate fails the required marked response pair (3*19+d,3*21). This proves only the stated failure of that recipe on that marked suffix, not a general impossibility of iterated returns.

## Verifier correction

The initial `verify.py` used Python assertions, making optimized replay vacuous. The author replaced them by explicit runtime checks before the audited replay. The independent `audit_phases.py` uses explicit exceptions throughout; its normal and optimized outputs byte-match.

## Artifacts

- `audit_phases.py`: independent permutation and fixed-name birth controls.
- `phase_audit.json`: canonical independent receipt.
- The authoritative full census and verifier are in the sibling `phases/` directory. Normal and optimized replay outputs matched; duplicate stdout snapshots are omitted from publication.

No finite-field or 81/27 observer identity follows from this prefix/suffix audit. No phase equality, unrooted isomorphism, or scalar count was used to replace a physical source.
