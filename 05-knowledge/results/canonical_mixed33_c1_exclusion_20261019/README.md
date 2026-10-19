# Canonical mixed (3,3,1) exclusion

**Transport setup:** The oversized gzip certificate is stored losslessly as raw-byte parts. Before running the audit or replay commands below, run `python3 reconstruct_transport.py` from this directory. See [the transport instructions](TRANSPORT_README.md) for original hashes, part verification, and safe reconstruction.

The theorem concerns roots selected by maximum strict-majority span and then maximum common-incomparable count. It does not claim an any-root exclusion or the full seven-thin conjecture.

Read `canonical_mixed33_c1_exclusion_20261019.md` for the complete argument, external six-thin dependency, and the two relaxation obstacles that motivated the refinements.

## Replay

Only Python's standard library is needed for all verification commands. The numerical search files require NumPy/SciPy, but they are not needed to verify the result.

From this directory:

```
python audit_terminal.py
python canonical331-cap-continuation-stage/replay_complete.py
python canonical331-continuation-stage/test_replay_rejections.py
python canonical331-cap-continuation-stage/test_cap_replay.py
python canonical331-stage/verify_simple_obstacle.py
python canonical331-continuation-stage/verify_simple_obstacle.py
```

The complete replay independently generates the core catalogue, checks every analytic exclusion, reconstructs every extension and boundary generator, and verifies the entire three-phase certificate union. It invokes each phase's original verifier. The directory names are retained verbatim to keep source locks and relative dependencies unchanged.

The first two progress records intentionally say `EXACT_CONE_OBSTACLE`: those phases stopped at documented relaxation failures. Their committed prefixes are reused without alteration. The third phase finishes all remaining cores. Only the final three-phase union is the complete certificate.

`first_obstacle.json` files are relaxation witnesses, not actual poset counterexamples. The first is excluded by the shared slack-one split; the second by quantitative boundary ancestor caps. The proof cites the separate general lemmas.

Exact integer weights can exceed machine-integer range. Replay uses arbitrary-precision Python integers throughout. Corrupted-record tests check omitted split children, wrong degree/slack, incorrect masks, negative weights, false margins, wrong coefficient convention, duplicate caps and wrong ancestors.

No minimum-support or novelty claim is made. General comparison, exchange and interleaving methods are attributed in the proof.
