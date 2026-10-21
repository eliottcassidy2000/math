# Fixed symmetry slack: canonical word-module rationality

Main deliverable: SYMMETRY_SLACK_RATIONALITY.md.

Status: PROVED in the stated scopes, with independent adversarial audit PASS and FINITE-EXACT controls. Its hypotheses use the proved bounded exceptional-symmetry theorem in continuation_width_symmetry_20261020/README.md. No Git changes have been made by this worker.

Main result: for each fixed d=3f(T)−s(Aut(T)), the ordinary vertex/feedback bivariate series is a finite positive sum z^a y^b L_d(z,y)^ell, ell≤2d+1, where L_d is the tail of1/(1−z−yz³) above vertex degree d. The self-converse channel uses paired words and palindromes and has denominator dividing(1−z²−y²z⁶)^(d+1). A finite exact compiler is specified; it has not been run at every d, and no efficient runtime is asserted.

Reproduce:

python3 verify_symmetry_slack.py 7
python3 -O verify_symmetry_slack.py 7

The order-seven transcripts byte-match and pass21,749 exact checks on all532 nonempty unlabelled tournaments through order seven, plus long-word replacement and synchronized-group hostiles. The smaller default order-six outputs are retained separately. Computations are finite controls, not the all-order proof.

Inherited contrast: fixed-feedback rationality was rediscovered, checked against current references, and explicitly not claimed as new. The genuinely new proved step closes the previously open all-fixed-symmetry-slack rationality question by canonicalizing the word modules and proving all long slots are fixed individually by the actual stabilizer. It retains cyclic quotient synchronization rather than replacing it by independent triangle rotations.

Independent review: session23-slack-review/AUDIT.md. The separate checker reports4,576 exact assertions on139 actual tournaments and548 proposed partitions, including nonmaximal partitions and long/converse hostiles through order15. This is separate from the exhaustive order-seven producer controls.
