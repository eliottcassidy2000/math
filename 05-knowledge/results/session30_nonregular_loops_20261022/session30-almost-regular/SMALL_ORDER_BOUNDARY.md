# Small-order failure boundary for the candidate truncations

FINITE-EXACT universe selection and H targets; floating-point engine predictions. Retrospective audit, not a holdout. All almost-regular isomorphism classes at n2,4,6 were tested.

Second order improves every one of the five n6 classes, but not the unique n2 or n4 classes. In particular:
- n2: exact y0.8, first0.291358025, second−0.645968602.
- n4: exact y0.959395656, first0.964974247, second1.029164505.

Hence neither positivity of the predicted y nor monotone improvement with truncation order holds at every finite order. These counterexamples do not refute an asymptotic expansion. At n6 the first-order predicted κ also ranks one H45 class below the H43 classes, while second order restores the ordering of distinct H levels in this complete universe. No general ranking guarantee follows; the independent n8 isospectral pair still has wrong ordering at second order.

Source data and reproduction: small_second_order.py, small_second_order.json, small_second_order.out. The script imports candidate predictions from the sibling decorated engine but computes exact H targets from the independently certified class packet.
