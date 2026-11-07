# Session 209: merge Sessions 14, 15 and 27 through leaf spectra

The inherited paired-black leaf classification now gives an exact joint
order/Hamilton-count law, with a precise distinction between class and
labelled measures. No all-tournament count spectrum or global minimum-order
formula is claimed.

**Main results**

- The Q-leaf count spectrum is exactly powers of3 and5 times even powers
  of3: H=3^k (k>=1) or H=5*9^k (k>=0).
- Their minimum orders among leaf realizations are3k+1 and6k+4,
  respectively. No15^k, k>=1, is a Q-leaf count, including Session208's
  strong order14 witness with H=15^6.
- The mixed S4 family is exponentially rare among leaves but maximizes H
  at every order n=0 or4 modulo6. Rarity does not remove extremizers.
- A uniform Q-leaf has typical log(H)/n about0.213410, whereas a uniform
  labelled tournament conditioned on its Q-image being a leaf has rate
  about0.146747. Exact algebraic constants, moment asymptotics and a normal
  limit are proved.

These concern the paired-black converse-merged MULTIGRAPH, with multiplicity
degree and loops counted twice. They do not concern leaves of the labelled
native graph or single-arc-flip graph. The two Q-leaves at order4 share one
edge; leaf counts are not automatically bridge counts.

**Reproduce**

    python3 verify.py
    python3 -O verify.py

Both runs regenerate the same RESULT.json and should match verify.out.
The verifier uses no third-party packages, performs literal Hamilton-path
orbit checks on all tournament classes represented by all labelled codes
through order5, and has explicit checks that survive optimization.

PROOF.md contains the complete derivations and scope. DEPENDENCIES.json
pins the primary Session14/15/27 source texts to a verified remote commit;
the all-order classification is inherited, not re-established by the small
census. The independent review and checker are included under audit/.

The strongest connection to Session27 is its automorphism-weighted SCC
formula, not a numerical analogy involving repunits or lattices. The bare
(n,H) spectrum forgets word order and labels. Session27's marked contexts,
and the companion one-port-blindness/two-port-decoding packet, address that
different information loss.
