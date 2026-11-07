# Numeric 14/15/27: an exact, degree-resolved comparison

See [proof and boundaries](PROOF.md). At order14, raw H15 has110 ordinary classes and55 converse-merged nodes; H27 has210 and105. Their exact fundamental paired-black degree histograms are respectively {3:5,5:41,15:9} and {1:28,9:65,27:12}. Thus H27 has28 genuine leaves, only4/15 of its stratum, while H15 has none. This is not a native-labelled degree statement or an adjacency classification.

The same standalone checker reconfirms m(15)=5, m(27)=6 and the inherited strong order14 witness with H=15^6, which supplies the7/3 construction-rate upper bound. No optimality is claimed.

Run `python3 verify.py` and `python3 -O verify.py`; both match RESULT.json. Inputs and source pins are bundled. Exact SCC consequences use an inherited complete order7 atlas and Busch's cited strong-tournament bound. The proof distinguishes these dependencies from newly recomputed controls and from the prior-session interpretation of14,15,27.
