# Intrinsic generalized-pair structure: reviewed packet

A strong C2 graph belongs to the generalized pair family Gamma exactly when it has a start-two root whose deletion begins with a C1 strong component. PROOF.md proves the criterion and a canonical rooted N/M/Q inverse; structure.py returns the actual source and full vertex map. P is literally Q after a singleton J move, including its singleton-source boundary.

The packet also proves the exact rooted presentation generating function and gives explicit graph-level obstructions at H77 and H289. Neither orientation of the H289 example belongs to Gamma, yet another graph in Gamma has the same H and order. These are graph exclusions, not scalar-spectrum omissions.

Source replay, in a copy because it writes its receipt:

    python verify.py
    python -O verify.py

Both must match the included validation receipt. The independent check uses separate endpoint DP, Kosaraju SCCs and literal births on a complete865-graph C2 input slice:

    python audit/check.py
    python -O audit/check.py

The independent checker imports no author code. Its CSV input and matching receipts are included. RELEASE_FILES.sha256 is the unchanged reviewed source manifest; PUBLICATION_FILES.sha256 additionally binds the independent artifacts and this wrapper. No earlier theorem or navigation file is changed.
