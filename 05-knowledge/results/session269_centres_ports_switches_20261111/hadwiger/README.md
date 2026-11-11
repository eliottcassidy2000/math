# Hadwiger–Nelson: exact core suspension and connector support

Read PROOF.md for the induced Euclidean construction and its full scope. `replay.py` reproduces the 14-point interacting square-core graph and its exact response. `equality.py` reproduces a 7-point connector-sensitive 3-colour equality gadget at squared port distance 1/18. Results are receipt.json and equality_receipt.json. No new plane chromatic bound is claimed.

All correctness checks use explicit require/raise calls and also run under python3 -O. Both normal and optimized runs passed with byte-identical receipts and stdout; see OPTIMIZATION_CHECKS.json. The main receipt records the balanced-triangle Gram control (24,24).
