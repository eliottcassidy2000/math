# Bounded exploratory probes, not all-order conclusions

`cut_probe.py` enumerates source words through length10 and records the largest ratio of best cut height squared to target minus one, for full and balanced seeds. Its observations led to the exact adjacent-cut theorem in ../PROOF.md; the script is not its proof.

`probe.py` checks a finite family of pairs of affine source templates in X=2^a. Each has a prefix and suffix with total length at most4, one O or Z run of length a+c with -3<=c<=2, and the appropriate t or b seed. It records no matching balanced polynomial 4X^2-4X+3. The fixed run length is always positive at a=4, avoiding the zero-length endpoint discontinuity. This finite failure is not an obstruction to arbitrary physical words or to other templates.

`pump.py` lists all balanced precursor words for227 and stretches existing O/Z runs with rates0..6 and total rate4..8. None supplies the checked sequence227,3971,65027,1046531. It does not insert initially empty runs, change the fixed J pattern, or establish any universal failure.

Run the scripts individually from any working directory. They import `probe.py` from their own directory. NumPy is required by the exploratory scripts only. The three text files contain their finite outputs.
