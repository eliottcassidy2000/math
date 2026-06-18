"""
Final rigor check of the AP-family two-witness theorem and the S1 placement claim.
"""
from fractions import Fraction as F
from math import gcd

def nrm(x):
    r=x-int(x); r=r+1 if r<0 else r
    return r if r<=F(1,2) else 1-r

R=[]
def out(s): R.append(s); print(s)

# Witness 1 residue table: tau=2/27, ||2 r /27|| for r=0..26
out("WITNESS 1 (tau=2/27): residues r where ||2r/27|| < 2/27")
bad=[r for r in range(27) if nrm(F(2*r,27))<F(2,27)]
out(f"  unsafe residues: {bad}  (claim: {{0,13,14}})  match={set(bad)=={0,13,14}}")
out(f"  small-part level min_{{v=1..12}} ||2v/27|| = {min(nrm(F(2*v,27)) for v in range(1,13))} (claim 2/27)")

# Witness 2: tau=14k/(182k+1), level = 14k/(182k+1); no-wraparound 12*(m/13)<m+1
out("")
out("WITNESS 2 (tau=14k/(182k+1)=(m/13)/(m+1)): level and no-wraparound")
allok=True
for k in range(1,200):
    m=182*k; t=F(m//13, m+1)
    lvl=min(nrm(F(v)*t) for v in list(range(1,13))+[m])
    claim=F(14*k,182*k+1)
    if lvl!=claim: allok=False; out(f"  MISMATCH k={k}: {lvl} vs {claim}")
out(f"  level == 14k/(182k+1) for k=1..199: {allok}")
out(f"  no-wraparound 12*(m/13) < m+1: 12 m/13 = m - m/13 < m < m+1 always (algebraic). "
    f"check k=1: 12*{182//13}={12*14}=168 < {182+1}=183: {12*14<183}")
out(f"  min combined floor = min(2/27, 28/365) = {min(F(2,27),F(28,365))} = {float(min(F(2,27),F(28,365))):.5f}")
out(f"  2/27 - 1/14 = {F(2,27)-F(1,14)} > 0 ; 28/365 - 1/14 = {F(28,365)-F(1,14)} > 0")

# S1 placement: {1..12,m} has exactly one speed >13
out("")
out("S1 PLACEMENT: {1..12,m}, m>=14 has k=#{v>13}=1 => case S1 (already PROVED).")
out("  => AP theorem rigorously re-closes an S1 slice; does NOT touch open S3.")

with open(r"C:\Users\Eliott\Documents\GitHub\math\05-knowledge\results\lrc14_ap_witness_check_kps-S4-wf.out","w") as f:
    f.write("\n".join(R))
out("\n[written]")
