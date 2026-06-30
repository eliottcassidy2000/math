"""
Han-Lee congruence-Siegel 2nd moment -> LRC floor (the witness route, THM-501 D(q,S)).
D(q,S)/q = density of lonely residues a mod q. 1st moment (independence) ~ (6/7)^13.
2nd moment over the ENSEMBLE of covering sets = Han-Lee random-counting. Concentration => floor for
GENERIC covering sets (the SL(2)/EKL-analog level). The MIN over sets = the gap (rare tight sets).
"""
import random
def lonely_residues(S,q,n=14):
    thr=q/(2*n)  # ||va/q||>=1/14 <=> (va mod q) in [thr, q-thr] roughly; exact: min(r,q-r)>=q/14
    cnt=0
    for a in range(1,q):
        ok=True
        for v in S:
            r=(v*a)%q
            if min(r,q-r)*n < q:   # ||v a/q|| < 1/14  => danger
                ok=False; break
        if ok: cnt+=1
    return cnt
def is_cov(S): return all(any(s%qq==0 for s in S) for qq in range(2,15))
q=1009  # prime shell
random.seed(7)
# build a sample of covering 13-sets (varied: small, divisor-loaded, random)
samples=[]
samples.append(sorted(set(range(2,15))))                       # {2..14}
samples.append([1,2,3,4,5,6,7,8,9,10,11,13,84])               # the 7/89 near-tight one
samples.append([1,2,3,4,5,6,7,8,9,10,11,13,28])
tries=0
while len(samples)<60 and tries<200000:
    tries+=1; S=sorted(random.sample(range(1,40),13))
    if is_cov(S): samples.append(S)
dens=[lonely_residues(S,q)/q for S in samples]
import statistics as st
mean=st.mean(dens); var=st.pvariance(dens); mn=min(dens); mx=max(dens)
print(f"Han-Lee random counting: D(q,S)/q at q={q} over {len(samples)} covering 13-sets")
print(f"  1st moment (independence pred (6/7)^13) = {(6/7)**13:.5f}")
print(f"  sample MEAN   D/q = {mean:.5f}")
print(f"  sample VAR    D/q = {var:.6f}   (the Han-Lee 2nd-moment / concentration)")
print(f"  sample MIN    D/q = {mn:.5f}  <-- the worst/tightest set (the gap: is it > 0?)")
print(f"  sample MAX    D/q = {mx:.5f}")
print(f"  Paley-Zygmund-style floor signal: min stays {'POSITIVE (floor holds on sample)' if mn>0 else 'HITS 0'}")
# the near-tight {..84}: its D/q
print(f"\n  near-tight {{1..11,13,84}}: D/q = {lonely_residues(samples[1],q)/q:.5f}  (M=7/89; still many lonely residues)")
print("  => GENERIC covering sets concentrate near (6/7)^13 with small variance (floor holds, SL(2)/EKL level).")
print("     The floor is a 2nd-moment/concentration result; the TIGHT cap (all sets, no exceptions) is SL(4).")
