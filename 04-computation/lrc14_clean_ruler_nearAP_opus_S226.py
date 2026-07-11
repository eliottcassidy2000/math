"""
opus-2026-07-11-S226: systematic clean-ruler certificates for the near-AP residual class (feeds kps THM-707 hB5).

kps THM-707: hB5 <= "every residual covering family v (13 speeds) has a CLEAN RULER q:
   liveCount(q) >= 1  AND  maxBand(v,q) <= 5",  where
   bandCount(v,q,p) = #{i : (v_i p mod q) NOT in [ceil(q/14), floor(13q/14)]}   (runners in the 1/14 danger arc),
   liveCount(q) = #{p in [1,q-1] : bandCount = 0},  maxBand = max_p bandCount.
Then B5 = liveCount - Penalty, Penalty=0 when maxBand<=5, so B5 = liveCount > 0 => live => lonely.

kps discharged the binding {1..12,26} via the pair-sum ruler q=27. This script tests SYSTEMATICALLY:
does every near-AP 13-speed family have a clean ruler, and is a PAIR-SUM ruler q=v_a+v_b always clean?
"""
def band_stats(v, q):
    lo=-(-q//14); hi=(13*q)//14   # ceil(q/14), floor(13q/14): safe band [lo,hi]
    live=0; maxband=0
    for p in range(1,q):
        bc=sum(1 for vi in v if not (lo <= (vi*p) % q <= hi))
        if bc==0: live+=1
        if bc>maxband: maxband=bc
    return live, maxband
def clean_rulers(v, qmax=60):
    out=[]
    for q in range(8, qmax+1):
        live,mb=band_stats(v,q)
        if live>=1 and mb<=5: out.append((q,live,mb))
    return out
def pair_sums(v): return sorted(set(a+b for i,a in enumerate(v) for b in v[i+1:]))

fams = {
  "AP {1..13}          ": list(range(1,14)),
  "AP dilated 2*       ": [2*i for i in range(1,14)],
  "AP dilated 3*       ": [3*i for i in range(1,14)],
  "near-AP {1..12,14}  ": list(range(1,13))+[14],
  "near-AP {1..12,26}  ": list(range(1,13))+[26],   # kps binding case, q=27
  "near-AP {1..12,20}  ": list(range(1,13))+[20],
  "near-AP {2..14}     ": list(range(2,15)),
  "perturbed {1..11,13,14}":[1,2,3,4,5,6,7,8,9,10,11,13,14],
  "2-block {1..7,20..25}":[1,2,3,4,5,6,7,20,21,22,23,24,25],
}
print(f"{'family':>24} {'#clean rulers<=60':>18} {'smallest clean (q,live,maxBand)':>34} {'pair-sum clean?':>16}")
for name,v in fams.items():
    cr=clean_rulers(v)
    ps=pair_sums(v)
    ps_clean=[(q,)+band_stats(v,q)[:2] for q in ps if q<=200 and band_stats(v,q)[0]>=1 and band_stats(v,q)[1]<=5]
    first=cr[0] if cr else None
    pscl = f"YES q={ps_clean[0][0]}" if ps_clean else "none<=200"
    print(f"{name:>24} {len(cr):>18} {str(first):>34} {pscl:>16}")

print("\n=== detail on the binding {1..12,26} (kps q=27) ===")
v=list(range(1,13))+[26]
for q in [26,27,28,39,52]:
    live,mb=band_stats(v,q); tag="CLEAN" if live>=1 and mb<=5 else ""
    print(f"  q={q}: liveCount={live}, maxBand={mb}  {tag}")
print("\n=== does the AP {1..13} itself (the extremal) have a clean ruler? ===")
v=list(range(1,14)); cr=clean_rulers(v,80)
print(f"  {v}: clean rulers <=80 = {cr[:6]}{'...' if len(cr)>6 else ''}  (total {len(cr)})")
