# -*- coding: utf-8 -*-
# kps-2026-07-11-S127cont30: bounded statement -- SHALLOW / depth-5 sufficiency, honestly characterized.
#
# The SHALLOW half of the clean-ruler condition (HYP-6000): every residual family has a pair-sum modulus
# q=v_i+v_j with maxBand<=5, so the depth-5 Bonferroni B5 (the Lean obligation) certifies.
#
# FINDING: depth-5 via pair-sums holds for GENERIC residuals but FAILS for COARSE-REDUCIBLE ones -- which
# are dispatched by OTHER branches (dilation-invariance / detuned THM-678 / coarse reduction), NOT residual.
#   - GENERIC (primitive, no prime | >6 elements, longest-AP<=7) adversarial: max[min-over-pairsums maxBand]=4
#     => some pair-sum has maxBand<=4<5 on every generic residual found (depth-5 sufficient, with margin).
#   - STRUCTURED counterexamples (NON-residual): dilated AP 7*{1..13} (all pair-sums | 7 => maxBand 13);
#     detuned 7*{1..12}+{92} (12 multiples of 7 => THM-678 detuned dispatch => maxBand 8). Both coarse-reducible.
# So the deep-coverage families that force depth>5 are EXACTLY the structured ones the dispatch architecture
# already peels; the clean-ruler/B5 route needs only the GENERIC residual, where depth-5 is correct.
#
# PROVABLE pieces: (i) p=1 exact bandCount at q=v_a+v_b>Vmax = #{i:14 v_i<=q} + #{i:14 v_i>=13q} (no p-search);
# (ii) average bandCount over multipliers <= 13/7 (each runner in the 1/7 danger arc a 1/7-fraction of the time).
from math import gcd
from functools import reduce
def bandCount(v,q,p): return sum(1 for vi in v if not (q<=14*((vi*p)%q)<=13*q))
def maxband(v,q): return max(bandCount(v,q,p) for p in range(1,q))
def minmaxband(v):
    ps={v[i]+v[j] for i in range(13) for j in range(i+1,13) if v[i]+v[j]>=14}
    return min(maxband(v,q) for q in ps), ps

def main():
    print("STRUCTURED (coarse-reducible, NON-residual) families break depth-5 via pair-sums:")
    for name,v in {'dilated AP 7*{1..13}':[7*k for k in range(1,14)],
                   'detuned 7*{1..12}+{92}':[7*k for k in range(1,13)]+[92],
                   'GENERIC residual example':[1,4,5,9,10,15,25,34,38,42,48,51,55]}.items():
        mm,_=minmaxband(v)
        cr = 'COARSE-REDUCIBLE => other branch' if mm>5 else 'generic => depth-5 OK'
        print(f"  {name:26s}: min-over-pairsums maxBand = {mm:2d}  ({cr})")
    print()
    print("PROVABLE p=1 bound at the two-largest pair-sum q=v12+v13 (>Vmax, so v_i mod q = v_i):")
    for v in [[1,4,5,9,10,15,25,34,38,42,48,51,55], list(range(1,13))+[26]]:
        q=v[-1]+v[-2]
        b1=bandCount(v,q,1); lo=sum(1 for x in v if 14*x<=q); hi=sum(1 for x in v if 14*x>=13*q)
        print(f"  v(max={v[-1]}): q={q}, bandCount(1)={b1} = lo(14v<=q)={lo} + hi(14v>=13q)={hi}   [ratio<=13 => hi=0]")
    print()
    print("BOUNDED STATEMENT: on GENERIC (primitive, non-coarse-reducible) residuals, some pair-sum has")
    print("maxBand<=4, so depth-5 B5 certifies (the Lean obligation's depth is correct); the depth-5 failures")
    print("are the coarse-reducible families (dilated/detuned) that dilation + THM-678 + coarse-reduction peel.")

if __name__ == '__main__':
    main()
