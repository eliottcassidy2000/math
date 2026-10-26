# Terras-shift-one merges of x = 2^(2L+1)+1 and y = 2^(2L+2)+1: first n with T^(n+1)(y) = T^n(x) before x reaches 1.
# Equivalent (THM-4581) to absorption of the pair chain (k,e) = (1,-1) for (T(y), x) = (3x-1, x).
import json, sys
from collections import defaultdict
def T(x): return x>>1 if x%2==0 else (3*x+1)>>1
def shift_one_merge(x, y):
    """Return (n, odd_steps_of_x, value) for the first n>=0 with T^(n+1)(y) == T^n(x) and T^n(x) > 2 ; None if x reaches 1 first."""
    u=T(y); v=x; n=0; odd=0
    while v>2:
        if u==v: return (n, odd, v)
        if v%2: odd+=1
        v=T(v); u=T(u); n+=1
    return None
Lmax=int(sys.argv[1]) if len(sys.argv)>1 else 1500
res={}
for L in range(1, Lmax+1):
    x=2**(2*L+1)+1; y=2**(2*L+2)+1
    r=shift_one_merge(x,y)
    res[L]=r
# (a) control: even L -> Syracuse index L+2 ; L = 3 mod 4 -> L+3
ctrl_even=all(res[L] is not None and res[L][1]==L+2 for L in range(2,Lmax+1,2))
ctrl_3=all(res[L] is not None and res[L][1]==L+3 for L in range(3,Lmax+1,4))
# (b) L = 1 mod 4: merge fraction and depths (Syracuse index = odd steps)
ones=[L for L in range(1,Lmax+1,4)]
merged={L:res[L] for L in ones if res[L] is not None}
frac_windows={}
for lo in range(1,Lmax+1,250):
    w=[L for L in ones if lo<=L<lo+250]
    frac_windows[f'{lo}-{lo+249}']=round(sum(1 for L in w if res[L] is not None)/len(w),3)
# (c) 2-adic periodicity test: if L = L' mod 2^s and the merge Terras time n satisfies n <= s+? then equal.
# The first N Terras steps of p=3^(L+1)+2 depend on L mod 2^(N-2) (order of 3 mod 2^N); merge at Terras time n after 2L+1
# means total Terras steps n-(2L+1) after the forced prefix. Test: group by L mod 2^s, compare merge times n' = n-(2L+1) when n' <= s-... (we test consistency for n' <= s+2 - 2 = s).
viol=0; tests=0
for s in range(2,12):
    groups=defaultdict(list)
    for L in ones:
        groups[L%(2**s)].append(L)
    for cls,Ls in groups.items():
        times=[(res[L][0]-(2*L+1)) if res[L] is not None else None for L in Ls]
        short=[t for t in times if t is not None and t<=s]
        if short:
            tests+=1
            # every L in the class must merge at the same n'
            if any(t!=short[0] for t in times): viol+=1
# (d) the merge tree: minimal residue classes (mod 2^s) forcing a merge at post-prefix Terras time n' <= s
tree={}
for s in range(2,12):
    groups=defaultdict(list)
    for L in ones: groups[L%(2**s)].append(L)
    for cls,Ls in sorted(groups.items()):
        times=set((res[L][0]-(2*L+1)) if res[L] is not None else None for L in Ls)
        if len(times)==1 and None not in times and list(times)[0]<=s:
            tree.setdefault(s,{})[cls]=list(times)[0]
summary={'Lmax':Lmax,'control_even_L_index_L+2':ctrl_even,'control_L=3mod4_index_L+3':ctrl_3,
 'L=1mod4_count':len(ones),'merged_before_1':len(merged),'merge_fraction_windows':frac_windows,
 'periodicity_tests':tests,'periodicity_violations':viol,
 'post_prefix_terras_merge_times_L=1mod4_first40':{L:(res[L][0]-(2*L+1), res[L][1]) if res[L] else None for L in ones[:40]},
 'forced_merge_classes_by_modulus':{str(2**s):{str(c):t for c,t in d.items()} for s,d in tree.items()}}
print(json.dumps(summary,indent=0))
