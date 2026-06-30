"""
Prediction 2: the AP (staircase) maximizes higher-arity additive resonances Sum t_v v = 0
(the SL(k)/zeta(k) richness). Count k-term signed (+-1) relations: subsets T,|T|=k with a
sign choice making sum +-v = 0. Arity k=2,3,4 = the SL(2),SL(3),SL(4) levels.
"""
from itertools import combinations, product
def count_pm1_relations(S, k):
    cnt=0
    for T in combinations(sorted(set(S)), k):
        found=False
        # fix first sign +, search others (avoid double count of global sign)
        for signs in product([1,-1], repeat=k-1):
            if T[0] + sum(s*v for s,v in zip(signs,T[1:]))==0:
                found=True; break
        if found: cnt+=1
    return cnt
sets={
 "AP {1..13} (staircase)": list(range(1,14)),
 "covering {1..11,13,84}": [1,2,3,4,5,6,7,8,9,10,11,13,84],
 "{2..14} covering": list(range(2,15)),
 "lacunary {1,2,4,..,4096}": [2**i for i in range(13)],
}
print("k-term signed (+-1) additive relations  (subsets T,|T|=k with +-v summing to 0):")
print(f"{'set':>30} {'k=2 (SL2/ze2)':>14} {'k=3 (SL3/ze3,Littlewood)':>25} {'k=4 (SL4/ze4,tight cap)':>24}")
for nm,S in sets.items():
    r2=count_pm1_relations(S,2); r3=count_pm1_relations(S,3); r4=count_pm1_relations(S,4)
    print(f"{nm:>30} {r2:>14} {r3:>25} {r4:>24}")
print()
print("Prediction 2: the AP/staircase should be RICHEST in higher-arity relations (the additive-relation")
print("-richest direction). k=4 = the SL(4)/degree-4 tight-cap level. If AP tops k=4, the tight cap's")
print("consec-extremality (HYP-2823) = AP maximizing quadruple resonances = the SL(4) Rogers term.")
