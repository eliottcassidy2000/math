import random, sys
from fractions import Fraction
exec(open('04-computation/lrc14_adversarial_verify_eliott.py').read().split('# ---- DECISIVE')[0].split('if __name__')[0])
random.seed(7); NR=20000; rbest=None; rcex=0; n=0
for _ in range(NR):
    miss=random.choice([1,2,4,8])
    pool=[x for x in range(1,50) if x!=miss]
    C=tuple(sorted(random.sample(pool,12)))
    if has_tower(C) or not primitive(C): continue
    n+=1
    L=lonely_measure(C)
    if rbest is None or L<rbest[0]: rbest=(L,C)
    if L<THR2: rcex+=1; print("RAND CEX",L,C)
print(f"random {n} valid tower-broken cores in [1,49]: min meas={rbest[0]}={float(rbest[0]):.6f}")
print(f"sub-THR2 counterexamples: {rcex}")
