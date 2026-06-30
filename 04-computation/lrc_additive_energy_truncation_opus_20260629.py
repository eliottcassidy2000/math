"""
|T|=3 additive-energy perturbation of the cyclotomic SOS. Inclusion-exclusion:
L(S)=Sum_T (-1)^|T| meas(cap_T D_i); truncation L^{<=k}=Sum_{j<=k}(-1)^j E_j, E_j=E_tau[C(X,j)] (factorial moments).
E_2,E_3 = the pair/triple additive ENERGY (set-dependent). |T|<=2 = the ~floor; E_3 = the perturbation.
TEST: does the AP maximize the energy E_3 (=> minimize L)?  Also true lonely measure L=P(X=0).
"""
from math import comb
def Xdist(S,Q,n=14):
    # danger count distribution X(tau)=#{i: ||v_i tau||<1/14}
    from collections import Counter
    c=Counter()
    for m in range(Q):
        x=0
        for v in S:
            r=(v*m)%Q
            if min(r,Q-r)*n < Q: x+=1   # ||v m/Q|| < 1/14
        c[x]+=1
    return c,Q
def moments(S,Q):
    c,Q=Xdist(S,Q)
    E=[sum(cnt*comb(x,j) for x,cnt in c.items())/Q for j in range(5)]
    L=c.get(0,0)/Q
    Lk=[sum((-1)**j*E[j] for j in range(k+1)) for k in range(5)]
    return E,L,Lk
Q=100800
sets={
 "AP {1..13} (M=1/14, tight)": list(range(1,14)),
 "{1..11,13,84} (covering, 7/89)": [1,2,3,4,5,6,7,8,9,10,11,13,84],
 "{2..14} (covering, loose 1/8)": list(range(2,15)),
 "{1,2,4,..,4096} (lacunary/Sidon)": [2**i for i in range(13)],
}
print(f"factorial moments E_j=E[C(X,j)] (=|T|=j additive energy), truncations L^<=k, true L=P(X=0). Q={Q}")
print(f"{'set':>34} {'E2(pair)':>9} {'E3(trip)':>9} {'E4':>8} {'L<=2':>8} {'L<=3':>8} {'L(true)':>9}")
for nm,S in sets.items():
    E,L,Lk=moments(S,Q)
    print(f"{nm:>34} {E[2]:>9.4f} {E[3]:>9.4f} {E[4]:>8.4f} {Lk[2]:>8.4f} {Lk[3]:>8.4f} {L:>9.5f}")
print()
print("E_1 = 13/7 = %.4f (set-independent, the union bound). E_2,E_3 set-DEPENDENT = the additive energy."%(13/7))
print("TEST: AP maximizes E_3 (most triple resonances) => drives L toward 0 (AP is the L-MINIMIZER, razor-thin).")
print("Covering sets have LESS energy => L>0 (the floor). Lacunary (Sidon, ~no relations) => largest L.")
