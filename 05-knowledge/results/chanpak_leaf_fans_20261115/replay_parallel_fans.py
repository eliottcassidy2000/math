"""Certificates for two unbounded finite-width windows, without the search."""
from replay_fan_channels import leaf, profiles

# (offset, v2(S-K), port, chronological parent word)
BANK = {
    1: [
        (1,0,1,"0"), (3,5,7,"00000404"), (5,8,8,"00000051"),
        (7,5,0,"00000045"), (9,4,8,"00000046"), (11,8,7,"00104050"),
        (13,7,7,"00002456"), (15,7,7,"00000202"), (17,7,6,"00104565"),
        (19,9,8,"00003051"), (21,7,6,"00002064"), (23,8,7,"00011155"),
        (25,8,7,"00000144"),
    ],
    17: [
        (23,5,8,"00000351"), (25,1,8,"00003356"), (27,2,8,"00000205"),
        (29,1,8,"00002067"), (31,1,8,"00004307"), (33,1,8,"00000020"),
        (35,3,8,"00022560"), (37,2,8,"00000435"), (39,4,7,"00011165"),
        (41,3,8,"00000254"), (43,3,8,"00023364"),
    ],
}

for slope, seeds in BANK.items():
    lo, hi = seeds[0][0], seeds[-1][0]
    if [c for c, _, _, _ in seeds] != list(range(lo,hi+1,2)):
        raise AssertionError(seeds)
    for c,e,p,word in seeds:
        R=(0,)
        for parent in word:
            R=leaf(R,int(parent))
        H, data=profiles(R)
        S,d,K,delta=data[p]
        H1, data1=profiles(leaf(R,p))
        eta=delta+2*(d-K)
        if (S-K,eta,H1,data1[p][3]) != (slope<<e,-c,(slope<<(e+1))+c,-c):
            raise AssertionError((slope,c,e,p,word,H,d,K,H1,data1[p]))
    print(f"Slope {slope}, odd offsets {lo}..{hi}, all heights j >= {1+max(e for _,e,_,_ in seeds)}")
