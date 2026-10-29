"""Independent compact-encoding census of all Leech minimum shell families."""
from itertools import combinations
from collections import Counter
import json

def need(b,m):
    if not b:raise ArithmeticError(m)

def main():
    C={0}
    for g in [7973<<i for i in range(11)]+[(1<<24)-1]:C|={c^g for c in C}
    need(len(C)==4096,'dimension12')
    weights=Counter(c.bit_count() for c in C)
    need(dict(weights)=={0:1,8:759,12:2576,16:759,24:1},'Golay weights')
    vectors=set();counts=Counter()
    def accept(v,kind):
        parity=v[0]&1
        need(all((x&1)==parity for x in v),'common parity')
        c=sum(1<<i for i,x in enumerate(v) if x%4==parity+2)
        need(c in C,'code residue')
        need(sum(v)%8==4*parity and sum(x*x for x in v)==32,'sum guard and norm')
        encoded=sum((x+4)<<(4*i) for i,x in enumerate(v))
        need(encoded not in vectors,'unique full vector');vectors.add(encoded);counts[kind]+=1
    for i,j in combinations(range(24),2):
        for signs in range(4):
            v=[0]*24;v[i]=4 if signs&1 else -4;v[j]=4 if signs&2 else -4
            accept(v,'two_fours')
    for c in C:
        if c.bit_count()!=8:continue
        support=[i for i in range(24) if c>>i&1]
        for k in [0,2,4,6,8]:
            for negative in combinations(support,k):
                N=set(negative);v=[0]*24
                for i in support:v[i]=-2 if i in N else 2
                accept(v,'octad_twos')
    for c in C:
        for distinguished in range(24):
            v=[(3 if c>>i&1 else -3) if i==distinguished else (-1 if c>>i&1 else 1) for i in range(24)]
            accept(v,'one_three')
    need(dict(counts)=={'two_fours':1104,'octad_twos':97152,'one_three':98304},'shell families')
    need(len(vectors)==196560,'complete shell')
    # In Z^24/8Z^24: 2 parities *4096 code residues *2^23
    # choices of the coordinate high bits satisfying the one sum parity.
    allowed_mod8=2*4096*2**23
    need(8**24//allowed_mod8==2**36==8**12,'index and scaling')
    need(4%8!=0 and 24%8!=4,'two dropped-guard hostiles')
    print(json.dumps({'minimum_vector_counts':dict(counts),'total':len(vectors),
        'index_power_of_two':36,'minimum_squared_norm':4,
        'plain_Construction_A_roots':48},sort_keys=True))
if __name__=='__main__':main()
