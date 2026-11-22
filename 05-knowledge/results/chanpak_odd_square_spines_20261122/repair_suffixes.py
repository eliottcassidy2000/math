"""Try an arbitrary high-word repair and up to two low-bit edits at a near witness."""
from argparse import ArgumentParser
from collections import Counter
from itertools import combinations
from math import comb

from rational_carries import TARGET, chart, joint_approximants, one_side


WITNESSES = [(735435556549430,1575929833923500), (4814468855604222,240731609083994)]


def search(width, radius, index):
    data = chart(4)
    x,y = WITNESSES[index]
    base = 1<<width
    seed = x%base + base*(y%base)
    stats = Counter(carries=0, divisor_pairs=0, face_reads=0)
    hits = []
    pairs = 0
    for size in range(radius+1):
        for bits in combinations(range(2*width),size):
            packed = seed ^ sum(1<<b for b in bits)
            j,i = divmod(packed,base)
            left,right = joint_approximants(width,i,j)
            hits.extend(one_side(data,TARGET,width,left,right,stats))
            pairs += 1
    assert pairs == sum(comb(2*width,r) for r in range(radius+1))
    print(index,width,radius,pairs,dict(stats),"hits",hits,flush=True)


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--width",type=int,default=40)
    parser.add_argument("--radius",choices=(0,1,2),type=int,default=2)
    parser.add_argument("--index",choices=(0,1),type=int)
    args = parser.parse_args()
    for index in ((args.index,) if args.index is not None else (0,1)):
        search(args.width,args.radius,index)
