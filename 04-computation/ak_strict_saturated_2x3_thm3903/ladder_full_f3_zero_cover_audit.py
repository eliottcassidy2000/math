#!/usr/bin/env python3
"""Test the live-cycle zero covers forced by a 3+3 first wave."""

from itertools import combinations, combinations_with_replacement
import sympy as sp

from ladder_all_cost10_kernel_search import make_topology
from ladder_full_f4_symbolic_audit import generator_matrix

W1_SUPPORT = None


def minimal_zero_covers(rest):
    b,c,g,h=sp.symbols("b c g h")
    vals=(b,c,g,h)
    w1=(b,-b,0,-g,g,0)
    w2=(0,c,-c,0,-h,h)
    minors=[sp.expand(w1[i]*w2[j]-w1[j]*w2[i]) for i,j in combinations(rest,2)]
    good=[]
    for size in range(5):
        for subset in combinations(range(4),size):
            if all(poly.subs({vals[i]:0 for i in subset})==0 for poly in minors):
                if not any(set(old).issubset(subset) for old in good):
                    good.append(subset)
    return tuple(good)


def main():
    cases=0; terminal_six=0; rank_counts={}; examples=[]
    for positions in combinations_with_replacement(range(6),3):
        topology=make_topology((0,1,2,3),positions)
        if topology is None: continue
        params,A=generator_matrix(topology)
        for fire in combinations(range(6),3):
            rest=tuple(v for v in range(6) if v not in fire)
            cols=tuple(2*v+j for v in rest for j in (0,1))
            for cover in minimal_zero_covers(rest):
                subs={params[0]:0}
                subs.update({params[1+i]:0 for i in cover})
                T=A[:,cols].subs(subs)
                rank=T.to_DM().to_field().rank()
                cases+=1;rank_counts[rank]=rank_counts.get(rank,0)+1
                if rank==6:
                    terminal_six+=1
                    if len(examples)<20:examples.append((positions,fire,cover))
    print("LADDER_FULL_F3_LIVE_ZERO_COVER_AUDIT")
    print(f"cases={cases};terminal_rank_six={terminal_six}")
    print(f"rank_counts={tuple(sorted(rank_counts.items()))}")
    print(f"examples={tuple(examples)}")


if __name__=="__main__":main()
