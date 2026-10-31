"""Exact common-dilation bank attachment and literal interval-clique controls."""
from fractions import Fraction as F
from itertools import combinations
from math import lcm, ceil
from collections import Counter
from pathlib import Path
import argparse
import json

def merge(intervals):
    answer=[]
    for left,right in sorted(intervals):
        if left > right:
            continue
        if answer and left <= answer[-1][1]:
            answer[-1]=(answer[-1][0],max(right,answer[-1][1]))
        else:
            answer.append((left,right))
    return answer

def intersection(left,right):
    result=[]
    i=j=0
    while i < len(left) and j < len(right):
        lo=max(left[i][0],right[j][0]); hi=min(left[i][1],right[j][1])
        if lo <= hi:
            result.append((lo,hi))
        if left[i][1] < right[j][1]:
            i+=1
        else:
            j+=1
    return merge(result)

def safe_set(speeds, radius):
    result=[(F(0),F(1))]
    for speed in sorted(set(speeds)):
        teeth=[((k+radius)/speed,(k+1-radius)/speed) for k in range(speed)]
        result=intersection(result, teeth)
    return result

def forward_gap(phase, bank):
    for left,right in bank:
        if left <= phase <= right:
            return F(0)
        if phase < left:
            return left-phase
    return 1+bank[0][0]-phase

def largest_gap(bank):
    if not bank:
        return None
    return max([bank[0][0]+1-bank[-1][1]]+
               [bank[i+1][0]-bank[i][1] for i in range(len(bank)-1)])

def attached(core, bank, multiplier):
    if not core or not bank:
        return None
    for left,right in core:
        gap=forward_gap((multiplier*left) % 1, bank)
        if gap <= multiplier*(right-left):
            return left+gap/multiplier
    return None

def complete_rejections(core, bank):
    if not core or not bank:
        return {'kind':'all'}
    width=max(right-left for left,right in core)
    if width:
        bound=max(1,ceil(largest_gap(bank)/width))
        return {'kind':'finite','tail':bound,
                'rejected':[d for d in range(1,bound) if attached(core,bank,d) is None]}
    period=lcm(*(left.denominator for left,_ in core))
    return {'kind':'periodic','period':period,
            'rejected_residues':[d for d in range(period) if attached(core,bank,d or period) is None]}

def clique_number(core, speeds, radius):
    """Closed interval graph, preserving actual tooth coordinates and colors."""
    vertices=[]
    for color,speed in enumerate(speeds):
        vertices.extend((left,right,color) for left,right in intersection(core,safe_set([speed],radius)))
    # A maximum overlap occurs at one of the closed left endpoints, including point components.
    return max((sum(left <= t <= right for left,right,_ in vertices)
                for t in {left for left,_,_ in vertices}), default=0), len(vertices)

def verify():
    checks=0
    def need(ok,name):
        nonlocal checks
        checks+=1
        if not ok:
            raise RuntimeError(name)
    native=0
    for core_speeds in combinations(range(1,7),2):
        for bank_speeds in combinations(range(1,7),2):
            for radius in (F(1,5),F(1,4),F(1,3)):
                core=safe_set(core_speeds,radius); bank=safe_set(bank_speeds,radius)
                for d in range(1,17):
                    witness=attached(core,bank,d)
                    truth=safe_set(tuple(core_speeds)+tuple(d*v for v in bank_speeds),radius)
                    need((witness is not None)==bool(truth),'bank compiler versus literal comb intersection')
                    if witness is not None:
                        need(all(min((v*witness)%1,1-(v*witness)%1)>=radius
                                 for v in tuple(core_speeds)+tuple(d*v for v in bank_speeds)),
                             'actual common time witness')
                    native+=1
    radius=F(1,14)
    banks={a:safe_set(a,radius) for a in combinations(range(1,9),2)}
    families=[]
    for c in combinations(range(1,14),11):
        core=safe_set(c,radius)
        for a,bank in banks.items():
            rejection=complete_rejections(core,bank)
            need(rejection['kind']=='finite' and not rejection['rejected'],
                 'actual thirteen-speed scaled-bank family')
            families.append({'core':c,'bank':a,**rejection})
    hostile_core=safe_set((1,2),F(1,4))
    need(attached(hostile_core,safe_set((3,),F(1,4)),1) is not None and
         attached(hostile_core,safe_set((12,),F(1,4)),1) is not None and
         attached(hostile_core,safe_set((3,12),F(1,4)),1) is None,
         'independent arrivals are not a common time')
    exceptional=complete_rejections(hostile_core,safe_set((3,12),F(1,4)))
    point=complete_rejections(safe_set((1,2),F(1,3)),safe_set((1,2),F(1,3)))
    need(point=={'kind':'periodic','period':3,'rejected_residues':[0]},'point components retain permanent residue exclusions')
    # A triangle in a graph on whole safe sets is insufficient; component intervals obey Helly.
    pairwise=[bool(safe_set(pair,F(1,3))) for pair in combinations((1,2,3),2)]
    overlap,vertices=clique_number([(F(0),F(1))],(1,2,3),F(1,3))
    need(all(pairwise) and not safe_set((1,2,3),F(1,3)) and overlap==2,
         'runner-level K3 versus literal interval clique')
    for speeds in combinations(range(1,8),3):
        for radius in (F(1,5),F(1,4),F(1,3),F(1,2)):
            cl,_=clique_number([(F(0),F(1))],speeds,radius)
            need((cl==len(speeds))==bool(safe_set(speeds,radius)),'interval clique exactly retains common phase')
    big=10**100+39
    witness=attached(safe_set((1,2),F(1,4)),safe_set((3,12),F(1,4)),big)
    need(witness is not None and all(min((v*witness)%1,1-(v*witness)%1)>=F(1,4)
             for v in (1,2,3*big,12*big)), 'hundred-digit multiplier without tooth enumeration')
    result={'native_comparisons':native,'checks':checks,'anchor_families':len(families),
            'anchor_core_count':78,'anchor_bank_count':28,
            'anchor_tail_histogram':dict(sorted(Counter(f['tail'] for f in families).items())),
            'exceptional_family_core_1_2_bank_3_12_radius_1_4':exceptional,
            'point_family':point,'runner_graph_hostile':{'pairwise':pairwise,'component_clique_number':overlap,
                                                       'component_vertices':vertices},
            'large_multiplier':str(big),'large_witness':str(witness),
            'scope':'All multipliers proved by the displayed stopping bound for each finite core/bank entry; no arbitrary-core LRC claim.'}
    return result,families

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--output',type=Path)
    ap.add_argument('--families',type=Path)
    args=ap.parse_args()
    result,families=verify()
    text=json.dumps(result,indent=2,sort_keys=True)+'\n'
    if args.output: args.output.write_text(text,encoding='utf-8',newline='\n')
    else: print(text,end='')
    if args.families: args.families.write_text(json.dumps(families,separators=(',', ':'))+'\n',encoding='utf-8',newline='\n')
