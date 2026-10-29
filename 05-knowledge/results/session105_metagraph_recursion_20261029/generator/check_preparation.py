#!/usr/bin/env python3
"""Exact transport-derived parent groups plus orbit-first full-port bank.
No child-sized canonicalization supplies production data. Guards survive -O.
"""
from pathlib import Path
from hashlib import sha256
from collections import Counter
from fractions import Fraction
from math import factorial
from statistics import median
from time import perf_counter
import argparse
import json
import sys

HERE = Path(__file__).resolve().parent
ap = argparse.ArgumentParser()
ap.add_argument('--baseline-root', type=Path)
ap.add_argument('--max-order', type=int, default=7)
ap.add_argument('--repeats', type=int, default=3)
args = ap.parse_args()
# Self-contained packet first; explicit override remains hash-verified.
roots = ([args.baseline_root] if args.baseline_root else [HERE.parent / 'vendor'] + list(HERE.parents))
rel = Path('05-knowledge/results/session88_symmetry_certificates_20261028/blockfill')
root = next((r for r in roots if (r / rel / 'check_adaptive.py').is_file()), None)
if root is None: raise RuntimeError('Need --baseline-root containing Session87 and Session88 published sources')
for name, expected in [('baseline_sources.py','684368faaf908d9d28b966bfa8f9642d0ba1f3848e1083e6531c86bd3cc8bdbc'), ('check_adaptive.py','103b1b0e9e1929f2b71885901ae2f5c84ea01fc3f63bd60d216fa9ba6042adf7'), ('check_blockfill.py','cb246f280786f6a7f8b1f829d9d9b9e0993cb2298eb399783809bedb27780bdb')]:
    if sha256((root/rel/name).read_bytes()).hexdigest() != expected: raise RuntimeError('Baseline source mismatch: '+name)
sys.path.insert(0,str(root/rel))
sys.argv=[sys.argv[0], '--baseline-root', str(root)]
import check_adaptive as adaptive
base=adaptive.base
old=adaptive.old
need=old.need
original_prepare=base.prepare


def prepare(parents,n):
    """Same exact parent_data and states, no raw-mask classification requests."""
    start=perf_counter()
    engine=old.ParentCanonicalizer(max_order=n)
    states={}; parent_data=[]; orbit_sizes=[]; seeded=0
    for raw in parents:
        code,incoming=engine.maps(raw)
        rows=old.matrix(n,code)
        inverse_first=base.inverse(incoming[0])
        group=tuple(sorted(base.compose(inverse_first,p) for p in incoming))
        need(group and all(old.encode(rows,p)==code for p in group), 'bad transported group')
        need(all(sorted(p)==list(range(n)) for p in group), 'bad group permutation')
        # The transport torsor gives ALL canonical-parent attaining maps.
        value=(code,group)
        if rows in engine.cache: need(engine.cache[rows]==value,'transport cache collision')
        else: engine.cache[rows]=value; seeded+=1
        parent_data.append((code,raw,incoming,group))
        seen=set()
        for mask in range(1<<n):
            if mask in seen: continue
            orbit={old.move_mask(mask,p) for p in group}
            need(mask==min(orbit),'orbit minimum traversal')
            stabilizer=sum(old.move_mask(mask,p)==mask for p in group)
            need(len(orbit)*stabilizer==len(group),'mask orbit-stabilizer')
            seen.update(orbit); orbit_sizes.append(len(orbit))
            states[(code,mask)]=old.extend(rows,mask)
        need(len(seen)==1<<n,'incomplete mask bank')
    return engine,states,tuple(sorted(parent_data)), {
        'parent_preparation_requests':len(parents), 'raw_state_classification_requests':0,
        'preparation_total_requests':engine.requests, 'preparation_searches':engine.searches,
        'preparation_orders_tested':engine.orders_tested,'preparation_seconds':perf_counter()-start,
        'transport_seeded_cache_entries':seeded, 'raw_masks_covered':sum(orbit_sizes),
        'mask_orbit_size_histogram':dict(sorted(Counter(orbit_sizes).items()))}


def build(parents,n,new):
    base.prepare=prepare if new else original_prepare
    try: return adaptive.build(parents,n)
    finally: base.prepare=original_prepare


def weighted(result,n):
    # Each pointed orbit has weight 1/|Aut(U,v)|; sum is N/|Aut(U)|.
    total=Fraction(0); block_total=Fraction(0)
    distinct={row for row in result['Q'].values()}
    for row in distinct:
        aut=None; lhs=Fraction(0)
        for key,m in row:
            G=len(result['automorphisms'][key])
            if aut is None: aut=G
            need(G==aut,'same child group order')
            need(G % m == 0,'nonintegral point stabilizer')
            lhs+=Fraction(1,G//m)
        need(lhs==Fraction(n,aut),'weighted deletion identity')
        block_total+=lhs
    for code,raw,incoming,group in result['parent_data']:
        total+=Fraction(1<<(n-1),len(group))
    need(total==block_total,'weighted extension/deletion identity')
    labelled=Fraction(factorial(n),n)*total
    need(labelled==1<<(n*(n-1)//2),'labelled mass')
    return {'weighted_extension_mass':str(total),'weighted_deletion_mass':str(block_total),'labelled_child_mass':str(labelled),'blocks':len(distinct)}


def main():
    need(1<=args.max_order<=8 and args.repeats>=1,'scope guard')
    parents=[()]; reports=[]
    for n in range(args.max_order):
        a=build(parents,n,False); b=build(parents,n,True)
        hashes=base.compare(a,b)
        independent=base.independent_audit(b,n+1)
        weights=weighted(b,n+1)
        timings={'baseline':[],'orbit_transport':[]}
        if n+1>=6:
            for k in range(args.repeats):
                for new in ([False,True] if k%2==0 else [True,False]):
                    t=perf_counter(); c=build(parents,n,new); elapsed=perf_counter()-t
                    base.compare(b,c)
                    timings['orbit_transport' if new else 'baseline'].append(elapsed)
        medians={key:median(v) for key,v in timings.items() if v}
        reports.append({'order':n+1,'baseline':a['counts'],'orbit_transport':b['counts'],
            'equal_complete_output_sha256':hashes,'independent_audit':independent,'weighted':weights,
            'timings':timings,'median_seconds':medians})
        parents=list(b['bank'])
    return {'status':'FINITE-EXACT','baseline':'Session88 adaptive with unchanged full-port output',
        'universe':'All ordinary tournament isomorphism classes; no converse merging or filters',
        'results':reports,'hostile_81_27':old.blowup_hostile(), 'row_guards':base.hostile_tests()}

if __name__=='__main__': print(json.dumps(main(),indent=2,sort_keys=True))
