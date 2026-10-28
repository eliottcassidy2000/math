"""Coarsest probabilistic congruence for a declared consecutive-pair observer."""
from array import array
from collections import Counter,defaultdict
from functools import reduce
from math import gcd
from pathlib import Path
import argparse,json,sys
HERE=Path(__file__).resolve().parent
def need(b,message):
    if not b: raise RuntimeError(message)
def intern(values):
    ids={};result=[]
    for value in values:
        if value not in ids:ids[value]=len(ids)
        result.append(ids[value])
    return result,len(ids)
def canonical_probability(counter):
    d=reduce(gcd,counter.values())
    return tuple((k,v//d) for k,v in sorted(counter.items()))
def main():
    ap=argparse.ArgumentParser();ap.add_argument('directory',type=Path);args=ap.parse_args()
    reports=[]
    for n in range(2,8):
        q=array('I');q.frombytes((args.directory/f'flag_states{n}.u32').read_bytes())
        if sys.byteorder!='little':q.byteswap()
        meta=json.loads((args.directory/f'flag_meta{n}.json').read_text(encoding='utf-8'))
        m=(n-1)*(n-2)//2;full=(1<<m)-1
        need(len(q)==n*(1<<m),'complete fixed-path rooted universe')
        tau=[((s//n)^full)*n+s%n for s in range(len(q))]
        root_groups=[[] for _ in meta]
        for s,t in enumerate(q):root_groups[t].append(s)
        # Only the ordinary child and its ordinary deletion class are outputs.
        observed=[(meta[t][2],meta[t][0]) for t in q]
        colors,count=intern(observed);trajectory=[count]
        while True:
            profiles,_=intern(canonical_probability(Counter(colors[s] for s in group))
                              for group in root_groups)
            upcoming,nextcount=intern((colors[s],colors[tau[s]],profiles[q[s]]) for s in range(len(q)))
            if nextcount==count:
                need(upcoming==colors,'stable exact refinement')
                break
            colors,count=upcoming,nextcount;trajectory.append(count)
        # A stronger rooted-type observer has the exact reset/involution pair quotient.
        pairs,pair_count=intern((q[s],q[tau[s]]) for s in range(len(q)))
        factor={}
        for s,p in enumerate(pairs):
            if p in factor:need(factor[p]==colors[s],'coarse carrier factors through rooted pair')
            else:factor[p]=colors[s]
        witnesses={}
        by_output=defaultdict(list);by_color=defaultdict(list)
        for s,o in enumerate(observed):by_output[o].append(s)
        for s,c in enumerate(colors):by_color[c].append(s)
        for group in by_output.values():
            s=group[0]
            t=next((t for t in group if colors[t]!=colors[s]),None)
            if t is not None:witnesses['equal_child_parent_classes_different_futures']=[s,t];break
        for group in by_color.values():
            s=group[0];t=next((t for t in group if q[t]!=q[s]),None)
            if t is not None:witnesses['distinct_pointed_types_same_declared_future']=[s,t];break
        out={'order':n,'literal_path_root_states':len(q),'pointed_types':len(meta),
             'child_parent_observations':trajectory[0],'coarsest_closed_carrier':count,
             'refinement_sizes':trajectory,'stronger_pointed_pair_carrier':pair_count,
             'witnesses':witnesses}
        reports.append(out);print(json.dumps(out),flush=True)
        # Replayable finite classifiers and decoding metadata, compressed later.
        import gzip
        payload={'order':n,'pointed_types':meta,'pointed_id_by_state':list(q),'closed_color_by_state':colors}
        (HERE/f'closed_flags{n}.json.gz').write_bytes(gzip.compress(
            (json.dumps(payload,separators=(',',':'))+'\n').encode(),mtime=0))
    (HERE/'closed_flags.json').write_text(json.dumps({'status':'FINITE-EXACT',
       'operations':['native fixed-path complement','uniform fresh re-marking preserving distinguished vertex'],
       'output':'ordinary child class and ordinary deletion-parent class',
       'results':reports},indent=2)+'\n',encoding='utf-8')
if __name__=='__main__':main()
