"""Exact binary inverse of linear-g responses and complete five-core charts.

No search cutoff is used in linear_inverse/full_inverse. max_short in
chart_inverse is explicitly a restricted optional probe, not completeness.
"""
from functools import lru_cache
from math import isqrt
from pathlib import Path
import importlib.util,json

HERE=Path(__file__).resolve().parent
RESULTS=HERE.parents[1]
CORES=((0,),(2,4,1),(2,4,9,3),(2,4,9,19,7),(2,4,25,3,11))


def require(ok,msg):
    if not ok:raise RuntimeError(msg)


def g(n):
    total=1
    while n>=2:
        q=n//2
        if n%2==0:total+=q+(q&1)
        n=q
    return total


def linear_inverse(y,d,c,e,stats=None):
    """All B>=0 with d*g(B)+c*B+e*(B%2)=y; d,c>0,0<=e<=d."""
    require(d>0 and c>0 and 0<=e<=d,'linear inverse coefficient domain')
    if stats is not None:stats['calls']=stats.get('calls',0)+1
    if y<d:return []
    lo=max(0,(y-d-e+c+d-1)//(c+d))
    hi=(y-d)//c
    if lo>hi:return []
    if hi-lo<=8:
        if stats is not None:stats['candidate_checks']=stats.get('candidate_checks',0)+hi-lo+1
        return [b for b in range(lo,hi+1) if d*g(b)+c*b+e*(b&1)==y]
    evens=linear_inverse(y,d,2*c+d,d,stats)
    odds=linear_inverse(y-c-e,d,2*c,0,stats)
    return [2*b for b in evens]+[2*b+1 for b in odds]


def free_axis_inverse(y,d,stats=None):
    """All even representatives B=2q of d*g(B)=y, trailing ones removed."""
    return [2*q for q in linear_inverse(y,d,d,d,stats)]


def chart_inverse(target,coeff,max_short=None,stats=None):
    """Return normalized(A,B) solutions for positive coeff=(h,d,L,i,o,s,e).

    If max_short is omitted the result is complete in the normalized grammar.
    If supplied, all solutions with min(A,B)<=max_short are returned.
    """
    h,d,L,i,o,s,e=coeff
    require(d>0 and L>0 and h>=d and min(i,o,s,e)>=0 and max(s,e)<=d,'actual positive coefficient domain')
    if target<h:return []
    K=isqrt((target-h)//L)
    if max_short is not None:
        require(max_short>=0,'nonnegative short-flank bound')
        K=min(K,max_short)
    found=set()
    # A is the bounded coordinate, then exchange the flanks. This duplicates
    # some interior pairs deliberately; the set removes them exactly.
    for transpose in (False,True):
        u,v,par_u,par_v=(o,i,e,s) if transpose else (i,o,s,e)
        for A in range(K+1):
            base=h+d*(g(A)-2)+u*A+par_u*(A&1)
            y=target-base
            c=L*A+v
            if c==0:
                require(A==0 and h==d and par_v==0,'uniform-axis guard')
                Bs=free_axis_inverse(y,d,stats)
            else:
                Bs=linear_inverse(y,d,c,par_v,stats)
            for B in Bs:
                found.add((B,A) if transpose else (A,B))
    return sorted(found)


def raw_records():
    return json.loads((RESULTS/'chanpak_five_core_bridges_20261110/cover/positive_raw_signatures.json').read_text())


def bank_inverse(target,max_short=None,first=False,stats=None):
    result=[]
    for r in raw_records():
        pairs=chart_inverse(target,r['positive_observer'],max_short,stats)
        for A,B in pairs:
            result.append({'core':r['core_index'],'cut':r['cut'],'A':A,'B':B,
                           'coeff':r['positive_observer']})
            if first:return result
    return result


def main():
    source=RESULTS/'session251_creative_completion_20261111/ternary/compiler.py'
    spec=importlib.util.spec_from_file_location('signed',source)
    oracle=importlib.util.module_from_spec(spec);spec.loader.exec_module(oracle)
    linear_tests=0
    for d in (1,3,5,9):
        for c in range(1,16):
            for e in range(d+1):
                for y in range(1,129):
                    expected=[b for b in range(max(0,(y-d)//c)+1) if d*g(b)+c*b+e*(b&1)==y]
                    require(sorted(linear_inverse(y,d,c,e))==expected,'linear inverse brute comparison')
                    linear_tests+=1
    chart_tests=0;pair_tests=0
    for record in raw_records():
        core=CORES[record['core_index']];cut=record['cut']
        sig=oracle.signature(core,cut)
        expected={n:set() for n in range(1,256,2)}
        for n,A,B in oracle.generate(sig,255):expected[n].add((A,B))
        for target in range(1,256,2):
            result=set(chart_inverse(target,record['positive_observer']))
            require(result==expected[target],f'complete chart inverse mismatch{record},{target}')
            chart_tests+=1;pair_tests+=len(result)
    large=[]
    for target in (425,33857,8820901,102*2**20+19,102*2**25+19):
        stats={}
        witness=bank_inverse(target,max_short=16,first=True,stats=stats)
        for w in witness:
            require(oracle.count(oracle.signature(CORES[w['core']],w['cut']),w['A'],w['B'])==target,'large signed witness')
        large.append({'target':target,'short_flank_probe':16,'witness':witness,'work':stats})
    report={'status':'PASS; exact inverse, no universal coverage claim',
        'linear_brute_cases':linear_tests,'complete_chart_target_cases':chart_tests,
        'matched_normalized_pairs':pair_tests,'large_short_flank_probes':large}
    out=json.dumps(report,indent=2,sort_keys=True)+'\n'
    (HERE/'chart_inverse.json').write_text(out,encoding='utf-8')
    print(out)


if __name__=='__main__':main()
