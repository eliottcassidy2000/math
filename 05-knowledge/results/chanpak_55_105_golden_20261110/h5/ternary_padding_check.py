#!/usr/bin/env python3
"""Small formula and graph controls for the one-active-singleton theorem."""
import hashlib
import importlib.util
import json
from collections import Counter
from pathlib import Path

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('fixed_core_compiler',HERE/'check.py')
c=importlib.util.module_from_spec(spec)
spec.loader.exec_module(c)


def strong(rows):
    full=(1<<len(rows))-1
    for i in range(len(rows)):
        reached=1<<i
        while True:
            more=reached
            for j in range(len(rows)):
                if reached>>j&1:
                    more|=rows[j]
            if more==reached:
                break
            reached=more
        if reached!=full:
            return False
    return True


def main():
    cores={'C3':(2,4,1),'Q4':(4,1,10,3),
           'H9a':(2,4,9,19,7),'H9b':(2,4,25,3,11)}
    checks=Counter()
    results={}
    for name,core in cores.items():
        c.require(c.tournament(core) and strong(core),'nontrivial strong core')
        cf=c.coefficients(core)
        h,q=cf['h'],cf['q']
        c.require(q>=(len(core)+1)*h,'full path-cut injection')
        c.require(min(cf['cover_start']+cf['cover_end'])>=h,'per-vertex path-cut injection')
        fullmask=(1<<len(core))-1
        for cut in cf['cuts']:
            outmask=cut['outmask']
            c.require(cut['e'] >= (len(core)-outmask.bit_count())*h,'end set lower bound')
            c.require(cut['s'] >= outmask.bit_count()*h,'start set lower bound')
            if outmask.bit_count()==len(core)-1:
                c.require(cut['t']+cut['left_parity']==h+cut['e'],'singleton exception left')
                checks['singleton_exception_identities']+=1
            if outmask.bit_count()==1:
                c.require(cut['t']+cut['right_parity']==h+cut['s'],'singleton exception right')
                checks['singleton_exception_identities']+=1
            for a in range(8):
                for b in range(8):
                    n=c.value(cf,cut,a,b)
                    if n in (3*h+2,3*h+4):
                        c.require(outmask not in (0,fullmask),'uniform cuts cannot branch')
                        c.require((a,b) in ((0,0),(1,0),(0,1)),'ternary active padding bound')
                        checks['small_box_branch_solutions']+=1
                    checks['formula_triples']+=1
            if outmask not in (0,fullmask):
                for k in (2,3,4):
                    c.require(c.value(cf,cut,(1<<k)-1,0)>=5*h,'long left flank fivefold threshold')
                    c.require(c.value(cf,cut,0,(1<<k)-1)>=5*h,'long right flank fivefold threshold')
                    checks['long_flank_bounds']+=2
        image,_=c.support(core,3*h+4)
        branch_rows={}
        for target in (3*h+2,3*h+4):
            candidates=[(cut['outmask'],a,b)
                        for cut in cf['cuts'] if cut['outmask'] not in (0,fullmask)
                        for a,b in ((0,0),(1,0),(0,1))
                        if c.value(cf,cut,a,b)==target]
            c.require(bool(candidates)==(target in image),'complete inverse equals three-state criterion')
            branch_rows[target]={'exists':bool(candidates),'candidate_count':len(candidates)}
            if candidates:
                outmask,a,b=candidates[0]
                witness=c.graph(core,a,b,outmask)
                c.require(c.hp(witness)==target and strong(witness),'strong minimal branch witness')
                # Restore unused neutral names on BOTH sides. Same active mask,
                # larger widths: old names are present, and only the scalar count is read.
                for widths in ((1,1),(2,2),(3,1)):
                    padded=c.graph(core,a,b,outmask,widths)
                    c.require(c.tournament(padded) and c.hp(padded)==target,'common-source padding reconstruction')
                    checks['padded_graph_witnesses']+=1
                branch_rows[target]['witness']={'outmask':outmask,'A':a,'B':b,'rows':witness}
            checks['branch_support_comparisons']+=1
        results[name]={'core':core,'h':h,'branches':branch_rows}
    equality=c.graph(cores['C3'],3,0,3)
    c.require(c.tournament(equality) and strong(equality),'strong sharp equality tournament')
    c.require(c.hp(equality)==c.literal_hp(equality)==15,'sharp fivefold equality')
    print(json.dumps({'status':'PASS','scope':'four named cores; finite controls of a proved all-order reduction',
                      'checks':dict(checks),'cores':results,'sharp_boundary':{'rows':equality,'H':15,'h':3,'active_singletons':2},
                      'dependencies_sha256':{'check.py':hashlib.sha256((HERE/'check.py').read_bytes()).hexdigest()},
                      'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2,sort_keys=True))


if __name__=='__main__':
    main()
