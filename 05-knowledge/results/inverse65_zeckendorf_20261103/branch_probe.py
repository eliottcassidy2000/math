"""Exact cone certificates and complete source-aware small obstructions."""
from pathlib import Path
from fractions import Fraction
import importlib.util,json,sys
from bank import *
HERE=Path(__file__).resolve().parent
RAYS=((2,1,1),(1,0,1),(0,1,1))  # physical A/B/C cone, with C>=H

def main():
    growth={};slack=[]
    for l,m in MATS.items():
        need(all(5*sum(m[i][j]*(6,9,5)[j] for j in range(3))>=9*(6,9,5)[i] for i in range(3)),('uniform norm',l))
        ratios=[Fraction(step(v,m)[0],v[0]) for v in RAYS if v[0]]
        need(step(RAYS[2],m)[0]>=0,'recession ray')
        growth[l]=str(min(ratios))
        if l!='a':need(min(ratios)>=5,('branch gap',l))
        if l=='a':continue
        for trailing in range(3):
            rat=[]
            for v in RAYS:
                w=step(v,m)
                for _ in range(trailing):w=step(w,MATS['a'])
                rat.append(Fraction(w[1]+w[2]-w[0],w[0]))
            slack.append((min(rat),l,trailing))
    need(min(r[0] for r in slack)==Fraction(1,96),'uniform long-word slack')
    # Prior graph DP is separate from our matrix decoder.
    prior=HERE.parent/'golden_hamilton_transfer_20261103';sys.path.insert(0,str(prior))
    spec=importlib.util.spec_from_file_location('graph_control',prior/'verify.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    queue=[((),(1,0,0))];sources=[]
    for word,v in queue:
        if v[0]==55:
            rows=actual_target(word);counts=[mod.h(mod.add(rows,cut)) for cut in range(1<<len(rows))]
            need(111 not in counts,('55->111 false positive',word))
            sources.append({'word':word,'state':v,'rows':rows,'attachment_H':counts})
        for l,m in MATS.items():
            w=step(v,m)
            if w[0]<=55:queue.append((word+(l,),w))
    need(len(queue)==334 and len(sources)==6,'all-order H55 census')
    # Exact two-singleton offset solutions: long words have H<=192; short words exhaustive.
    ss,pp,_=closure(192)
    exceptions={v:recover(i,pp) for i,v in enumerate(ss) if v[1]+v[2]-v[0] in (1,2)}
    short=[((),(1,0,0))];short_high=[]
    for depth in range(1,4):
        nxt=[]
        for word,v in short:
            for l,m in MATS.items():
                w=step(v,m);nw=word+(l,);nxt.append((nw,w))
                if w[1]+w[2]-w[0] in (1,2):
                    exceptions.setdefault(w,nw)
                    if w[0]>192:short_high.append((nw,w))
        short=nxt
    need(not short_high,'unexpected large short-word exception')
    need(len(exceptions)==41 and max(v[0] for v in exceptions)==155,'complete aa offset census')
    controls=json.loads((prior/'verification.json').read_text())
    positive=controls['positive55to111'];need(mod.h(positive['core'])==55 and mod.h(mod.add(positive['core'],positive['mask']))==111,'old positive control')
    prior_hostile=HERE.parent/'primegame_port_spectrum_20261103';sys.path.insert(0,str(prior_hostile))
    spc=importlib.util.spec_from_file_location('hostile_replay',prior_hostile/'verify.py');hm=importlib.util.module_from_spec(spc);spc.loader.exec_module(hm)
    hostile=hm.hostile()
    need(hostile=={'groups':[81,27],'common_H':[462915,462915],
                  'port8_after_sink':[2390013,1754919],'common_cyclic_triples':67},'literal hostile control')
    result={'status':'PASS','growth_factors':growth,'minimum_suffix_slack':str(min(x[0] for x in slack)),
            'minimum_slack_cases':[(l,t) for x,l,t in slack if x==Fraction(1,96)],
            'words_H_le55':len(queue),'H55_sources':sources,'positive55to111':positive,'hostile81_27':hostile,
            'two_singleton_positive_ternary_offsets':[
                {'state':v,'word':w,'H':v[0],'offset':2*(v[1]+v[2]-v[0])}
                for v,w in sorted(exceptions.items())],
            'short_words_tested':sum(65**j for j in range(1,4))}
    (HERE/'branches.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'PASS','H55_sources':len(sources),'minimum_suffix_slack':result['minimum_suffix_slack'],
                      'ternary_ss_exceptions':result['two_singleton_positive_ternary_offsets'],'hostile81_27':hostile}))
if __name__=='__main__':main()
