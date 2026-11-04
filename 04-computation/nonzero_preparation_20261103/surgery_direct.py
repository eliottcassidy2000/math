#!/usr/bin/env python3
"""Zero C2 followed by actual old-arc repair, keeping every old name.
All three declared programs are exact; no count-neutral preparation is assumed.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import surgery_native as native
import surgery_ready as ready

GATES=0

def need(ok,label):
    global GATES
    GATES+=1
    if not ok:raise RuntimeError(label)

def relabel(rows,p):
    need(sorted(p)==list(range(len(rows))),'physical permutation')
    out=[0]*len(rows)
    for i,row in enumerate(rows):
        for j in range(len(rows)):
            if row>>j&1:out[p[i]]|=1<<p[j]
    return tuple(out)

CONFIG = [{'H': 4289, 'gap': 2, 'core': [6, 4, 24, 19, 3], 'ports': [0], 'word': [1, 1, 0, 0, 0, 1], 'mapping': [7, 2, 9, 11, 10, 1, 4, 5, 6, 3, 0, 8, 12], 'aux': 12, 'root': 8, 'source': [7934, 7812, 7816, 7922, 7814, 8086, 7862, 7424, 7775, 4224, 6656, 512, 2048], 'target': [7934, 7812, 7296, 8182, 7814, 8086, 8118, 7424, 3607, 4228, 6656, 4608, 256], 'pairs': [[2, 3], [2, 9], [3, 8], [6, 8], [8, 12], [11, 12]], 'allstrong': True, 'bad': 0}, {'H': 4291, 'gap': 1, 'core': [2, 4, 25, 19, 3], 'ports': [0], 'word': [0, 1, 0, 1, 0, 1], 'mapping': [7, 2, 9, 11, 10, 1, 4, 5, 6, 3, 0, 8, 12], 'aux': 12, 'root': 8, 'source': [7934, 7812, 7816, 7922, 7814, 8086, 7862, 7424, 7775, 4224, 6656, 512, 2048], 'target': [7934, 8068, 7296, 8182, 7814, 8086, 7862, 7936, 3669, 4100, 6656, 4608, 256], 'pairs': [[1, 8], [2, 3], [2, 9], [3, 8], [7, 9], [8, 12], [11, 12]], 'allstrong': True, 'bad': 0}, {'H': 4291, 'gap': 2, 'core': [2, 12, 9, 17, 7], 'ports': [0], 'word': [0, 0, 0, 0, 0, 1], 'mapping': [7, 9, 10, 2, 12, 1, 4, 5, 6, 3, 0, 8, 11], 'aux': 11, 'root': 8, 'source': [7934, 7812, 7816, 7922, 7814, 8086, 7862, 7424, 7775, 4224, 6656, 512, 2048], 'target': [7934, 8068, 3584, 8182, 8070, 8086, 8118, 7428, 5637, 6272, 6656, 256, 2052], 'pairs': [[1, 8], [2, 3], [2, 7], [2, 12], [3, 8], [4, 8], [6, 8], [8, 11], [9, 11]], 'allstrong': True, 'bad': 0}]


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();base=Path(__file__).parent
    for name,pin in [('surgery_ready.py','024ac0b483d832865c6f12eabf0fe0f8676ac39188631e44910a0a4c9e66ca2f'),
                     ('surgery_native.py','66eb9d682e0128bdf087c79b24a44dcfe127daf810ec560de85daa9689e220e2')]:
        need(sha256((base/name).read_bytes()).hexdigest()==pin,'frozen source dependency '+name)
    old,r=native.prepared(3);first,source=native.c2(old,9)
    old_obs=ready.counts(old,(8,));zero_obs=ready.counts(source,(8,))
    need(r==8 and len(old)==11 and (old_obs['H'],old_obs['F'][8])==(1429,1563),'actual eleven-name source/root')
    need(len(source)==13 and (zero_obs['H'],zero_obs['F'][8])==(4287,4689),'literal zero-C2 source response')
    prefix=[]
    for graph in (old,first,source):
        f=ready.counts(graph,(8,));b=native.signature(graph,1<<8)
        need(native.strong(graph) and (f['H'],f['F'][8])==(b[0],b[2]),'common strong prefix and independent observers')
        prefix.append({'rows':graph,'H':b[0],'F':b[2],'defect':b[2]-b[0]})
    cases=[];uniformity=[]
    for cfg,expected_F in zip(CONFIG,(7653,6403,8371)):
        T,canonical_root,canonical_aux=ready.ready_target(tuple(cfg['core']),tuple(cfg['ports']),tuple(cfg['word']))
        target=relabel(T,cfg['mapping'])
        need(cfg['mapping'][canonical_root]==8 and cfg['mapping'][canonical_aux]==cfg['aux'],'retained root and marked auxiliary')
        need(target==tuple(cfg['target']) and source==tuple(cfg['source']),'literal full source/target matrices')
        names=[f'v{i}'for i in range(13)];names[8]='r'
        rec=native.execute(names,source,target,[tuple(p)for p in cfg['pairs']])
        need(all(rec['native_stage_strongness']),'all native intermediate states strong')
        obs=native.joint_readout(target,8,cfg['aux'])
        need(obs['H']==cfg['H'] and obs['H'] in(3*1429+2,3*1429+4),'actual direct nonzero target count')
        need(obs['F']==expected_F and obs['sigma']==cfg['gap'] and obs['tau']==0,'actual next joint guard and changed F')
        next_first,next_child=native.c2(target,cfg['aux'])
        child_obs=native.signature(next_child,1<<8)
        need(child_obs[0]==3*obs['H']+2*cfg['gap'] and child_obs[2]==3*obs['F'],'immediate strong next nonzero step')
        for graph in(target,next_first,next_child):
            f=ready.counts(graph,(8,));b=native.signature(graph,1<<8)
            need((f['H'],f['F'][8])==(b[0],b[2]),'independent actual root observers')
        rec.update({'construction':cfg,'joint':obs,'first_digit':obs['H']-3*1429,
                    'next_digit':2*cfg['gap'],'root_index':8,'aux_index':cfg['aux'],
                    'original_rows':old,'zero_first_rows':first,'zero_child_rows':source,
                    'permanent_births_to_ready':2,'ready_order':13,'repair_peak_order':14,
                    'next_first_rows':next_first,'next_child_rows':next_child,
                    'next_H':child_obs[0],'next_F':child_obs[2],'next_defect':child_obs[2]-child_obs[0],
                    'permanent_births_to_next_child':4,'next_child_order':15})
        cases.append(rec)
        # Old0 is the whole original prefix. Its exact suffix signature pays
        # the complete module-inflation formula, not a sample interpolation.
        vertices=[i for i in range(1,13) if i!=8]
        need((target[8]&1)==1 and all((target[0]>>i)&1 for i in vertices),
             'literal old0 prefix interface')
        suffix=native.induced(target,vertices)
        ports=sum(((target[i]>>8)&1)<<j for j,i in enumerate(vertices))
        sig=native.signature(suffix,ports)
        coefficient,intercept=sig[2],sig[0]-sig[4]
        need(coefficient in (3821,3197,4179) and coefficient*2+intercept==cfg['H'],
             'exact finite suffix supplier')
        inflated=native.attachment(suffix,ports,2)
        ir=2+len(suffix);ia=2+vertices.index(cfg['aux'])
        io=native.joint_readout(inflated,ir,ia)
        expected=3*(756*4-83)+rec['first_digit']
        need(native.strong(inflated) and io['H']==4*coefficient+intercept
             and io['sigma']==cfg['gap'] and io['tau']==0,'first inflated-clock actual guard')
        need(io['H']!=expected and io['H']-expected==2*(coefficient-2268),
             'first-clock failure and exact count-map error')
        uniformity.append({'target_H':cfg['H'],'gap':cfg['gap'],'suffix_rows':suffix,
                           'suffix_ports':ports,'signature':sig,'H_coefficient':coefficient,
                           'H_intercept':intercept,'F_intercept':sig[0],
                           'error_coefficient':coefficient-2268,'clock2_observers':io,
                           'clock2_requested_H':expected})
    total=GATES+native.GATES+ready.GATES
    result={'scope':'three actual all-strong zero-first repair programs, each returns a declared positive C2 guard; finite scope',
            'root_index':8,'common_prefix':prefix,'cases':cases,'uniformity':uniformity,'gates':total,
            'native_complements':sum(c['native_complements']for c in cases),
            'helper_only_deletions':sum(c['helper_deletions']for c in cases),
            'strong_native_states':sum(len(c['native_stage_strongness'])for c in cases)}
    raw=(json.dumps(result,indent=2,sort_keys=True)+'\n').encode();args.output.write_bytes(raw)
    print(json.dumps({'gates':total,'cases':len(cases),'native_complements':result['native_complements'],
                      'strong_native_states':result['strong_native_states'],'receipt_bytes':len(raw),
                      'receipt_sha256':sha256(raw).hexdigest()},sort_keys=True))


if __name__=='__main__':main()
