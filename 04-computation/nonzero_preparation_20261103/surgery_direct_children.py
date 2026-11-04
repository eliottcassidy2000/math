#!/usr/bin/env python3
"""Two strong second-round repairs from the direct4289 physical target.
No additional permanent preparation names; both outputs carry positive guards.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import surgery_native as native
import surgery_ready as ready
import surgery_direct as direct

GATES=0

def need(ok,label):
    global GATES
    GATES+=1
    if not ok:raise RuntimeError(label)

CONFIG = [{'H': 12869, 'gap': 1, 'core': [2, 4, 1], 'ports': [0], 'word': [1, 0, 0, 1, 1, 0, 0, 0, 1, 1], 'mapping': [12, 13, 11, 10, 7, 2, 9, 1, 4, 5, 6, 3, 0, 8, 14], 'aux': 14, 'root': 8, 'source': [32510, 32388, 31872, 32758, 32390, 32662, 32694, 32000, 28183, 28804, 31232, 29184, 16640, 4096, 8192], 'target': [32510, 32388, 32128, 32502, 32646, 32662, 32694, 32000, 11787, 31876, 30720, 24576, 18688, 20480, 256], 'pairs': [[2, 8], [3, 8], [4, 8], [8, 14], [9, 10], [9, 11], [11, 12], [13, 14]], 'allstrong': True, 'bad': 0}, {'H': 12871, 'gap': 2, 'core': [2, 4, 9, 3], 'ports': [0, 1], 'word': [1, 1, 1, 0, 0, 1, 0, 0, 1], 'mapping': [11, 7, 13, 14, 2, 9, 10, 1, 5, 4, 6, 3, 0, 8, 12], 'aux': 12, 'root': 8, 'source': [32510, 32388, 31872, 32758, 32390, 32662, 32694, 32000, 28183, 28804, 31232, 29184, 16640, 4096, 8192], 'target': [32510, 32644, 30848, 32758, 32422, 32646, 32694, 22784, 26133, 30852, 31364, 28928, 256, 4224, 12288], 'pairs': [[1, 8], [2, 10], [4, 5], [7, 10], [7, 13], [8, 11], [9, 11], [12, 14]], 'allstrong': True, 'bad': 0}]


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();base=Path(__file__).parent
    need(sha256((base/'surgery_direct.py').read_bytes()).hexdigest()
         =='10c5e39d6439057818622607de83d97d666eddf187c3ff30099e051702538db9','frozen direct physical supplier')
    original_cfg=direct.CONFIG[0]
    canonical,r,a=ready.ready_target(tuple(original_cfg['core']),tuple(original_cfg['ports']),tuple(original_cfg['word']))
    old=direct.relabel(canonical,original_cfg['mapping'])
    need(old==tuple(original_cfg['target']) and original_cfg['mapping'][r]==8,'actual direct4289 physical input')
    aux=original_cfg['mapping'][a]
    first,source=native.c2(old,aux)
    prefix=[]
    for rows in(old,first,source):
        f=ready.counts(rows,(8,));b=native.signature(rows,1<<8)
        need(native.strong(rows) and (f['H'],f['F'][8])==(b[0],b[2]),'common actual strong C2 prefix')
        prefix.append({'rows':rows,'H':b[0],'F':b[2],'defect':b[2]-b[0]})
    need((prefix[0]['H'],prefix[0]['F'])==(4289,7653)
         and(prefix[-1]['H'],prefix[-1]['F'])==(12871,22959),'physical source distinguished by full root response')
    cases=[]
    for cfg,expected in zip(CONFIG,((16599,4311,12288),(22291,5907,16384))):
        T,r,a=ready.ready_target(tuple(cfg['core']),tuple(cfg['ports']),tuple(cfg['word']))
        target=direct.relabel(T,cfg['mapping'])
        need(cfg['mapping'][r]==8 and cfg['mapping'][a]==cfg['aux'],'full retained port map')
        need(target==tuple(cfg['target']) and source==tuple(cfg['source']),'literal full named matrices')
        names=[f'v{i}'for i in range(15)];names[8]='r'
        rec=native.execute(names,source,target,[tuple(p)for p in cfg['pairs']])
        need(all(rec['native_stage_strongness']),'all native intermediate states strong')
        obs=native.joint_readout(target,8,cfg['aux'])
        f=ready.counts(target,(8,))
        need((f['H'],f['F'][8])==(obs['H'],obs['F']),'independent actual target H/F')
        need(obs['H']==cfg['H'] and obs['sigma']==cfg['gap'] and obs['tau']==0,
             'actual positive output joint guard')
        need((obs['F'],obs['L'],obs['K'])==expected,'exact joint readout control')
        need(obs['H'] in(3*4289+2,3*4289+4),'two actual nonzero count children')
        rec.update({'construction':cfg,'joint':obs,'root_index':8,'aux_index':cfg['aux'],
                    'digit_from4289':obs['H']-3*4289,'ready_order':15,
                    'permanent_births_from4289':2,'peak_order':16})
        cases.append(rec)
    # A single physical first+4 child supplies both branches: use its gap1
    # immediately, or keep H fixed while rearming to the literal gap2 target.
    B=tuple(direct.CONFIG[1]['target']);C=tuple(direct.CONFIG[2]['target'])
    switch_pairs=[(2,7),(2,9),(2,12),(4,8),(6,8),(7,9),(8,11),(8,12),(9,11),(11,12)]
    names=[f'v{i}'for i in range(13)];names[8]='r'
    switch=native.execute(names,B,C,switch_pairs)
    need(all(switch['native_stage_strongness']),'every same-H rearm native state strong')
    before=native.joint_readout(B,8,direct.CONFIG[1]['aux'])
    after=native.joint_readout(C,8,direct.CONFIG[2]['aux'])
    need((before['H'],before['F'],before['sigma'],before['tau'])==(4291,6403,1,0)
         and(after['H'],after['F'],after['sigma'],after['tau'])==(4291,8371,2,0),
         'same physical H with changed positive guard and full root response')
    switch.update({'root_index':8,'source_joint':before,'target_joint':after,
                   'source_aux':direct.CONFIG[1]['aux'],'target_aux':direct.CONFIG[2]['aux']})
    leaves=[{'H':c['joint']['H'],'F':c['joint']['F'],'defect':c['joint']['root_defect'],
             'rows':c['target_rows'],'positive_guard':c['joint']['sigma']}for c in cases]
    for rows,aux,expected in((B,direct.CONFIG[1]['aux'],(12875,19209)),
                            (C,direct.CONFIG[2]['aux'],(12877,25113))):
        first,child=native.c2(rows,aux)
        f=ready.counts(child,(8,));b=native.signature(child,1<<8)
        need((f['H'],f['F'][8])==expected==(b[0],b[2]),'actual same-child +4 side leaves')
        leaves.append({'H':b[0],'F':b[2],'defect':b[2]-b[0],'rows':child,'positive_guard':None})
    need({leaf['H']for leaf in leaves}=={12869,12871,12875,12877}
         and all(len(leaf['rows'])==15 and native.strong(tuple(leaf['rows']))for leaf in leaves),
         'four actual strong leaves, exactly fifteen retained names each')
    total=GATES+native.GATES+ready.GATES+direct.GATES
    result={'scope':'two actual all-strong direct second-round repairs; both outputs have declared positive joint guards, no unbounded closure',
            'root_index':8,'common_prefix':prefix,'cases':cases,'same_H_rearm':switch,
            'four_leaves':leaves,'gates':total,
            'native_complements':sum(c['native_complements']for c in cases)+switch['native_complements'],
            'helper_only_deletions':sum(c['helper_deletions']for c in cases)+switch['helper_deletions'],
            'strong_native_states':sum(len(c['native_stage_strongness'])for c in cases)+len(switch['native_stage_strongness'])}
    raw=(json.dumps(result,indent=2,sort_keys=True)+'\n').encode();args.output.write_bytes(raw)
    print(json.dumps({'gates':total,'programs':len(cases)+1,'native_complements':result['native_complements'],
                      'strong_native_states':result['strong_native_states'],'receipt_bytes':len(raw),
                      'receipt_sha256':sha256(raw).hexdigest()},sort_keys=True))


if __name__=='__main__':main()
