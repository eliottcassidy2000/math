#!/usr/bin/env python3
"""Literal second-round rearmings: three strong branches and one weak padding.

Imports only the explicitly pinned first-round verifier/selector. All targets
are reconstructed from core, word, padding, and an explicit physical permutation.
No stochastic search runs here and no optimal edit claim is made.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import surgery_native as native
import surgery_ready as ready
from surgery_native import (strong,flip,compile_pair,to_sets,to_bits,native_sets,signature)

GATES=0

def require(ok,label):
    global GATES
    GATES+=1
    if not ok: raise RuntimeError(label)


def relabel(rows,mapping):
    n=len(rows)
    require(sorted(mapping)==list(range(n)), 'literal permutation of every old name')
    out=[0]*n
    for i,row in enumerate(rows):
        for j in range(n):
            if row>>j&1: out[mapping[i]] |= 1<<mapping[j]
    return tuple(out)


def prepend_chain(rows,length):
    n=len(rows)+length
    return (tuple(((1<<n)-1)^((1<<(i+1))-1) for i in range(length))
            +tuple(row<<length for row in rows))

def execute(names,source,target,pairs):
    require(len(names) == len(source) == len(target), 'same physical name set')
    current = source
    program = []
    stage_strong = [strong(source)]
    native_strong = []
    root_index = names.index('r')
    def observer(rows):
        value = signature(rows,1 << root_index)
        require(value[0] % 2 == 1 and value[2] >= value[0], 'actual intermediate root observer')
        return {'H':value[0],'F_root':value[2],'defect':value[2]-value[0]}
    stage_observers = [observer(source)]
    for index,(a,b) in enumerate(pairs):
        U,V,W,p,q,case,inverse,block_word = compile_pair(current,a,b)
        require(all(0 < row[-1].bit_count() < len(names) and
                    min(row[-1].bit_count(),len(names)-row[-1].bit_count()) <= 2
                    for row in (U,W)), 'both helper endpoint cuts are nontrivial and sparse')
        helper = f'helper{index}'
        local = names+[helper]
        graph = to_sets(current,names)
        out_names = {names[j] for j in range(len(names)) if (U[-1] >> j) & 1}
        graph[helper] = set(out_names)
        for name in names:
            if name not in out_names: graph[name].add(helper)
        require(to_bits(graph,local) == U, 'independent helper-cut reconstruction')
        first = [local[i] for i in p]; second = [local[i] for i in q]
        graph = native_sets(graph,first)
        require(to_bits(graph,local) == V, 'independent first native state')
        graph = native_sets(graph,second)
        require(to_bits(graph,local) == W, 'independent second native state')
        del graph[helper]
        for name in names: graph[name].discard(helper)
        after = to_bits(graph,names)
        require(after == flip(current,a,b), 'independent helper-only deletion endpoint')
        require(set(graph) == set(names), 'every old physical name survives')
        flags = [strong(U),strong(V),strong(W)]
        native_observers = [observer(state) for state in (U,V,W)]
        native_strong.extend(flags)
        stage_strong.append(strong(after))
        stage_observers.append(observer(after))
        program.append({'old_pair':[names[a],names[b]],'helper':helper,
            'helper_out_names':sorted(out_names),'case':case,'inverse':inverse,
            'pair_block_word':[names[i] for i in block_word],
            'guards':[first,second],'birth_rows':U,'after_first_rows':V,
            'after_second_rows':W,'after_helper_deletion_rows':after,
            'native_stage_strongness':flags,'native_root_observers':native_observers})
        # A legal guard's reverse and a repeated name are explicit hostile tests.
        initial = to_sets(U,local)
        for bad in (list(reversed(first)),first[:-1]+[first[0]]):
            try:
                native_sets(initial,bad)
            except ValueError:
                pass
            else:
                raise RuntimeError('invalid native guard was accepted')
        current = after
    require(current == target, 'entire literal target matrix')
    root = root_index
    source_obs = signature(source,1 << root)
    target_obs = signature(target,1 << root)
    return {'names':names,
            'source_rows':source,'target_rows':target,'source_H':source_obs[0],
            'target_H':target_obs[0],'source_root_F':source_obs[2],
            'target_root_F':target_obs[2],'program':program,
            'completed_adjacency_strongness':stage_strong,
            'completed_adjacency_root_observers':stage_observers,
            'native_stage_strongness':native_strong,
            'native_complements':2*len(pairs),'temporary_births':len(pairs),'helper_deletions':len(pairs),
            'permanent_births':0,'peak_order':len(names)+1}


CONFIG = [{'sourceH': 4291, 'gap': 1, 'mapping': [9, 11, 12, 10, 2, 1, 5, 4, 6, 3, 0, 8, 7], 'padding': 0, 'core': [2, 4, 1], 'ports': [0], 'word': [0, 0, 0, 0, 1, 1, 1, 1], 'pairs': [[1, 9], [2, 8], [7, 12], [8, 10], [9, 11]], 'aux': 7, 'target': [7934, 8068, 8064, 7926, 7846, 8070, 7862, 256, 6233, 4480, 7040, 640, 2176]}, {'sourceH': 4291, 'gap': 2, 'mapping': [9, 10, 2, 11, 12, 1, 5, 4, 6, 3, 0, 8, 7], 'padding': 0, 'core': [2, 12, 9, 17, 7], 'ports': [0], 'word': [0, 0, 0, 0, 0, 1], 'pairs': [[1, 9], [2, 9], [2, 11], [3, 8], [4, 8], [6, 8], [7, 12], [10, 11]], 'aux': 7, 'target': [7934, 8068, 5248, 8182, 8102, 8070, 8118, 256, 7173, 6532, 4736, 1156, 2176]}, {'sourceH': 4289, 'gap': 1, 'mapping': [11, 12, 10, 7, 2, 9, 1, 4, 5, 6, 3, 0, 8, 13], 'padding': 0, 'core': [0], 'ports': [0], 'word': [1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1], 'pairs': [[2, 9], [6, 8], [8, 13], [12, 13]], 'aux': 13, 'target': [16126, 16004, 15488, 16374, 16006, 16022, 16310, 15616, 4151, 15748, 14592, 8448, 10240, 256]}, {'sourceH': 4289, 'gap': 2, 'mapping': [3, 13, 10, 7, 12, 9, 2, 1, 4, 5, 6, 0, 8, 11], 'padding': 1, 'core': [2, 4, 9, 3], 'ports': [1, 2], 'word': [0, 1, 0, 1, 1, 0, 1], 'pairs': [[1, 8], [6, 8], [7, 12], [7, 13], [11, 13], [0, 3]], 'aux': 11, 'target': [16118, 16260, 16000, 16375, 16006, 16022, 16310, 3328, 12341, 15744, 14592, 256, 2176, 6272]}]


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    base=Path(__file__).parent
    for name,pin in [('surgery_ready.py','024ac0b483d832865c6f12eabf0fe0f8676ac39188631e44910a0a4c9e66ca2f'),
                     ('surgery_native.py','66eb9d682e0128bdf087c79b24a44dcfe127daf810ec560de85daa9689e220e2')]:
        require(sha256((base/name).read_bytes()).hexdigest()==pin,'frozen first-round dependency '+name)
    sources={}
    for gap in (1,2):
        old,source,target,r,a,pairs=native.surgery(3,gap)
        child=native.c2(target,a)[1]
        H=signature(child,1<<r)[0]
        sources[H]=child
    cases=[]
    for cfg in CONFIG:
        source=sources[cfg['sourceH']]
        original,r,a=ready.ready_target(tuple(cfg['core']),tuple(cfg['ports']),tuple(cfg['word']))
        padded=prepend_chain(original,cfg['padding'])
        r+=cfg['padding'];a+=cfg['padding']
        require(cfg['mapping'][r]==8 and cfg['mapping'][a]==cfg['aux'], 'physical retained root and auxiliary map')
        target=relabel(padded,cfg['mapping'])
        require(target==tuple(cfg['target']), 'independent literal target reconstruction')
        require(len(target)==len(source), 'no old-name deletion and no preparation birth')
        names=[f'v{i}' for i in range(len(source))];names[8]='r'
        rec=execute(names,source,target,[tuple(e) for e in cfg['pairs']])
        require(all(rec['native_stage_strongness']), 'every temporary-helper native state strong')
        expected_strong=cfg['padding']==0
        require(strong(target)==expected_strong, 'strong or explicitly weak target scope')
        if expected_strong:
            require(all(rec['completed_adjacency_strongness']), 'strong branch every adjacency stage')
        else:
            require(all(rec['completed_adjacency_strongness'][:-1]), 'weak branch loses strongness only at final helper deletion')
            prefix=cfg['mapping'][0]
            require(target[prefix]==((1<<len(target))-1)^(1<<prefix), 'literal global source prefix')
        obs=native.joint_readout(target,8,cfg['aux'])
        require(obs['H']==cfg['sourceH'] and obs['sigma']==cfg['gap'] and obs['tau']==0,
                'actual neutral second-round ready response')
        n=len(target)
        first=native.birth(target,((1<<n)-1)^(1<<cfg['aux']))
        child=native.birth(first,(1<<n)-1)
        child_sig=signature(child,1<<8)
        require(child_sig[0]==3*cfg['sourceH']+2*cfg['gap'] and child_sig[2]==3*obs['F'],
                'actual named child with its changed physical root response')
        require(strong(first)==expected_strong and strong(child)==expected_strong,
                'no hidden strongness promotion after C2')
        for graph in (source,target,first,child):
            f=ready.counts(graph,(8,));b=signature(graph,1<<8)
            require((f['H'],f['F'][8])==(b[0],b[2]), 'independent forward/backward observer control')
        rec.update({'source_H':cfg['sourceH'],'gap':cfg['gap'],'padding':cfg['padding'],
                    'construction':cfg,'root_index':8,'aux_index':cfg['aux'],'joint':obs,
                    'target_strong':expected_strong,'C2_first_rows':first,'child_rows':child,
                    'child_H':child_sig[0],'child_F':child_sig[2],
                    'child_defect':child_sig[2]-child_sig[0],
                    'child_strong':expected_strong,'total_permanent_births':2,'full_peak_order':len(child)})
        cases.append(rec)
    total=GATES+native.GATES+ready.GATES
    result={'scope':'actual finite two-round binary tree; three second-round strong leaves and one explicitly weak leaf',
            'cases':cases,'gates':total,'new_checks':GATES,'compiler_checks':native.GATES,
            'forward_checks':ready.GATES,
            'native_complements':sum(c['native_complements'] for c in cases),
            'helper_only_deletions':sum(c['helper_deletions'] for c in cases),
            'strong_native_states':sum(len(c['native_stage_strongness']) for c in cases)}
    raw=(json.dumps(result,indent=2,sort_keys=True)+'\n').encode();args.output.write_bytes(raw)
    print(json.dumps({'gates':total,'cases':len(cases),'native_complements':result['native_complements'],
                      'strong_native_states':result['strong_native_states'],'receipt_bytes':len(raw),
                      'receipt_sha256':sha256(raw).hexdigest()},sort_keys=True))


if __name__=='__main__': main()
