"""Compile the actual fourth +2 repair, then replay every native guard."""
from pathlib import Path
import hashlib, json, sys

HERE = Path(__file__).resolve().parent
PACKET = HERE.parent
RESULTS = PACKET.parent
OLD = RESULTS / 'zero_ready_component_returns_20261103'
sys.path.insert(0, str(OLD/'continuation'))
from native_tree import compile_both, digest, replay
import native
from reader import need


def main():
    path = PACKET/'selector/named_repair.json'
    data = json.loads(path.read_text())
    source, target = data['source'], data['target']
    names = source['names']; rows = tuple(source['rows'])
    need(names == target['names'] and not target['new_names'], 'No permanent name changes')
    saved = HERE/'native.json'
    if saved.exists():
        program = json.loads(saved.read_text())['program']
    else:
        changes = []
        for j, (a, b) in enumerate(target['old_arc_differences']):
            ai, bi = names.index(a), names.index(b)
            U,V,W,p,q,case,inverse = compile_both(rows, ai, bi)
            need(native.old.native(U,p)==V and native.old.native(V,q)==W, 'Literal guards')
            after = tuple(x & ((1 << len(names))-1) for x in W[:-1])
            need(after == native.old.flip(rows,ai,bi), 'Exactly the intended old arc')
            helper = 30000+j; local = names+[helper]
            changes.append(dict(changed_pair=[a,b], helper=helper, helper_out_mask=U[-1],
                paths=[[local[i] for i in word] for word in (p,q)], case=case, inverse=inverse,
                after_delete_sha256=digest(after)))
            rows = after
        need(rows == tuple(target['rows']), 'Actual full target matrix')
        program = dict(parent_H=source['H'], target_H=target['H'], branch='2222',
                       permanent_births=[], changes=changes)
    receipt = replay(source, program, target)
    pins = [path, Path(__file__), OLD/'continuation/native_tree.py', OLD/'inverse/native.py',
            RESULTS/'adaptive_metagraph_reformulation_20261103/metagraph_probe.py']
    report = dict(status='FINITE-EXACT', program=program, receipt=receipt,
        scope='One actual named fourth +2 edge. All 23 old names and the four-vertex prefix survive. No intermediate H or strongness invariant claimed.',
        source_sha256={str(p.relative_to(RESULTS)): hashlib.sha256(p.read_bytes()).hexdigest() for p in pins})
    saved.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(receipt,sort_keys=True))


if __name__ == '__main__':
    main()
