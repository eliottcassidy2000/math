"""Actual named component replacements for the first two return branches."""
from pathlib import Path
from functools import lru_cache
import json
from reader import need, birth, strong

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
OLD = ROOT/'05-knowledge/results/triangle_return_170807_20261103'


def arrange(old_rows, old_names, candidates, mark, word, fresh):
    """Minimum old arc differences in this specified ordered-tail template."""
    loc = {v: i for i, v in enumerate(old_names)}
    count = len(candidates)
    full = (1 << count)-1
    need(count+len(fresh) == len(word), 'Template occupancy')
    @lru_cache(None)
    def dp(mask, used_new):
        position = mask.bit_count()+used_new
        if position == len(word):
            return 0, ()
        best = (10**9, ())
        if used_new < len(fresh):
            cost, suffix = dp(mask, used_new+1)
            best = cost, (fresh[used_new],)+suffix
        for j, v in enumerate(candidates):
            if mask >> j & 1:
                continue
            backwards = sum(bool(mask >> k & 1) and bool(old_rows[loc[v]] >> loc[u] & 1)
                            for k, u in enumerate(candidates))
            mismatch = bool(old_rows[loc[mark]] >> loc[v] & 1) != (word[position] == '1')
            cost, suffix = dp(mask | (1 << j), used_new)
            option = backwards+int(mismatch)+cost, (v,)+suffix
            if option < best:
                best = option
        return best
    return dp(0, 0)


def wrap(component, component_names, source_names):
    names = list(source_names)+[x for x in component_names if x not in source_names]
    edge = set()
    for i, a in enumerate(component_names):
        for j, b in enumerate(component_names):
            if component[i] >> j & 1:
                edge.add((a, b))
    edge.add((50, 51))
    for a in (50, 51):
        for b in component_names:
            edge.add((a, b))
    for v in names:
        if v != 39:
            edge.add((31, 39) if v == 31 else (39, v))
    rows = tuple(sum(1 << j for j, b in enumerate(names) if (a, b) in edge) for a in names)
    need(strong(rows) and strong(component), 'Strong component and wrapper')
    for i in range(len(rows)):
        for j in range(i):
            need(((rows[i] >> j) ^ (rows[j] >> i)) & 1, 'Tournament relation')
    return names, rows


def main():
    source = json.loads((OLD/'interface.json').read_text())
    old = json.loads((OLD/'zero_ready.json').read_text())[1]
    B, labels = tuple(old['suffix_rows']), old['suffix_names']
    outputs = []
    for offset in (2, 4):
        if offset == 4:
            word = '10101001010010100'
            cost, order = arrange(B, labels, [x for x in labels if x != 31], 31, word, [55,56,57,58])
            names = [31]+list(order)
            component = tuple([sum(1 << (i+1) for i, b in enumerate(word) if b == '1')]
                              + [sum(1 << j for j in range(i+2, len(names)))
                                 | (1 if bit == '0' else 0) for i, bit in enumerate(word)])
            feedback = [39,31]
            seed = {'rows':[2,0], 'root':0, 'word':word[1:]}
        else:
            word = '10111010010110'
            cost, order = arrange(B, labels, [x for x in labels if x not in (31,6,8)], 31, word, [56,57,58])
            names = [6,8,55,31]
            component = (2,12,9,1)  # 6->8->55->6; mark31 beats only6.
            seed = {'rows':list(component), 'root':3, 'word':word}
            for v, bit in zip(order, word):
                component = birth(component, 0 if bit == '1' else 1 << 3)
                names.append(v)
            feedback = [39,31,55]
        full_names, rows = wrap(component, names, source['target_names'])
        source_loc = {v:i for i,v in enumerate(source['target_names'])}
        target_loc = {v:i for i,v in enumerate(full_names)}
        differences = [[a,b] for i,a in enumerate(source['target_names'])
                       for b in source['target_names'][i+1:]
                       if ((source['target'][source_loc[a]] >> source_loc[b]) & 1)
                       != ((rows[target_loc[a]] >> target_loc[b]) & 1)]
        need(all(a in labels and b in labels for a,b in differences), 'Only old component arcs change')
        outputs.append(dict(offset=offset,H=3*170807+offset,names=full_names,rows=rows,
                            component_names=names,component_rows=component,mark=31,root=39,
                            spare=6,prefix=[50,51],structural_feedback_names=feedback,
                            new_names=[55,56,57,58],old_arc_differences=differences,
                            seed_reader=seed,template_tail_cost=cost,
                            cost_scope='Tail-order optimum inside the explicitly fixed seed template only.'))
    record=dict(status='FINITE_EXACT_NAMED_COMPONENT_RETURNS',
                source_names=source['target_names'],source_rows=source['target'],targets=outputs)
    (HERE/'named_targets.json').write_text(json.dumps(record,indent=2)+'\n')
    print([(x['H'],len(x['names']),len(x['old_arc_differences'])) for x in outputs])


if __name__ == '__main__':
    main()
