#!/usr/bin/env python3
"""Literal old-name-preserving native compiler for the SCC6 +4 repair.

The general compiler is inherited (citations in component/NATIVE_SIX.md).
This self-contained implementation and separate set-valued interpreter record
all guards at clocks1,2,3,4,5, including two finite boundary patches. No producer
module is imported.
"""
import argparse
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
import json

GATES = 0


def require(ok, message):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(message)


def flip(rows,a,b):
    out = list(rows)
    out[a] ^= 1 << b
    out[b] ^= 1 << a
    return tuple(out)


def strong(rows):
    reach = [row | (1 << v) for v,row in enumerate(rows)]
    for k in range(len(rows)):
        for v in range(len(rows)):
            if (reach[v] >> k) & 1:
                reach[v] |= reach[k]
    return all(r == (1 << len(rows))-1 for r in reach)


def native_bits(rows,path):
    n = len(rows)
    require(len(path) == n and set(path) == set(range(n)), 'compiler word names')
    require(all((rows[a] >> b) & 1 for a,b in zip(path,path[1:])), 'compiler word directed')
    keep = {frozenset((a,b)) for a,b in zip(path,path[1:])}
    out = rows
    for a in range(n):
        for b in range(a+1,n):
            if frozenset((a,b)) not in keep:
                out = flip(out,a,b)
    return out


def pair_word(rows,a,b):
    # The inherited polynomial pair-block construction, not Hamilton search.
    word = []
    for v in range(len(rows)):
        if v in (a,b): continue
        pos = next((i for i,w in enumerate(word) if (rows[v] >> w) & 1),len(word))
        word.insert(pos,v)
    pos = next((i for i,w in enumerate(word)
                if ((rows[a] | rows[b]) >> w) & 1),len(word))
    if pos == len(word) or (rows[b] >> word[pos]) & 1:
        block = [a,b]
    else:
        block = [b,a]
    word[pos:pos] = block
    u,v = block
    oriented = rows if (rows[u] >> v) & 1 else flip(rows,a,b)
    require(all((oriented[x] >> y) & 1 for x,y in zip(word,word[1:])),
            'pair-block path in one of the two orientations')
    return tuple(word),oriented


def compile_pair(rows,a,b):
    n = len(rows)
    require(n >= 3 and a != b, 'traced recipe domain')
    path,base = pair_word(rows,a,b)
    at = next(i for i in range(n-1) if {path[i],path[i+1]} == {a,b})
    u,v = path[at:at+2]
    if v != path[-1]:
        helper_out = 1 << u
        first = path+(n,)
        case = 'append'
    else:
        helper_out = ((1 << n)-1) ^ (1 << v)
        first = (n,)+path
        case = 'prepend'
    U = tuple(row | (0 if (helper_out >> i) & 1 else 1 << n)
              for i,row in enumerate(base))+(helper_out,)
    second = path[:at+1]+(n,)+path[at+1:]
    V = native_bits(U,first)
    W = native_bits(V,second)
    old_mask = (1 << n)-1
    require(tuple(row & old_mask for row in W[:-1]) == flip(base,a,b),
            'traced single old-pair identity')
    inverse = base != rows
    if inverse:
        U,W = W,U
        first,second = second,first
    require(tuple(row & old_mask for row in U[:-1]) == rows,
            'actual helper birth preserves every old pair')
    require(native_bits(U,first) == V and native_bits(V,second) == W,
            'actual orientation guarded replay')
    return U,V,W,first,second,case,inverse,path


def to_sets(rows,names):
    return {name:{names[j] for j in range(len(names)) if (rows[i] >> j) & 1}
            for i,name in enumerate(names)}


def to_bits(graph,names):
    return tuple(sum(1 << j for j,name in enumerate(names) if name in graph[v])
                 for v in names)


def native_sets(graph,path):
    # Independent elementary interpreter: labeled directed adjacency sets.
    if len(path) != len(graph) or len(set(path)) != len(path) or set(path) != set(graph):
        raise ValueError('invalid physical name list')
    if any(b not in graph[a] for a,b in zip(path,path[1:])):
        raise ValueError('non-directed guard')
    protected = {frozenset((a,b)) for a,b in zip(path,path[1:])}
    result = {v:set() for v in graph}
    for a in graph:
        for b in graph[a]:
            if frozenset((a,b)) in protected:
                result[a].add(b)
            else:
                result[b].add(a)
    return result


@lru_cache(None)
def signature(rows,ports):
    n = len(rows)
    pred = tuple(sum(1 << u for u in range(n) if (rows[u] >> v) & 1)
                 for v in range(n))
    def end_table(incoming):
        out = [[0]*n for _ in range(1 << n)]
        for v in range(n): out[1 << v][v] = 1
        for mask in range(1,1 << n):
            if mask & (mask-1) == 0: continue
            vertices = mask
            while vertices:
                bit = vertices & -vertices; vertices -= bit
                v = bit.bit_length()-1
                old = mask ^ bit
                choices = old & incoming[v]
                while choices:
                    u = choices & -choices; choices -= u
                    out[mask][v] += out[old][u.bit_length()-1]
        return out
    end = end_table(pred)
    reverse = end_table(rows)
    h = [sum(v) for v in end]; h[0] = 1
    e = [sum(v[i] for i in range(n) if (ports >> i) & 1) for v in end]
    s = [sum(v[i] for i in range(n) if (ports >> i) & 1) for v in reverse]
    full = (1 << n)-1
    return (h[full],sum(h[m]*h[full^m] for m in range(full+1)),
            sum(e[m]*h[full^m] for m in range(full+1)),
            sum(h[m]*s[full^m] for m in range(full+1)),
            sum(e[m]*s[full^m] for m in range(full+1)),e[full],s[full])


def transitive(n):
    return tuple(sum(1 << j for j in range(i+1,n)) for i in range(n))


def attachment(suffix,ports,t):
    n = len(suffix)
    root = t+n
    full_suffix = ((1 << n)-1) << t
    rows = [v | full_suffix for v in transitive(t)]
    rows += [(row << t) | ((1 << root) if (ports >> i) & 1 else 0)
             for i,row in enumerate(suffix)]
    rows += [((1 << t)-1) | (full_suffix ^ (ports << t))]
    return tuple(rows)


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
    require(all(stage_strong), 'all completed adjacency stages strong')
    root = len(names)-1
    source_obs = signature(source,1 << root)
    target_obs = signature(target,1 << root)
    require(target_obs[0] == 3*source_obs[0]+4, 'actual uniform +4 count')
    return {'names':names,
            'source_rows':source,'target_rows':target,'source_H':source_obs[0],
            'target_H':target_obs[0],'source_root_F':source_obs[2],
            'target_root_F':target_obs[2],'program':program,
            'completed_adjacency_strongness':stage_strong,
            'completed_adjacency_root_observers':stage_observers,
            'native_stage_strongness':native_strong,
            'native_complements':2*len(pairs),'temporary_births':len(pairs),'helper_deletions':len(pairs),
            'permanent_births':0,'peak_order':len(names)+1}


def one_clock(t,target_suffix):
    names = [f'p{i}' for i in range(t)]+[f'b{i}' for i in range(5)]+['r']
    source = attachment(transitive(5),21,t)
    target = attachment(target_suffix,80,t-2)
    require(len(names) == t+6, 'bulk old order')
    require(attachment(transitive(7),84,t-2) == source,
            'two old prefix names only change frame')
    changes = [(i,j) for i in range(6) for j in range(i+1,6)
               if not ((target_suffix[i] >> j) & 1)]
    require(len(changes) == 10, 'ten suffix old-pair discrepancies')
    pairs = [(t-2+i,t-2+j) for i,j in changes]+[(t,t+5)]
    require(names[pairs[-1][0]] == 'b0' and names[pairs[-1][1]] == 'r', 'root-port pair')
    out = execute(names,source,target,pairs)
    require(out['source_H'] == 21*(1 << t)-11, 'actual bulk source count')
    require(out['source_root_F']-out['source_H'] == 12 and
            out['target_root_F']-out['target_H'] == 46, 'actual bulk root defect')
    out.update(clock=t,case='bulk',source_prefix=names[:t],
               target_prefix=names[:t-2],target_suffix_names=names[t-2:-1])
    return out


def boundary(t):
    require(t in (1,2), 'finite boundary parameter')
    suffix,ports,expected = ((18,8,3,5,14),19,(13,104,55,58,26,8,9)) if t == 1 else \
                           ((56,33,35,38,46,0),33,(11,200,113,37,14,11,1))
    require(signature(suffix,ports) == expected, 'inherited scalar boundary actual suffix')
    names = [f'p{i}' for i in range(t)]+[f'b{i}' for i in range(5)]+['r']
    source = attachment(transitive(5),21,t)
    target = attachment(suffix,ports,1)
    pairs = [(i,j) for i in range(len(source)) for j in range(i+1,len(source))
             if ((source[i] ^ target[i]) >> j) & 1]
    require(len(pairs) == (9 if t == 1 else 11), 'finite boundary discrepancy count')
    out = execute(names,source,target,pairs)
    require((out['source_H'],out['target_H']) == ((31,97) if t == 1 else (73,223)),
            'finite boundary actual counts')
    require(out['source_root_F']-out['source_H'] == 12 and
            out['target_root_F']-out['target_H'] == expected[4], 'finite boundary root defects')
    out.update(clock=t,case='finite_patch',source_prefix=names[:t],
               target_prefix=names[:1],target_suffix_names=names[1:-1])
    return out


def compressed_patterns(case):
    """Finite controls for the analytically proved prefix-block invariance."""
    prefix = case['target_prefix']
    k = len(prefix)
    mapping = {v:f's{i}' for i,v in enumerate(case['target_suffix_names'])}
    mapping['r'] = 'r'
    patterns = []
    for step in case['program']:
        mapping[step['helper']] = 'helper'
        words = []
        for word in step['guards']:
            at = word.index(prefix[0])
            require(word[at:at+k] == prefix, 'literal consecutive prefix guard block')
            require(at > 0 and at+k < len(word), 'both guard endpoints outside prefix')
            words.append([mapping[v] for v in word[:at]]+['PREFIX']+
                         [mapping[v] for v in word[at+k:]])
        out = set(step['helper_out_names'])
        require(out.intersection(prefix) in (set(),set(prefix)), 'uniform helper cut on prefix')
        V = step['after_first_rows']
        rmask = (1 << k)-1
        for i in range(k):
            expected = (1 << (i+1) if i+1 < k else 0) + sum(1 << j for j in range(i-1))
            require(V[i] & rmask == expected, 'first-native reverse-nonconsecutive prefix')
        outside = tuple(row >> k for row in V[k:])
        pattern = {'guards':words,'case':step['case'],'inverse':step['inverse'],
                   'helper_beats_prefix':bool(out.intersection(prefix)),
                   'helper_outside_out_names':sorted(mapping[v] for v in out-set(prefix)),
                   'first_native_outside_rows':outside}
        if k == 3:
            # Contract the strong R_3 block; the quotient is a directed graph,
            # not necessarily a tournament because both directions can occur.
            q = [0]*10
            for i in range(9):
                q[i+1] = outside[i] << 1
                if V[k+i] & rmask: q[i+1] |= 1
                if any((V[j] >> (k+i)) & 1 for j in range(k)): q[0] |= 1 << (i+1)
            require(strong(tuple(q)), 'contracted k3 first-native directed quotient')
        patterns.append(pattern)
    return patterns


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=Path(__file__).with_suffix('.json'))
    args = parser.parse_args()
    core = (2,4,57,3,11,27)
    target = tuple(v | 64 for v in core)+(0,)
    require(strong(core), 'actual SCC6 core')
    require(target == (66,68,121,67,75,91,0), 'literal target suffix')
    require(signature(transitive(5),21) == (1,32,21,21,12,1,1), 'source suffix signature')
    require(signature(transitive(7),84) == (1,128,84,21,12,1,0), 'reframed source signature')
    require(signature(target,80) == (17,448,252,105,46,17,4), 'target suffix signature')
    cases = [boundary(1),boundary(2)]+[one_clock(t,target) for t in (3,4,5)]
    patterns = [compressed_patterns(case) for case in cases if case['case'] == 'bulk']
    require(patterns[0] == patterns[1] == patterns[2], 'same compressed bulk compiler at k1,k2,k3')
    empty_source = attachment(transitive(5),21,2)
    empty_target = attachment(target,80,0)
    hs = signature(empty_source,1 << (len(empty_source)-1))[0]
    ht = signature(empty_target,1 << (len(empty_target)-1))[0]
    require(hs == 73 and ht == 219 and ht != 3*hs+4, 'empty-prefix hostile')
    receipt = {'scope':'FINITE-EXACT native programs at t=1..5; all-t piecewise proof uses inherited compiler',
        'source_formula':'21*2^t-11','target_formula':'63*2^t-29',
        'bulk_graph_domain':'integers t>=3; two separate finite patches at t=1,2',
        'uniform_algorithm_domain':'all integers t>=1','core_rows':core,'suffix_rows':target,
        'target_port_mask':80,'cases':cases,'bulk_compressed_compiler_patterns':patterns[0],
        'empty_prefix_hostile':{'clock':2,'source_H':hs,'target_H':ht,'required_3H_plus4':3*hs+4},
        'always_active_gates':GATES,'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest()}
    content = json.dumps(receipt,sort_keys=True,indent=2)+'\n'
    args.output.write_text(content)
    print(json.dumps({'gates':GATES,'cases':[{'t':c['clock'],'H':(c['source_H'],c['target_H']),
        'strong_native_stages':sum(c['native_stage_strongness']),
        'native_stages':len(c['native_stage_strongness'])} for c in cases],
        'receipt_sha256':sha256(content.encode()).hexdigest()},sort_keys=True))


if __name__ == '__main__':
    main()
