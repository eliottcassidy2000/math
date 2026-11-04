#!/usr/bin/env python3
"""Literal old-name-preserving six-flip rearming of actual H1429.

The pair-block compiler and set-valued native interpreter are copied from
its audited predecessor native_six.py. This file constructs NEW actual source
and target graphs. No predecessor receipt or runtime producer import is used.
Both physical targets retain root8 and all eleven old names.
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


def birth(rows,incoming):
    n = len(rows)
    return tuple(row | (((incoming >> i)&1) << n)
                 for i,row in enumerate(rows))+(((1 << n)-1)^incoming,)


def c2(rows,aux):
    n = len(rows)
    first = birth(rows,((1 << n)-1)^(1 << aux))
    final = birth(first,(1 << n)-1)
    require(strong(first) and strong(final), 'both C2 permanent births strong')
    return first,final


def prepared(clock):
    k = clock-2
    require(k>=1, 'legal original clock')
    D = (2,4,57,3,11,27)
    suffix = tuple(row|64 for row in D)+(0,)
    parent = attachment(suffix,80,k)
    root,aux = k+7,k+6
    out = c2(parent,aux)[1]
    require(strong(out), 'literal PREPARED_FOUR output')
    return out,root


def surgery(clock,gap):
    old,root = prepared(clock)
    k = clock-2
    if gap==2:
        source = old
        aux = k+6
        pairs = [(k,k+7),(k,k+8),(k+1,k+2),(k+3,k+4),(k+6,k+9),(k+7,k+8)]
    else:
        require(gap==1, 'two declared surgery branches')
        aux = len(old)
        source = birth(old,((1 << len(old))-1)^(1 << root))
        pairs = [(k+1,k+2),(k+2,k+7),(k+4,k+7),(k+7,k+8),(k+7,k+9),(k+8,k+9)]
        require(strong(source), 'one permanent preparatory birth strong')
    target = source
    for a,b in pairs:
        target = flip(target,a,b)
    return old,source,target,root,aux,pairs


def induced(rows,vertices):
    return tuple(sum(((rows[i] >> j)&1) << b for b,j in enumerate(vertices))
                 for i in vertices)


def joint_readout(rows,root,aux):
    require(root!=aux and rows[aux]==1 << root, 'sole-out-neighbour actual joint gate')
    whole = signature(rows,1 << root)
    start = signature(rows,1 << aux)[3]
    keep = [i for i in range(len(rows)) if i!=aux]
    deleted = induced(rows,keep)
    L = signature(deleted,1 << keep.index(root))[2]
    keep2 = [i for i in keep if i!=root]
    K = signature(induced(rows,keep2),0)[1]
    require(L+K==whole[2], 'independent deleted-graph L+K=F identity')
    return {'H':whole[0],'F':whole[2],'S':start,'L':L,'K':K,'G':L+K,
            'sigma':start-whole[0],'tau':L+K-whole[2],
            'root_defect':whole[2]-whole[0]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    cases=[]
    for gap in (2,1):
        old,source,target,root,aux,pairs = surgery(3,gap)
        names = [f'v{i}' for i in range(len(source))]
        names[root]='r'
        rec = execute(names,source,target,pairs)
        require(all(rec['native_stage_strongness']), 'every actual native intermediate strong')
        initial = signature(old,1 << root)
        require(initial[0]==1429 and initial[2]==1563 and root==8,
                'old actual output and physical root, never helper')
        obs = joint_readout(target,root,aux)
        require(obs['H']==1429 and obs['sigma']==gap and obs['tau']==0,
                'count-neutral ready target with exact joint response')
        require(obs['F']==(1493 if gap==2 else 1923), 'changed physical root F')
        first,child=c2(target,aux)
        child_sig=signature(child,1 << root)
        require(child_sig[0]==3*1429+2*gap and child_sig[2]==3*obs['F'],
                'actual nonzero ternary child and transported root')
        # Full forward DP is an independent observer implementation.
        import surgery_ready as sr
        for graph in (old,source,target,first,child):
            forward=sr.counts(graph,(root,))
            backward=signature(graph,1 << root)
            require((forward['H'],forward['F'][root])==(backward[0],backward[2]),
                    'forward/backward actual named observer agreement')
        rec.update({'gap':gap,'root_index':root,'aux_index':aux,
                    'original_rows':old,'original_H':initial[0],'original_F':initial[2],
                    'permanent_preparation_births':int(gap==1),
                    'joint':obs,'C2_first_rows':first,'child_rows':child,
                    'child_H':child_sig[0],'child_F':child_sig[2],
                    'child_defect':child_sig[2]-child_sig[0],
                    'total_permanent_births':2+int(gap==1),
                    'full_peak_order':len(child)})
        cases.append(rec)
    clock_controls=[]
    for clock in range(3,7):
        for gap in (2,1):
            old,source,target,root,aux,pairs=surgery(clock,gap)
            k=clock-2
            obs=joint_readout(target,root,aux)
            expected=(744*2**k-59,744*2**k+5) if gap==2 else (961*2**k-493,961*2**k+1)
            require((obs['H'],obs['F'])==expected and obs['sigma']==gap,
                    'actual clock control of analytic response formula')
            current=source
            for a,b in pairs:
                current=flip(current,a,b)
                require(strong(current), 'actual clock macrostep strong')
            # Recover the finite suffix and its actual incoming-root mask.
            suffix_vertices=[i for i in range(k,len(target)) if i!=root]
            suffix=induced(target,suffix_vertices)
            ports=sum(((target[i]>>root)&1)<<j for j,i in enumerate(suffix_vertices))
            sig=signature(suffix,ports)
            require((sig[0],sig[2],sig[4])==((5,744,64) if gap==2 else (1,961,494)),
                    'finite suffix h,f,j supplier, not formula fitting')
            clock_controls.append({'clock':clock,'gap':gap,'joint':obs,
                                   'suffix_rows':suffix,'suffix_ports':ports,'signature':sig})
    result={'scope':'two actual H1429 all-strong native preparations; exact all-clock adjacency response theorem',
            'cases':cases,'clock_controls':clock_controls,'gates':GATES,
            'native_complements':sum(c['native_complements'] for c in cases),
            'helper_only_deletions':sum(c['helper_deletions'] for c in cases),
            'native_strong_states':sum(len(c['native_stage_strongness']) for c in cases)}
    raw=(json.dumps(result,indent=2,sort_keys=True)+'\n').encode()
    args.output.write_bytes(raw)
    print(json.dumps({'gates':GATES,'cases':len(cases),'native_complements':result['native_complements'],
                      'strong_native_states':result['native_strong_states'],
                      'receipt_bytes':len(raw),'receipt_sha256':sha256(raw).hexdigest()},sort_keys=True))


if __name__=='__main__':
    main()

