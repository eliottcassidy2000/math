"""Compile named component changes into legal native complements.

The old one-helper theorem is used algorithmically; the saved program is
checked by replaying every guard and every changed old arc. H is not held
fixed during a source-changing excursion.
"""
from pathlib import Path
from functools import lru_cache
from itertools import combinations
import hashlib
import importlib.util
import json
from reader import need, birth

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
SOURCE = ROOT/'05-knowledge/results/adaptive_metagraph_reformulation_20261103/metagraph_probe.py'
spec = importlib.util.spec_from_file_location('inherited_native', SOURCE)
old = importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)


def traced_path(rows, u, v):
    """Subset search with the requested u->v adjacency enforced."""
    full = (1 << len(rows))-1
    @lru_cache(None)
    def walk(mask, last):
        if mask == full:
            return ()
        if last == u:
            available = (1 << v) & ~mask
        else:
            available = (full if last < 0 else rows[last]) & ~mask & ~(1 << v)
        while available:
            bit = available & -available
            available -= bit
            w = bit.bit_length()-1
            # Do not strand the forced successor when entering u.
            if w == u and mask >> v & 1:
                continue
            tail = walk(mask | bit, w)
            if tail is not None:
                return (w,)+tail
        return None
    return walk(0, -1)


def compile_arc(rows, a, b):
    u, v = (a, b) if rows[a] >> b & 1 else (b, a)
    path = traced_path(rows, u, v)
    if path is not None:
        U, V, W, p, q, case = old.traced_gadget(rows, u, v, path)
        return U, V, W, p, q, case, False
    flipped = old.flip(rows, u, v)
    path = traced_path(flipped, v, u)
    need(path is not None, 'One orientation has a traced Hamilton word')
    U, V, W, p, q, case = old.traced_gadget(flipped, v, u, path)
    return W, V, U, q, p, case, True


def main():
    data = json.loads((HERE/'named_targets.json').read_text())
    programs = []
    for target in data['targets']:
        names = list(data['source_names'])
        rows = tuple(data['source_rows'])
        final_names, final_rows = target['names'], target['rows']
        loc = {v:i for i,v in enumerate(final_names)}
        permanent_births = []
        for new in target['new_names']:
            out = sum(1 << i for i,v in enumerate(names) if final_rows[loc[new]] >> loc[v] & 1)
            rows = birth(rows, out)
            names.append(new)
            permanent_births.append(dict(name=new,out_mask=out,rows=rows))
        need(names == final_names, 'Permanent physical-name order')
        moves = []
        for number, (a,b) in enumerate(target['old_arc_differences']):
            ai, bi = names.index(a), names.index(b)
            U,V,W,p,q,case,inverse = compile_arc(rows, ai, bi)
            helper = 1000+100*target['offset']+number
            local = names+[helper]
            need(old.native(U,p) == V and old.native(V,q) == W, 'Legal literal replay')
            after = tuple(x & ((1 << len(rows))-1) for x in W[:-1])
            need(after == old.flip(rows,ai,bi), 'Only requested old arc differs')
            moves.append(dict(changed_pair=[a,b],helper=helper,
                              helper_out_mask=U[-1],birth_rows=U,
                              paths=[[local[v] for v in path] for path in (p,q)],
                              native_rows=[V,W],after_delete=after,
                              case=case,inverse=inverse))
            rows = after
        need(rows == tuple(final_rows), 'Exact named endpoint')
        programs.append(dict(target_H=target['H'],permanent_births=permanent_births,
                             native_moves=2*len(moves),temporary_births=len(moves),
                             deleted_old_names=[],peak_order=len(names)+1,
                             changes=moves,final_names=names,final_rows=rows,
                             preserved_at_macro_endpoints={'prefix':[50,51],
                                'root':39,'component_mark':31,'spare':6}))
        print('Compiled',target['H'],2*len(moves),'native moves;',len(names),'endpoint vertices')
    report=dict(status='FINITE_EXACT_NATIVE_PROGRAMS',programs=programs,
                source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
                target_graph_sha256=hashlib.sha256((HERE/'named_targets.json').read_bytes()).hexdigest(),
                scope='Every old name survives. Old arcs outside B agree at each completed helper gadget. Native intermediates need not preserve H, strongness or port stars.')
    (HERE/'native.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__ == '__main__':
    main()
