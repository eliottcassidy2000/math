"""Bounded source-reselection tests; no global absence inference."""
from itertools import combinations
from pathlib import Path
import csv
import hashlib
import json
import sys
import check as f

HERE = Path(__file__).resolve().parent
CATALOG = f.ROOT/'05-knowledge/results/actual_hamilton_multiplicities_20261021/session24-multiplicities/n8_classes.csv'
SHIFTED = HERE.parent/'frames/search_receipt.json'
R77 = (22,124,88,113,96,69,1)
R79 = (64,17,67,7,13,31,58)
R79B = (4,5,80,7,43,15,187,63)
Q_SOURCE = (16,5,9,35,46,7)


def difference(a,b):
    return {k:a.get(k,0)-b.get(k,0) for k in a.keys() | b.keys()
            if a.get(k,0) != b.get(k,0)}


def two_birth(rows,root):
    n = len(rows)
    return tuple(r | (1 << n) | ((1 << (n+1)) if i != root else 0)
                 for i,r in enumerate(rows)) + (1 << (n+1),1 << root)


def ready_charts():
    result = []
    for rows,root,port,h,b in ((R77,6,1,77,45),(R79,0,5,79,41),(R79B,2,0,79,40)):
        ordinary = f.response(rows,(port,))
        marked = f.marked_response(rows,port,root)
        f.need(ordinary == {(1,):h,(2,):b},'ready ordinary response')
        f.need(difference(marked,ordinary) == {(1,):2},'ready marked response')
        after = two_birth(rows,root)
        after_poly = f.response(after,(port,))
        f.need(after_poly == {(1,):3*h+4,(2,):3*b},'ready generic extension')
        root_defects = [sum(f.marked_parity(after,r)[1:])-f.counts(after)[-1]
                        for r in range(len(after))]
        f.need(2 not in root_defects,'this output must not silently rearm')
        plus2_cuts = []
        for cut in range(1 << len(rows)):
            child = tuple(row | ((1 << len(rows)) if not (cut >> i) & 1 else 0)
                          for i,row in enumerate(rows)) + (cut,)
            if f.counts(child)[-1] != 3*h+2:
                continue
            if f.response(child,(port,)) == {(1,):3*h+2,(2,):3*b}:
                plus2_cuts.append(cut)
        result.append({'rows':rows,'root':root,'module_port':port,'H_coefficients':[h,b],
                       'G_coefficients':[h+2,b],'plus4_rows':after,'plus4_root_defects':root_defects,
                       'one_birth_plus2_cuts':plus2_cuts})
    return result


def catalogue_audit():
    wanted = {37:27,91:54,199:108,149:72,151:72,293:144,295:144,581:288,583:288}
    include = set(wanted) | {77,79,233,235,241}
    classes = {h:0 for h in include}
    scalar_ready = []
    matched = []
    ready_two_run = []
    for item in csv.DictReader(CATALOG.open(encoding='utf-8')):
        h = int(item['H'])
        if h not in include:
            continue
        classes[h] += 1
        code = int(item['code'])
        rows = tuple((code >> (8*i)) & 255 for i in range(8))
        f.need(f.counts(rows)[-1] == h,'catalogue H control')
        for root in range(8):
            _,even,odd = f.marked_parity(rows,root)
            if even+odd-h != 2:
                continue
            scalar_ready.append({'H':h,'root':root,'rows':rows})
            for port in range(8):
                ordinary = f.response(rows,(port,))
                if h in wanted and ordinary == {(1,):h,(2,):wanted[h]}:
                    # Include same-root module slots in the initial necessary
                    # response test. No such slot survives in this audit.
                    matched.append({'H':h,'root':root,'port':port,'rows':rows})
                if port == root or any(k[0] > 2 for k in ordinary):
                    continue
                marked = f.marked_response(rows,port,root)
                if difference(marked,ordinary) == {(1,):2}:
                    ready_two_run.append({'H':h,'a2':ordinary.get((2,),0),'root':root,'port':port})
    f.need(not matched,'unexpected shifted-ready candidate')
    return {'scope':'complete inherited ordinary order8 bank, selected H values; lower frames embed by neutral prefix',
            'class_counts':classes,'scalar_ready_roots':scalar_ready,
            'required_slopes':wanted,'matched_shifted_ready_candidates':matched,
            'all_two_run_ready_slots_in_selected_heads':ready_two_run}


def local_241_repairs():
    target = two_birth(R79,0)
    port = 5
    ordinary = f.response(target,(port,))
    f.need(ordinary == {(1,):241,(2,):123},'H241 target response')
    edges = [(i,j) for j in range(9) for i in range(j) if port not in (i,j)]
    result = []
    for radius in (1,2,3):
        scalar = 0
        marked_candidates = 0
        matched = []
        for flips in combinations(edges,radius):
            rows = list(target)
            for i,j in flips:
                rows[i] ^= 1 << j
                rows[j] ^= 1 << i
            if f.counts(rows)[-1] != 241:
                continue
            scalar += 1
            roots = []
            for root in range(9):
                if root == port:
                    continue
                h,e,o = f.marked_parity(rows,root)
                if e+o-h == 2:
                    roots.append(root)
            if not roots or f.response(tuple(rows),(port,)) != ordinary:
                continue
            marked_candidates += 1
            for root in roots:
                if difference(f.marked_response(rows,port,root),ordinary) == {(1,):2}:
                    matched.append({'flips':flips,'root':root,'rows':rows})
        f.need(not matched,'unexpected local rearming success')
        result.append({'radius':radius,'tested':sum(1 for _ in combinations(edges,radius)),
                       'scalar_H241_candidates':scalar,'full_response_scalar_ready_candidates':marked_candidates,
                       'universal_rearms':matched})
    return {'target_rows':target,'module_port':port,'allowed_edges':edges,'radii':result}


def marked_embedding(source,target,source_port,target_port,omitted):
    pool = [v for v in range(len(target)) if v != omitted]
    mask = sum(1 << v for v in pool)
    degree = [r.bit_count() for r in source]
    choices = {u:[v for v in pool if (target[v]&mask).bit_count() == degree[u]
                   and (u != source_port or v == target_port)] for u in range(len(source))}
    if any(not v for v in choices.values()):
        return None
    order = sorted(choices,key=lambda u:len(choices[u]))
    mapping,used = {},set()
    def visit(position):
        if position == len(order):
            return dict(mapping)
        u = order[position]
        for v in choices[u]:
            if v in used or any(((source[u] >> w) & 1) != ((target[v] >> z) & 1)
                                for w,z in mapping.items()):
                continue
            mapping[u] = v
            used.add(v)
            answer = visit(position+1)
            if answer is not None:
                return answer
            used.remove(v)
            del mapping[u]
        return None
    return visit(0)


def shifted_audit():
    hits = json.loads(SHIFTED.read_text(encoding='utf-8'))['hits']
    source,_ = f.substitute(Q_SOURCE,{0:f.transitive(3)})
    matches = []
    trials = 0
    for index,hit in enumerate(hits):
        for omitted in range(9):
            if omitted == hit['root']:
                continue
            for clone in range(3):
                trials += 1
                mapping = marked_embedding(source,hit['rows'],clone,hit['root'],omitted)
                if mapping is not None:
                    matches.append({'presentation':index,'variable_clone':clone,'omitted':omitted,'map':mapping})
    f.need(len(hits)==18 and trials==432 and not matches,'marked shifted-target embedding census')
    literal = []
    for k in range(1,7):
        rows,_ = f.substitute(Q_SOURCE,{0:f.transitive(k)})
        h = f.counts(rows)[-1]
        defects = [sum(f.marked_parity(rows,r)[1:])-h for r in range(len(rows))]
        f.need(2 not in defects,'Q literal source root check')
        literal.append({'k':k,'G_minus_H':defects})
    split,_ = f.substitute(Q_SOURCE,{0:f.transitive(2)})
    f.need(difference(f.marked_response(Q_SOURCE,0,4),f.response(Q_SOURCE,(0,))) == {(1,):22},
           'fixed exterior root22 identity')
    f.need(difference(f.marked_response(split,0,1),f.response(split,(0,))) == {(1,):5},
           'last clone5 identity')
    return {'source_rows':source,'presentations':len(hits),'marked_cases':trials,'marked_embeddings':matches,
            'scope':'variable source slot is one of three clones; target marked slot must be retained',
            'Q_literal_root_defects':literal,'fixed_old_root4_defect_coefficient':22,
            'last_clone_defect_coefficient_nonempty_remainder':5}


def main():
    report = {'status':'FINITE-EXACT','ready_charts':ready_charts(),'catalogue':catalogue_audit(),
              'local_H241_repairs':local_241_repairs(),'shifted_Q_source':shifted_audit(),
              'source_pins':[{'path':str(p.relative_to(f.ROOT)).replace('\\','/'),
                              'sha256_lf':hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()}
                             for p in (CATALOG,SHIFTED)]}
    suffix = '.optimized' if sys.flags.optimize else ''
    (HERE/('rearming_receipt'+suffix+'.json')).write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':report['status'],'ready_charts':3,'catalogue_classes':sum(report['catalogue']['class_counts'].values()),
                      'Q_marked_embedding_cases':report['shifted_Q_source']['marked_cases'],
                      'H241_local_radii':report['local_H241_repairs']['radii']}))


if __name__ == '__main__':
    main()
