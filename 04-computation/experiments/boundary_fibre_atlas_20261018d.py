#!/usr/bin/env python3
"""Source-marked boundary response cells from bipartite ternary fibres.

This compresses existence witnesses in a declared bank, not all colourings.
All checks stay enabled with python -O. Standard library only.
"""
import argparse
from collections import deque
from itertools import combinations, product
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CHECKS = 0


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(message)


def nat(x, bound=None):
    need(type(x) is int and x >= 0 and (bound is None or x < bound), 'exact index')


def canonical(word):
    names = {}
    return tuple(names.setdefault(c, len(names)) for c in word)


def matchings(items):
    if not items:
        return [()]
    a = items[0]
    return [((a,b),)+tail for b in items[1:]
            for tail in matchings(tuple(x for x in items if x not in (a,b)))]


def audit_graph(n, edges, ports):
    nat(n); need(n > 0, 'nonempty graph')
    need(type(edges) is tuple and len(edges) == len(set(edges)), 'distinct edges')
    for edge in edges:
        need(type(edge) is tuple and len(edge) == 2, 'edge tuple')
        a,b = edge; nat(a,n); nat(b,n); need(a < b, 'ordered edge')
    need(type(ports) is tuple and len(set(ports)) == len(ports), 'distinct ports')
    for v in ports: nat(v,n)


def fibres(n, edges, labels):
    """Return every component and parity, or None for a monochromatic odd cycle."""
    need(type(labels) is tuple and len(labels) == n, 'complete ternary word')
    for a in labels: nat(a,3)
    adj = [[] for _ in range(n)]
    for a,b in edges:
        if labels[a] == labels[b]:
            adj[a].append(b); adj[b].append(a)
    component = [-1]*n; parity = [0]*n; number = 0
    for seed in range(n):
        if component[seed] >= 0: continue
        component[seed] = number; todo = deque([seed])
        while todo:
            v = todo.popleft()
            for w in adj[v]:
                if component[w] < 0:
                    component[w] = number; parity[w] = parity[v]^1; todo.append(w)
                elif parity[w] != parity[v]^1:
                    return None
        number += 1
    return tuple(component), tuple(parity), number


def boundary_cell(labels, data, ports):
    comp, parity, number = data
    visible = tuple(dict.fromkeys(comp[v] for v in ports))
    blocks = tuple(tuple(i for i,v in enumerate(ports) if comp[v] == c) for c in visible)
    offsets = tuple(tuple(parity[ports[i]]^parity[ports[block[0]]] for i in block) for block in blocks)
    signature = (tuple(labels[v] for v in ports), blocks, offsets)
    return signature, visible, number-len(visible)


def named_boundary_words(labels, data, ports):
    comp, parity, number = data
    visible = tuple(dict.fromkeys(comp[v] for v in ports))
    for bits in product(range(2), repeat=len(visible)):
        phases = dict(zip(visible,bits))
        word = tuple((3*(parity[v]^phases[comp[v]])+4*labels[v])%6 for v in ports)
        yield word, bits


def lift_named(labels, data, ports, boundary):
    """Same fixed ternary word; both coordinates of every pin are compulsory."""
    comp, parity, number = data
    need(type(boundary) is tuple and len(boundary) == len(ports), 'literal boundary')
    phases = {}
    for v,c in zip(ports,boundary):
        nat(c,6); need(c%3 == labels[v], 'pinned ternary label')
        phase = c%2 ^ parity[v]
        need(comp[v] not in phases or phases[comp[v]] == phase, 'pinned path parity')
        phases[comp[v]] = phase
    return tuple((3*(parity[v]^phases.get(comp[v],0))+4*labels[v])%6 for v in range(len(labels)))


def lift_orbit(labels, data, ports, target):
    """A global palette permutation is retained; it may change ternary labels."""
    for raw,bits in named_boundary_words(labels,data,ports):
        if canonical(raw) != canonical(target): continue
        mapping = dict(zip(raw,target))
        need(len(set(mapping.values())) == len(mapping), 'injective palette map')
        remaining = iter(c for c in range(6) if c not in mapping.values())
        for c in range(6):
            if c not in mapping: mapping[c] = next(remaining)
        fixed = lift_named(labels,data,ports,raw)
        return tuple(mapping[c] for c in fixed), bits, tuple(mapping[c] for c in range(6))
    raise ValueError('target orbit absent from cell')


def validate_colouring(n, edges, ports, boundary, colours):
    need(len(colours) == n and all(type(c) is int and 0 <= c < 6 for c in colours), 'six-colour word')
    need(tuple(colours[v] for v in ports) == boundary, 'exact source')
    need(all(colours[a] != colours[b] for a,b in edges), 'all actual edges')


def consume_receipts(result, source):
    """Reconstruct saved receipts without trusting stored component diagrams.

    The actual inherited graph is authoritative; all parities are recomputed.
    This is a consumer of a declared atlas, not a complete-colouring census.
    """
    n = source['point_count']
    edges = tuple(tuple(e) for e in source['edges']); ports = tuple(source['ports'])
    audit_graph(n, edges, ports)
    need(result['graph_sha256'] == source['graph_sha256'], 'bound graph hash')
    need(result['coordinate_sha256'] == source['coordinate_sha256'], 'bound coordinates')
    patterns = tuple(tuple(row['pattern']) for row in source['rows'])
    need(len(result['receipts']) == len(patterns), 'one receipt per source pattern')
    cell_data = []
    for record in result['selected_cells']:
        raw_labels = record['labels']
        need(type(raw_labels) is str and len(raw_labels) == n and set(raw_labels) <= set('012'),
             'literal complete ternary string')
        labels = tuple(map(int, raw_labels))
        data = fibres(n, edges, labels); need(data is not None, 'consumer odd-cycle guard')
        signature, visible, hidden = boundary_cell(labels, data, ports)
        # JSON-normalization prevents a tuple/list difference from obscuring equality.
        need(json.dumps(signature) == json.dumps(record['signature']), 'recomputed diagram')
        nat(record['hidden_components'])
        need(hidden == record['hidden_components'], 'recomputed hidden multiplicity')
        support = {canonical(w) for w, _ in named_boundary_words(labels, data, ports)}
        need(type(record['covered_patterns']) in (tuple,list), 'support index sequence')
        for i in record['covered_patterns']: nat(i, len(patterns))
        need(len(set(record['covered_patterns'])) == len(record['covered_patterns']), 'distinct support indices')
        need(support == {patterns[i] for i in record['covered_patterns']}, 'recomputed support')
        cell_data.append((labels, data, visible))
    for target, receipt in zip(patterns, result['receipts']):
        need(type(receipt['pattern']) in (tuple,list), 'literal source pattern')
        for c in receipt['pattern']: nat(c, 6)
        need(tuple(receipt['pattern']) == target, 'receipt source order')
        cid = receipt['cell']; nat(cid, len(cell_data))
        labels, data, visible = cell_data[cid]
        bits = tuple(receipt['component_phases'])
        need(len(bits) == len(visible), 'all visible phases')
        for bit in bits: nat(bit, 2)
        palette = tuple(receipt['palette'])
        need(len(palette) == 6 and set(palette) == set(range(6)), 'full palette permutation')
        for c in palette: nat(c, 6)
        comp, parity, _ = data; phases = dict(zip(visible, bits))
        colours = tuple(palette[(3*(parity[v]^phases.get(comp[v],0))+4*labels[v])%6]
                        for v in range(n))
        validate_colouring(n, edges, ports, target, colours)
        need(hashlib.sha256(bytes(colours)).hexdigest() == receipt['colour_sha256'], 'receipt digest')
    return len(result['receipts'])


def equality_criterion(n, edges, s, t):
    """Unpinned graph: every feasible ternary fibre forces an even terminal path."""
    exists = False
    for labels in product(range(3), repeat=n):
        data = fibres(n, edges, labels)
        if data is None: continue
        exists = True
        comp, parity, _ = data
        if labels[s] != labels[t] or comp[s] != comp[t] or parity[s] != parity[t]:
            return False
    return exists


def controls():
    # Independent all-bit enumeration; includes disconnected and empty-boundary cases.
    graphs = 0; label_cases = 0
    for n in range(1,5):
        pairs = tuple(combinations(range(n),2))
        for mask in range(1<<len(pairs)):
            edges = tuple(e for i,e in enumerate(pairs) if mask>>i&1)
            ports = tuple(range(min(n,2)))
            audit_graph(n,edges,ports); graphs += 1
            for labels in product(range(3), repeat=n):
                label_cases += 1; data = fibres(n,edges,labels)
                brute = {}
                for bits in product(range(2), repeat=n):
                    colours = tuple((3*s+4*t)%6 for s,t in zip(bits,labels))
                    if all(colours[a] != colours[b] for a,b in edges):
                        word = tuple(colours[v] for v in ports)
                        brute[word] = brute.get(word,0)+1
                need((data is not None) == bool(brute), 'odd-cycle equivalence')
                if data is None: continue
                signature,visible,hidden = boundary_cell(labels,data,ports)
                generated = dict(named_boundary_words(labels,data,ports))
                need(set(generated) == set(brute), 'full named response')
                need(set(brute.values()) == {2**hidden}, 'hidden-component multiplicity')
                for boundary in brute:
                    validate_colouring(n,edges,ports,boundary,lift_named(labels,data,ports,boundary))
    # Same boundary label is not enough: retain path parity and connectivity.
    for edges,bad in [(((0,1),),(0,0)), (((0,1),(1,2)),(0,3))]:
        n = max(max(e) for e in edges)+1; labels=(0,)*n; ports=(0,n-1)
        try: lift_named(labels,fibres(n,edges,labels),ports,bad)
        except ValueError: pass
        else: raise ValueError('bad pinned path accepted')
    try: lift_named((1,),fibres(1,(),(1,)),(0,),(0,))
    except ValueError: pass
    else: raise ValueError('wrong isolated pin label accepted')
    empty = boundary_cell((0,1),fibres(2,(),(0,1)),())
    need(empty[2] == 2, 'empty boundary keeps hidden component count')
    # Complete-minus-an-edge is a positive logical control, not planar geometry.
    k7e = tuple(e for e in combinations(range(7),2) if e != (0,1))
    need(equality_criterion(7,k7e,0,1), 'K7 minus edge forces equal terminals')
    need(not equality_criterion(7,tuple(combinations(range(7),2)),0,1), 'exclude vacuous K7')
    need(not equality_criterion(3,((0,1),(1,2)),0,2), 'even path alone does not force six-colour equality')
    return {'simple_graphs': graphs, 'ternary_label_cases': label_cases,
            'terminal_controls': ['K7-e', 'K7 vacuity', 'unrestricted even path']}


def mask_indices(mask):
    while mask:
        bit = mask & -mask
        yield bit.bit_length()-1
        mask ^= bit


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--verify-only', action='store_true')
    args = ap.parse_args(); small = controls()
    input_path = ROOT/'05-knowledge/results/coupled_support_branches_20261018c.json'
    raw = input_path.read_bytes().replace(b'\r\n',b'\n')
    source = json.loads(raw)['geometric_candidates']['lifted_candidate']
    if args.verify_only:
        result = json.loads(args.out.read_text(encoding='utf-8'))
        need(result['input_lf_sha256'] == hashlib.sha256(raw).hexdigest(), 'bound input file')
        count = consume_receipts(result, source)
        print(json.dumps({'receipt_count':count, 'checks':CHECKS, 'status':'VERIFIED saved receipts'}))
        return
    n = source['point_count']; edges = tuple(tuple(e) for e in source['edges']); ports=tuple(source['ports'])
    audit_graph(n,edges,ports)
    patterns = tuple(tuple(row['pattern']) for row in source['rows'])
    need(len(patterns)==len(set(patterns))==570, 'declared complete source bank')
    index = {p:i for i,p in enumerate(patterns)}
    # Rebuild source partitions independently of the inherited generator.
    rebuilt = set()
    for c in product(range(6),repeat=7):
        word=(0,)+c
        if canonical(word)==word and all(word[a]!=word[b] for a,b in edges if b<8): rebuilt.add(word)
    need(rebuilt==set(patterns), 'independent source partition census')
    pairings=matchings(tuple(range(6))); need(len(pairings)==15,'all palette matchings')
    by_support={}; signatures=set(); histogram={}; component_hist={}
    for row_id,row in enumerate(source['rows']):
        colours=tuple(map(int,row['colours'])); validate_colouring(n,edges,ports,patterns[row_id],colours)
        for chart_id,matching in enumerate(pairings):
            lookup={c:i for i,pair in enumerate(matching) for c in pair}
            labels=tuple(lookup[c] for c in colours)
            data=fibres(n,edges,labels); need(data is not None, 'colouring supplies a bipartite fibre')
            signature,visible,hidden=boundary_cell(labels,data,ports); signatures.add(signature)
            cover={canonical(word) for word,bits in named_boundary_words(labels,data,ports)}
            need(cover<=set(patterns), 'cell support belongs to actual source')
            mask=sum(1<<index[p] for p in cover)
            need(mask>>row_id&1, 'generator remains represented')
            histogram[len(cover)]=histogram.get(len(cover),0)+1
            component_hist[len(visible)]=component_hist.get(len(visible),0)+1
            if mask not in by_support:
                by_support[mask]={'origin_row':row_id,'matching':matching,'labels':labels,
                                 'signature':signature,'hidden_components':hidden}
    universe=(1<<len(patterns))-1; uncovered=universe; chosen=[]
    while uncovered:
        mask=max(by_support,key=lambda m:((m&uncovered).bit_count(),m.bit_count(),-m))
        need(mask&uncovered != 0,'atlas covers next source')
        chosen.append(mask); uncovered &= ~mask
    # Delete redundant cells; this is inclusion-minimal, never claimed minimum.
    for mask in chosen[:]:
        others=0
        for m in chosen:
            if m!=mask: others |= m
        if others==universe: chosen.remove(mask)
    # A pairwise-incompatible set is a certified lower bound within this bank.
    co=[0]*570
    for mask in by_support:
        for i in mask_indices(mask): co[i] |= mask
    available=universe; packing=[]
    while available:
        i=min(mask_indices(available),key=lambda j:((co[j]&available).bit_count(),j))
        packing.append(i); available &= ~co[i]
    need(all(sum(bool(mask>>i&1) for i in packing)<=1 for mask in by_support),'bank-specific packing bound')
    cells=[]; owner={}
    for cid,mask in enumerate(chosen):
        record=by_support[mask]; labels=record['labels']; data=fibres(n,edges,labels)
        for i in mask_indices(mask): owner.setdefault(i,cid)
        cells.append(dict(record,labels=''.join(map(str,labels)),covered_patterns=list(mask_indices(mask))))
    receipts=[]
    for i,target in enumerate(patterns):
        cid=owner[i]; labels=tuple(map(int,cells[cid]['labels'])); data=fibres(n,edges,labels)
        colours,bits,palette=lift_orbit(labels,data,ports,target)
        validate_colouring(n,edges,ports,target,colours)
        receipts.append({'pattern':target,'cell':cid,'component_phases':bits,'palette':palette,
                         'colour_sha256':hashlib.sha256(bytes(colours)).hexdigest()})
    result={'status':'PROVED finite-interface identities; FINITE-EXACT declared witness atlas',
            'input':input_path.name,'input_lf_sha256':hashlib.sha256(raw).hexdigest(),
            'coordinate_sha256':source['coordinate_sha256'],'graph_sha256':source['graph_sha256'],
            'points':n,'ports':ports,'boundary_patterns':570,'input_words':570,
            'pairing_charts':15,'candidate_cells':8550,'distinct_boundary_diagrams':len(signatures),
            'distinct_supports':len(by_support),'support_size_histogram':histogram,
            'boundary_component_histogram':component_hist,'selected_cells':cells,'receipts':receipts,
            'bank_cover_bounds':{'packing_lower':len(packing),'greedy_upper':len(cells),'packing_patterns':packing},
            'controls':small,
            'limitations':['Not a minimum over all possible ternary labelings.',
                           'The atlas preserves existence, not the count of all six-colourings.',
                           'Palette permutations are retained and may change arithmetic label names.',
                           'A reduced boundary can be glued only along its marked vertices; otherwise enlarge it.',
                           'No new unit-distance edge or equality fragment is supplied.']}
    result['consumed_receipts'] = consume_receipts(result, source)
    result['checks'] = CHECKS
    args.out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({k:result[k] for k in ('candidate_cells','distinct_boundary_diagrams','distinct_supports',
                                         'support_size_histogram','bank_cover_bounds','checks')},indent=2))


if __name__ == '__main__': main()
