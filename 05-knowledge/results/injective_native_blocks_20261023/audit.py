"""Independent final audit: literal ports and a rowwise image certificate.

No sibling checker imports and no discovery-cache inputs. The inherited owner
map supplies the finite isomorphism atlas; no new isomorphism algorithm is claimed.
"""
from array import array
from collections import Counter
import csv
import gzip
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent
NATIVE = RESULTS/'metagraph_order8_20261023'
CHECKS = 0


def need(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(why)


def main():
    doc = json.loads((HERE/'search/selected_blocks.json').read_text())
    images = json.loads((HERE/'electrical/image_witnesses.json').read_text())
    native = json.loads((RESULTS/'metagraph_order8_electrical_20261023/receipt.json').read_text())['order8']
    nodes = list(csv.DictReader((NATIVE/'merged_nodes.tsv').open(), delimiter='\t'))
    merged = {int(row[k]): int(row['merged_id']) for row in nodes for k in ('class_a', 'class_b')}
    owner = array('H')
    owner.frombytes(gzip.decompress((NATIVE/'owner_map.u16.gz').read_bytes()))
    if sys.byteorder != 'little':
        owner.byteswap()
    edges = [tuple(int(row[k]) for k in ('source', 'target', 'blue', 'black'))
             for row in csv.DictReader(gzip.decompress((NATIVE/'merged_edges.tsv.gz').read_bytes()).decode().splitlines(), delimiter='\t')]
    adj = [[] for _ in nodes]
    for i, j, blue, black in edges:
        adj[i].append((j, blue+black))
        adj[j].append((i, blue+black))
    def lap(x):
        # Independent row-sum implementation; the producer uses edge accumulation.
        return [sum(w*(x[i]-x[j]) for j, w in row) for i, row in enumerate(adj)]
    certs = []
    for record in images['solutions']:
        d = [0]*len(nodes)
        for i, c in record['charge']:
            d[i] = 2*c
        z = record['potential_mod256']
        correction = dict(record['integer_correction_div256_nonzero'])
        need(len(z) == len(nodes) and sum(d) == 0, 'full grounded charge carrier')
        lhs = lap(z)
        need(all(lhs[i] == d[i]+256*correction.get(i, 0) for i in range(len(nodes))), 'exact integer image with retained correction')
        certs.append((dict(record['charge']), z))
    u, v = doc['targets']
    def get_certificate(charge):
        return next(z for d, z in certs if d == charge)
    z1 = get_certificate({1425:-1,1477:1,1921:-1,2837:1})
    zv = get_certificate({u:1,v:-1})
    z2 = get_certificate({1394:-1,1477:1,1841:1,2185:-1})
    x = native['higher_layers'][-1]['lifts'][0]
    ys = native['radical_lifts_mod4'][1:]
    need(all(a % 256 == 0 for a in lap(x)), 'original integer lift is harmonic mod256')
    chords = list((i,j) for i in range(8) for j in range(i+2,8))
    pairs = list(combinations(range(8),2))
    full = (1 << len(chords))-1
    def word(a, p):
        return sum(a[p[i]][p[j]] << k for k,(i,j) in enumerate(chords))
    state_words = [[], []]
    for slot in doc['slots']:
        canon = slot['canon']
        a = [[0]*8 for _ in range(8)]
        for k,(i,j) in enumerate(pairs):
            a[i][j] = (canon >> (27-k)) & 1
            a[j][i] = 1-a[i][j]
        for bit, key in enumerate(('before','after')):
            p = slot[key]['path']
            need(sorted(p) == list(range(8)) and all(a[i][j] for i,j in zip(p,p[1:])), 'literal path legal')
            m = word(a,p)
            need(m == slot[key]['mask'], 'actual indexed forward word')
            need(merged[owner[m]] == slot['source'] and merged[owner[full^m]] == slot[key]['target'], 'actual source and complement target')
            changed = [[1-a[i][j] if i != j else 0 for j in range(8)] for i in range(8)]
            for i,j in zip(p,p[1:]):
                changed[i][j], changed[j][i] = a[i][j], a[j][i]
            need(word(changed,p) == (full^m), 'literal K_P leaves exactly path arcs')
            cm = word(list(map(list, zip(*a))), list(reversed(p)))
            need(owner[m] != owner[cm] and merged[owner[cm]] == slot['source'], 'actual distinct converse source owners')
            need(merged[owner[full^cm]] == slot[key]['target'], 'actual converse target owner')
            state_words[bit].extend((m,cm))
    for words in state_words:
        opposite = [full^m for m in words]
        need(len(set(words)) == len(set(opposite)) == 16, 'injective source AND target normalized ports')
        need(not set(words).intersection(opposite), 'source and target endpoints disjoint')
    # Derive actual selected conductance changes directly, without a rectangle formula.
    changes = Counter()
    for slot in doc['slots']:
        changes[tuple(sorted((slot['source'],slot['before']['target'])))] -= 2
        changes[tuple(sorted((slot['source'],slot['after']['target'])))] += 2
    delta = {}
    for (i,j), coefficient in changes.items():
        if not coefficient:
            continue
        delta.setdefault(i,Counter())[i] += coefficient
        delta.setdefault(i,Counter())[j] -= coefficient
        delta.setdefault(j,Counter())[j] += coefficient
        delta.setdefault(j,Counter())[i] -= coefficient
    charges = Counter()
    for slot in doc['slots']:
        charges[slot['source']] += slot['sign']
    for node, row in delta.items():
        if node == u:
            z = [a-b for a,b in zip(z1,z2)]
        elif node == v:
            z = [b-a for a,b in zip(z1,z2)]
        else:
            z = [charges[node]*a for a in zv]
        lhs = lap(z)
        need(all((lhs[j]-row.get(j,0)) % 256 == 0 for j in range(len(nodes))), 'physical DeltaL row lies in imL modulo256')
        need(sum(row.values()) == 0 and sum(c*x[j] for j,c in row.items()) == 0, 'zero row sum and exact original integer current')
    # Repeat the complete restricted binary-state count by building actual edge changes.
    accepted = []
    for bits in product(range(2), repeat=8):
        records = [slot['after' if bit else 'before'] for slot,bit in zip(doc['slots'],bits)]
        if len({(slot['source'],tuple(row['path'])) for slot,row in zip(doc['slots'],records)}) != 8:
            continue
        if any(records[a]['path'][i] != records[b]['path'][j] for (a,i),(b,j) in doc['shared_owner_equalities']):
            continue
        if Counter(row['target'] for row in records) != Counter(slot['before']['target'] for slot in doc['slots']):
            continue
        d = Counter()
        for slot,bit in zip(doc['slots'],bits):
            if bit:
                d[slot['source']] += slot['sign']
        if sum(c*x[s] for s,c in d.items()) % 128:
            continue
        if any(sum(c*y[s] for s,c in d.items()) % 2 for y in ys):
            continue
        accepted.append(''.join(map(str,bits)))
    need(accepted == ['00000000','01110111','10001000','11111111'], 'exact complete restricted state fibre')
    blocks = ((0,4),(1,2,3,5,6,7))
    factored = [''.join(map(str,bits)) for bits in product(range(2),repeat=8)
                if all(len({bits[i] for i in block}) == 1 for block in blocks)]
    need(factored == accepted, 'authenticated native fibre has exactly two implication blocks')
    output = {'status':'FINITE-EXACT independent audit', 'checks':CHECKS,
              'literal_source_ports_per_state':16, 'literal_target_ports_per_state':16,
              'image_certificates':len(certs), 'nonzero_Laplacian_rows':len(delta),
              'accepted_eight_slot_states':accepted,
              'exact_implication_blocks':[list(b) for b in blocks],
              'all_dyadic_extension_dependency':'original certified 2-primary exponent256',
              'claim':'Actual selected port block; universal original dyadic aggregate guards; not all-modulus or full native rewiring.'}
    path = HERE/'audit_receipt.json'
    if '--write' in sys.argv:
        path.write_text(json.dumps(output,indent=2)+'\n',encoding='utf-8')
    else:
        need(json.loads(path.read_text()) == output, 'frozen independent audit receipt')
    print(json.dumps(output,sort_keys=True))


if __name__ == '__main__':
    main()
