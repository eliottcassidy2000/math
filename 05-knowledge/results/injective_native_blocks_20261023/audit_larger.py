"""Independent literal-word and subset-sum audit of the sixteen-source circuit."""
from array import array
import csv
import gzip
from itertools import combinations
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent
CHECKS = 0


def need(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(why)


def main():
    doc = json.loads((HERE/'search/larger16_selected.json').read_text())
    block = doc['circuits'][0]
    electrical = json.loads((RESULTS/'metagraph_order8_electrical_20261023/receipt.json').read_text())['order8']
    x = electrical['higher_layers'][-1]['lifts'][0]
    ys = electrical['radical_lifts_mod4'][1:]
    native = RESULTS/'metagraph_order8_20261023'
    nodes = list(csv.DictReader((native/'merged_nodes.tsv').open(),delimiter='\t'))
    merged = {int(r[k]):int(r['merged_id']) for r in nodes for k in ('class_a','class_b')}
    owner = array('H'); owner.frombytes(gzip.decompress((native/'owner_map.u16.gz').read_bytes()))
    if sys.byteorder != 'little':
        owner.byteswap()
    pairs = list(combinations(range(8),2))
    chords = [(i,j) for i in range(8) for j in range(i+2,8)]
    full = (1 << len(chords))-1
    def word(a,p):
        return sum(a[p[i]][p[j]] << k for k,(i,j) in enumerate(chords))
    inventories = [set(),set()]
    sources = [slot['source'] for slot in block['slots']]
    need(len(sources) == len(set(sources)) == 16, 'distinct actual sources')
    for slot in block['slots']:
        a = [[0]*8 for _ in range(8)]
        canon = int(nodes[slot['source']]['canon_a'])
        for bit,(i,j) in enumerate(pairs):
            a[i][j] = (canon >> (27-bit)) & 1
            a[j][i] = 1-a[i][j]
        for state,key in enumerate(('before','after')):
            p = slot[key]['path']
            need(sorted(p) == list(range(8)) and all(a[i][j] for i,j in zip(p,p[1:])), 'actual path')
            m = word(a,p)
            cm = word(list(map(list,zip(*a))), list(reversed(p)))
            need(m == slot[key]['mask'], 'exact indexed word')
            for q in (m,cm):
                need(merged[owner[q]] == slot['source'], 'actual merged source')
                need(merged[owner[full^q]] == slot[key]['target'], 'actual inverse target')
                need(q not in inventories[state] and (full^q) not in inventories[state], 'no reused endpoint port')
                inventories[state].update((q,full^q))
    need([len(s) for s in inventories] == [64,64], 'all endpoint capacities retained')
    need(len(inventories[0] | inventories[1]) == 128, 'entire fixed-background reservoir')
    u,v = doc['targets']
    need(x[u] == x[v] and all((y[u]-y[v]) % 2 == 0 for y in ys), 'all target-row debts vanish')
    # Doubling a list of partial sums is independent of the producer's full bit-word loop.
    partial = [(0,0,0,0)]  # signed margin, exact X pairing, four binary sums, subset mask
    for i,slot in enumerate(block['slots']):
        sign,source = slot['sign'],slot['source']
        binary = sum((y[source] % 2) << j for j,y in enumerate(ys))
        partial += [(a+sign,b+sign*x[source],c^binary,mask|(1 << i)) for a,b,c,mask in partial]
    counts = {}
    for k in range(2,9):
        accepted = [mask for a,b,c,mask in partial if a == 0 and c == 0 and b % (1 << (k-1)) == 0]
        counts[str(1 << k)] = len(accepted)
    need(list(counts.values()) == [378,194,86,38,18,4,2], 'complete dyadic acceptance ladder')
    need(accepted == [0,65535], 'no nontrivial proper full-bank submove')
    need(partial[-1][0:3] == (0,0,0), 'full circuit pays exact integer X and all binary debts')
    need(sum(x[s] for s in block['positive']) == sum(x[s] for s in block['negative']) == 1016, 'exact positive and negative source sums')
    out = {'status':'FINITE-EXACT independent sixteen-source audit', 'checks':CHECKS,
           'subset_states':len(partial),'accepted_by_original_dyadic_modulus':counts,
           'full_bank_accepted_masks':accepted,'per_state_endpoint_ports':64,
           'fixed_background_reservoir':128,'exact_X_sum_each_side':1016,
           'scope':'Sharp sixteen-slot budget in supplied binary menu; aggregate guards; no global wiring claim.'}
    path = HERE/'audit_larger_receipt.json'
    if '--write' in sys.argv:
        path.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    else:
        need(json.loads(path.read_text()) == out, 'saved independent larger receipt')
    print(json.dumps(out,sort_keys=True))


if __name__ == '__main__':
    main()
