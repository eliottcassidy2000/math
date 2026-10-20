"""Complete labelled three-K search on the bounded probe's residuals only."""
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import json
import probe as P


def path_masks(t, pair_bits):
    return [(word, sum(pair_bits[P.edge(u,v)] for u,v in zip(word,word[1:])))
            for word in P.paths(t) if not P.blue(t,word)]


def recover_path(mask, pairs, n):
    if mask.bit_count() != n-1:
        return None
    neighbors = [[] for _ in range(n)]
    while mask:
        bit = mask & -mask
        mask -= bit
        u,v = pairs[bit.bit_length()-1]
        neighbors[u].append(v)
        neighbors[v].append(u)
        if len(neighbors[u]) > 2 or len(neighbors[v]) > 2:
            return None
    ends = [v for v in range(n) if len(neighbors[v]) == 1]
    if len(ends) != 2 or any(not row for row in neighbors):
        return None
    result = [ends[0]]
    previous = None
    while len(result) < n:
        options = [v for v in neighbors[result[-1]] if v != previous]
        if len(options) != 1:
            return None
        previous = result[-1]
        result.append(options[0])
        if len(set(result)) != len(result):
            return None
    return tuple(result)


def runs(word):
    result = []
    i = 0
    while i < len(word):
        j = i+1
        direction = word[j]-word[i] if j < len(word) else 0
        if abs(direction) == 1:
            while j < len(word) and word[j]-word[j-1] == direction:
                j += 1
        else:
            j = i+1
        result.append(list(word[i:j]))
        i = j
    return result


def search(t):
    n = len(t)
    partner = P.complement(t,tuple(range(n)))
    pairs = list(combinations(range(n),2))
    pair_bits = {pair:1<<i for i,pair in enumerate(pairs)}
    original = sum(pair_bits[(i,i+1)] for i in range(n-1))
    first = path_masks(t,pair_bits)
    last = path_masks(partner,pair_bits)
    successes = 0
    best = []
    for q,qmask in first:
        source_middle = P.complement(t,q)
        for s,smask in last:
            r = recover_path(original ^ qmask ^ smask,pairs,n)
            if r is None:
                continue
            if not P.path(source_middle,r):
                r = r[::-1]
            if not P.path(source_middle,r) or P.blue(source_middle,r):
                continue
            words = (q,r,s)
            evaluation = P.evaluate(t,words)
            P.need(evaluation['success'],'three-K search verification')
            successes += 1
            block_runs = [runs(word) for word in words]
            cost = (sum(len(row) for row in block_runs),
                    sum(item['endpoint_excess'] == 0 for item in evaluation['steps']),
                    max(len(row) for row in block_runs),words)
            cert = {'words':words,'runs':block_runs,'evaluation':evaluation,
                    'compiler':P.compile_ledger(n,words)}
            best.append((cost,cert))
            best.sort(key=lambda item:item[0])
            del best[12:]
    return {'first_black_paths':len(first),'last_black_paths':len(last),
            'pairs_examined':len(first)*len(last),'three_K_certificates':successes,
            'best': [certificate for _,certificate in best]}


def main():
    path = Path(__file__).with_name('probe.json')
    source = json.loads(path.read_text())
    records = []
    for item in source['residuals']:
        result = search(tuple(item['rows']))
        record = {key:item[key] for key in ['n','rows','partner_rows','x','y','branch']}
        record.update(result)
        records.append(record)
        print(json.dumps({key:record[key] for key in ['n','x','y','three_K_certificates']}),flush=True)
        if result['best']:
            print('best',result['best'][0]['words'],[s['endpoint_excess'] for s in result['best'][0]['evaluation']['steps']],flush=True)
    receipt = {'status':'FINITE-EXACT complete labelled three-K search for the nine specified residual presentations',
               'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
               'probe_source_sha256':sha256(Path(P.__file__).read_bytes()).hexdigest(),
               'input_receipt_sha256':sha256(path.read_bytes()).hexdigest(),
               'scope':'Every black first and last word on each exact labelled endpoint; no isomorphism or converse endpoint allowance',
               'records':records}
    Path(__file__).with_suffix('.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    return receipt


if __name__ == '__main__':
    main()
