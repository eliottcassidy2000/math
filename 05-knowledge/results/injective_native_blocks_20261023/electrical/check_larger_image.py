"""Independent literal image and full subset checks for the 12/16-source circuits."""
from collections import Counter
import csv
import gzip
from hashlib import sha256
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
REPO = next(p for p in HERE.parents if (p/'01-canon').is_dir())
CHECKS = 0


def need(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok: raise RuntimeError(why)


def main():
    pins = json.loads((HERE/'SOURCE_PINS.json').read_text())
    blobs = {}
    for key in ('edges','electrical_receipt'):
        pin = pins[key];raw = (REPO/pin['path']).read_bytes()
        normalized = raw if pin['hash_mode'] == 'raw' else raw.replace(b'\r\n',b'\n')
        need(sha256(normalized).hexdigest() == pin['sha256'], 'Original input hash '+key)
        blobs[key] = raw
    edges = [tuple(int(row[k]) for k in ('source','target','blue','black'))
             for row in csv.DictReader(gzip.decompress(blobs['edges']).decode().splitlines(),delimiter='\t')]
    original = json.loads(blobs['electrical_receipt'])['order8']
    xs = original['higher_layers'][-1]['lifts'][0]
    ys = original['radical_lifts_mod4'][1:]
    data = json.loads((HERE/'larger_image.json').read_text())
    need(data['modulus'] == 256 and data['input_edges_sha256'] == pins['edges']['sha256'],
         'Larger certificate original edge hash and modulus')
    need(len(data['solutions']) == 2, 'Both larger source circuits present')
    expected = [dict([(i,1) for i in (1480,1691,2071,2185,2581,3000)] +
                     [(i,-1) for i in (1638,1923,2176,2553,2835,2945)]),
                dict([(i,1) for i in (176,1477,1691,1974,2071,2087,2185,2905)] +
                     [(i,-1) for i in (1394,1425,1457,1638,1841,1921,2553,2581)])]
    reports = []
    for record,charge in zip(data['solutions'],expected):
        need(dict(record['charge']) == charge and record['coefficient'] == 2, 'Exact signed owner charge')
        potential = record['potential_mod256']
        need(len(potential) == 3528 and potential[3527] == 0, 'Full grounded potential')
        current = [0]*3528
        for i,j,b,k in edges:
            value = (b+k)*(potential[i]-potential[j]);current[i] += value;current[j] -= value
        rhs = [2*charge.get(i,0) for i in range(3528)]
        for i,(a,b) in enumerate(zip(current,rhs)):
            need((a-b) % 256 == 0, 'Literal source image at node '+str(i))
        correction = [[i,(a-b)//256] for i,(a,b) in enumerate(zip(current,rhs)) if a != b]
        need(correction == record['integer_correction_div256_nonzero'], 'Full integer correction ledger')
        positive = sum(xs[i] for i,c in charge.items() if c == 1)
        negative = sum(xs[i] for i,c in charge.items() if c == -1)
        need(positive == negative, 'Exact primitive-potential cancellation')
        terms = sorted(charge.items());size = len(terms)
        signs = [c for i,c in terms]
        xvalues = [c*xs[i] for i,c in terms]
        yvalues = [sum((y[i] % 2) << j for j,y in enumerate(ys)) for i,c in terms]
        states = [(0,0,0)]*(1 << size)
        counts = Counter()
        accepted_full = []
        for mask in range(1 << size):
            if mask:
                bit = mask & -mask;j = bit.bit_length()-1
                balance,cyclic,binary = states[mask ^ bit]
                states[mask] = (balance+signs[j],cyclic+xvalues[j],binary ^ yvalues[j])
            balance,cyclic,binary = states[mask]
            if balance: continue
            counts[2] += 1  # Every paired charge is zero modulo two.
            for k in range(2,9):
                if binary == 0 and cyclic % (1 << (k-1)) == 0:
                    counts[1 << k] += 1
                    if k == 8: accepted_full.append(mask)
        need(accepted_full == [0,(1 << size)-1], 'No proper balanced full-dyadic source subcircuit')
        need(counts[256] == 2, 'Exactly two full-dyadic accepted endpoint subsets')
        if size == 16:
            need([counts[1 << k] for k in range(2,9)] == [378,194,86,38,18,4,2],
                 'Exact sixteen-slot dyadic pruning profile')
        reports.append({'support':size,'charge':record['charge'],'primitive_positive_sum':positive,
                        'primitive_negative_sum':negative,'all_binary_subsets_tested':1 << size,
                        'degree_and_dyadic_accepted_counts':dict(sorted(counts.items())),
                        'accepted_full_dyadic_subsets':accepted_full,
                        'scope':'Selected signed source circuit; native path/capacity authentication is in search/'})
    result = {'status':'FINITE-EXACT direct original image and source-circuit controls',
              'larger_image_sha256_lf':sha256((HERE/'larger_image.json').read_bytes().replace(b'\r\n',b'\n')).hexdigest(),
              'original_source_pins':{key:pins[key] for key in ('edges','electrical_receipt')},
              'circuits':reports,'checks_before_replay':CHECKS}
    result = json.loads(json.dumps(result))
    target = HERE/'larger_image_receipt.json'
    if '--write' in sys.argv: target.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    else: need(result == json.loads(target.read_text()), 'Frozen larger image receipt')
    print(json.dumps({'status':'PASS','checks':CHECKS,'circuits':reports}))


if __name__ == '__main__': main()
