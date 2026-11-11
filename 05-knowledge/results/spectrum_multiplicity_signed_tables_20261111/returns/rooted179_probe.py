"""Decide the exact rooted179 target in the inherited complete strong fibre."""
from pathlib import Path
import hashlib
import json
import subprocess

from native_return import converse, ends, strong

HERE = Path(__file__).resolve().parent
REV = '63b86513f8be'
BANK = '04-computation/experiments/h_spectrum_multiplicities_20261021/audit_S50/inventory_375.json'


def need(condition, label):
    if not condition:
        raise RuntimeError(label)


def decode(n, code):
    # Exact convention of the inherited audit_S50/alib.py:decode.
    T = [0] * n
    bit = 0
    for i in range(n):
        for j in range(i + 1, n):
            u, v = (i, j) if code >> bit & 1 else (j, i)
            T[u] |= 1 << v
            bit += 1
    return T


def main():
    raw = subprocess.check_output(['git', 'show', REV + ':' + BANK], cwd=HERE)
    inventory = json.loads(raw)
    records = []
    candidates = []
    for order, rows in inventory['orders'].items():
        n = int(order)
        for code, H, self_converse, aut in rows:
            if H != 179:
                continue
            T = decode(n, int(code))
            e = ends(T)
            h = [sum(row) for row in e]
            h[0] = 1
            full = (1 << n) - 1
            need(h[full] == 179 and strong(T), 'literal strong 179 class')
            record = {'order': n, 'code': code, 'endpoint_counts': e[full]}
            records.append(record)
            for root, E in enumerate(e[full]):
                if E == 46:
                    s = ends(converse(T))
                    U = sum(s[m][root] * h[full ^ m] for m in range(1, full + 1))
                    candidates.append({'order': n, 'code': code, 'root': root,
                                       'H_E_U': [179, E, U], 'rows': T})
    counts = {n: sum(row['order'] == n for row in records) for n in (8, 9, 10)}
    need(counts == {8: 26, 9: 9, 10: 4}, 'global strong fibre universe')
    need([(r['order'], r['code'], r['root'], r['H_E_U'][2]) for r in candidates] == [
        (8, '167805088', 1, 232), (9, '18254201090', 1, 276)], 'exact endpoint46 candidates')
    receipt = {
        'status': 'FINITE-EXACT; global conclusion uses inherited complete strong fibre and prime-prefix reduction in ONE_EXIT.md',
        'target_H_E_U': [179, 46, 278], 'hits': 0,
        'bank_commit': REV, 'bank_path': BANK, 'bank_sha256': hashlib.sha256(raw).hexdigest(),
        'ordinary_class_counts_by_order': counts, 'all_class_endpoint_vectors': records,
        'endpoint46_candidates': candidates,
    }
    (HERE / 'rooted179_receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'status': 'pass', 'strong_classes': len(records),
                      'roots': sum(r['order'] for r in records),
                      'endpoint46_U_values': [r['H_E_U'][2] for r in candidates],
                      'target_hits': 0}))


if __name__ == '__main__':
    main()
