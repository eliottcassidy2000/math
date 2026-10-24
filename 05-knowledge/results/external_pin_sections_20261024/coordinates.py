"""Authenticate the recovered code-basis adapter for incoming graded sections."""
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent
COUNT = 0


def need(ok, message):
    global COUNT
    COUNT += 1
    if not ok:
        raise RuntimeError(message)


def add(a, b):
    return ((((a >> 4) + (b >> 4)) % 128) << 4) | ((a ^ b) & 15)


def forward(q):
    x, y = q >> 4, q & 15
    return (x << 4) | (y >> 1) | (((x & 1) ^ (y & 1)) << 3)


def inverse(q):
    x, y = q >> 4, q & 15
    return (x << 4) | ((y & 7) << 1) | ((x & 1) ^ (y >> 3))


def main():
    electric = json.loads((RESULTS/'metagraph_order8_electrical_20261023/receipt.json').read_text())['order8']
    x, y = electric['higher_layers'][-1]['lifts'][0], electric['radical_lifts_mod4']
    for q in range(2048):
        need(inverse(forward(q)) == q and forward(inverse(q)) == q, 'two-sided full-group bijection')
        for g in (16,1,2,4,8):
            need(forward(add(q,g)) == add(forward(q),forward(g)), 'additivity at all group generators')
    for s in range(len(x)):
        incoming = ((x[s]%128)<<4) | sum((y[j][s]%2)<<j for j in range(4))
        native = ((x[s]%128)<<4) | sum((y[j+1][s]%2)<<j for j in range(4))
        need(forward(incoming) == native, 'actual owner-code adapter from canonical electrical relation')
    bank = json.loads((RESULTS/'three_target570_20261024/bank.json').read_text())
    ids = bank['reservoir_sources']
    labels = [((x[s]%128)<<4) | sum((y[j+1][s]%2)<<j for j in range(4)) for s in ids]
    witnesses = 0
    for k, section in enumerate(bank['sections']):
        for q, mask in section.items():
            value = 0
            for i, label in enumerate(labels):
                if mask >> i & 1:
                    value = add(value,label)
            need(mask.bit_count() == k and value == forward(int(q)),
                 'every incoming full subset witness retains grade and physical code after adaptation')
            witnesses += 1
    print(f'PASS {COUNT}; all2048 code values, {len(x)} actual owners, {witnesses} incoming grade witnesses')


if __name__ == '__main__':
    main()
