"""Independent verifier: no imports from check.py or inherited research code.

Verifies supplied XOR proof vectors directly, without Gaussian elimination.
Polynomial bit 12 is the constant term; lower bits are the 12 free bits.
"""
from itertools import combinations, product
from pathlib import Path
import argparse
import hashlib
import json

HERE = Path(__file__).resolve().parent
C = 1 << 12
P = tuple(range(10))
PAIRS = list(combinations(P, 2))
CUTS = [(1, 8), (2, 8), (2, 7)]


def require(test, message):
    if not test:
        raise RuntimeError(message)


def pvalue(polynomial, value):
    return ((polynomial & value).bit_count() ^ (polynomial >> 12)) & 1


def source(incoming):
    orbits = sorted({min((i,j),(9-j,9-i)) for i,j in combinations(range(1,9),2) if j>i+1})
    require(len(orbits) == 12, 'dimension')
    position = {p:i for i,p in enumerate(orbits)}
    matrix = [[None]*10 for _ in P]
    for i,j in PAIRS:
        if i == 0:
            f = C if j != 9 and j not in incoming else 0
        elif j == 9:
            f = 0 if 9-i in incoming else C
        elif j == i+1:
            f = C
        else:
            f = 1 << position[min((i,j),(9-j,9-i))]
        matrix[i][j], matrix[j][i] = f, f ^ C
    require(all(matrix[i][i+1] == C for i in range(9)), 'standard path')
    require(all(matrix[i][j] == matrix[9-j][9-i] for i,j in PAIRS), 'blue source')
    return matrix, orbits


def k_move(matrix, word):
    require(tuple(sorted(word)) == P, 'word is permutation')
    retained = {frozenset(p) for p in zip(word, word[1:])}
    return [[None if i==j else matrix[i][j] ^ (0 if frozenset((i,j)) in retained else C)
             for j in P] for i in P]


def xor_proof(equations, indices):
    require(len(indices) == len(set(indices)), 'duplicate proof indices')
    require(all(isinstance(i,int) and 0 <= i < len(equations) for i in indices), 'proof range')
    value = 0
    for index in indices:
        value ^= equations[index]
    return value


def rows_of(matrix, value):
    return [sum(1<<j for j in P if i != j and pvalue(matrix[i][j],value)) for i in P]


def direct_replay(rows, words):
    now = list(rows)
    witnesses = []
    for q in words:
        require(tuple(sorted(q)) == P, 'numeric permutation')
        require(all(now[u] >> v & 1 for u,v in zip(q,q[1:])), 'numeric directed')
        rev = {q[i]:q[9-i] for i in P}
        failures = [(u,v,rev[v],rev[u]) for u in P for v in P if u != v
                    and bool(now[u]>>v&1) != bool(now[rev[v]]>>rev[u]&1)]
        require(failures, 'numeric actual anti failure')
        witnesses.append(failures[0])
        edges = {frozenset(e) for e in zip(q,q[1:])}
        nxt = [0]*10
        for u,v in PAIRS:
            forward = bool(now[u] >> v & 1) ^ (frozenset((u,v)) not in edges)
            a,b = (u,v) if forward else (v,u)
            nxt[a] |= 1<<b
        now = nxt
    target = list(rows)
    for u,v in PAIRS:
        if v != u+1:
            target[u] ^= 1<<v
            target[v] ^= 1<<u
    require(now == target, 'numeric exact labelled target')
    return witnesses


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--certificates', type=Path, default=HERE/'certificates.json')
    parser.add_argument('--output', type=Path, default=HERE/'independent_receipt.json')
    args = parser.parse_args()
    data = json.loads(args.certificates.read_text())
    require(data['order'] == 10 and data['dimension'] == 12 and data['cut_pairs'] == [list(p) for p in CUTS], 'header')
    expected = {(incoming,bits) for incoming in combinations(range(2,9),2) for bits in product(range(2),repeat=3)}
    observed, singleton, contradictions, direction_checks = set(), set(), 0, 0
    for cell in data['cells']:
        incoming,bits = tuple(cell['incoming']),tuple(cell['cell'])
        key = incoming,bits
        require(key in expected and key not in observed, 'unique valid cell')
        observed.add(key)
        original, orbits = source(incoming)
        require(data['variables'] == [sorted([list(p) for p in {(i,j),(9-j,9-i)}]) for i,j in orbits], 'named variable ordering')
        guards = [original[i][j] ^ (bit*C) for (i,j),bit in zip(CUTS,bits)]
        guard_values = {}
        for equation in guards:
            coefficient = equation & (C-1)
            require(coefficient and coefficient & (coefficient-1) == 0 and coefficient not in guard_values, 'independent cuts')
            guard_values[coefficient] = equation & C
        require(len(cell['words']) == len(cell['proofs']) == 3, 'three steps')
        now = original
        for word, proof in zip(cell['words'],cell['proofs']):
            require(tuple(sorted(word)) == P, 'symbolic permutation')
            for u,v in zip(word,word[1:]):
                equation = now[u][v] ^ C
                coefficient = equation & (C-1)
                if coefficient:
                    require(coefficient in guard_values, 'direction depends on an unfixed bit')
                    equation ^= coefficient ^ guard_values[coefficient]
                require(equation == 0, 'symbolic direction under guards')
                direction_checks += 1
            rev = {word[i]:word[9-i] for i in P}
            equations = guards + [now[u][v] ^ now[rev[v]][rev[u]] for u,v in PAIRS]
            if 'contradiction' in proof:
                require(set(proof) == {'contradiction'}, 'proof shape')
                require(xor_proof(equations,proof['contradiction']) == C, 'XOR is not 0=1')
                contradictions += 1
            else:
                require(set(proof) == {'assignment','units'}, 'solution shape')
                value = proof['assignment']
                require(0 <= value < C and len(proof['units']) == 12, 'solution range')
                require(all(pvalue(e,value) == 0 for e in equations), 'solution violates equation')
                for k,indices in enumerate(proof['units']):
                    require(xor_proof(equations,indices) == (1<<k)^(((value>>k)&1)*C), 'rank proof')
                singleton.add((incoming,value))
            now = k_move(now,word)
        require(now == k_move(original,P), 'symbolic target')
    require(observed == expected, 'complete source partition')
    fixes = json.loads((HERE/'exceptions.json').read_text())
    require({(tuple(r['incoming']),r['mask']) for r in fixes} == singleton, 'all and only residuals repaired')
    witnesses = []
    for fix in fixes:
        incoming = tuple(fix['incoming']); matrix,_ = source(incoming)
        require(rows_of(matrix,fix['mask']) == fix['rows'], 'literal source')
        require(rows_of(k_move(matrix,P),fix['mask']) == fix['partner_rows'], 'literal target')
        require(fix['outgoing'] == sorted(9-i for i in incoming), 'outgoing ports')
        witnesses.append({'incoming':incoming, 'anti_failures':direct_replay(fix['rows'],fix['words'])})
    require(contradictions == 501 and len(singleton) == 3 and direction_checks == 4536, 'summary counts')
    receipt = {'status':'PASS: independent direct-XOR certificate verification', 'cells':len(observed),
        'inconsistent_blue_systems':contradictions,'singleton_residuals':len(singleton),
        'symbolic_direction_checks':direction_checks,'numeric_replays':witnesses,
        'source_hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),args.certificates,HERE/'exceptions.json']}}
    args.output.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))


if __name__ == '__main__':
    main()
