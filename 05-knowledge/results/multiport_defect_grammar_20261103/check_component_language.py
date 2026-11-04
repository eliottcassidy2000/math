"""Exact singleton/C3 component-language inverse and independent finite audit."""
from functools import lru_cache
from itertools import combinations, product
from math import comb
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def atom_signature(atom):
    radix, k = atom
    return ((1, 2, k, k, 0, k, k) if radix == 2 else
            (3, 12, 3*k, 3*k, k*(k-1), k, k))


def ordered(left, right):
    h,p,f,g,j,e,s = left
    H,P,F,G,J,E,S = right
    return (h*H, p*P, f*H+p*F, g*P+h*G,
            j*H+h*J+e*S+g*F, h*E, s*H)


def signature(word):
    current = atom_signature(word[0])
    for atom in word[1:]:
        current = ordered(current, atom_signature(atom))
    return current


@lru_cache(None)
def inverse_radices(q, triangles):
    """Low component first; no artificial empty or high-zero terminal atom."""
    if q <= 0 or triangles < 0:
        return ()
    output = []
    for radix in (2, 4):
        remaining = triangles-(radix == 4)
        if remaining < 0:
            continue
        digit, quotient = q % radix, q // radix
        if quotient == 0:
            # A leading C3 digit1 is legal even though q has only one bit.
            if digit > 0 and remaining == 0:
                output.append(((radix, digit),))
        else:
            output.extend(((radix, digit),)+word
                          for word in inverse_radices(quotient, remaining))
    return tuple(output)


def domino_words(q, triangles):
    """Independent binary positions: dominoes may cover the final zero pad."""
    length = q.bit_length()
    output = set()
    for starts in combinations(range(length), triangles):
        if any(y == x+1 for x,y in zip(starts, starts[1:])):
            continue
        starts = set(starts)
        word = []
        i = 0
        while i < length:
            if i in starts:
                word.append((4, (q >> i) & 3))
                i += 2
            else:
                word.append((2, (q >> i) & 1))
                i += 1
        output.add(tuple(word))
    return output


def inverse(f):
    need(f > 0, 'positive target marginal required')
    output = []
    triangles = 0
    q = f
    while True:
        words = set(inverse_radices(q, triangles))
        need(words == domino_words(q, triangles), ('independent inverse', f, triangles))
        top = q.bit_length()-triangles+1
        expected = comb(top, triangles) if top >= triangles else 0
        need(len(words) == expected, ('word count', f, triangles, len(words), expected))
        for word in sorted(words):
            need(word[-1][1] > 0, ('terminal port', word))
            need(sum(r == 4 for r,k in word) == triangles, ('triangle count', word))
            need(signature(word)[2] == f, ('target marginal', word))
            output.append(word)
        if q % 3:
            break
        q //= 3
        triangles += 1
    need(len(output) == len(set(output)), 'inverse duplicate')
    return output


def literal_word(word):
    rows = []
    ports = 0
    for radix, k in word:
        first = len(rows)
        size = 1 if radix == 2 else 3
        fresh = ((1 << size)-1) << first
        rows = [r | fresh for r in rows] + [0]*size
        if size == 3:
            for u,v in ((0,1),(1,2),(2,0)):
                rows[first+u] |= 1 << (first+v)
        for offset in range(k):
            ports |= 1 << (first+offset)
    return rows, ports


def literal_signature(rows, ports):
    """Subset path DP, independent of the ordered-component transfer."""
    n = len(rows)
    full = (1 << n)-1
    dp = [[0]*n for _ in range(full+1)]
    starts = [[0]*n for _ in range(full+1)]
    for v in range(n):
        dp[1 << v][v] = 1
        starts[1 << v][v] = (ports >> v) & 1
    for mask in range(1, full+1):
        for end in range(n):
            total, marked = dp[mask][end], starts[mask][end]
            if not total:
                continue
            possible = rows[end] & ~mask
            while possible:
                bit = possible & -possible
                possible ^= bit
                vertex = bit.bit_length()-1
                dp[mask | bit][vertex] += total
                starts[mask | bit][vertex] += marked
    paths = [sum(row) for row in dp]
    paths[0] = 1
    endings = [sum(count for v,count in enumerate(row) if ports >> v & 1) for row in dp]
    beginnings = [sum(row) for row in starts]
    return (paths[full],
            sum(paths[m]*paths[full ^ m] for m in range(full+1)),
            sum(endings[m]*paths[full ^ m] for m in range(full+1)),
            sum(paths[m]*beginnings[full ^ m] for m in range(full+1)),
            sum(endings[m]*beginnings[full ^ m] for m in range(full+1)),
            endings[full], beginnings[full])


def full_attachment(rows, ports, clock):
    n = len(rows)
    root = clock+n
    result = [0]*(root+1)
    for i in range(clock):
        result[i] = ((1 << (root-i-1))-1) << (i+1)
        result[root] |= 1 << i
    for i,row in enumerate(rows):
        result[clock+i] = row << clock
        if ports >> i & 1:
            result[clock+i] |= 1 << root
        else:
            result[root] |= 1 << (clock+i)
    return result


def strong(rows):
    full = (1 << len(rows))-1
    for source in range(len(rows)):
        reached = todo = 1 << source
        while todo:
            bit = todo & -todo
            todo ^= bit
            new = rows[bit.bit_length()-1] & ~reached
            todo |= new
            reached |= new
        if reached != full:
            return False
    return True


def main():
    need(inverse_radices(1, 1) == (((4, 1),),), 'top C3 digit1 was lost')
    need(((4, 2),) in inverse_radices(2, 1), 'single C3 digit2 was lost')
    need(((4, 3),) in inverse_radices(3, 1), 'single C3 digit3 was lost')
    need(inverse_radices(1, 2) == (), 'impossible extra triangle accepted')
    c4_rows = [8,1,3,6]
    c4_signatures = [literal_signature(c4_rows, 1 << v) for v in range(4)]
    need(strong(c4_rows), 'richer atom strongness')
    need(all(sig[:2] == (5,32) for sig in c4_signatures), 'richer atom h,p')
    need([sig[2] for sig in c4_signatures] == [7,9,6,5], 'richer atom one-port marginals')
    expected_hits = {1:['2','4'], 2:['22','42'], 3:['422']}
    expected_counts = {0:1, 1:4, 2:7, 3:13}
    levels = []
    literal_checks = attachment_checks = 0
    for level in range(4):
        f = 7*3**level
        bank = []
        for word in inverse(f):
            sig = signature(word)
            rows, ports = literal_word(word)
            need(literal_signature(rows, ports) == sig, ('literal signature', word))
            literal_checks += 1
            h,p,ff,g,j,e,s = sig
            for clock in (1,2):
                full = full_attachment(rows, ports, clock)
                need(strong(full), ('actual strongness', word, clock))
                actual_h = literal_signature(full, 0)[0]
                need(actual_h == h-j+(1 << clock)*ff, ('actual attachment count', word, clock))
                attachment_checks += 1
            bank.append(dict(word=[('T' if r == 2 else 'C', k) for r,k in word],
                             signature=sig, intercept=h-j, rows=rows,
                             ports=[i for i in range(len(rows)) if ports >> i & 1]))
        need(len(bank) == expected_counts[level], ('candidate count', level))
        targets = []
        for letters in product((2,4), repeat=level):
            intercept = -5
            for digit in letters:
                intercept = 3*intercept+digit
            accepted = [i for i,record in enumerate(bank) if record['intercept'] == intercept]
            targets.append(dict(letters=''.join(map(str, letters)),
                                intercept=intercept, candidate_indices=accepted))
        if level:
            hits = [t['letters'] for t in targets if t['candidate_indices']]
            need(hits == expected_hits[level], ('exact target classification', level, hits))
            print('LEVEL %s: f=%s; %s component words; accepted=%s' % (level, f, len(bank), ','.join(hits)))
        levels.append(dict(level=level, f=f, candidates=bank, targets=targets,
                           distinct_intercepts=sorted(set(b['intercept'] for b in bank))))
    need((literal_checks, attachment_checks) == (25,50), 'audit universe changed')
    report = dict(
        status='FINITE_EXACT_COMPLETE_FIXED_MARGINAL_ATOM_LANGUAGE_ACCEPTANCE',
        alphabet=dict(T='singleton with port bit0 or1', C='directed triangle with0..3 ports'),
        scope='Ordered sums, nonempty terminal port set, common transitive prefix clock; no native or old-name transformation claim.',
        source=dict(f=7, intercept=-5, word=[['T',1],['T',1],['T',1]]),
        inverse='For every c with3^c dividing f, decode f/3^c in radices2/4 with exactlyc radix4 positions and nonzero top digit.',
        independent_controls=dict(binary_domino_inverse=True,
            literal_seven_coordinate_signatures=literal_checks,
            literal_full_attachment_counts_and_strongness=attachment_checks,
            final_C3_single_bit_boundary=True),
        richer_atom_control=dict(rows=c4_rows, h=5, p=32,
            one_port_marginals=[7,9,6,5],
            scope='Known strong C4 breaks integral normalized atom ratios; no extended-language search performed.'),
        levels=levels,
        reproduction=[
            'python 05-knowledge/results/multiport_defect_grammar_20261103/check_component_language.py',
            'python -O 05-knowledge/results/multiport_defect_grammar_20261103/check_component_language.py'],
        guards='Explicit RuntimeError checks; no Python assertions.')
    (HERE/'component_language.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8', newline='\n')
    print('PASS: complete inverse equals binary-domino inverse; 25 literal signatures and 50 strong full attachments')


if __name__ == '__main__':
    main()
