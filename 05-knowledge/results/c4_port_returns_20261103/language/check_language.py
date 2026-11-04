"""Complete finite-f acceptance for ordered singleton/C3/C4 atoms.

The prefix-product bound is proved in PROOF.md. No chosen order cutoff is
used. All guards survive optimization. Writes only this owned directory.
"""
from pathlib import Path
from itertools import permutations, product
from math import gcd
import hashlib
import importlib.util
import json

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1]
OLD_PATH = RESULTS/'multiport_defect_grammar_20261103'/'check_component_language.py'
spec = importlib.util.spec_from_file_location('inherited_component_checker', OLD_PATH)
old = importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def permutation_signature(rows, ports):
    """Independent small-atom audit by every ordered vertex permutation."""
    n = len(rows)
    paths = [[] for _ in range(1 << n)]
    paths[0] = [()]
    for mask in range(1, 1 << n):
        vertices = [i for i in range(n) if mask >> i & 1]
        for word in permutations(vertices):
            if all(rows[u] >> v & 1 for u,v in zip(word, word[1:])):
                paths[mask].append(word)
    full = (1 << n)-1
    p = f = g = j = 0
    for first in range(full+1):
        for left in paths[first]:
            for right in paths[full ^ first]:
                end = bool(left and ports >> left[-1] & 1)
                start = bool(right and ports >> right[0] & 1)
                p += 1
                f += end
                g += start
                j += end and start
    return (len(paths[full]), p, f, g, j,
            sum(ports >> word[-1] & 1 for word in paths[full]),
            sum(ports >> word[0] & 1 for word in paths[full]))


def make_atoms():
    atoms = []
    for kind,rows,masks in [('T',[0],range(2)), ('C',[2,4,1],[0,1,3,7]),
                            ('Q',[8,1,3,6],range(16))]:
        need(old.strong(rows), ('strong atom', kind))
        for mask in masks:
            sig = old.literal_signature(rows, mask)
            need(sig == permutation_signature(rows, mask), ('atom independent count', kind, mask))
            atoms.append(dict(kind=kind, mask=mask, rows=rows, signature=sig))
    return atoms


ATOMS = make_atoms()
BY_KIND = {kind:tuple(i for i,a in enumerate(ATOMS) if a['kind'] == kind)
           for kind in ('T','C','Q')}
KIND_HP = {kind:ATOMS[ids[0]]['signature'][:2] for kind,ids in BY_KIND.items()}
MIN_F = {kind:min(ATOMS[i]['signature'][2] for i in ids if ATOMS[i]['mask'])
         for kind,ids in BY_KIND.items()}


def prefix_inverse(target):
    """Exact reduction: p(prefix) <= target for every valid final atom."""
    words = []
    states = 0
    def visit(word, sig):
        nonlocal states
        states += 1
        for i,atom in enumerate(ATOMS):
            a = atom['signature']
            result = a if sig is None else old.ordered(sig, a)
            newword = word+(i,)
            if atom['mask'] and result[2] == target:
                words.append(newword)
            # A future nonempty terminal port gives a strictly positive
            # additional term, so f(prefix) must be strictly below target.
            if result[1] <= target and result[2] < target:
                visit(newword, result)
    visit((), None)
    need(len(words) == len(set(words)), 'prefix inverse duplicate')
    return set(words), states


def skeleton_inverse(target):
    """Independent inverse via the nonnegative coefficient equation."""
    found = set()
    skeleton_count = 0
    coefficient_nodes = 0
    gcd_rejections = 0

    def solve(kinds):
        nonlocal skeleton_count,coefficient_nodes,gcd_rejections
        skeleton_count += 1
        weights = []
        for i in range(len(kinds)):
            weight = 1
            for kind in kinds[:i]:
                weight *= KIND_HP[kind][1]
            for kind in kinds[i+1:]:
                weight *= KIND_HP[kind][0]
            weights.append(weight)
        options = [BY_KIND[kind] for kind in kinds]
        options[-1] = tuple(i for i in options[-1] if ATOMS[i]['mask'])
        suffix_gcd = [0]*(len(kinds)+1)
        suffix_min = [0]*(len(kinds)+1)
        suffix_max = [0]*(len(kinds)+1)
        for index in range(len(kinds)-1,-1,-1):
            values = [weights[index]*ATOMS[i]['signature'][2] for i in options[index]]
            g = suffix_gcd[index+1]
            for v in values:
                g = gcd(g,v)
            suffix_gcd[index] = g
            suffix_min[index] = min(values)+suffix_min[index+1]
            suffix_max[index] = max(values)+suffix_max[index+1]

        def choose(index, remaining, word):
            nonlocal coefficient_nodes,gcd_rejections
            coefficient_nodes += 1
            if not suffix_min[index] <= remaining <= suffix_max[index]:
                return
            divisor = suffix_gcd[index]
            if divisor and remaining % divisor:
                gcd_rejections += 1
                return
            if index == len(kinds):
                if remaining == 0:
                    found.add(word)
                return
            for i in options[index]:
                choose(index+1, remaining-weights[index]*ATOMS[i]['signature'][2], word+(i,))
        choose(0,target,())

    def skeletons(prefix, p):
        for kind in ('T','C','Q'):
            if p*MIN_F[kind] <= target:
                solve(prefix+(kind,))
            next_p = p*KIND_HP[kind][1]
            if next_p <= target:
                skeletons(prefix+(kind,), next_p)
    skeletons((),1)
    return found,dict(skeletons=skeleton_count, coefficient_nodes=coefficient_nodes,
                      gcd_rejections=gcd_rejections)


def word_signature(word):
    sig = ATOMS[word[0]]['signature']
    for i in word[1:]:
        sig = old.ordered(sig, ATOMS[i]['signature'])
    return sig


def literal_word(word):
    rows = []
    ports = 0
    for i in word:
        atom = ATOMS[i]
        offset = len(rows)
        fresh = ((1 << len(atom['rows']))-1) << offset
        rows = [r | fresh for r in rows] + [r << offset for r in atom['rows']]
        ports |= atom['mask'] << offset
    return rows, ports


def main():
    expected = {21:(5,35),63:(9,147),189:(18,739),
                126:(12,403),252:(20,1059),504:(26,2819),1008:(42,7555)}
    target_intercepts = {21:[-13,-11], 63:[-37,-35,-31,-29],
                        189:[-109,-107,-103,-101,-91,-89,-85,-83]}
    literal_checks = full_checks = 0
    reports = []
    for target in (21,63,189,126,252,504,1008):
        direct,states = prefix_inverse(target)
        independent,stats = skeleton_inverse(target)
        need(direct == independent, ('two complete inverse implementations', target))
        need((len(direct),states) == expected[target], ('declared universe', target))
        records = []
        for word in sorted(direct):
            sig = word_signature(word)
            rows,ports = literal_word(word)
            need(old.literal_signature(rows, ports) == sig, ('literal full signature', target, word))
            literal_checks += 1
            need(5*sig[1] <= 32*target, ('full product bound', target, word))
            need((1 << len(rows)) <= sig[1], ('vertex bound', target, word))
            if target in (21,63,189):
                for clock in (1,2):
                    full = old.full_attachment(rows, ports, clock)
                    need(old.strong(full), ('literal attachment strongness', target, word))
                    need(old.literal_signature(full,0)[0] == sig[0]-sig[4]+(1 << clock)*target,
                         ('literal full H', target, word))
                    full_checks += 1
            records.append(dict(word=[dict(kind=ATOMS[i]['kind'],mask=ATOMS[i]['mask']) for i in word],
                                signature=sig, intercept=sig[0]-sig[4], rows=rows,
                                port_mask=ports, contains_C4=any(ATOMS[i]['kind']=='Q' for i in word)))
        requested = target_intercepts.get(target,[-35,-29])
        acceptance = {str(c):[i for i,r in enumerate(records) if r['intercept']==c] for c in requested}
        if target in (63,126,252,504,1008):
            need(not acceptance['-35'] and not acceptance['-29'], ('24/44 obstruction', target))
        reports.append(dict(f=target,prefix_states=states,independent_inverse_stats=stats,
                            candidate_count=len(records),candidates=records,
                            distinct_intercepts=sorted(set(r['intercept'] for r in records)),
                            requested_intercepts=acceptance))
        print('f=%d: %d words; %d prefix states; accepted intercepts=%s' %
              (target,len(records),states,[int(c) for c,w in acceptance.items() if w]))
    need((literal_checks,full_checks)==(132,64),'audit universe changed')
    report=dict(status='FINITE_EXACT_COMPLETE_FIXED_MARGINAL_C4_ATOM_LANGUAGE',
        source_pin='1b50a1d763',
        scope='All ordered words over T, C3, and strong C4 with positive terminal ports; fixed-f acceptance is complete without a chosen order cutoff. No physical conversion or native legality inferred.',
        inherited_checker=dict(path=str(OLD_PATH.relative_to(RESULTS)).replace('\\','/'),
            normalized_sha256=hashlib.sha256(OLD_PATH.read_bytes().replace(b'\r\n',b'\n')).hexdigest()),
        atom_bank=ATOMS,
        independent_controls=dict(atom_permutation_signatures=len(ATOMS),
            literal_word_signatures=literal_checks,full_attachment_count_and_strongness=full_checks,
            prefix_inverse_equals_coefficient_inverse=True),
        common_clock_f=[21,63,189],
        shifted_clock=dict(source_f=63,shifts=[1,2,3,4],output_clock='t-shift, valid for t>=shift+1',
                           target_intercepts=[-35,-29],physical_prefix_absorption_proved=False),
        results=reports,
        reproduction=['python 05-knowledge/results/c4_port_returns_20261103/language/check_language.py',
                      'python -O 05-knowledge/results/c4_port_returns_20261103/language/check_language.py'])
    (HERE/'receipt.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS: 22 atom permutation checks; 132 independent word signatures; 64 strong full attachments')


if __name__=='__main__':
    main()
