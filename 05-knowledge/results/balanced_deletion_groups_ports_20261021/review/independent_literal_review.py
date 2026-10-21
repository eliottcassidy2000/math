#!/usr/bin/env python3
"""Independent review: enumerate literal directed words, not producer states.

No import of a producer or repository file. The cut defect is measured by
XORing two independently computed target tournaments, then checked against
gap endpoints. All validations survive python -O.
"""
from collections import Counter, defaultdict
from itertools import combinations, product
import json


def require(value, message):
    if not value:
        raise RuntimeError(message)


def matrix(k):
    n = 2*k+1
    out = [[a < b for b in range(n)] for a in range(n)]
    for a,b in ((0,k),(k,2*k)):
        out[a][b], out[b][a] = False, True
    return out


def directed_words(out):
    n = len(out)
    followers = [sum(1 << b for b in range(n) if out[a][b]) for a in range(n)]
    all_bits = (1 << n)-1
    def recurse(word, remaining):
        if not remaining:
            yield tuple(word)
            return
        allowed = remaining if not word else remaining & followers[word[-1]]
        while allowed:
            b = allowed & -allowed
            allowed -= b
            v = b.bit_length()-1
            word.append(v)
            yield from recurse(word, remaining-b)
            word.pop()
    yield from recurse([], all_bits)


def edge_set(word):
    return {tuple(sorted((a,b))) for a,b in zip(word,word[1:])}


def reverse_path(out, word):
    result = [row.copy() for row in out]
    for a,b in zip(word,word[1:]):
        result[a][b], result[b][a] = result[b][a], result[a][b]
    return result


def restricted(out, survivors):
    return [[out[a][b] for b in survivors] for a in survivors]


def gap_ports(word, survivors):
    rank = {v:i for i,v in enumerate(survivors)}
    selected = [(pos,v) for pos,v in enumerate(word) if v in rank]
    return {tuple(sorted((rank[a],rank[b])))
            for (i,a),(j,b) in zip(selected,selected[1:]) if j>i+1}


def blue(out, word):
    reverse = dict(zip(word,reversed(word)))
    return all(out[a][b] == out[reverse[b]][reverse[a]]
               for a in range(len(out)) for b in range(len(out)) if a != b)


def toggle(word, k, axis):
    """Move a literal inserted vertex between its two sorted shore runs."""
    centre = k+1
    cut = word.index(centre)
    before, after = word[:cut], word[cut+1:]
    la = [v for v in before if v < centre]
    lb = [v for v in before if v > centre]
    ra = [v for v in after if v < centre]
    rb = [v for v in after if v > centre]
    if axis == 0:
        if not ra:  # The entire left shore is before the centre.
            return word
        old,new = (la,ra) if k in la else (ra,la)
        old.remove(k)
        new.append(k)
        new.sort()
    else:
        if not lb:  # The entire right shore is after the centre.
            return word
        old,new = (lb,rb) if k+2 in lb else (rb,lb)
        old.remove(k+2)
        new.append(k+2)
        new.sort()
    return tuple(la+lb+[centre]+ra+rb)


def g_action(word,k,g):
    for a in range(2):
        if g & (1 << a):
            word = toggle(word,k,a)
    return word


def reflected(word):
    n = len(word)
    return tuple(n-1-v for v in reversed(word))


def balanced_checks():
    reports = []
    for k in range(2,7):
        out = matrix(k+1)
        n = len(out)
        survivors = [v for v in range(n) if v not in (k,k+2)]
        rank = {v:i for i,v in enumerate(survivors)}
        parent = restricted(out,survivors)
        require(parent == matrix(k), 'restricted tournament mismatch')
        words = tuple(directed_words(out))
        word_set = set(words)
        q = 2**(k-1)+1
        require(len(words) == (2*q-1)**2, 'literal DFS count mismatch')
        fibres = defaultdict(set)
        defects = Counter()
        e_a_c = Counter()
        stabilizers = Counter()
        blue_total = 0
        targets = {}
        for word in words:
            smaller = tuple(rank[v] for v in word if v in rank)
            require(all(parent[a][b] for a,b in zip(smaller,smaller[1:])),
                    'deletion not a directed path')
            fibres[smaller].add(word)
            child_target = restricted(reverse_path(out,word),survivors)
            parent_target = reverse_path(parent,smaller)
            defect = { (a,b) for a,b in combinations(range(len(survivors)),2)
                       if child_target[a][b] != parent_target[a][b] }
            ports = gap_ports(word,survivors)
            require(defect == ports, 'literal target defect not equal gap ports')
            require(ports <= edge_set(smaller), 'port not a parent edge')
            # Derive the two candidate seams solely from the retained word.
            pivot = smaller.index(k)
            ia = [v for v in smaller[:pivot] if v < k]
            jb = [v for v in smaller[:pivot] if v > k]
            ac = [v for v in smaller[pivot+1:] if v < k]
            bc = [v for v in smaller[pivot+1:] if v > k]
            pre = (ia[-1] if ia else None, jb[0] if jb else k)
            post = (ac[-1] if ac else k, bc[0] if bc else None)
            epsilon = word.index(k) < word.index(k+1)
            eta = word.index(k+2) > word.index(k+1)
            seams = {tuple(sorted(pair)) for pair,occupied in
                     ((pre,epsilon or not eta),(post,not epsilon or eta))
                     if occupied and None not in pair}
            require(seams == defect, 'literal two-seam reconstruction')
            defects[len(defect)] += 1
            e = int(word[0] in (k,k+2))+int(word[-1] in (k,k+2))
            a = int(abs(word.index(k)-word.index(k+2)) == 1)
            require(len(defect) == 2-e-a, 'endpoint adjacency census equation')
            e_a_c[e,a,len(defect)] += 1
            orbit = {g_action(word,k,g) for g in range(4)}
            require(orbit <= word_set, 'literal toggle is not a legal path')
            require(all(tuple(rank[v] for v in w if v in rank) == smaller for w in orbit),
                    'toggle changed deletion')
            stab = tuple(g for g in range(4) if g_action(word,k,g) == word)
            stabilizers[stab] += 1
            left_inf, right_inf = word[0] == 0, word[-1] == n-1
            expected = tuple(g for g in range(4)
                             if (not (g&1) or left_inf) and (not (g&2) or right_inf))
            require(stab == expected, 'wrong stabilizer line')
            for g,h in product(range(4),repeat=2):
                require(g_action(g_action(word,k,g),k,h) == g_action(word,k,g^h),
                        'V4 group law')
            require(reflected(toggle(reflected(word),k,0)) == toggle(word,k,1),
                    'reflection fails to swap axes')
            is_blue = blue(out,word)
            require(is_blue == (word == reflected(word)), 'blue reflection discrepancy')
            blue_total += is_blue
            targets[word] = child_target
        require(blue_total == 2*q-1, 'blue count')
        require(len(fibres) == q*q, 'deletion surjectivity cardinality')
        for members in fibres.values():
            word = next(iter(members))
            require(members == {g_action(word,k,g) for g in range(4)}, 'fibre is not orbit')
        require(defects == Counter({0:2*q,1:2*q*q-4,2:2*q*q-6*q+5}),
                'all-order defect formula differs from literal DFS')
        s = q-1
        require(sum(count*a for (e,a,c),count in e_a_c.items()) == 2*s*(s+1),
                'adjacency count')
        require(sum(count*e for (e,a,c),count in e_a_c.items()) == 6*s+2,
                'endpoint incidence count')
        require(sum(count for (e,a,c),count in e_a_c.items() if e == 2) == 2,
                'double endpoint count')
        require(sum(count for (e,a,c),count in e_a_c.items() if e and a) == 2*s,
                'endpoint adjacency overlap count')
        group_permutations = {tuple(reflected(g_action(w,k,g)) if r else g_action(w,k,g)
                                    for w in words) for g,r in product(range(4),range(2))}
        require(len(group_permutations) == 8, 'D8 not faithful')
        standard = tuple(range(n))
        require(gap_ports(standard,survivors) == {(k-1,k),(k,k+1)}, 'standard defect')
        report = dict(parent_k=k, source_order=n, literal_DFS_paths=len(words),
                      defect_census=dict(sorted(defects.items())),
                      stabilizer_path_counts={str(a):b for a,b in sorted(stabilizers.items())},
                      endpoint_adjacency_defect_counts={str(a):b for a,b in sorted(e_a_c.items())},
                      dihedral_permutations=len(group_permutations))
        if k == 2:
            hostile = ((4,6,3,0,1,2,5),(2,4,6,3,0,1,5))
            require(toggle(hostile[0],2,0) == hostile[1], 'hostile toggle')
            require(all(not blue(out,w) for w in hostile), 'hostile not both black')
            spectra = [sorted(map(sum,targets[w])) for w in hostile]
            require(spectra == [[0,1,2,3,4],[0,2,2,3,3]], 'hostile score vectors')
            require(spectra[0] != sorted(4-a for a in spectra[1]), 'converse collision')
            report['hostile_scores'] = spectra
            # Do not promote a labelled non-neighbour to an unlabelled one.
            parent_words = tuple(directed_words(parent))
            transitive_parent_paths = []
            for w in parent_words:
                r_target = reverse_path(parent,w)
                k_target = [[False if i == j else not r_target[i][j]
                             for j in range(5)] for i in range(5)]
                require(r_target != targets[hostile[0]], 'hostile unexpectedly labelled R-neighbour')
                require(k_target != targets[hostile[0]], 'hostile unexpectedly labelled K-neighbour')
                if sorted(map(sum,r_target)) == list(range(5)):
                    transitive_parent_paths.append(w)
            require((3,4,2,0,1) in transitive_parent_paths,
                    'lost unlabelled neighbour countercontrol')
            report['hostile_labelled_R_K_nonneighbour'] = True
            report['unlabelled_transitive_neighbour_parent_paths'] = transitive_parent_paths
        reports.append(report)
    return reports


def general_defect_checks():
    """All n=5 labelled tournaments/HPs and every nontrivial survivor set."""
    n = 5
    pairs = list(combinations(range(n),2))
    cases = 0
    histogram = Counter()
    for mask in range(1 << len(pairs)):
        out = [[False]*n for _ in range(n)]
        for bit,(a,b) in enumerate(pairs):
            out[a][b],out[b][a] = bool(mask & (1<<bit)),not bool(mask & (1<<bit))
        for word in directed_words(out):
            child_target = reverse_path(out,word)
            for survivors_mask in range(1,1<<n):
                survivors = [a for a in range(n) if survivors_mask & (1<<a)]
                rank = {v:i for i,v in enumerate(survivors)}
                smaller = tuple(rank[v] for v in word if v in rank)
                base = restricted(out,survivors)
                if not all(base[a][b] for a,b in zip(smaller,smaller[1:])):
                    continue
                lhs = restricted(child_target,survivors)
                rhs = reverse_path(base,smaller)
                difference = {(a,b) for a,b in combinations(range(len(survivors)),2)
                              if lhs[a][b] != rhs[a][b]}
                ports = gap_ports(word,survivors)
                require(difference == ports, 'general exact defect identity')
                require(len(ports) <= n-len(survivors), 'general r-port bound')
                # The complement/converse has exactly the same discrepancy.
                require(all((not lhs[a][b]) != (not rhs[a][b])
                            for a,b in difference), 'K converse discrepancy')
                cases += 1
                histogram[len(ports)] += 1
    return dict(order=n,tournaments=1<<len(pairs),valid_path_deletions=cases,
                defect_histogram=dict(sorted(histogram.items())))


def all_depth_checks():
    reports = []
    for k,d in product(range(2,5),range(1,4)):
        c = k+d
        out = matrix(c)
        left_new = set(range(k,c))
        right_new = set(range(c+1,k+2*d+1))
        survivors = [v for v in range(2*c+1) if v not in left_new|right_new]
        rank = {v:i for i,v in enumerate(survivors)}
        require(restricted(out,survivors) == matrix(k), 'all-depth underlying deletion')
        defects = Counter()
        boundary_types = defaultdict(Counter)
        paths = 0
        for word in directed_words(out):
            pivot = word.index(c)
            E = set(word[:pivot]) & left_new
            H = {2*c-v for v in word[pivot+1:] if v in right_new}
            smaller = tuple(rank[v] for v in word if v in rank)
            pc = smaller.index(k)
            ia = [v for v in smaller[:pc] if v < k]
            jb = [v for v in smaller[:pc] if v > k]
            ac = [v for v in smaller[pc+1:] if v < k]
            bc = [v for v in smaller[pc+1:] if v > k]
            pre = (ia[-1] if ia else None,jb[0] if jb else k)
            post = (ac[-1] if ac else k,bc[0] if bc else None)
            seams = {tuple(sorted(pair)) for pair,occupied in
                     ((pre,bool(E) or H != left_new),(post,E != left_new or bool(H)))
                     if occupied and None not in pair}
            require(seams == gap_ports(word,survivors), 'all-depth seam criterion')
            require(len(seams) <= 2, 'all-depth two-port bound')
            lhs = restricted(reverse_path(out,word),survivors)
            rhs = reverse_path(matrix(k),smaller)
            difference = {(a,b) for a,b in combinations(range(2*k+1),2)
                          if lhs[a][b] != rhs[a][b]}
            require(difference == seams, 'all-depth literal target correction')
            paths += 1
            defects[len(seams)] += 1
            left_type = 'empty' if not ia else 'infinity' if 0 in ia else 'nonempty'
            right_type = 'empty' if not bc else 'infinity' if not jb else 'nonempty'
            boundary_types[left_type,right_type][len(seams)] += 1
        q,D = 2**(k-1)+1,2**d
        expected = Counter({0:D*D+2*(q-2),
                            1:2*(q-2)*D*D+2*D+2*(q-2)**2,
                            2:(q-2)**2*(D*D-2)+2*(q-2)*(D-1)+1})
        require(defects == expected, 'all-depth census A1')
        rows = {
            ('empty','empty'): {0:D*D},
            ('empty','nonempty'): {0:q-2,1:(q-2)*(D*D-1)},
            ('empty','infinity'): {1:D},
            ('nonempty','nonempty'): {1:2*(q-2)**2,2:(q-2)**2*(D*D-2)},
            ('nonempty','infinity'): {1:q-2,2:(q-2)*(D-1)},
            ('infinity','infinity'): {2:1}}
        for types,counts in rows.items():
            require(boundary_types[types] == Counter(counts), 'all-depth boundary row')
            require(boundary_types[types[::-1]] == Counter(counts), 'all-depth reflected row')
        reports.append(dict(parent_k=k,depth=d,paths=paths,
                            defects=dict(sorted(defects.items())),
                            boundary_strata_checked=len(boundary_types)))
    return reports


if __name__ == '__main__':
    result = dict(status='PASS independent literal DFS and target-XOR audit',
                  balanced=balanced_checks(), general=general_defect_checks(),
                  all_depth=all_depth_checks())
    print(json.dumps(result,indent=2,sort_keys=True))
