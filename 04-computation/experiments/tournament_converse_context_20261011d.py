"""Converse, confluent atom fusion, and reversible marked star fans.

The proof is all-order; the inherited rank atlas is explicitly limited to n<=6.
No labelled equivariant bijection is inferred from the unlabelled type codec.
"""
from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations, permutations, product
import json

import tournament_even_bijection_20261011c as old

CHECKS = 0


def need(ok, message):
    if not ok:
        raise ValueError(message)


def check(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(message)


def graph(value):
    need(type(value) is tuple and len(value) == 2, 'exact (order,mask) graph')
    old.source(*value)
    need(old.automorphism_even(*value), 'automorphism-even graph')
    return value


def strong_dual(letter):
    old.word_input((letter,))
    n, mask = letter
    return n, old.canonical(n, mask ^ ((1 << len(old.pairs(n))) - 1), tournament=True)


def dual_word(word):
    old.word_input(word)
    return tuple(strong_dual(letter) for letter in reversed(word))


def factors(value):
    n, mask = graph(value)
    inverse = old._atlas(n)[1]
    return tuple(sorted((inverse[p] for p in old.components(n, mask, tournament=False)),
                        reverse=True))


def blocks_input(blocks):
    need(type(blocks) is tuple, 'exact ordered block tuple')
    for block in blocks:
        old.word_input(block)
        need(old.lyndon(block), 'each block is a nonempty Lyndon word')
    old.word_input(tuple(letter for block in blocks for letter in block))


def normalize(blocks, side='left'):
    """Confluent conditional fusion; each step removes one block boundary."""
    blocks_input(blocks)
    need(type(side) is str and side in ('left', 'right'), 'declared rewrite schedule')
    current = blocks
    while True:
        choices = tuple(i for i in range(len(current)-1) if current[i] < current[i+1])
        if not choices:
            return current
        i = choices[0] if side == 'left' else choices[-1]
        merged = current[i] + current[i+1]
        need(old.lyndon(merged), 'inherited Lyndon concatenation lemma')
        current = current[:i] + (merged,) + current[i+2:]


def from_factors(blocks):
    blocks_input(blocks)
    need(all(x >= y for x, y in zip(blocks, blocks[1:])), 'complete normal form')
    n = sum(old.weight(w) for w in blocks)
    atlas = old._atlas(n)[0]
    return old.assemble(tuple(atlas[w] for w in blocks), tournament=False)


def context_product(left, right):
    graph(left); graph(right)
    return from_factors(normalize(factors(left) + factors(right)))


def union_compatible(left, right):
    a, b = factors(left), factors(right)
    need(left[0]+right[0] <= old.LIMIT, 'combined order in declared atlas')
    return not a or not b or a[-1] >= b[0]


def converse_even(value):
    """Read components, reverse/converse SCC letters, then fuse to normal form."""
    word = tuple(letter for block in factors(value) for letter in block)
    return from_factors(normalize(tuple((letter,) for letter in dual_word(word))))


def self_converse(value):
    word = tuple(letter for block in factors(value) for letter in block)
    return word == dual_word(word)


@dataclass(frozen=True)
class FanFrame:
    n: int
    centre: int
    keep: tuple


def audit_frame(frame):
    need(type(frame) is FanFrame, 'exact fan frame')
    need(type(frame.n) is int and 4 <= frame.n <= old.LIMIT, 'marked fan order4..6')
    need(type(frame.centre) is int and 0 <= frame.centre < frame.n, 'exact centre label')
    need(type(frame.keep) is tuple and len(frame.keep) == 2
         and all(type(v) is int and 0 <= v < frame.n for v in frame.keep),
         'two exact retained leaf labels')
    need(frame.keep[0] < frame.keep[1] and frame.centre not in frame.keep,
         'unordered retained pair in increasing order, separate from centre')


def fan_graph(frame, joined):
    audit_frame(frame)
    need(type(joined) is int and joined in (0, 1), 'exact joining bit')
    leaves = tuple(v for v in range(frame.n) if v != frame.centre) if joined else frame.keep
    edges = {tuple(sorted((frame.centre, v))) for v in leaves}
    return frame.n, sum(1 << i for i, e in enumerate(old.pairs(frame.n)) if e in edges)


def fan_read(frame, value):
    audit_frame(frame); graph(value)
    need(value[0] == frame.n, 'same labelled source order')
    choices = tuple(bit for bit in (0, 1) if fan_graph(frame, bit) == value)
    need(len(choices) == 1, 'graph must match the same marked fan frame')
    return choices[0]


def fan_toggle(frame, value):
    return fan_graph(frame, 1-fan_read(frame, value))


def fan_mark(value, keep=None):
    """Path+isolates determines its frame; the star requires its retained pair."""
    n, mask = graph(value)
    need(4 <= n <= old.LIMIT, 'marked fan order4..6')
    matrix = old.adjacency(n, mask, False)
    degrees = tuple(sum(row) for row in matrix)
    if sorted(degrees) == [0]*(n-3) + [1, 1, 2]:
        centre = degrees.index(2)
        actual = tuple(v for v in range(n) if matrix[centre][v])
        chosen = actual if keep is None else keep
    else:
        need(sorted(degrees) == [1]*(n-1) + [n-1], 'star or P3 plus isolates')
        centre = degrees.index(n-1)
        need(keep is not None, 'a star must retain its two-leaf sidecar')
        chosen = keep
    frame = FanFrame(n, centre, chosen)
    audit_frame(frame); fan_read(frame, value)
    return frame


def relabel_frame(frame, labels):
    audit_frame(frame)
    need(type(labels) is tuple and len(labels) == frame.n
         and all(type(v) is int for v in labels) and sorted(labels) == list(range(frame.n)),
         'exact frame transporter')
    return FanFrame(frame.n, labels[frame.centre], tuple(sorted(labels[v] for v in frame.keep)))


def relabel_graph(value, labels):
    n, mask = graph(value)
    need(type(labels) is tuple and len(labels) == n
         and all(type(v) is int for v in labels) and sorted(labels) == list(range(n)),
         'exact graph transporter')
    index = {p:i for i,p in enumerate(old.pairs(n))}
    return n, sum(1 << index[tuple(sorted((labels[u], labels[v])))]
                  for i,(u,v) in enumerate(old.pairs(n)) if mask >> i & 1)


@lru_cache(None)
def all_normal_forms(blocks):
    answers = set()
    for i in range(len(blocks)-1):
        if blocks[i] < blocks[i+1]:
            answers.update(all_normal_forms(blocks[:i] + (blocks[i]+blocks[i+1],) + blocks[i+2:]))
    return frozenset(answers) if answers else frozenset((blocks,))


def all_lyndon_cuts(word):
    if not word:
        return ((),)
    out = []
    for mask in range(1 << (len(word)-1)):
        cut, start = [], 0
        for i in range(1,len(word)):
            if mask >> (i-1) & 1:
                cut.append(word[start:i]); start=i
        cut.append(word[start:]); cut=tuple(cut)
        if all(old.lyndon(b) for b in cut):
            out.append(cut)
    return tuple(out)


def reject(fn):
    try:
        fn()
    except ValueError:
        check(True, 'typed false-domain input rejected')
    else:
        check(False, 'required rejection')


def main():
    catalogue, counts, pairs_by_order = {}, [], {}
    schedules = 0
    for n in range(old.LIMIT+1):
        catalogue[n] = tuple((n,m) for m in old._catalogue(n,False)[0]
                             if old.automorphism_even(n,m))
        fixed, strong_fixed, exchanged = 0, 0, []
        for value in catalogue[n]:
            dual = converse_even(value)
            tournament = old.even_to_tournament(*value)
            direct = old.tournament_to_even(n, tournament[1] ^ ((1 << len(old.pairs(n)))-1))
            check(dual == direct, 'independent literal-arc reversal matches rewrite carrier')
            check(converse_even(dual) == value, 'converse is an involution on graph types')
            check(self_converse(value) == (dual == value), 'exact SCC palindrome criterion')
            fixed += value == dual
            if value < dual:
                exchanged.append((value[1],dual[1]))
            word=tuple(letter for block in factors(value) for letter in block)
            for cut in all_lyndon_cuts(word):
                nf=old.factor(word)
                check(all_normal_forms(cut) == frozenset((nf,)), 'all rewrite branches converge')
                check(normalize(cut,'left') == normalize(cut,'right') == nf, 'opposite schedules agree')
                schedules += 1
        for m in old._atlas(n)[2].get(n,()):
            strong_fixed += strong_dual((n,m)) == (n,m)
        counts.append({'n':n,'types':len(catalogue[n]),'fixed':fixed,'strong_fixed':strong_fixed})
        pairs_by_order[n] = exchanged
    # Independent generating-series convolution for palindromic SCC words.
    for n,row in enumerate(counts):
        expected = len(catalogue[n//2]) if n%2==0 else 0
        for centre in range(1,n+1):
            if (n-centre)%2==0:
                expected += counts[centre]['strong_fixed'] * len(catalogue[(n-centre)//2])
        check(expected == row['fixed'], 'central fixed letter and arbitrary half-word count')
    pair_controls=triple_controls=0
    for n in range(old.LIMIT+1):
        for i in range(n+1):
            for g,h in product(catalogue[i],catalogue[n-i]):
                joined=context_product(g,h)
                union=old.assemble(tuple(x for x in (g,h) if x[0]),tournament=False)
                check((joined==union) == union_compatible(g,h), 'sharp one-seam union criterion')
                check(converse_even(joined) == context_product(converse_even(h),converse_even(g)),
                      'converse reverses ordered context product')
                pair_controls += 1
        for i in range(n+1):
            for j in range(n-i+1):
                for g,h,k in product(catalogue[i],catalogue[j],catalogue[n-i-j]):
                    check(context_product(context_product(g,h),k)
                          == context_product(g,context_product(h,k)), 'context product associative')
                    triple_controls += 1
    # Exact all-order star formula, controlled in the available atlas.
    a=(1,0); b=(3,old._atlas(3)[2][3][0])
    for n in range(3,old.LIMIT+1):
        r=n-3
        star=(n,(1<<(n-1))-1)
        separate=old.assemble(((3,3),)+((1,0),)*r,tournament=False)
        word=(a,)*r+(b,)
        check(old._atlas(n)[0][word] == star, 'first Lyndon atom is exactly the star')
        check(converse_even(star)==separate and converse_even(separate)==star,
              'direct all-order fan operation is the actual chosen converse')
    check(converse_even((4,3)) == (4,7), 'minimal componentwise-converse hostile')
    check(converse_even((3,3)) == (3,3) and converse_even((1,0)) == (1,0),
          'both components separately fixed but their union is not')
    check(strong_dual((6,18)) == (6,22), 'reversal must also converse each strong letter')
    # The fan coincidence does not extend to chromatic-number preservation.
    source_edges={e for i,e in enumerate(old.pairs(6)) if 110 >> i & 1}
    target_edges={e for i,e in enumerate(old.pairs(6)) if 123 >> i & 1}
    two=(0,0,1,1,1,0); three=(0,1,2,0,1,1)
    check(converse_even((6,110))==(6,123) and bool(source_edges)
          and all(two[u]!=two[v] for u,v in source_edges)
          and {(0,1),(0,2),(1,2)} <= target_edges
          and all(three[u]!=three[v] for u,v in target_edges),
          'explicit converse pair has chromatic numbers2 and3')
    # Complete labelled fan-frame orbit, all source transports.
    fan_counts=[]; covariances=0
    for n in range(4,old.LIMIT+1):
        fibres=Counter(); paths=set()
        for centre in range(n):
            for keep in combinations(tuple(v for v in range(n) if v!=centre),2):
                frame=FanFrame(n,centre,keep)
                path,star=fan_graph(frame,0),fan_graph(frame,1)
                paths.add(path); fibres[star]+=1
                check(fan_mark(path)==frame and fan_mark(star,keep)==frame, 'frame reconstruction')
                check(fan_toggle(frame,path)==star and fan_toggle(frame,star)==path, 'labelled inverse pair')
                check(old.canonical(n,star[1],tournament=False)==converse_even(path)[1],
                      'marked toggle implements the unmarked type involution')
                for p in permutations(range(n)):
                    moved=relabel_frame(frame,p)
                    check(fan_graph(moved,0)==relabel_graph(path,p)
                          and fan_graph(moved,1)==relabel_graph(star,p),
                          'marked fan operation is relabelling-equivariant')
                    covariances+=1
        multiplicity=(n-1)*(n-2)//2
        check(len(paths)==n*multiplicity and len(fibres)==n
              and set(fibres.values())=={multiplicity}, 'exact labelled sidecar fibre size')
        # Every unmarked star automorphism acts transitively on retained pairs.
        target=FanFrame(n,0,(1,2))
        orbit={relabel_frame(target,(0,)+p).keep for p in permutations(range(1,n))}
        check(len(orbit)==multiplicity, 'no invariant retained pair on an unmarked star')
        fan_counts.append({'n':n,'paths_with_isolates':len(paths),'stars':len(fibres),
                           'retained_pair_fibre':multiplicity})
    # State-dependent edge toggling requires the same frame, not just graph type.
    frame=FanFrame(4,0,(1,2)); path=fan_graph(frame,0); star=fan_graph(frame,1)
    reject(lambda: fan_mark(star))
    reject(lambda: fan_read(FanFrame(4,0,(1,3)),path))
    reject(lambda: fan_graph(FanFrame(4,False,(1,2)),0))
    reject(lambda: fan_graph(FanFrame(4,0,(True,2)),0))
    reject(lambda: fan_graph(frame,True))
    reject(lambda: fan_mark((4,3.0)))
    reject(lambda: converse_even((True,0)))
    reject(lambda: converse_even((3,7)))
    reject(lambda: strong_dual((3,0)))
    reject(lambda: normalize(((a,a),)))
    reject(lambda: normalize(((a,),),'greedy'))
    reject(lambda: from_factors(((a,),(b,))))
    reject(lambda: context_product((6,0),(1,0)))
    reject(lambda: relabel_frame(frame,(0,1,2,3.0)))
    print(json.dumps({'scope':'type converse via confluent factors; labelled marked fan lift',
                      'checks':CHECKS,'complete_orders':'0..6','counts':counts,
                      'exchanged_graph_masks':pairs_by_order,
                      'lyndon_cut_schedule_controls':schedules,
                      'context_pairs':pair_controls,'context_triples':triple_controls,
                      'chromatic_hostile':{'source_mask':110,'target_mask':123,'chi':[2,3]},
                      'labelled_fans':fan_counts,'frame_covariance_controls':covariances},indent=2))


if __name__=='__main__':
    main()
