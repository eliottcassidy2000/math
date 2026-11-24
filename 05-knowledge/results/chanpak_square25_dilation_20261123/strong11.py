"""Three apex parents at the sharp order for strong transitive precedence.

--census verifies the catalog of 911 ten-vertex parents and tests every apex mask.
"""
from collections import Counter
from functools import lru_cache
from itertools import combinations
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/"chanpak_tiling_viewpoint_20261123"))
from precedence_bank import count_and_order


PARENTS = {
    153: (990, 1020, 1016, 1008, 992, 705, 896, 768, 32, 256),
    1563: (542, 1020, 1016, 1008, 864, 705, 385, 785, 33, 320),
    1243: (654, 1020, 1016, 1008, 737, 321, 897, 800, 529, 32),
}
# Each tuple specifies the old height, the new apex's outneighbor mask, and the
# exact new height. Pick one witness for each odd congruence class modulo 24.
SELECTOR = [
    (153, 256, 289), (1243, 1008, 2547), (153, 8, 1349), (1243, 8, 7783),
    (153, 258, 2721), (1563, 1000, 5363), (153, 510, 2797), (1243, 1000, 3807),
    (153, 760, 929), (1563, 1008, 3187), (153, 260, 1821), (1563, 1010, 6191),
]


def children(rows, wins):
    apex = 1 << len(rows)
    return tuple(r if wins>>i&1 else r|apex for i,r in enumerate(rows)) + (wins,)


def strong(rows):
    full = (1 << len(rows))-1
    reach = [r | (1<<i) for i,r in enumerate(rows)]
    for v in range(len(rows)):
        for w in range(len(rows)):
            if reach[w]>>v&1:
                reach[w] |= reach[v]
    return all(r == full for r in reach)


def decode(encoded):
    encoded = int(encoded,16)
    size = 10
    bits = f"{encoded:045b}"
    assert len(bits) == 45
    rows = [0]*size
    index = 0
    for i in range(size):
        for j in range(i+1,size):
            if bits[index] == "1":
                rows[i] |= 1 << j
            else:
                rows[j] |= 1 << i
            index += 1
    return tuple(rows)


def positions(bits):
    while bits:
        bit = bits & -bits
        bits ^= bit
        yield bit.bit_length()-1


def end_table(rows):
    size = len(rows)
    full = (1 << size)-1
    incoming = [full ^ row ^ (1 << i) for i,row in enumerate(rows)]
    table = [[0]*size for _ in range(full + 1)]
    for mask in range(1, full+1):
        for v in positions(mask):
            rem = mask ^ (1 << v)
            table[mask][v] = sum(table[rem][j] for j in positions(rem & incoming[v])) if rem else 1
    return table


def apex_heights(rows):
    # Split a word at the apex. For a nonempty prefix and suffix their one
    # meeting pair (i,j) requires the apex to lose to i and beat j.
    size=len(rows)
    full=(1 << size)-1
    ends=end_table(rows)
    starts=end_table([full ^ row ^ (1 << i) for i,row in enumerate(rows)])
    cross=[[0]*size for _ in range(size)]
    for mask in range(1,full):
        other=full ^ mask
        for i in positions(mask):
            left=ends[mask][i]
            if left:
                for j in positions(other):
                    cross[i][j]+=left*starts[other][j]
    linear=[starts[full][j]-ends[full][j]+sum(cross[i][j] for i in range(size))
            for j in range(size)]
    out=[sum(ends[full])]*(full+1)
    for mask in range(1,full+1):
        bit=mask & -mask
        j=bit.bit_length()-1
        rem=mask ^ bit
        out[mask]=out[rem]+linear[j]-sum(cross[i][j]+cross[j][i] for i in positions(rem))
    return out


def transitive_children(rows):
    # One Python integer tests the 1024 apex masks simultaneously.
    n=len(rows)
    lanes=1 << n
    all_lanes=(1 << lanes)-1
    size=n+1
    arcs=[[0]*size for _ in range(size)]
    for i,row in enumerate(rows):
        for j in positions(row):
            arcs[i][j]=all_lanes
        wins=sum(1 << mask for mask in range(lanes) if mask>>i&1)
        arcs[n][i]=wins
        arcs[i][n]=all_lanes ^ wins
    cycles=[0]*(1 << size)
    ends=[[0]*size for _ in cycles]
    for vertices in range(1,len(cycles)):
        root_bit=vertices & -vertices
        root=root_bit.bit_length()-1
        others=vertices ^ root_bit
        for j in positions(others):
            rem=others ^ (1 << j)
            subtotal=0 if rem else arcs[root][j]
            smaller=ends[rem | root_bit]
            for v in positions(rem):
                subtotal ^= smaller[v] & arcs[v][j]
            ends[vertices][j]=subtotal
            cycles[vertices] ^= subtotal & arcs[j][root]
    # Superset transform sums the cyclic supports containing a specified pair.
    for i in range(size):
        bit=1 << i
        for start in range(0,len(cycles),2*bit):
            for s in range(start,start+bit):
                cycles[s] ^= cycles[s+bit]
    precedence=[[arcs[i][j] ^ cycles[(1<<i)|(1<<j)] for j in range(size)] for i in range(size)]
    survivors=all_lanes
    for i in range(size):
        for j in range(i+1,size):
            for k in range(j+1,size):
                a,b,c=precedence[i][j],precedence[j][k],precedence[i][k]
                survivors &= ~(a^c) | (a^b)
    return positions(survivors)


@lru_cache(None)
def canonical_ancestor(rows):
    if not strong(rows):
        return None
    n=len(rows)
    full=(1 << n)-1
    precede=[0]*n

    def walks(used, choices):
        if used == full:
            return 1
        odd=0
        for v in positions(choices):
            updated=used | (1<<v)
            ways=walks(updated,rows[v] & ~updated)
            if ways:
                precede[v] ^= used
                odd ^= 1
        return odd

    assert walks(0,full) == 1
    order=sorted(range(n),key=lambda v:precede[v].bit_count())
    if [precede[v].bit_count() for v in order] != list(range(n)):
        return None
    return tuple(sum(((rows[v]>>w)&1)<<j for j,w in enumerate(order)) for v in order)


def ancestry(current):
    all_odd=set(range(1,24,2))
    seven={1,5,9,13,17,19,21}
    for smaller,expected in ((9,126),(8,29),(7,12),(6,3),(5,1),(4,1)):
        below={}
        for rows,support in current.items():
            for v in range(len(rows)):
                low=(1<<v)-1
                deleted=tuple((r & low) | ((r>>(v+1))<<v) for j,r in enumerate(rows) if j!=v)
                canonical=canonical_ancestor(deleted)
                if canonical is not None:
                    below.setdefault(canonical,set()).update(support)
        total=set().union(*below.values())
        assert len(below) == expected and total == (all_odd if smaller==9 else seven)
        if smaller == 9:
            counts=Counter(tuple(sorted(vals)) for vals in below.values())
            split=(1,5,9,13,17,21)
            mixed=(1,5,9,13,17,19,21)
            rare_a=(11,15,19,23)
            rare_b=(3,7,15)
            assert all(len(set(x+y)) < 12 for x,y in combinations(counts,2))
            covers={frozenset((x,y,z)) for x,y,z in combinations(counts,3)
                    if len(set(x+y+z))==12}
            assert covers == {frozenset((split,rare_a,rare_b)),frozenset((mixed,rare_a,rare_b))}
            assert [counts[key] for key in (split,mixed,rare_a,rare_b)] == [78,8,2,2]
            for key,expected_h in ((rare_a,751),(rare_b,591)):
                assert {count_and_order(t)[0] for t,res in below.items()
                        if tuple(sorted(res)) == key} == {expected_h}
            print("optimal two-step nine-vertex parent triples:", (78+8)*2*2)
        print("hereditary size",smaller,"ancestors",len(below),"union",sorted(total))
        current=below


def census():
    words=(HERE/"strong11_parents.txt").read_text().splitlines()
    assert len(set(words)) == len(words) == 911
    old_bins=Counter()
    breadth=Counter()
    six=set()
    profiles=Counter()
    strong_profiles=Counter()
    heights_by_profile={}
    strong_parents=0
    marks=0
    ancestry_roots={}
    for word in words:
        rows=decode(word)
        h,order=count_and_order(rows)
        assert order == list(range(10))
        old_bins[h%24]+=1
        parent_is_strong=strong(rows)
        strong_parents+=parent_is_strong
        heights=apex_heights(rows)
        support=set()
        for mask in transitive_children(rows):
            if strong(children(rows,mask)):
                support.add(heights[mask]%24)
                marks+=1
        signature=tuple(sorted(support))
        profiles[signature]+=1
        strong_profiles[signature]+=int(parent_is_strong)
        heights_by_profile.setdefault(signature,set()).add(h)
        if parent_is_strong:
            ancestry_roots[rows]=support
        breadth[len(support)]+=1
        if len(support)==6:
            six.add(tuple(sorted(support)))
    assert [old_bins[r] for r in range(1,24,2)] == [186,2,123,8,170,0,131,4,134,2,141,10]
    assert breadth == {1:112,2:156,3:167,4:175,5:164,6:137}
    assert six == {(1,5,9,13,17,21)}
    assert (strong_parents, marks) == (468,11203)
    split=(1,5,9,13,17,21)
    four=(11,15,19,23)
    rare_a=(3,7,15)
    rare_b=(3,7,11)
    assert max(len(set(x+y)) for x,y in combinations(profiles,2)) == 10
    covers=[(x,y,z) for x,y,z in combinations(profiles,3) if len(set(x+y+z))==12]
    assert {frozenset(c) for c in covers} == {frozenset((split,four,rare_a)),
                                               frozenset((split,four,rare_b))}
    assert [(profiles[x], strong_profiles[x]) for x in (split,four,rare_a,rare_b)] == [
        (137,70), (2,2), (2,2), (2,0)]
    assert [heights_by_profile[x] for x in (four,rare_a,rare_b)] == [{1563},{1243},{591}]
    print("optimal unordered parent triples, arbitrary and strong:",137*2*4,70*2*2)
    print("parents by residue support size", sorted(breadth.items()), "strong masks", marks)
    ancestry(ancestry_roots)


if __name__ == "__main__":
    for old, rows in PARENTS.items():
        assert count_and_order(rows)[0] == old and strong(rows)
    for r, (old,mask,new) in enumerate(SELECTOR):
        rows=children(PARENTS[old],mask)
        assert count_and_order(rows)[0] == new and new%24 == 2*r+1 and strong(rows)
    print("strong orders 10 -> 11; all twelve certified heights", [h for _,_,h in SELECTOR])
    if sys.argv[1:] == ["--census"]:
        census()
