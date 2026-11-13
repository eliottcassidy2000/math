"""Literal all-module core handoff; small named controls, not a spectrum scan."""
from pathlib import Path
import json
from itertools import combinations

INHERITED_REF = "06847e4c6e6709585bf1b4cfa6f84713081034ea"
S = (2, 20, 9, 35, 45, 7)


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def tt(n):
    return tuple(((1 << n) - 1) ^ ((1 << (i+1)) - 1) for i in range(n))


def ordinal(*blocks):
    n = sum(map(len, blocks))
    rows, offset = [], 0
    for block in blocks:
        later = ((1 << n) - 1) ^ ((1 << (offset+len(block))) - 1)
        rows.extend((row << offset) | later for row in block)
        offset += len(block)
    return tuple(rows)


def substitute(frame, blocks):
    slots, offset = [], 0
    for block in blocks:
        slots.append(tuple(range(offset, offset+len(block))))
        offset += len(block)
    rows = [0] * offset
    for i, block in enumerate(blocks):
        for local, vertex in enumerate(slots[i]):
            rows[vertex] = block[local] << slots[i][0]
            for j, targets in enumerate(slots):
                if frame[i] >> j & 1:
                    rows[vertex] |= sum(1 << v for v in targets)
    return tuple(rows), slots


def pi_frame():
    rows = list(tt(9))
    for start, end in ((4, 0), (5, 0), (7, 2), (8, 2)):
        rows[end] ^= 1 << start
        rows[start] ^= 1 << end
    return tuple(rows)


def inventory(rows):
    n = len(rows)
    dp = [[0]*n for _ in range(1 << n)]
    H = [1] + [0]*((1 << n)-1)
    for mask in range(1, 1 << n):
        for v in range(n):
            if mask >> v & 1:
                before = mask ^ (1 << v)
                dp[mask][v] = 1 if not before else sum(
                    dp[before][u] for u in range(n) if rows[u] >> v & 1)
        H[mask] = sum(dp[mask])
    return H, dp


def response(rows):
    H, _ = inventory(rows)
    full = len(H)-1
    return H[-1], sum(H[m]*H[full ^ m] for m in range(full+1))


def third_bin(rows):
    H, _ = inventory(rows)
    full, total = len(H)-1, 0
    for first in range(full+1):
        rest, second = full ^ first, full ^ first
        while True:
            total += H[first]*H[second]*H[rest ^ second]
            if second == 0:
                break
            second = (second-1) & rest
    return total


def run_table(frame, module, maximum):
    """Exact exterior words with repeated, nonadjacent module placeholders."""
    ext = [v for v in range(len(frame)) if v != module]
    n, symbol = len(ext), len(ext)
    end = [[[0]*(n+1) for _ in range(maximum+1)] for _ in range(1 << n)]
    total = [[0]*(maximum+1) for _ in range(1 << n)]
    total[0][0] = 1
    for mask in range(1 << n):
        for k in range(maximum+1):
            for i, vertex in enumerate(ext):
                if mask >> i & 1:
                    previous = mask ^ (1 << i)
                    value = int(previous == 0 and k == 0)
                    value += sum(end[previous][k][j] for j, other in enumerate(ext)
                                 if frame[other] >> vertex & 1)
                    if frame[module] >> vertex & 1:
                        value += end[previous][k][symbol]
                    end[mask][k][i] = value
            if k:
                end[mask][k][symbol] = int(mask == 0 and k == 1) + sum(
                    end[mask][k-1][i] for i, vertex in enumerate(ext)
                    if frame[vertex] >> module & 1)
            if mask or k:
                total[mask][k] = sum(end[mask][k])
    return ext, total, end


def universal_u_certificate():
    frame = pi_frame()
    ext, total, end = run_table(frame, 1, 3)
    full, u = (1 << len(ext))-1, ext.index(3)
    opposite = tuple(((1 << 9)-1) ^ (1 << i) ^ row for i, row in enumerate(frame))
    _, _, opposite_end = run_table(opposite, 1, 3)
    coefficients = dict(H=[], s=[], q=[], F=[])
    for runs in range(1, 4):
        coefficients['H'].append(total[full][runs])
        coefficients['s'].append(opposite_end[full][runs][u])
        coefficients['q'].append(total[full ^ (1 << u)][runs])
        coefficients['F'].append(sum(end[mask][first_runs][u]*total[full ^ mask][runs-first_runs]
                                      for mask in range(full+1) if mask >> u & 1
                                      for first_runs in range(runs+1)))
    need(coefficients == dict(H=[361,180,0], s=[57,0,0], q=[121,60,0], F=[601,606,186]),
         'complete universal clock-root run certificate')
    return coefficients


def state(rows, root):
    H, E = inventory(rows)
    full = (1 << len(rows))-1
    opposite = tuple(full ^ (1 << i) ^ row for i, row in enumerate(rows))
    _, opposite_E = inventory(opposite)
    return (H[-1], opposite_E[full][root], H[full ^ (1 << root)],
            sum(E[m][root]*H[full ^ m] for m in range(1, full+1)))


def source(K):
    return substitute(S, (ordinal(tt(2), K),) + ((0,),)*5)


def target(K, j=1):
    blocks = [(0,)]*9
    blocks[1], blocks[3] = tuple(K), tt(j)
    return substitute(pi_frame(), blocks)


def handoff_certificate(K):
    old, old_slots = source(K)
    new, new_slots = target(K)
    prefix_a, prefix_b = old_slots[0][:2]
    old_K = old_slots[0][2:]
    fixed = {i: old_slots[i][0] for i in range(1, 6)}
    # Target positions a,c,u,d0,d1 are retained old names.
    mapping = {prefix_a: new_slots[0][0], fixed[1]: new_slots[2][0],
               fixed[4]: new_slots[3][0], fixed[3]: new_slots[4][0],
               fixed[2]: new_slots[5][0]}
    mapping.update(zip(old_K, new_slots[1]))
    deleted = (prefix_b, fixed[5])
    inserted = tuple(new_slots[i][0] for i in (6, 7, 8))
    flips = {frozenset((fixed[2], fixed[3])), frozenset((fixed[1], fixed[3])),
             frozenset((fixed[4], prefix_a))}
    flips |= {frozenset((fixed[i], v)) for i in (2, 3, 4) for v in old_K}
    observed = set()
    for x in mapping:
        for y in mapping:
            if x < y and ((old[x] >> y) & 1) != ((new[mapping[x]] >> mapping[y]) & 1):
                observed.add(frozenset((x, y)))
    need(observed == flips, "exact retained-old arc changes")
    need(len(new) == len(old)+1, "net one vertex")
    need(set(range(len(old))) - set(mapping) == set(deleted), "exact deleted vertices")
    need(set(range(len(new))) - set(mapping.values()) == set(inserted), "exact inserted vertices")
    h, p = response(K)
    source_h, target_h = inventory(old)[0][-1], inventory(new)[0][-1]
    need(source_h == 60*p-h and target_h == 180*p+h == 3*source_h+4*h,
         "universal scalar handoff")
    source_state = state(old, fixed[4])
    target_state = state(new, new_slots[3][0])
    w = third_bin(K)
    need(source_state == (60*p-h,12*p+2*h,12*p+3*h,60*p), 'source ready profile')
    need(target_state == (180*p+h,57*h,60*p+h,186*w+48*p-53*h), 'target clock profile')
    return dict(module=K, h=h, p=p, source_rows=old, target_rows=new,
                source_h=source_h, target_h=target_h, old_to_new=mapping,
                deleted=deleted, inserted=inserted,
                retained_arc_flips=sorted(sorted(pair) for pair in flips),
                clock_u=new_slots[3][0], clock_v=new_slots[1],
                third_bin=w, source_state=source_state, target_state=target_state)


def clock_step(R, k, j):
    """R may be empty, but both before/after dyadic modules must be nonempty."""
    need(k >= 1 and (R or k >= 2), "nonempty dyadic module after transfer")
    before, old_slots = target(ordinal(R, tt(k)), j)
    after, new_slots = target(ordinal(R, tt(k-1)), j+1)
    moved = old_slots[1][-1]
    mapping = {}
    for block in range(9):
        if block == 1:
            mapping.update(zip(old_slots[1][:-1], new_slots[1]))
        elif block == 3:
            mapping[moved] = new_slots[3][0]
            mapping.update(zip(old_slots[3], new_slots[3][1:]))
        else:
            mapping.update(zip(old_slots[block], new_slots[block]))
    changes = []
    for x in range(len(before)):
        for y in range(x+1, len(before)):
            if ((before[x] >> y) & 1) != ((after[mapping[x]] >> mapping[y]) & 1):
                changes.append((x, y))
    need(changes == [tuple(sorted((moved, old_slots[2][0])))], "one physical clock arc")
    old_h, new_h = inventory(before)[0][-1], inventory(after)[0][-1]
    anchor = response(R)[0] if R else 1
    need(2*new_h == 3*old_h-anchor, "anchored transfer law")
    return dict(R=R, k=k, j=j, before_h=old_h, after_h=new_h, anchor=anchor,
                changed_arc=changes[0], old_to_new=mapping)


def select_target(Y):
    """Exact finite inverse for all Π[TT_j,K], j>=1; no runtime claim."""
    need(Y > 0 and Y % 2 == 1, 'positive odd target')
    rejected, accepted = [], []
    j = 1
    while 120*3**j+1 <= Y:
        amplitude = 60*3**j
        maximum_h = Y//(2*amplitude+1)
        for h in range(Y % (2*amplitude), maximum_h+1, 2*amplitude):
            if h <= 0:
                continue
            p = (Y-h)//amplitude
            nmax = min(p//h-1, p.bit_length()-1)
            witness = None
            if h == 1:
                # The inherited H=1 iff transitive theorem makes this exact.
                if p >= 2 and p & (p-1) == 0:
                    witness = tt(p.bit_length()-1)
            else:
                for n in range(1, nmax+1):
                    pairs = tuple(combinations(range(n), 2))
                    for mask in range(1 << len(pairs)):
                        rows = [0]*n
                        for bit, (a,b) in enumerate(pairs):
                            if mask >> bit & 1:
                                rows[a] |= 1 << b
                            else:
                                rows[b] |= 1 << a
                        if response(rows) == (h,p):
                            witness = tuple(rows)
                            break
                    if witness is not None:
                        break
            item = dict(j=j, h=h, p=p, maximum_module_order=nmax)
            if witness is None:
                rejected.append(item)
            else:
                rows, slots = target(witness, j)
                actual = inventory(rows)[0][-1]
                need(actual == Y, 'inverse-selected literal target')
                item.update(module=witness, source_h=60*p-h,
                            first_handoff_h=180*p+h, target_rows=rows,
                            clock_u=slots[3], clock_v=slots[1])
                accepted.append(item)
        j += 1
    return dict(target=Y, accepted=accepted, rejected_physical_candidates=rejected)


def main():
    certificate = universal_u_certificate()
    modules = [tt(1), tt(2), tt(3), (2, 4, 1), (2, 4, 9, 3)]
    controls = [handoff_certificate(K) for K in modules]
    clocks = [clock_step((), 2, 1), clock_step((), 3, 1),
              clock_step((2, 4, 1), 1, 1)]
    # The active clock root has a separate interface from readiness.
    u_states = []
    for m in range(1, 5):
        graph, slots = target(tt(m))
        u_states.append(dict(m=m, state=state(graph, slots[3][0])))
    output = dict(inherited_ref=INHERITED_REF,
                  scope="Five named module handoffs; three retained-vertex clock transfers; four active-u state observations",
                  pi_frame=pi_frame(), universal_u_coefficients=certificate,
                  controls=controls, clock_controls=clocks,
                  u_states=u_states,
                  exact_inverse_controls=[select_target(Y) for Y in (361,721,1081)])
    need(output['exact_inverse_controls'][2]['rejected_physical_candidates'] ==
         [dict(j=1,h=1,p=6,maximum_module_order=2)], 'arithmetic versus physical inverse hostile')
    Path(__file__).with_name("receipt.json").write_text(json.dumps(output, indent=2)+"\n")
    print("PASS: five universal-move controls and three one-arc clock transfers")


if __name__ == "__main__":
    main()
