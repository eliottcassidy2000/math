"""Exact finite-program compiler into a marked Hamiltonian tournament."""
from itertools import product, permutations
import json


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def encode_bits(machine):
    s = len(machine)//2
    need(s >= 1 and len(machine) == 2*s, 'positive binary-state machine')
    width = s.bit_length()
    bits = [1]*s + [0]
    for write, direction, target in machine:
        need(write in (0,1) and direction in (-1,1) and 0 <= target <= s,
             'transition field domain')
        bits.extend((write, int(direction == 1)))
        bits.extend((target >> j)&1 for j in range(width-1,-1,-1))
    return bits


def encode(machine):
    bits = encode_bits(machine)
    n = 2
    while (n-1)*(n-2)//2 < len(bits):
        n += 1
    rows = [0]*n
    at = 0
    for u in range(n):
        for v in range(u+1,n):
            if v == u+1:
                bit = 1
            else:
                bit = bits[at] if at < len(bits) else 0
                at += 1
            if bit:
                rows[u] |= 1 << v
            else:
                rows[v] |= 1 << u
    return tuple(rows), tuple(range(n))


def arc(rows,u,v):
    return (rows[u] >> v)&1


def decode(rows, path):
    n = len(rows)
    need(sorted(path) == list(range(n)), 'marked vertex order')
    need(all(not arc(rows,u,u) for u in range(n)), 'no loops')
    need(all(arc(rows,u,v)+arc(rows,v,u)==1 for u in range(n) for v in range(u+1,n)),
         'one direction per unordered pair')
    need(all(arc(rows,u,v) for u,v in zip(path,path[1:])), 'directed marked spine')
    bits = [arc(rows,path[i],path[j]) for i in range(n) for j in range(i+2,n)]
    need(0 in bits, 'terminated unary header')
    s = bits.index(0)
    need(s >= 1, 'state count is positive')
    width = s.bit_length()
    length = s+1+2*s*(2+width)
    need(length <= len(bits), 'complete finite table')
    need(all(b==0 for b in bits[length:]), 'canonical zero padding')
    machine = []
    at = s+1
    for _ in range(2*s):
        write, right = bits[at:at+2]
        at += 2
        target = 0
        for b in bits[at:at+width]:
            target = 2*target+b
        at += width
        need(target <= s, 'valid target code')
        machine.append((write, 1 if right else -1, target))
    return tuple(machine)


def relabel(rows,path,p):
    out = [0]*len(rows)
    for u in range(len(rows)):
        for v in range(len(rows)):
            if arc(rows,u,v):
                out[p[u]] |= 1 << p[v]
    return tuple(out), tuple(p[u] for u in path)


def k_move(rows,path):
    edges = {frozenset((u,v)) for u,v in zip(path,path[1:])}
    out = list(rows)
    for u in range(len(rows)):
        for v in range(u+1,len(rows)):
            if frozenset((u,v)) not in edges:
                out[u] ^= 1 << v
                out[v] ^= 1 << u
    return tuple(out)


def trajectory(machine, bound=100):
    s = len(machine)//2
    state = head = 0
    ones = set()
    seen = {}
    for time in range(bound+1):
        if state == s:
            return {'kind':'halt','steps':time}
        key = state,head,tuple(sorted(ones))
        if key in seen:
            return {'kind':'repeat','first':seen[key],'steps':time,'period':time-seen[key]}
        seen[key] = time
        if time == bound:
            return {'kind':'bounded','steps':time}
        write,direction,state = machine[2*state+int(head in ones)]
        if write:
            ones.add(head)
        else:
            ones.discard(head)
        head += direction


def pad_states(machine, new_count):
    """Preserve every old execution, renaming only the halt target."""
    old_count = len(machine)//2
    need(new_count >= old_count >= 1, 'padding state count')
    out = [(w,d,new_count if t == old_count else t) for w,d,t in machine]
    for state in range(old_count,new_count):
        out.extend(((0,1,state),(1,1,state)))
    return tuple(out)


def prefix(machine, steps, initial=(0,0,())):
    """Literal configurations with a common terminal marker for comparison."""
    state,head,tape = initial
    ones = set(tape)
    s = len(machine)//2
    out = []
    for _ in range(steps+1):
        out.append(('H' if state == s else state,head,tuple(sorted(ones))))
        if state == s:
            break
        write,direction,state = machine[2*state+int(head in ones)]
        if write:
            ones.add(head)
        else:
            ones.discard(head)
        head += direction
    return out


def main():
    cases = relabel_cases = rejected = 0
    for s in (1,2):
        options = tuple(product((0,1),(-1,1),range(s+1)))
        for machine in product(options,repeat=2*s):
            rows,path = encode(machine)
            need(decode(rows,path)==machine, 'exact syntax round trip')
            need(len(rows)==2 or (len(rows)-2)*(len(rows)-3)//2 < len(encode_bits(machine)),
                 'least vertex order for this serialization capacity')
            try:
                decode(k_move(rows,path),path)
            except RuntimeError as err:
                need(str(err)=='state count is positive','K boundary caused by initial header bit')
                rejected += 1
            else:
                raise RuntimeError('canonical K image must leave positive-header code language')
            cases += 1
    loop = ((0,1,1),(0,1,2),(0,-1,0),(0,-1,0))
    halt = ((1,1,1),(0,1,2),(0,-1,0),(0,-1,0))
    a,p = encode(loop)
    b,_ = encode(halt)
    changed = [(u,v) for u in range(len(a)) for v in range(u+1,len(a))
               if arc(a,u,v)!=arc(b,u,v)]
    need(len(changed)==1,'one exact program bit gives one tournament arc')
    need([(d,t) for _,d,t in loop]==[(d,t) for _,d,t in halt],
         'same read-labelled state/movement control graph')
    need(trajectory(loop)=={'kind':'repeat','first':0,'steps':2,'period':2},'literal loop proof')
    need(trajectory(halt)=={'kind':'halt','steps':3},'literal halt proof')
    for perm in permutations(range(len(a))):
        ar,pr = relabel(a,p,perm)
        need(decode(ar,pr)==loop, 'joint vertex/frame equivariance')
        relabel_cases += 1
    padding_cases = 0
    one_options = tuple(product((0,1),(-1,1),range(2)))
    for machine in product(one_options,repeat=2):
        padded = pad_states(machine,25)
        need(decode(*encode(padded))==padded,'padded syntax round trip')
        for tape in ((),(0,),(-1,1),(-2,0,3)):
            need(prefix(machine,80,(0,0,tape))==prefix(padded,80,(0,0,tape)),
                 'padding preserves literal tape/head/native clock')
            padding_cases += 1
    for machine in (loop,halt):
        for size in (5,15,25):
            padded = pad_states(machine,size)
            need(trajectory(padded)==trajectory(machine),'padding preserves hostile outcomes')
            need(prefix(machine,20)==prefix(padded,20),'padding preserves hostile native trace')
            padding_cases += 1
    print(json.dumps({'status':'PASS','complete_transition_tables':cases,
        'canonical_K_images_outside_code':rejected,'joint_frame_relabellings':relabel_cases,
        'padding_trace_controls':padding_cases,
        'hostile':{'order':len(a),'changed_pair':changed,'loop':trajectory(loop),'halt':trajectory(halt)},
        'scope':'finite program syntax with marked Hamiltonian frame; no halting equivalence for unframed graph moves'},
        indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
