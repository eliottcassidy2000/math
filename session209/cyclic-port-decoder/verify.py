#!/usr/bin/env python3
"""Exact, standard-library checks. Run with python3 verify.py or python3 -O verify.py."""
from itertools import permutations
from collections import Counter
import json
from pathlib import Path


def require(value, message):
    if not value:
        raise RuntimeError(message)


def cyclic(n):
    return [sum(1 << ((i+j) % n) for j in range(1, (n+1)//2)) for i in range(n)]


def converse(a):
    return [((1 << len(a))-1) ^ (1 << i) ^ row for i, row in enumerate(a)]


def join(blocks):
    n = sum(map(len, blocks)); out = []; off = 0
    for block in blocks:
        later = ((1 << n)-1) ^ ((1 << (off+len(block)))-1)
        out.extend((row << off) | later for row in block)
        off += len(block)
    return out


def extension(a, masks):
    k = len(a); out = list(a)+[0]*len(masks)
    for i,b in enumerate(masks):
        for u in range(k):
            if b >> u & 1:
                out[u] |= 1 << (k+i)
            else:
                out[k+i] |= 1 << u
        for j in range(i+1,len(masks)):
            out[k+i] |= 1 << (k+j)
    return out


def hp(a):
    n = len(a)
    if not n:
        return 1
    d = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        d[1 << v][v] = 1
    for s in range(1,1 << n):
        for v in range(n):
            if not (s >> v & 1):
                continue
            p = s ^ (1 << v)
            if p:
                d[s][v] = sum(d[p][u] for u in range(n) if p >> u & 1 and a[u] >> v & 1)
    return sum(d[-1])


def hp_permutations(a):
    return sum(all(a[p[i]] >> p[i+1] & 1 for i in range(len(a)-1)) for p in permutations(range(len(a))))


def is_isomorphism(a,b,p):
    return all(((a[i] >> j)&1) == ((b[p[i]] >> p[j])&1) for i in range(len(a)) for j in range(len(a)))


def reflection(n,mask):
    for c in range(n):
        p = [(c-i) % n for i in range(n)]
        if all(((mask >> i)&1)==((mask >> p[i])&1) for i in range(n)):
            return p
    return None


def quotient_module(n,r,code):
    """Substitute C_n for vertex0 in an arbitrary labelled (r+1)-tournament."""
    out = cyclic(n)+[0]*r; bit=0
    blocks = [list(range(n))]+[[n+i] for i in range(r)]
    for i in range(r+1):
        for j in range(i+1,r+1):
            first,second=(blocks[i],blocks[j]) if code >> bit & 1 else (blocks[j],blocks[i])
            for u in first:
                for v in second:
                    out[u] |= 1 << v
            bit += 1
    return out


def flip_first(a,n):
    out=list(a)
    for i in range(n):
        out[i] ^= ((1 << n)-1) ^ (1 << i)
    return out


def run():
    result = {"status":"PASS", "module_controls":{}, "local_two_port":{}, "chain_controls":[]}
    for n in (3,5):
        a=cyclic(n); b=converse(a)
        require(hp(a)==hp_permutations(a),"local count paths disagree")
        checks=0
        # Every binary color of the module has a color-preserving reflection.
        for mask in range(1 << n):
            p=reflection(n,mask)
            require(p is not None,"binary color has no reflection")
            require(is_isomorphism(extension(a,[mask]),extension(b,[mask]),p+[n]),"local reflection fails")
            require(hp(extension(a,[mask]))==hp(extension(b,[mask])),"one-port local count mismatch")
        # Exhaust all quotient orientations through four exterior vertices.
        # Exterior incidence of the added port is arbitrary: independently
        # test all its masks, not merely uniform ones.
        for r in range(5):
            for code in range(1 << ((r+1)*r//2)):
                a=quotient_module(n,r,code); b=flip_first(a,n)
                for mask in range(1 << (n+r)):
                    p=reflection(n,mask & ((1 << n)-1))+list(range(n,n+r+1))
                    require(is_isomorphism(extension(a,[mask]),extension(b,[mask]),p),"embedded module reflection fails")
                    checks+=1
        result["module_controls"][str(n)]={"max_exterior_order":4,"isomorphism_checks":checks}
        a=extension(cyclic(n),[1,2]); b=extension(converse(cyclic(n)),[1,2])
        x,y=hp(a),hp(b)
        require((x,y)==((11,13) if n==3 else (89,91)),"wrong local contrast")
        require((x,y)==(hp_permutations(a),hp_permutations(b)),"local permutation contrast mismatch")
        result["local_two_port"][str(n)]={"masks":[1,2],"counts":[x,y],"mod4":[x%4,y%4],"permutation_counts":True}
    # Smallest nonmodule hostile: the second port can expose a triangle reversal.
    a=extension(cyclic(3),[1]); b=extension(converse(cyclic(3)),[1])
    x,y=hp(extension(a,[10])),hp(extension(b,[10]))
    require((x,y)==(11,13),"nonmodule hostile fails")
    result["nonmodule_hostile"]={"core_order":4,"core_counts":[hp(a),hp(b)],"port_mask":10,"responses":[x,y]}
    # Non-palindromic rigid word; palindromic S4 mixed leaf; nonleaf C5 boundary.
    for sizes in ([1,3,1,3,3],[3,1,4,1,3],[1,5,3]):
        blocks=[cyclic(t) if t in (3,5) else ([2,4,9,3] if t==4 else [0]) for t in sizes]
        slots=[j for j,t in enumerate(sizes) if t in (3,5)]
        base_h=1
        for t in sizes:
            base_h *= {1:1,3:3,4:5,5:15}[t]
        recovered=[]
        # All 2^s fixed-label chirality variants, s local binary observations each.
        for variant in range(1 << len(slots)):
            curr=[list(b) for b in blocks]
            for t,j in enumerate(slots):
                if variant >> t & 1:
                    curr[j]=converse(curr[j])
            core=join(curr)
            require(hp(core)==base_h,"SCC product mismatch")
            observations=[]
            for t,j in enumerate(slots):
                off=sum(sizes[:j]); left=(1 << off)-1
                masks=[left | (1 << off),left | (1 << (off+1))]
                observed=hp(extension(core,masks))
                q=base_h//({3:3,5:15}[sizes[j]])
                local=({3:11,5:89}[sizes[j]])+2*((variant >> t)&1)
                require(observed==q*local,"localized two-port factor mismatch")
                observations.append(((observed-q*({3:11,5:89}[sizes[j]]))//2)%2)
            require(sum(bit << t for t,bit in enumerate(observations))==variant,"mod4 decoder fails")
            recovered.append(observations)
        result["chain_controls"].append({"block_orders":list(sizes),"raw_H":base_h,"variants":1 << len(slots),"bits_recovered":len(slots),"all_localized_DPs_pass":True})
    # Exhaust full one-port response tables for two triangle SCCs, independently
    # of the explicit isomorphism checks above.
    tables=[]
    for bits in range(4):
        blocks=[converse(cyclic(3)) if bits >> i & 1 else cyclic(3) for i in range(2)]
        core=join(blocks)
        tables.append([hp(extension(core,[mask])) for mask in range(64)])
    require(all(t==tables[0] for t in tables),"two-triangle one-port fibre fails")
    result["two_triangle_full_one_port_tables"]={"variants":4,"masks_each":64,"all_equal":True,"histogram":dict(sorted(Counter(tables[0]).items()))}
    # The mixed-leaf center has no one-port ambiguity: exhaust every labelled
    # strong four-vertex core and its full sixteen-mask response table.
    strong4=[]
    for code in range(64):
        a=[0]*4; bit=0
        for i in range(4):
            for j in range(i+1,4):
                u,v=(i,j) if code >> bit & 1 else (j,i)
                a[u] |= 1 << v; bit+=1
        strong=True
        for start in range(4):
            seen=1 << start
            while True:
                new=seen
                for u in range(4):
                    if seen >> u & 1:
                        new |= a[u]
                if new==seen:
                    break
                seen=new
            if seen!=15:
                strong=False
                break
        if strong:
            require(hp(a)==5,"S4 count differs")
            strong4.append(tuple(hp(extension(a,[mask])) for mask in range(16)))
    require(len(strong4)==24 and len(set(strong4))==24,"S4 one-port collision")
    result["mixed_center_one_port_faithfulness"]={"labelled_strong_S4":24,"masks_each":16,"distinct_full_tables":24,"complete_order4_sweep":64}
    # Reflection argument cannot automatically extend to C7: a binary coloring
    # need not have a reflection, although its naked one-port counts may agree.
    missing=[b for b in range(128) if reflection(7,b) is None]
    require(bool(missing),"C7 hostile unexpectedly absent")
    result["color_symmetry_boundary"]={"C7_nonreflective_mask":missing[0],"all_nonreflective_masks":len(missing),"no_count_failure_inferred":True}
    return result

if __name__ == "__main__":
    data=run()
    text=json.dumps(data,indent=2,sort_keys=True)+"\n"
    Path(__file__).with_name("RESULT.json").write_text(text)
    print(text,end="")
