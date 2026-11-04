"""Exhaustive original-module XOR controls through five named vertices.

Modules and strong trees are built by definition. Every induced-subset
ordered-cover table is exact. Legal normal forms are replayed with a fresh
module check before each step. No removable assertions or third-party code.
"""
import hashlib
import itertools
import json


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def rows_from_code(n, code):
    rows = [0]*n
    for bit,(i,j) in enumerate(itertools.combinations(range(n),2)):
        u,v = (j,i) if code >> bit & 1 else (i,j)
        rows[u] |= 1 << v
    return tuple(rows)


def module(rows, mask):
    return all((row & mask) in (0,mask) for i,row in enumerate(rows) if not mask >> i & 1)


def all_modules(rows):
    return [mask for mask in range(1,1 << len(rows)) if module(rows,mask)]


def transitive(rows):
    return sorted(row.bit_count() for row in rows) == list(range(len(rows)))


def opposite(rows):
    full = (1 << len(rows))-1
    return tuple(full ^ row ^ (1 << v) for v,row in enumerate(rows))


def quotient(rows, children):
    representatives = [(mask & -mask).bit_length()-1 for mask in children]
    return tuple(sum(1 << j for j,v in enumerate(representatives) if rows[u] >> v & 1)
                 for u in representatives)


def tree(rows, modules):
    strong = [M for M in modules
              if all((M & X) in (0,M,X) for X in modules)]
    nodes = {}
    for M in strong:
        if M.bit_count() < 2:
            continue
        proper = [X for X in strong if X != M and X & M == X]
        children = tuple(X for X in proper
                         if not any(X != Y and X & Y == X for Y in proper))
        need(sum(children) == M and all(not (X & Y) for X,Y in itertools.combinations(children,2)),
             "strong children do not partition")
        Q = quotient(rows,children)
        linear = transitive(Q)
        if not linear:
            need(all(X.bit_count() in (1,len(Q)) for X in all_modules(Q)), "Gallai quotient not prime")
        nodes[M] = {"children":children,"quotient":Q,"linear":linear}
    return strong,nodes


def path_cover_table(rows):
    n = len(rows)
    endpoint = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        endpoint[1 << v][v] = 1
    for mask in range(1,1 << n):
        for v,value in enumerate(endpoint[mask]):
            unused = rows[v] & ~mask
            while value and unused:
                bit = unused & -unused
                unused -= bit
                endpoint[mask | bit][bit.bit_length()-1] += value
    H = [1]+[sum(row) for row in endpoint[1:]]
    covers = [[0]*(n+1) for _ in range(1 << n)]
    covers[0][0] = 1
    for mask in range(1,1 << n):
        part = mask
        while part:
            rest = mask ^ part
            for k in range(1,n+1):
                covers[mask][k] += H[part]*covers[rest][k-1]
            part = (part-1)&mask
    return H,tuple(value for row in covers for value in row)


def cycle_mask(rows):
    result = 0
    for bit,triple in enumerate(itertools.combinations(range(len(rows)),3)):
        mask = sum(1 << v for v in triple)
        if all((rows[v] & mask).bit_count() == 1 for v in triple):
            result |= 1 << bit
    return result


def internal_edges(n, mask):
    return sum(1 << bit for bit,(u,v) in enumerate(itertools.combinations(range(n),2))
               if mask >> u & 1 and mask >> v & 1)


def legal_normal_form(n, start, target, nodes):
    current = start
    trace = []
    target_rows = rows_from_code(n,target)

    def flip(mask):
        nonlocal current
        if mask.bit_count() < 2:
            return
        need(module(rows_from_code(n,current),mask), ("illegal intermediate module", n,current,mask))
        current ^= internal_edges(n,mask)
        trace.append(mask)

    def visit(M):
        if M not in nodes:
            return
        node = nodes[M]
        children = node["children"]
        current_q = quotient(rows_from_code(n,current),children)
        target_q = quotient(target_rows,children)
        if not node["linear"]:
            need(target_q in (current_q,opposite(current_q)), "prime target orientation")
            if target_q != current_q:
                flip(M)
                for child in children:
                    flip(child)
        else:
            need(transitive(current_q) and transitive(target_q), "linear normal form target")
            order_now = sorted(range(len(children)),key=lambda i:-current_q[i].bit_count())
            order_target = sorted(range(len(children)),key=lambda i:-target_q[i].bit_count())
            for position,wanted in enumerate(order_target):
                at = order_now.index(wanted)
                while at > position:
                    left,right = order_now[at-1],order_now[at]
                    flip(children[left] | children[right])
                    flip(children[left])
                    flip(children[right])
                    order_now[at-1],order_now[at] = right,left
                    at -= 1
        need(quotient(rows_from_code(n,current),children) == target_q, "quotient not arranged")
        for child in children:
            visit(child)

    visit((1 << n)-1)
    need(current == target, ("normal form final mismatch",n,start,target,current))
    return trace


def main():
    census = []
    totals = {"families":0,"scalar_equal_families":0,"scalar_strict_families":0,
              "induced_cover_coordinate_comparisons":0,"strong_node_actions":0,
              "legal_normal_forms":0,"legal_reversals":0,"max_legal_reversals":0,
              "descendant_module_restrictions":0}
    examples = []
    digest = hashlib.sha256()
    for n in range(1,6):
        tournament_count = 1 << (n*(n-1)//2)
        cache = []
        for code in range(tournament_count):
            rows = rows_from_code(n,code)
            modules = all_modules(rows)
            H,q = path_cover_table(rows)
            strong,nodes = tree(rows,modules)
            cache.append((rows,set(modules),H,q,cycle_mask(rows),strong,nodes))
        family_count,equal_count,strict_count = 0,0,0
        for start,(rows,modules,H,q,cycles,strong,nodes) in enumerate(cache):
            selected_modules = sorted(M for M in modules if M.bit_count() >= 2)
            edge_masks = [internal_edges(n,M) for M in selected_modules]
            # Restriction of every original operation to every strong node.
            for M in strong:
                vertices = [v for v in range(n) if M >> v & 1]
                local_rows = tuple(sum(1 << j for j,v in enumerate(vertices) if rows[u] >> v & 1)
                                   for u in vertices)
                for X in selected_modules:
                    intersection = sum(1 << j for j,v in enumerate(vertices) if X >> v & 1)
                    if intersection:
                        need(module(local_rows,intersection), "invalid descendant operation")
                        totals["descendant_module_restrictions"] += 1
            current,previous_gray = start,0
            checked_legal = set()
            for family in range(1 << len(selected_modules)):
                gray = family ^ (family >> 1)
                if family:
                    changed = gray ^ previous_gray
                    current ^= edge_masks[changed.bit_length()-1]
                previous_gray = gray
                ur,um,uH,uq,uc,_,_ = cache[current]
                need(all(M in um for M in strong), "original strong module lost")
                need(cycles & uc == cycles, "original directed triangle destroyed")
                need(all(a <= b for a,b in zip(q,uq)), ("hereditary cover decrease",n,start,current))
                totals["induced_cover_coordinate_comparisons"] += len(q)
                bad = False
                for M,node in nodes.items():
                    K = quotient(ur,node["children"])
                    if node["linear"]:
                        bad |= not transitive(K)
                    else:
                        need(K in (node["quotient"],opposite(node["quotient"])), "prime quotient changed partially")
                    totals["strong_node_actions"] += 1
                same = H[-1] == uH[-1]
                need(same == (not bad) == (cycles == uc) == (H == uH) == (q == uq),
                     ("equality criterion",n,start,current))
                need(same or H[-1] < uH[-1], "scalar decrease")
                family_count += 1
                equal_count += same
                strict_count += not same
                if same and current not in checked_legal:
                    checked_legal.add(current)
                    trace = legal_normal_form(n,start,current,nodes)
                    totals["legal_normal_forms"] += 1
                    totals["legal_reversals"] += len(trace)
                    totals["max_legal_reversals"] = max(totals["max_legal_reversals"],len(trace))
                    digest.update(json.dumps([n,start,current,trace],separators=(',',':')).encode())
                    if current != start and current != (tournament_count-1-start) and len(examples) < n:
                        examples.append({"n":n,"start_rows":rows,"target_rows":ur,
                                         "legal_module_masks":trace,"H":H[-1]})
        census.append({"n":n,"tournaments":tournament_count,"original_module_families":family_count,
                       "scalar_equal":equal_count,"scalar_strict":strict_count})
        totals["families"] += family_count
        totals["scalar_equal_families"] += equal_count
        totals["scalar_strict_families"] += strict_count
    return {"status":"PASS","census":census,"totals":totals,
            "legal_normal_form_sha256":digest.hexdigest(),"legal_examples":examples,
            "scope":"Original-module XOR monotonicity and scalar/C3/legal-class reduction. No full named-response collision or injectivity theorem is asserted."}


if __name__ == '__main__':
    print(json.dumps(main(),sort_keys=True,indent=2))
