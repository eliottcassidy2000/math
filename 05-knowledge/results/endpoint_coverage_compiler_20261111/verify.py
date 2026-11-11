"""Independent physical-graph recount, with no endpoint-cover DP imported.

Recover W from four ordinary one-vertex H counts and U from endpoint counts.
Invert the triangular P_j(a) table to certify every coefficient of each named
infinite family, rather than trusting a short initial list of matching sizes.
"""
from pathlib import Path
from itertools import combinations
import json
import construct as con

HERE = Path(__file__).resolve().parent


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def count(rows):
    n = len(rows)
    for i in range(n):
        require(not (rows[i] >> i & 1), 'loop')
        for j in range(i):
            require(bool(rows[i]>>j&1) != bool(rows[j]>>i&1), 'not a tournament')
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1,1 << n):
        for v in range(n):
            if dp[mask][v]:
                choices = rows[v] & ~mask
                while choices:
                    bit = choices & -choices
                    choices -= bit
                    dp[mask|bit][bit.bit_length()-1] += dp[mask][v]
    return sum(dp[-1]), dp[-1]


def recursive_count(rows):
    n = len(rows)
    def walk(last,unused):
        if not unused:
            return 1
        choices = unused & rows[last]
        total = 0
        while choices:
            bit = choices & -choices
            choices -= bit
            total += walk(bit.bit_length()-1,unused^bit)
        return total
    return sum(walk(v,((1<<n)-1)^(1<<v)) for v in range(n))


def ear(rows, selected):
    n = len(rows)
    out = list(rows)+[sum(1<<v for v in selected)]
    for v in range(n):
        if v not in selected:
            out[v] |= 1 << n
    return out


def endpoint_state(rows, selected):
    h, end = count(rows)
    singles = {v:count(ear(rows,[v]))[0] for v in selected}
    w = sum(singles[a]+singles[b]-count(ear(rows,[a,b]))[0]-h
            for a,b in combinations(selected,2))
    u = sum(singles[v]-h+end[v] for v in selected)
    return h,w,u


def stirling_weights(a, limit):
    # Standard Stirling recurrence, separate from ordered-partition recurrence.
    s = [[0]*(limit+1) for _ in range(a+1)]
    s[0][0] = 1
    for i in range(1,a+1):
        for j in range(1,limit+1):
            s[i][j] = s[i-1][j-1]+j*s[i-1][j]
    factorial, out = 1,[]
    for j in range(1,limit+1):
        factorial *= j
        out.append(factorial*s[a][j])
    return out


def invert(values):
    out = []
    for a,value in enumerate(values,1):
        weights = stirling_weights(a,a)
        numerator = value-sum(x*y for x,y in zip(out,weights))
        require(numerator%weights[-1] == 0,'integer cover coefficient')
        out.append(numerator//weights[-1])
    return out


def substitute(rows, marked, module):
    labels = [v for v in range(len(rows)) if v != marked]
    out = [0]*(len(labels)+len(module))
    for i,j in combinations(range(len(out)),2):
        if j < len(labels):
            forward = bool(rows[labels[i]] >> labels[j] & 1)
        elif i < len(labels):
            forward = bool(rows[labels[i]] >> marked & 1)
        else:
            forward = bool(module[i-len(labels)] >> (j-len(labels)) & 1)
        out[i if forward else j] |= 1 << (j if forward else i)
    return out,labels+[marked]*len(module)


def main():
    bank = json.loads((HERE / 'inflation_bank_7.json').read_text())['records']
    named, tails, recursive, nontransitive = [],0,0,0
    for family,key in con.FAMILIES.items():
        item = next(r for r in bank if (r['n'],r['code'],r['selected'],r['inflated_vertex']) == key)
        defects, slopes = [],[]
        for size in range(1,item['n']+1):
            rows,selected,metadata = con.seed(family,size)
            h,w,u = endpoint_state(rows,selected)
            defects.append(h-w)
            slopes.append(u)
            require(h-w == 1 and u == metadata['U'],'actual seed state')
        recovered_d, recovered_u = invert(defects),invert(slopes)
        require(recovered_d == item['defect_coefficients'],'all defect coefficients')
        expected_u = item['slope_cover_coefficients']+[0]*(item['n']-len(item['slope_cover_coefficients']))
        require(recovered_u == expected_u,'all slope coefficients')
        for size in range(1,4):
            for length in (1,2):
                out = con.construct(family,size,length)
                h,_ = count(out['rows'])
                require(h == out['expected_H'],'actual tail')
                tails += 1
                if out['order'] <= 9:
                    require(recursive_count(out['rows']) == h,'independent path backtracking')
                    recursive += 1
        expanded, labels = substitute(item['rows'],item['inflated_vertex'],[2,4,1])
        selected = [i for i,v in enumerate(labels) if item['selected'] >> v & 1]
        h,w,u = endpoint_state(expanded,selected)
        require(h-w == 3,'C3 substitution must change the offset to 3')
        nontransitive += 1
        named.append(dict(family=family,base=key,defect_coefficients=recovered_d,slope_coefficients=recovered_u))
    ray = []
    for k in range(1,7):
        out = con.construct('ray99',length=k)
        h,_ = count(out['rows'])
        require(h == 99*2**k+23,'99 ray actual count')
        ray.append(dict(k=k,n=out['order'],H=h))
        if out['order'] <= 9:
            require(recursive_count(out['rows']) == h,'99 ray backtracking')
            recursive += 1
    base = [40,17,3,6,13,30]
    reverse = [sum(1<<b for b in range(6) if base[b] >> a & 1) for a in range(6)]
    twin = [endpoint_state(rows,[0,3]) for rows in (base,reverse)]
    require(twin == [(37,14,99),(37,14,117)],'same count and defect do not fix slope')
    left, ls, _ = con.seed('ternary30',2)
    right, rs, _ = con.seed('ternary90',1)
    permutation = [0,2,3,4,5,1,6]
    require({permutation[v] for v in ls} == set(rs), 'selected-set isomorphism')
    require(all(bool(left[i]>>j&1) == bool(right[permutation[i]]>>permutation[j]&1)
                for i in range(7) for j in range(7)), 'duplicate ternary family isomorphism')
    sporadic = []
    for code,mask,marked in [(9281,17,2),(131713,68,5)]:
        rows = con.ib.cc.decode(7,code)
        ds = []
        for size in range(1,7):
            labels,expanded = con.ib.inflate(rows,marked,size)
            selected = [i for i,v in enumerate(labels) if mask >> v & 1]
            h,w,u = endpoint_state(expanded,selected)
            ds.append(h-w)
        expected = [1-4*stirling_weights(a,3)[2] if code == 9281 else
                    1+2*stirling_weights(a,3)[1]-2*stirling_weights(a,3)[2] for a in range(1,7)]
        require(ds == expected,'sporadic size controls')
        sporadic.append(dict(code=code,selected=mask,marked=marked,defects_sizes_1_to_6=ds))
    report = dict(named_families=named,actual_tail_controls=tails,backtracking_controls=recursive,
                  C3_substitution_controls=nontransitive,ray99=ray,converse_states=twin,
                  ternary_duplicate_isomorphism=permutation,sporadic=sporadic)
    (HERE / 'verification.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    print(json.dumps({k:v for k,v in report.items() if k!='named_families'},indent=2))


if __name__ == '__main__':
    main()
