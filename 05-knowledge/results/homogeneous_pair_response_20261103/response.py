"""Full named seed response polynomials, cubic contrasts, and five module runs."""
import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import subprocess
import tempfile

import check

HERE = Path(__file__).resolve().parent
SEEDS = [[16,5,33,39,46,3], [16,5,49,39,42,3], [18,4,49,39,42,3]]
need = check.need


def seed(j, r, reverse=False):
    q, n = SEEDS[j], r+5
    rows = [q[u] & 31 for u in range(5)] + [0]*r
    for u in range(5):
        for v in range(5, n):
            if q[u] >> 5 & 1:
                rows[u] |= 1 << v
            else:
                rows[v] |= 1 << u
    for u in range(5, n):
        for v in range(u+1, n):
            rows[v if reverse else u] |= 1 << (u if reverse else v)
    return rows


def run_histogram(rows, module):
    n = len(rows)
    dp = [dict() for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][(v, int(v in module))] = 1
    for mask in range(1, 1 << n):
        for (v, k), count in dp[mask].items():
            todo = rows[v] & ~mask
            while todo:
                bit = todo & -todo
                todo -= bit
                w = bit.bit_length()-1
                state = (w, k+int(w in module and v not in module))
                dest = dp[mask | bit]
                dest[state] = dest.get(state, 0)+count
    totals = Counter()
    for (_, k), count in dp[-1].items():
        totals[k] += count
    return dict(sorted(totals.items()))


def masks(n):
    return [0]+[1 << u for u in range(n)]+[
        (1 << u) | (1 << v) for u, v in combinations(range(n), 2)]


def subs(mask):
    answer = [mask]
    while mask:
        mask = (mask-1) & answer[0]
        answer.append(mask)
    return answer


def coefficients(table, ms):
    values = {(a,c):table[i*len(ms)+j] for i,a in enumerate(ms) for j,c in enumerate(ms)}
    out = []
    for a in ms:
        for c in ms:
            coef = sum((-1)**((a.bit_count()+c.bit_count())-(b.bit_count()+d.bit_count()))*
                       values[b,d] for b in subs(a) for d in subs(c))
            out.append(coef)
    return out


def evaluate(coef, ms, a, c):
    return sum(value for k,value in enumerate(coef)
               if ms[k//len(ms)] & a == ms[k//len(ms)]
               and ms[k%len(ms)] & c == ms[k%len(ms)])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=HERE/'response.json')
    parser.add_argument('--tables', type=Path, default=HERE/'response_tables.json')
    args = parser.parse_args()
    pins = {name:sha256((HERE/name).read_bytes()).hexdigest()
            for name in ('response.py', 'response.cpp', 'check.py')}
    blocks, specs = [], []
    for r in (2, 7):
        for j in range(3):
            rows = seed(j,r)
            n, ms = len(rows), masks(len(rows))
            full = (1 << n)-1
            high = [(full,full),(full,0),(0,full),(full^1,full^2),
                    (full//3,2*(full//3)),(full^32,full^64)]
            queries = [(a,c) for a in ms for c in ms]+high
            for graph in (rows,check.opposite(rows)):
                blocks.append((graph,queries))
            specs.append((j,r,rows,ms,high))
    # Actual partial-module queries and a quartic hostile are separate from whole converse.
    extras = []
    for r in range(2,8):
        for j in range(3):
            for reverse in (False,True):
                graph = seed(j,r,reverse)
                extras.append((graph,[(32,0),(64,0)]))
    blocks += extras
    hostile = [24,17,19,6,8]
    reversed_hostile = [8,17,19,6,9]
    support = [1,3,5+1,5+4]
    hostile_queries = []
    for m in range(16):
        mask = sum(1 << v for i,v in enumerate(support) if m>>i & 1)
        hostile_queries.append((mask & 31,mask >> 5))
    blocks += [(hostile,hostile_queries),(reversed_hostile,hostile_queries)]
    text = [str(len(blocks))]
    for rows,queries in blocks:
        text += [str(len(rows))+' '+str(len(queries)), ' '.join(map(str,rows))]
        text += [str(a)+' '+str(c) for a,c in queries]
    with tempfile.TemporaryDirectory(prefix='pair-response-') as temp:
        exe = Path(temp)/'response.exe'
        subprocess.run(['g++','-std=c++17','-O3',str(HERE/'response.cpp'),'-o',str(exe)],check=True)
        result = subprocess.run([str(exe)],input='\n'.join(text)+'\n',text=True,
                                capture_output=True,check=True)
    stream = iter(map(int,result.stdout.split()))
    values = [[next(stream) for _ in queries] for rows,queries in blocks]
    need(next(stream,None) is None,'unexpected counter output')
    cases, summaries = [], []
    independent_queries = 0
    for idx,(j,r,rows,ms,high) in enumerate(specs):
        n, m = len(rows), len(ms)**2
        first, second = values[2*idx],values[2*idx+1]
        f,g = coefficients(first[:m],ms),coefficients(second[:m],ms)
        delta = [x-y for x,y in zip(f,g)]
        degrees = Counter()
        for k,value in enumerate(delta):
            a,c = ms[k//len(ms)],ms[k%len(ms)]
            degree = a.bit_count()+c.bit_count()
            need(not value or degree <= 3,'whole-converse quartic')
            need(not value or a == 0 or c == 0 or degree >= 2,'degree type')
            if degree == 2 and (not a or not c):
                need(value == 0,'pure quadratic contrast')
            if value:
                degrees[degree] += 1
        for k,(a,c) in enumerate(high):
            need(evaluate(f,ms,a,c)==first[m+k], 'full polynomial high query')
            need(evaluate(g,ms,a,c)==second[m+k], 'op polynomial high query')
        battery = [k for k in range(m)
                   if ms[k//len(ms)].bit_count()+ms[k%len(ms)].bit_count() <= 1
                   or (ms[k//len(ms)] and ms[k%len(ms)]
                       and ms[k//len(ms)].bit_count()+ms[k%len(ms)].bit_count()<=3)]
        need(len(battery)==1+2*n+n**3, 'cubic battery size')
        # Rebuild the difference from only the allowed raw queries; pure pair values inferred.
        only = {}
        for k in battery:
            only[ms[k//len(ms)],ms[k%len(ms)]] = first[k]-second[k]
        for a in ms:
            for c in ms:
                if not a and c.bit_count()==2:
                    only[a,c] = sum(only[0,1<<v] for v in range(n) if c>>v&1)
                if not c and a.bit_count()==2:
                    only[a,c] = sum(only[1<<v,0] for v in range(n) if a>>v&1)
        reconstructed = [0]*m
        for k in battery:
            a,c = ms[k//len(ms)],ms[k%len(ms)]
            reconstructed[k] = sum((-1)**(a.bit_count()+c.bit_count()-b.bit_count()-d.bit_count())*
                                   only[b,d] for b in subs(a) for d in subs(c))
        need(reconstructed==delta,'full named cubic from restricted queries')
        # Independent Python forward DP on fixed low and high masks at each order.
        controls = [(0,0),(32,0),(64,0),high[3]]
        for a,c in controls:
            expected = evaluate(f,ms,a,c)
            need(check.h(check.extend(rows,a,c))==expected,'independent forward DP')
            independent_queries += 1
        if r == 2:
            # Literal word enumeration at n+2=9 on one nonuniform query per state.
            need(check.literal_h(check.extend(rows,32,0))==evaluate(f,ms,32,0),'literal seed query')
        pair = [evaluate(f,ms,1<<v,0) for v in (5,6)]
        op_pair = [evaluate(g,ms,1<<v,0) for v in (5,6)]
        need(pair[1]>pair[0] and op_pair[1]<op_pair[0], 'pair sign in full table')
        summaries.append(dict(state=j,r=r,n=n,H=first[0],full_battery=m,
                              cubic_battery=len(battery),nonzero_contrast_by_degree=dict(degrees),
                              pair_queries=pair,converse_pair_queries=op_pair))
        cases.append(dict(state=j,r=r,rows=rows,masks=ms,query_order='row-major (a,c)',
                          values=first[:m],converse_values=second[:m],
                          coefficients=f,converse_coefficients=g,contrast_coefficients=delta,
                          high_queries=[dict(a=a,c=c,F=first[m+k],F_op=second[m+k])
                                        for k,(a,c) in enumerate(high)]))
    pair_cases = []
    offset = 2*len(specs)
    for k in range(0,len(extras),2):
        f,g = values[offset+k:offset+k+2]
        need(f[1]>f[0] and g[1]<g[0], 'partial reversal sign')
        pair_cases.append(dict(r=2+k//6,state=(k//2)%3,F=f,F_module_op=g))
    hv,hg = values[-2:]
    quartic = sum((-1)**(4-m.bit_count())*(hv[m]-hg[m]) for m in range(16))
    need(quartic==2, 'incoming partial-reversal quartic')
    word = [5,10,6,11,7,0,4,3,8,1,2,9]
    histograms = []
    for j in range(3):
        graph = seed(j,5)
        ext = check.extend(graph,32,64)
        need(sorted(word)==list(range(12)), 'five-run word universe')
        need(all(ext[u]>>v&1 for u,v in zip(word,word[1:])), 'five-run word arcs')
        hist = run_histogram(ext,set(range(5,10)))
        need(max(hist)==5 and hist[5]>0, 'sharp run cap')
        need(sum(hist.values())==check.h(ext), 'histogram/ordinary DP')
        histograms.append(dict(state=j,run_counts=hist,query=[32,64],witness=word))
    # Independent recheck of the incoming scalar law, r1..7.
    for r in range(1,8):
        for j in range(3):
            need(check.h(seed(j,r))==4*3**r+9*2**r-3+2*j,'scalar seed law')
    artifact = dict(status='FINITE-EXACT complete named Boolean response polynomials for six fixed cores',
                    names='old0..4 exterior; module5..r+4; x=n,y=n+1; x->y',
                    cases=cases)
    tables_text = json.dumps(artifact,sort_keys=True,separators=(',',':'))+'\n'
    args.tables.write_text(tables_text,encoding='utf8',newline='\n')
    receipt = dict(status='PASS: finite named controls; all-order claims proved separately',
                   cases=summaries,actual_compiled_queries=sum(len(q) for _,q in blocks),
                   independent_forward_queries=independent_queries,literal_seed_queries=3,
                   partial_module_pairs=pair_cases,partial_module_quartic=dict(coefficient=quartic,
                   monomial='a1*a3*c1*c4',rows=hostile,reversed_rows=reversed_hostile),
                   five_run_controls=histograms,
                   tables_sha256=sha256(tables_text.encode()).hexdigest(),source_sha256=pins)
    need(all(sha256((HERE/name).read_bytes()).hexdigest()==digest for name,digest in pins.items()),
         'source changed during run')
    args.output.write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n',encoding='utf8',newline='\n')


if __name__=='__main__':
    main()

