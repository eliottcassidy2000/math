#!/usr/bin/env python3
"""Native reset continuation; retains source weights and full integer clocks.

No convergence oracle. Snapshot audits cover every observed edge, including
exterior excursions. Old stage-19 hashes are a separate reproduction control.
"""
import argparse
import hashlib
import json
import sys
import time
from pathlib import Path
import numpy as np
from collatz_compact_resets_20261015f import Graph
from collatz_clock_tail_20261017b import audit, selection, tail, controls
from collatz_interval_atlas_20261016e import norms


def rooted_bracket(X, graph, values):
    # Only sources already connected to stopped ROOT contribute to R. The
    # signed nonROOT terminal l1 mass supplies a rigorous uncertainty radius.
    value = np.full(graph.M+1, -1, dtype=np.int64)
    value[2::3] = 2
    value[(graph.clock & 1) != 0] *= -1
    value[graph.endpoint != 1] = 0
    value[0] = 0
    prefix = np.cumsum(value)
    row = np.arange(X, 2*X, dtype=np.int64)
    rooted = prefix[4*row]-prefix[row]
    upper = np.array(values,dtype=np.int64)
    radius = upper-np.abs(rooted)
    if not np.all(radius>=0):
        raise RuntimeError('nonROOT radius is negative')
    lower = np.maximum(0,np.abs(rooted)-radius)
    if not np.all(lower<=upper):raise RuntimeError('invalid rooted bracket')
    detail = dict(certified_ROOT_sources=int(np.count_nonzero(graph.endpoint[1:]==1)),
                  total_sources=graph.M,nonROOT_terminal_count=int(len(np.unique(graph.endpoint[graph.endpoint!=1]))-1),
                  lower_tail=tail(list(map(int,lower))),rooted_tail=tail(list(map(int,np.abs(rooted)))),
                  radius_tail=tail(list(map(int,radius))),rows_with_positive_lower=int(np.count_nonzero(lower)),
                  rows_exact=int(np.count_nonzero(radius==0)),
                  bracket_endpoints_sha256=hashlib.sha256(np.column_stack((rooted-radius,rooted+radius)).tobytes()).hexdigest())
    return detail


def run(k, stages, output, root_stats=False, audit_receipt=None):
    X = 1 << k
    graph = Graph(8*X)
    rows = []
    before = None
    start = time.monotonic()
    for stage in range(1, max(stages)+1):
        cols, inverse = np.unique(graph.endpoint[1:], return_inverse=True)
        graph.extend(cols)
        graph.close(cols, inverse)
        del cols, inverse
        print(f'k={k} stage={stage} edges={graph.edge_count} elapsed={time.monotonic()-start:.1f}s', file=sys.stderr, flush=True)
        if stage not in stages:
            continue
        if audit_receipt is None:
            receipt = audit(graph)
        else:
            raw=Path(audit_receipt).read_bytes();saved=json.loads(raw)
            match=[s for r in saved['runs'] if r['k']==k for s in r['stages'] if s['stage']==stage]
            if len(match)!=1:raise RuntimeError('no unique inherited full clock audit')
            receipt=match[0]['audit']
            if any(hashlib.sha256(v.tobytes()).hexdigest()!=receipt[key] for v,key in
                   ((graph.endpoint,'endpoint_sha256'),(graph.clock,'clock_sha256'))):
                raise RuntimeError('independent rerun does not reproduce the fully audited source arrays')
            receipt=dict(receipt,inherited_receipt_sha256=hashlib.sha256(raw).hexdigest())
        if stage==max(stages):
            # The completed audit is the edge receipt. The final coefficient
            # calculations use only authenticated endpoint and clock arrays.
            graph.outside.clear()
        values = norms(X, graph)
        if before is not None and not all(b <= a for a, b in zip(before, values)):
            raise RuntimeError('literal native row monotonicity failed')
        before = values
        record = dict(stage=stage, audit=receipt, fine_tail=tail(values),
                      selections=[selection(X, values, max(1, (1 << (k//2))//4), stage)],
                      max_clock=int(graph.clock.max()), observed_edges=graph.edge_count)
        if root_stats:record['rooted_bracket']=rooted_bracket(X,graph,values)
        rows.append(record)
        payload = dict(status='FINITE-EXACT native reset continuation; uniform debt OPEN',
                       runs=[dict(k=k, X=X, stages=rows)], completed_through=stage,
                       requested_stages=stages)
        Path(output).write_text(json.dumps(payload, sort_keys=True, separators=(',', ':'))+'\n')
        print(f'snapshot k={k} stage={stage} saved', file=sys.stderr, flush=True)
    return payload


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--k', type=int, required=True)
    p.add_argument('--stages', type=int, nargs='+', required=True)
    p.add_argument('--output', required=True)
    p.add_argument('--root-stats', action='store_true')
    p.add_argument('--audit-receipt', help='Reuse a full audit only after exact endpoint/clock hash reproduction')
    a = p.parse_args()
    if a.k < 4 or any(s < 1 for s in a.stages):
        raise ValueError('positive stages and k>=4 required')
    controls()
    run(a.k, sorted(set(a.stages)), a.output,a.root_stats,a.audit_receipt)
