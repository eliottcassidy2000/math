#!/usr/bin/env python3
"""DEVELOPMENT compositor for projected k=3 heights 270 through 247.

The final theorem is not canonical until the upstream source/output and this
file's semantic digest are non-None and normal/-O transcripts agree.  The
mathematical terminal is the exact Cayley bound alpha(G_d)=d/R, where edges
have difference-order at most seven and R is the largest divisor of d at
most seven.  Thus |S|>d/R is a sharp cardinality-only torsion certificate.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import multiprocessing as mp
import re
from collections import Counter, defaultdict
from fractions import Fraction as Q
from math import gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ENGINE = ROOT / "04-computation/lrc14_j7_k3_z297_ray_status_torsion_closure_thm2941.py"
ATLAS_SOURCE = ROOT / "04-computation/lrc14_j7_k3_projected_scalar_body_atlas_thm2941.py"
ATLAS = ROOT / "05-knowledge/results/lrc14_j7_k3_projected_scalar_body_atlas_thm2941.out"
UPSTREAM = ROOT / "04-computation/lrc14_j7_k3_z275_to_z272_septimal_torsion_descent_thm2941.py"
UPSTREAM_OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_k3_z275_to_z272_septimal_torsion_descent_thm2941.out"
OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_k3_z270_to_z247_cardinality_torsion_descent_thm2941.out"

ENGINE_SHA256 = "d062c7ac8ebf6a433c8fb1543293e941c85625e2eb40b82fcf05fc2404539b0a"
ATLAS_SOURCE_SHA256 = "2af6d96882f336a409a8657070ed76a75c09a53b3789101b83103b051e864ded"
ATLAS_SHA256 = "cee82237ce1f51729813b9c916edd3353204c18172abe1d71278dee2c5562eda"
UPSTREAM_SHA256 = None
UPSTREAM_OUTPUT_SHA256 = None
SEMANTIC_SHA256 = None

LEVELS = (270, 268, 265, 260, 259, 257, 256, 255, 254, 253, 252, 251, 250, 249, 247)
ROW_COUNTS = {270:26,268:27,265:3,260:140,259:16,257:4,256:1,255:3,254:1,253:4,252:1,251:3,250:176,249:10,247:8}
LEVEL_TOTALS = {
    270:(5123,1577,3397,149), 268:(2801,1581,1161,59),
    265:(3,1,2,0), 260:(126049,53003,69850,3196),
    259:(407,124,274,9), 257:(271,123,124,24), 256:(2,0,2,0),
    255:(17,3,14,0), 254:(0,0,0,0), 253:(240,40,165,35),
    252:(2,2,0,0), 251:(3,2,1,0), 250:(66653,27629,38343,681),
    249:(391,171,156,64), 247:(845,227,538,80),
}
# zero-high hostile passes, one-high cases, body-distinct low pairs,
# high-ray recurrence checks, least cardinality-forced r histogram.
TERMINAL = {
    270:(146,159,18,300073,((2,115),(3,23),(4,11),(5,3),(6,6),(7,1))),
    268:(57,60,7,117819,((2,55),(3,2),(4,2),(6,1))),
    260:(3016,6236,606,1292385,((2,5276),(3,507),(4,200),(5,126),(6,127))),
    259:(8,9,1,31081,((2,8),(4,1))),
    257:(21,24,2,13487,((2,22),(3,2))),
    253:(32,35,2,92111,((2,30),(3,2),(4,2),(6,1))),
    250:(645,886,149,1475580,((2,703),(3,74),(4,51),(5,23),(6,32),(7,3))),
    249:(58,64,3,41118,((2,58),(3,5),(5,1))),
    247:(74,80,7,168058,((2,68),(3,5),(4,4),(5,2),(6,1))),
}
TOTALS = (202807,84483,114027,4297)
TERMINAL_TOTALS = (4057,7553,795,3531712)
WALL = Q(13,132)
ALIGNED_CAP = Q(36,91)
INHERITED_LEDGER = 375674
FINAL_LEDGER = 375251
NEXT_HEIGHT = 246
NEXT_COUNT = 194


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


require(sha(ENGINE) == ENGINE_SHA256, "engine changed")
require(sha(ATLAS_SOURCE) == ATLAS_SOURCE_SHA256, "atlas source changed")
require(sha(ATLAS) == ATLAS_SHA256, "atlas changed")
eng = load("z270_z247_engine", ENGINE)


def upstream_pins(development):
    pinned = UPSTREAM_SHA256 is not None and UPSTREAM_OUTPUT_SHA256 is not None
    require(pinned or development, "upstream is unpinned; pass --development only for a noncanonical replay")
    if not pinned:
        return ("UNPINNED-DEVELOPMENT", "UNPINNED-DEVELOPMENT")
    require(sha(UPSTREAM) == UPSTREAM_SHA256, "upstream source changed")
    require(sha(UPSTREAM_OUTPUT) == UPSTREAM_OUTPUT_SHA256, "upstream output changed")
    require("projected k3 cap<=270;next exact frontier=270" in UPSTREAM_OUTPUT.read_text(), "upstream conclusion changed")
    return UPSTREAM_SHA256, UPSTREAM_OUTPUT_SHA256


def atlas_tasks():
    pattern = re.compile(r"^row=E=([0-9,]+);.*;L=([0-9]+);high=([0-9]+);z1=([0-9]+);.*;suffix=(.*);lower=")
    counts = Counter()
    tasks = []
    for line in ATLAS.read_text().splitlines():
        match = pattern.match(line)
        if match is None:
            continue
        body = tuple(map(int,match.group(1).split(",")))
        L, high, first = map(int,match.group(2,3,4))
        counts[first] += 1
        if first in LEVELS:
            # The ray quotient's projected high obligation is active exactly
            # when the fixed first label lies below the wall.  The atlas's
            # printed ``HIGH-TAIL`` token describes its older slotwise
            # maximizer; it is not the logical gate.
            gate = first < high
            require(high == WALL.numerator*L//WALL.denominator+1, (first,body,"high"))
            tasks.append((first,body,L,high,gate))
    require({z:counts[z] for z in LEVELS} == ROW_COUNTS, "row counts changed")
    for z in range(NEXT_HEIGHT+1,max(LEVELS)+1):
        require(counts[z] == ROW_COUNTS.get(z,0), (z,counts[z]))
    require(counts[NEXT_HEIGHT] == NEXT_COUNT, "next frontier changed")
    require(len(tasks) == 423, len(tasks))
    return tuple(tasks), tuple((z,counts[z]) for z in range(NEXT_HEIGHT,max(LEVELS)+1))


def frozen_state(state):
    return tuple(sorted(state.items()))


def evaluate(task):
    first,body,atlas_L,atlas_high,gate = task
    eng.FIRST = first
    eng.ray.FIRST = first
    stream = eng.ray.Stream(body)
    require((stream.L,stream.high_floor)==(atlas_L,atlas_high),(first,body,"atlas"))
    trials,states,checks,signs = eng.ray.ray_quotient_states(stream)
    crude,status,residual = eng.exact_common_status_screen(stream,states)
    require(set(states)==set(crude)|set(status)|set(residual),(first,body,"partition"))
    require(not(set(crude)&set(status)),(first,body,"overlap"))
    require(not residual or gate,(first,body,"residual without high gate"))
    for ds,witness in crude.items():
        gap,q,M,target,capacity = witness
        D=lcm(*ds)
        require(M==D//q and target==dict(stream.target_data(D))[q],(first,body,ds,"crude data"))
        require(capacity==sum(eng.ray.local.fibre_cap(D,d,q) for d in ds),(first,body,ds,"capacity"))
        require(gap==target-capacity and gap>0,(first,body,ds,"crude gap"))
    instance_digest=hashlib.sha256()
    representative=None
    arcs_cache={}
    histogram_cache={}
    verified=0
    for ds,witness in sorted(status.items()):
        q,M,marginals,cap_set,histogram,certificate=witness
        D=lcm(*ds)
        require(M==D//q,(first,body,ds,"status cofactor"))
        rebuilt,capacities=eng.ray.local.hunter_status_data(D,ds,q)
        require(rebuilt==marginals and tuple(sorted(set(capacities)))==cap_set,(first,body,ds,"status data"))
        if D not in arcs_cache:
            arcs_cache[D]=eng.ray.fibre.projected_support_arcs(D,stream.ranges)
        key=(D,q)
        if key not in histogram_cache:
            histogram_cache[key]=eng.ray.fibre.residue_load_histogram(arcs_cache[D],q)
        require(histogram_cache[key]==histogram,(first,body,ds,"histogram"))
        eng.independent_farkas_check(q,marginals,capacities,histogram,certificate)
        instance=witness[:-1]
        if representative is None:
            representative=(first,body,ds,instance)
        instance_digest.update(f"{first}|{body}|{ds}|{instance}\n".encode())
        verified+=1
    stage_digest=hashlib.sha256()
    for ds in sorted(states):
        if ds in crude:
            stage,witness="crude",crude[ds]
        elif ds in status:
            stage,witness="status",status[ds][:-1] # solver basis is not semantic
        else:
            stage,witness="residual",None
        stage_digest.update(repr((ds,frozen_state(states[ds]),stage,witness)).encode()+b"\n")
    for ds in residual:
        labels=states[ds]["labels"]
        require(labels[0]==first and len(labels)==len(set(labels))==4,(first,body,ds,"labels"))
        require(tuple(sorted(stream.L//gcd(stream.L,z) for z in labels))==ds,(first,body,ds,"orders"))
        require(sum(z>=stream.high_floor for z in labels[1:])==1,(first,body,ds,"maximizer grammar"))
    sign_totals={s:sum(n for (_d,t),n in signs.items() if t==s) for s in (-1,0,1)}
    require(sign_totals[-1]==sign_totals[1],(first,body,"sign symmetry"))
    return (first,body,stream.L,stream.high_floor,stream.first_d,gate,trials,checks,
            tuple(sorted(sign_totals.items())),len(states),len(crude),len(status),len(residual),
            tuple(residual),tuple(sorted(Counter(w[1] for w in status.values()).items())),
            instance_digest.hexdigest(),verified,representative,stage_digest.hexdigest())


_ALPHA={}


def cayley_alpha(d):
    if d in _ALPHA:
        return _ALPHA[d]
    divisors=tuple(r for r in range(2,8) if d%r==0)
    R=max(divisors,default=1)
    alpha=d//R
    # alpha residue classes modulo alpha are R-cliques; [0,alpha) is independent.
    require(all(R//gcd(R,j)<=7 for j in range(1,R)),(d,R,"cliques"))
    require(all(d//gcd(d,j)>7 for j in range(1,alpha)),(d,R,"equality control"))
    _ALPHA[d]=(R,alpha,divisors)
    return _ALPHA[d]


def torsion_certificate(cells,d):
    cell_for={}
    for cell in cells:
        cell_for.setdefault(cell%d,cell)
    residues=tuple(sorted(cell_for))
    R,alpha,divisors=cayley_alpha(d)
    require(len(residues)>alpha,(d,len(residues),alpha,"sharp cardinality gate"))
    least=next(r for r in divisors if len(residues)>d//r)
    quotient=d//least
    buckets=defaultdict(list)
    for residue in residues:
        buckets[residue%quotient].append(residue)
    crowded=next(row for row in buckets.values() if len(row)>=2)
    a,b=crowded[:2]
    shift=(b-a)%d
    effective=d//gcd(d,shift)
    require(2<=effective<=least<=7,(d,least,effective))
    phase=tuple(min(u,effective-u) for u in range(1,effective) if gcd(u,effective)==1)
    require(phase and 7*min(phase)>=effective,(d,effective,"phase"))
    ca,cb=cell_for[a],cell_for[b]
    require((cb-ca)%d==shift,(d,"cell shift"))
    return (least,effective,R,alpha,len(residues),len(residues)-alpha,quotient,
            ca,cb,a,b,shift,min(phase),len(cells))


def terminal(task):
    first,body,residual=task
    eng.FIRST=first
    eng.ray.FIRST=first
    stream=eng.ray.Stream(body)
    needed={d for ds in residual for d in eng.suffix_slots(ds,stream.first_d)}
    low,high,low_signs,ray_checks=eng.build_literal_tables(stream,needed)
    require(dict(low_signs).get(-1,0)>0,(first,body,"negative lows lost"))
    gap,gap_witness=eng.duplicate_two_high_gap(stream,residual,low,high)
    require(gap>0,(first,body,"two-high gap",gap))
    zero=eng.zero_high_scalar_passes(stream,residual,low)
    cases=eng.one_high_cases(stream,residual,low,high)
    low_pairs={tuple(sorted(z for _d,z in rows)) for _ds,_hd,rows,_e in cases}
    clean_cache={}
    witness_cache={}
    least_hist=Counter(); effective_hist=Counter(); R_hist=Counter()
    minimum_slack=None
    digest=hashlib.sha256()
    for ds,high_d,low_rows,excess in cases:
        labels=tuple(sorted(z for _d,z in low_rows))
        if labels not in clean_cache:
            clean_cache[labels]=eng.fixed_safe_cells(stream,labels)
        key=(labels,high_d)
        if key not in witness_cache:
            witness_cache[key]=torsion_certificate(clean_cache[labels],high_d)
        witness=witness_cache[key]
        least_hist[witness[0]]+=1; effective_hist[witness[1]]+=1; R_hist[witness[2]]+=1
        minimum_slack=witness[5] if minimum_slack is None else min(minimum_slack,witness[5])
        require(all(eng.cell_clean(c,z,stream.L) for c in (witness[7],witness[8]) for z in (first,*labels)),(first,body,ds,"cell"))
        digest.update(repr((ds,high_d,low_rows,excess,witness)).encode()+b"\n")
    require(least_hist==effective_hist,(first,body,"effective order"))
    return (first,body,stream.L,stream.high_floor,stream.lower-stream.first_delta,len(residual),
            gap,gap_witness,len(zero),len(cases),len(low_pairs),low_signs,ray_checks,
            tuple(sorted(least_hist.items())),tuple(sorted(effective_hist.items())),
            tuple(sorted(R_hist.items())),minimum_slack,len(witness_cache),digest.hexdigest())


def ft(q):
    return "NONE" if q is None else f"{q.numerator}/{q.denominator}"


def render(records,terminal_records,atlas_counts,pins):
    totals=tuple(sum(row[i] for row in records) for i in (9,10,11,12))
    require(totals==TOTALS,totals)
    level_totals={z:tuple(sum(row[i] for row in records if row[0]==z) for i in (9,10,11,12)) for z in LEVELS}
    require(level_totals==LEVEL_TOTALS,level_totals)
    grouped=defaultdict(list)
    for row in terminal_records: grouped[row[0]].append(row)
    require(set(grouped)==set(TERMINAL),tuple(sorted(grouped)))
    summaries={}; global_least=Counter(); global_effective=Counter(); global_R=Counter()
    aggregate=[0,0,0,0]
    for z,expected in TERMINAL.items():
        rows=grouped[z]
        values=(sum(r[8] for r in rows),sum(r[9] for r in rows),sum(r[10] for r in rows),sum(r[12] for r in rows))
        require(values==expected[:4],(z,values))
        least=Counter(); effective=Counter(); Rh=Counter()
        for row in rows:
            least.update(dict(row[13])); effective.update(dict(row[14])); Rh.update(dict(row[15]))
        require(tuple(sorted(least.items()))==expected[4],(z,least))
        require(least==effective,(z,least,effective))
        global_least.update(least); global_effective.update(effective); global_R.update(Rh)
        summaries[z]=(values,tuple(sorted(least.items())),tuple(sorted(Rh.items())),min(r[16] for r in rows))
        aggregate=[a+b for a,b in zip(aggregate,values)]
    require(tuple(aggregate)==TERMINAL_TOTALS,aggregate)
    require(global_least==global_effective and sum(global_least.values())==TERMINAL_TOTALS[1],global_least)
    require(INHERITED_LEDGER-len(records)==FINAL_LEDGER,"ledger")
    semantic_payload=(LEVELS,WALL,ALIGNED_CAP,records,terminal_records,atlas_counts,level_totals,summaries,
                      tuple(sorted(global_least.items())),tuple(sorted(global_R.items())),pins,
                      INHERITED_LEDGER,FINAL_LEDGER,NEXT_HEIGHT,NEXT_COUNT)
    semantic=hashlib.sha256(repr(semantic_payload).encode()).hexdigest()
    if SEMANTIC_SHA256 is not None:
        require(semantic==SEMANTIC_SHA256,"semantic digest changed")
    lines=[
        "LRC14 projected k=3 z1=270 through z1=247 cardinality-torsion descent",
        f"engine_sha256={sha(ENGINE)}",f"atlas_source_sha256={sha(ATLAS_SOURCE)}",f"atlas_output_sha256={sha(ATLAS)}",
        f"upstream_source_sha256={pins[0]}",f"upstream_output_sha256={pins[1]}",
        "scope=423 exact atlas body rows at fifteen occupied heights270..247;all intervening heights checked empty;no finite high-label horizon",
        f"frontier_totals=states:{totals[0]};crude:{totals[1]};status:{totals[2]};residual:{totals[3]}",
        f"independent_exact_farkas_checks={totals[2]}/{totals[2]}:PASS;solver_basis_not_frozen",
        "logical_split=every residual row has first<high_floor,which is exactly the ray quotient's projected high-label gate;the atlas's printed HIGH-TAIL token is not used;duplicate-permitting >=2-high upper has strict positive exact gap on every residual body;therefore exactly one high",
        f"terminal_reduction=zero_high_hostile_passes:{aggregate[0]};one_high_cases:{aggregate[1]};body_distinct_low_pairs:{aggregate[2]};unit_ray_checks:{aggregate[3]}",
        "cayley_alpha=G_d edges have nonzero difference-order<=7;R=max({r:r|d,2<=r<=7} union {1});cosets mod d/R are R-cliques and [0,d/R) is independent;therefore alpha(G_d)=d/R",
        "cardinality_boundary=every fixed-safe residue set has |S|>d/R;this forces torsion order<=7;the independent equality set proves strictness is optimal for cardinality-only arguments",
        f"least_cardinality_forced_r_histogram={dict(sorted(global_least.items()))};optimal_R_histogram={dict(sorted(global_R.items()))}",
        "strict_seam=primitive high units preserve effective order s<=7;phase gap>=1/s>=1/7;strict-open radius1/14 dangers are disjoint,with only excluded endpoints meeting at s=7",
    ]
    for z in LEVELS:
        lines.append(f"LEVEL;z1={z};rows={ROW_COUNTS[z]};states={level_totals[z][0]};crude={level_totals[z][1]};status={level_totals[z][2]};residual={level_totals[z][3]};terminal={summaries.get(z)}")
    for row in records:
        first,body,L,high,d1,gate,trials,checks,signs,states,crude,status,residual,residual_ds,mhist,instance_digest,verified,representative,stage_digest=row
        require(verified==status,(first,body,"verified"))
        lines.append(f"BODY;z1={first};E={body};L={L};high={high};d1={d1};high_gate={gate};trials={trials};checks={checks};signs={dict(signs)};states={states};crude={crude};status={status};residual={residual};M={dict(mhist)};instance_sha256={instance_digest};verified={verified};representative={representative};stage_sha256={stage_digest};residual_sha256={hashlib.sha256(repr(residual_ds).encode()).hexdigest()}")
    for row in terminal_records:
        first,body,L,high,required,residual,gap,witness,zero,cases,pairs,low_signs,checks,least,effective,Rh,slack,certs,digest=row
        lines.append(f"TERMINAL;z1={first};E={body};L={L};high={high};required={ft(required)};residual={residual};two_high_gap={ft(gap)};gap_witness={witness};zero_high={zero};cases={cases};low_pairs={pairs};low_signs={dict(low_signs)};ray_checks={checks};least_r={dict(least)};effective={dict(effective)};optimal_R={dict(Rh)};minimum_optimal_slack={slack};certificates={certs};case_sha256={digest}")
    lines += [f"atlas_exact_check=selected_rows:423;all gaps270..247 empty;next occupied z1={NEXT_HEIGHT} with {NEXT_COUNT} rows",
              f"ledger_decrement={INHERITED_LEDGER}-423={FINAL_LEDGER};decrement counts body rows,not quotient states",
              "nonconsequence=projected scalar-atlas descent only;does not close z1=246,k<=1,the rung,or LRC(14)",
              "next_sidecar=THM-2984 primitive-unit reach is related but unused here",
              f"conclusion=all projected k3 rows at occupied heights270..247 are empty;with pinned z275..z272 package,projected k3 cap<={NEXT_HEIGHT};next exact frontier={NEXT_HEIGHT}",
              f"semantic_sha256={semantic}","all_exact_controls=PASS"]
    return "\n".join(lines)+"\n"


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--processes",type=int,default=8)
    parser.add_argument("--output",type=Path,default=OUTPUT)
    parser.add_argument("--development",action="store_true")
    args=parser.parse_args()
    require(args.processes>=1,"positive process count required")
    pins=upstream_pins(args.development)
    tasks,atlas_counts=atlas_tasks()
    if args.processes==1:
        records=tuple(evaluate(t) for t in tasks)
    else:
        with mp.get_context("spawn").Pool(min(args.processes,len(tasks))) as pool:
            records=tuple(pool.map(evaluate,tasks,chunksize=1))
    terminal_tasks=tuple((row[0],row[1],row[13]) for row in records if row[12])
    if args.processes==1:
        terminal_records=tuple(terminal(t) for t in terminal_tasks)
    else:
        with mp.get_context("spawn").Pool(min(args.processes,len(terminal_tasks))) as pool:
            terminal_records=tuple(pool.map(terminal,terminal_tasks,chunksize=1))
    payload=render(records,terminal_records,atlas_counts,pins)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(payload)
    print(payload,end="")


if __name__=="__main__":
    main()
