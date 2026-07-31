#!/usr/bin/env python3
"""Canonical exact compositor for projected k=3 heights 270 through 247.

The final theorem is not canonical until the upstream source/output and this
file's semantic digest are non-None and normal/-O transcripts agree.  Every
terminal already has |S|>alpha=d/R.  THM-2984 makes primitive-unit phases
height-free, while the local translated high-danger band has length d/7 and
contains at most ceil(d/7)<=alpha residue classes.  Thus one clean cell
escapes every local band without selecting a torsion pair.  The smaller
centered-band number beta(d) is deliberately not used.
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
UPSTREAM_SHA256 = "4137ab250def3ad6a66b4c75a5e1b5b1a82ba4100b00ea5f8616faa46fb501a9"
UPSTREAM_OUTPUT_SHA256 = "eea98955f91371d38d95cdeeb88b60a2305d34d0bddd2ea26570af8eede1b8e3"
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
# high-ray recurrence checks.
TERMINAL = {
    270:(146,159,18,300073), 268:(57,60,7,117819),
    260:(3016,6236,606,1292385), 259:(8,9,1,31081),
    257:(21,24,2,13487), 253:(32,35,2,92111),
    250:(645,886,149,1475580), 249:(58,64,3,41118),
    247:(74,80,7,168058),
}
TOTALS = (202807,84483,114027,4297)
TERMINAL_TOTALS = (4057,7553,795,3531712)
WALL = Q(13,132)
ALIGNED_CAP = Q(36,91)
INHERITED_LEDGER = 375674
FINAL_LEDGER = 375251
NEXT_HEIGHT = 246
NEXT_COUNT = 194
TOKEN_ABSENT_CONTROL = (250,(1,4,9,10,12,14))
TRANSLATED_BAND_CONTROL = (28,(0,1,2,3),Q(-1,2),Q(7,2))
TERMINAL_BODY_COUNT = 69
MIN_TWO_HIGH_GAP = Q(1438897,5584915336)
GLOBAL_R_HISTOGRAM = ((2,177),(3,93),(4,358),(5,933),(6,1065),(7,4927))
MIN_ALPHA_SLACK = 1


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    payload=path.read_bytes()
    require(b"\r" not in payload.replace(b"\r\n",b""),f"bare CR in {path}")
    return hashlib.sha256(payload.replace(b"\r\n",b"\n")).hexdigest()


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
    token_absent_control_seen = False
    for line in ATLAS.read_text().splitlines():
        match = pattern.match(line)
        if match is None:
            continue
        body = tuple(map(int,match.group(1).split(",")))
        L, high, first = map(int,match.group(2,3,4))
        counts[first] += 1
        if first in LEVELS:
            suffix_text = match.group(5)
            # The ray quotient needs an explicit later-high slot exactly when
            # the fixed first label lies below the wall.  Above it, strict
            # label ordering makes every actual later label high, while this
            # quotient harmlessly relaxes order.  The printed ``HIGH-TAIL``
            # token describes an older slotwise maximizer, not this predicate.
            gate = first < high
            require(high == WALL.numerator*L//WALL.denominator+1, (first,body,"high"))
            if (first,body) == TOKEN_ABSENT_CONTROL:
                require(gate and "HIGH-TAIL:" not in suffix_text and "1810:" in suffix_text,
                        (first,body,"token-absent control changed"))
                token_absent_control_seen = True
            tasks.append((first,body,L,high,gate))
    require({z:counts[z] for z in LEVELS} == ROW_COUNTS, "row counts changed")
    for z in range(NEXT_HEIGHT+1,max(LEVELS)+1):
        require(counts[z] == ROW_COUNTS.get(z,0), (z,counts[z]))
    require(counts[NEXT_HEIGHT] == NEXT_COUNT, "next frontier changed")
    require(len(tasks) == 423, len(tasks))
    require(token_absent_control_seen,"token-absent hostile control disappeared")
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


def translated_band_certificate(cells,d):
    """Freeze the deterministic local-band data, with no pair choice."""
    residues=len({cell%d for cell in cells})
    R,alpha,_divisors=cayley_alpha(d)
    ceiling=(d+6)//7
    require(ceiling<=alpha,(d,ceiling,alpha,"threshold hierarchy"))
    require(residues>alpha,(d,residues,alpha,"translated-band gate"))
    return (d,R,alpha,residues,ceiling,residues-ceiling,residues-alpha,len(cells))


def translated_band_hostile_control():
    """Refute use of THM-2984's smaller centered beta after translation."""
    d,residues,left,right=TRANSLATED_BAND_CONTROL
    beta=2*((d-1)//14)+1
    ceiling=(d+6)//7
    require(right-left==Q(d,7),(d,left,right,"band length"))
    require(all(left<r<right for r in residues),(d,residues,left,right))
    require((beta,ceiling,len(residues))==(3,4,4),(beta,ceiling,len(residues)))
    return (d,residues,left,right,beta,ceiling)


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
    certificate_cache={}
    R_hist=Counter()
    minimum_alpha_slack=None; minimum_band_slack=None; minimum_alpha_gain=None
    digest=hashlib.sha256()
    for ds,high_d,low_rows,excess in cases:
        labels=tuple(sorted(z for _d,z in low_rows))
        if labels not in clean_cache:
            clean_cache[labels]=eng.fixed_safe_cells(stream,labels)
        key=(labels,high_d)
        if key not in certificate_cache:
            certificate_cache[key]=translated_band_certificate(clean_cache[labels],high_d)
        certificate=certificate_cache[key]
        R_hist[certificate[1]]+=1
        alpha_slack=certificate[6]; band_slack=certificate[5]
        alpha_gain=certificate[2]-certificate[4]
        minimum_alpha_slack=alpha_slack if minimum_alpha_slack is None else min(minimum_alpha_slack,alpha_slack)
        minimum_band_slack=band_slack if minimum_band_slack is None else min(minimum_band_slack,band_slack)
        minimum_alpha_gain=alpha_gain if minimum_alpha_gain is None else min(minimum_alpha_gain,alpha_gain)
        digest.update(repr((ds,high_d,low_rows,excess,certificate)).encode()+b"\n")
    return (first,body,stream.L,stream.high_floor,stream.lower-stream.first_delta,len(residual),
            gap,gap_witness,len(zero),len(cases),len(low_pairs),low_signs,ray_checks,
            tuple(sorted(R_hist.items())),minimum_alpha_slack,minimum_band_slack,
            minimum_alpha_gain,len(certificate_cache),digest.hexdigest())


def ft(q):
    return "NONE" if q is None else f"{q.numerator}/{q.denominator}"


def render(records,terminal_records,atlas_counts,pins):
    totals=tuple(sum(row[i] for row in records) for i in (9,10,11,12))
    require(totals==TOTALS,totals)
    level_totals={z:tuple(sum(row[i] for row in records if row[0]==z) for i in (9,10,11,12)) for z in LEVELS}
    require(level_totals==LEVEL_TOTALS,level_totals)
    grouped=defaultdict(list)
    for row in terminal_records: grouped[row[0]].append(row)
    require(len(terminal_records)==TERMINAL_BODY_COUNT,len(terminal_records))
    require(set(grouped)==set(TERMINAL),tuple(sorted(grouped)))
    summaries={}; global_R=Counter()
    aggregate=[0,0,0,0]
    for z,expected in TERMINAL.items():
        rows=grouped[z]
        values=(sum(r[8] for r in rows),sum(r[9] for r in rows),sum(r[10] for r in rows),sum(r[12] for r in rows))
        require(values==expected[:4],(z,values))
        Rh=Counter()
        for row in rows:
            Rh.update(dict(row[13]))
        global_R.update(Rh)
        summaries[z]=(values,tuple(sorted(Rh.items())),min(r[14] for r in rows),
                      min(r[15] for r in rows),min(r[16] for r in rows))
        aggregate=[a+b for a,b in zip(aggregate,values)]
    require(tuple(aggregate)==TERMINAL_TOTALS,aggregate)
    require(sum(global_R.values())==TERMINAL_TOTALS[1],global_R)
    require(tuple(sorted(global_R.items()))==GLOBAL_R_HISTOGRAM,global_R)
    require(min(row[6] for row in terminal_records)==MIN_TWO_HIGH_GAP,"two-high minimum")
    require(min(row[14] for row in terminal_records)==MIN_ALPHA_SLACK,"alpha slack")
    require(INHERITED_LEDGER-len(records)==FINAL_LEDGER,"ledger")
    control_record=next(row for row in records if (row[0],row[1])==TOKEN_ABSENT_CONTROL)
    require(control_record[5] and control_record[12]>0,
            (TOKEN_ABSENT_CONTROL,"projected-gate control lost"))
    control_terminal=next(row for row in terminal_records if (row[0],row[1])==TOKEN_ABSENT_CONTROL)
    require(control_terminal[6]>0,(TOKEN_ABSENT_CONTROL,"two-high control gap"))
    band_control=translated_band_hostile_control()
    semantic_payload=(LEVELS,WALL,ALIGNED_CAP,records,terminal_records,atlas_counts,level_totals,summaries,
                      tuple(sorted(global_R.items())),pins,INHERITED_LEDGER,FINAL_LEDGER,
                      NEXT_HEIGHT,NEXT_COUNT,band_control)
    semantic=hashlib.sha256(repr(semantic_payload).encode()).hexdigest()
    if SEMANTIC_SHA256 is not None:
        require(semantic==SEMANTIC_SHA256,"semantic digest changed")
    lines=[
        "LRC14 projected k=3 z1=270 through z1=247 cardinality-torsion descent",
        f"engine_sha256={sha(ENGINE)}",f"atlas_source_sha256={sha(ATLAS_SOURCE)}",f"atlas_output_sha256={sha(ATLAS)}",
        f"upstream_source_sha256={pins[0]}",f"upstream_output_sha256={pins[1]}",
        "dependency_hash_basis=SHA-256 after CRLF-to-LF normalization;bare CR rejected",
        "scope=423 exact atlas body rows at fifteen occupied heights270..247;all intervening heights checked empty;no finite high-label horizon",
        f"frontier_totals=states:{totals[0]};crude:{totals[1]};status:{totals[2]};residual:{totals[3]}",
        f"independent_exact_farkas_checks={totals[2]}/{totals[2]}:PASS;solver_basis_not_frozen",
        "logical_split=every residual row here has first<high_floor,so THM-2941(25f) supplies an explicit later-high slot;when first>=high_floor strict ordering instead makes all actual later labels high;the atlas's printed HIGH-TAIL token is not used;duplicate-permitting >=2-high upper has strict positive exact gap on every residual body;therefore exactly one high in this residual universe",
        f"token_absent_hostile_control=z1:{TOKEN_ABSENT_CONTROL[0]};E:{TOKEN_ABSENT_CONTROL[1]};first:{control_record[0]};high_floor:{control_record[3]};first_below_high_floor:{control_record[5]};printed_HIGH-TAIL:false;printed_exact_high:1810;residual_states:{control_record[12]};two_high_gap:{ft(control_terminal[6])};first_failed_implication=absence of the literal tail token does not negate the projected maximum wall",
        "repaired_gate=THM-2941(25f) gives max(Z)>13L/132;integer labels and first<floor(13L/132)+1 force some later label>=high_floor;the representative may be an exact high label or HIGH-TAIL",
        f"terminal_reduction=residual_bodies:{TERMINAL_BODY_COUNT};zero_high_hostile_passes:{aggregate[0]};one_high_cases:{aggregate[1]};body_distinct_low_pairs:{aggregate[2]};unit_ray_checks:{aggregate[3]};minimum_two_high_gap:{ft(MIN_TWO_HIGH_GAP)}",
        "cayley_alpha=G_d edges have nonzero difference-order<=7;R=max({r:r|d,2<=r<=7} union {1});cosets mod d/R are R-cliques and [0,d/R) is independent;therefore alpha(G_d)=d/R",
        f"pair_graph_boundary=every fixed-safe residue set has |S|>alpha=d/R;minimum_alpha_slack:{MIN_ALPHA_SLACK};the equality interval proves alpha sharp only for forcing a short-order pair;no pair is selected,hashed,or used",
        f"translated_band_gate=THM-2984 makes primitive-unit cell phases height-free;for every local translate the strict high-danger band has length d/7 and contains at most ceil(d/7)<=alpha residues;every terminal has |S|>alpha;minimum_band_slack:{min(row[15] for row in terminal_records)};minimum_alpha_minus_ceil7:{min(row[16] for row in terminal_records)};each primitive unit and local translate has a fixed-safe cell",
        f"centered_beta_hostile=d:{band_control[0]};S:{band_control[1]};open_interval:({ft(band_control[2])},{ft(band_control[3])});length:d/7;beta:{band_control[4]};ceil7:{band_control[5]};all_four_residues_fit:true;first_failed_implication=centered beta does not bound a translated local danger band",
        f"optimal_R_histogram={dict(sorted(global_R.items()))}",
        "strict_local_boundary=an open circular interval of length d/7 contains at most ceil(d/7) lattice residues;the bound is sharp and endpoint equality is safe",
    ]
    for z in LEVELS:
        lines.append(f"LEVEL;z1={z};rows={ROW_COUNTS[z]};states={level_totals[z][0]};crude={level_totals[z][1]};status={level_totals[z][2]};residual={level_totals[z][3]};terminal={summaries.get(z)}")
    for row in records:
        first,body,L,high,d1,gate,trials,checks,signs,states,crude,status,residual,residual_ds,mhist,instance_digest,verified,representative,stage_digest=row
        require(verified==status,(first,body,"verified"))
        lines.append(f"BODY;z1={first};E={body};L={L};high={high};d1={d1};high_gate={gate};trials={trials};checks={checks};signs={dict(signs)};states={states};crude={crude};status={status};residual={residual};M={dict(mhist)};instance_sha256={instance_digest};verified={verified};representative={representative};stage_sha256={stage_digest};residual_sha256={hashlib.sha256(repr(residual_ds).encode()).hexdigest()}")
    for row in terminal_records:
        first,body,L,high,required,residual,gap,witness,zero,cases,pairs,low_signs,checks,Rh,alpha_slack,band_slack,alpha_gain,certs,digest=row
        lines.append(f"TERMINAL;z1={first};E={body};L={L};high={high};required={ft(required)};residual={residual};two_high_gap={ft(gap)};gap_witness={witness};zero_high={zero};cases={cases};low_pairs={pairs};low_signs={dict(low_signs)};ray_checks={checks};optimal_R={dict(Rh)};minimum_alpha_slack={alpha_slack};minimum_translated_band_slack={band_slack};minimum_alpha_minus_ceil7={alpha_gain};certificates={certs};case_sha256={digest}")
    lines += [f"atlas_exact_check=selected_rows:423;all gaps270..247 empty;next occupied z1={NEXT_HEIGHT} with {NEXT_COUNT} rows",
              f"ledger_decrement={INHERITED_LEDGER}-423={FINAL_LEDGER};decrement counts body rows,not quotient states",
              "nonconsequence=projected scalar-atlas descent only;does not close z1=246,k<=1,the rung,or LRC(14)",
              "dependency=THM-2984 primitive-unit height-free phase law is used;its centered-beta,gcd-stratum,transporter,and nonflag-complex gates are not used for the translated local band",
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
