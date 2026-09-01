#!/usr/bin/env python3
"""Exact cover certificates for the frozen three-row rank-8/9 atlas."""
import argparse, csv
from pathlib import Path
from ortools.sat.python import cp_model

N=145
MIXED=['022083a5','00b0832c','22c0a124','0228832c','10e12602','1284812c','10923282','10550a81','106a600a','14050236']
RANK8=['30612208','00c22289','02209125','02508125','1038220a','00510a89','0280c125','028085d0','1210b208','126c8008','0246a088','28c42202','101c1016','121c0007','18472200','082a0a88','03308124','14a1020a','12658200','32209208','30b48100','08168281']
DUAL={4:2165,10:1007,12:917,13:185,16:3983,17:4283,26:1211,27:5,28:1103,32:2069,34:293,36:941,37:125,38:545,42:959,43:2693,50:1769,52:2315,53:1829,55:1001,58:713,64:1523,68:809,69:3893,73:959,81:1871,87:767,91:3317,94:185,97:3755,98:1823,100:3941,101:347,107:293,117:3893,118:2123,120:1253,121:809,123:917,124:1151,127:4373,128:581,129:767,130:47,132:1487,133:1055,136:3893,137:1871,138:2363,139:3893,141:4235,142:1823}

def require(c,m):
    if not c:raise RuntimeError(m)

def fnv(words):
    h=0xcbf29ce484222325
    for x in words:
        for b in range(8):h^=(x>>(8*b))&255;h=(h*0x100000001b3)&((1<<64)-1)
    return h

def load(path):
    out=[]
    for r in csv.DictReader(path.open(newline='')):
        v=int(r['q96_lo'],16)|(int(r['q96_hi'],16)<<64)|(int(r['q100'],16)<<116)|(int(r['q210'],16)<<129)
        r['v']=v;r['bits']=tuple(i for i in range(N) if v>>i&1);r['count8']=int(r['count8']);r['count9']=int(r['count9']);out.append(r)
    require(len(out)==14619,'type count');return out

def maximal(rows,rank8):
    xs=[r for r in rows if not rank8 or r['count8']];xs.sort(key=lambda r:(-len(r['bits']),r['v']));out=[]
    for r in xs:
        if not any(r['v']&u['v']==r['v'] for u in out):out.append(r)
    return out

def selected(rows,masks,rank8):
    key='least8' if rank8 else 'least9';by={r[key]:r for r in rows if r[key]!='-'};require(all(m in by for m in masks),'cover mask missing');return [by[m] for m in masks]

def seg(r):
    return (sum(i<116 for i in r['bits']),sum(116<=i<129 for i in r['bits']),sum(i>=129 for i in r['bits']))

def verify_cover(rows,masks,rank8):
    chosen=selected(rows,masks,rank8);u=0
    for r in chosen:u|=r['v']
    require(u==(1<<N)-1,'cover incomplete');return chosen

def rank8_infeasible21(max8):
    m=cp_model.CpModel();x=[m.new_bool_var(f'x{j}') for j in range(len(max8))]
    for i in range(N):m.add(sum(x[j] for j,r in enumerate(max8) if r['v']>>i&1)>=1)
    m.add(sum(x)<=21);s=cp_model.CpSolver();s.parameters.num_search_workers=1;s.parameters.subsolvers.append('max_lp');s.parameters.random_seed=0;s.parameters.max_time_in_seconds=60
    st=s.solve(m);require(st==cp_model.INFEASIBLE,'rank8 <=21 not refuted');return s

def main():
    p=argparse.ArgumentParser();p.add_argument('atlas',type=Path);a=p.parse_args();rows=load(a.atlas);maxmix=maximal(rows,False);max8=maximal(rows,True)
    cm=verify_cover(rows,MIXED,False);c8=verify_cover(rows,RANK8,True)
    loads=[sum(DUAL.get(i,0) for i in r['bits']) for r in rows];require(max(loads)<=10000,'dual overload');require(sum(DUAL.values())==90128,'dual sum')
    s=rank8_infeasible21(max8)
    print('LRC595_EXACT_COVER_AUDIT_V2')
    print(f'TYPES ALL {len(rows)} MAXIMAL_MIXED {len(maxmix)} RANK8 {sum(r["count8"]>0 for r in rows)} MAXIMAL_RANK8 {len(max8)}')
    print('MIXED_COVER SIZE 10 FNV %016x'%fnv(int(m,16) for m in MIXED))
    for m,r in zip(MIXED,cm):print(f'MASK {m} RANK {int(m,16).bit_count()} TAGGED {seg(r)[0]}/{seg(r)[1]}/{seg(r)[2]} TICKS {r["least9_tick96"]}/{r["least9_tick100"]}/{r["least9_tick210"]}')
    print('MIXED_DUAL DENOM 10000 SUM 90128 VALUE 90128/10000 MAX_LOAD %d/10000 SATURATED_TYPES %d SUPPORT '%(max(loads),sum(x==max(loads) for x in loads))+','.join(f'{i}:{w}' for i,w in sorted(DUAL.items())))
    print('MIXED_EXACT_MINIMUM 10 REASON ceil(90128/10000)=10 AND EXPLICIT_SIZE10_COVER')
    print('RANK8_COVER SIZE 22 FNV %016x'%fnv(int(m,16) for m in RANK8))
    for m,r in zip(RANK8,c8):print(f'MASK {m} TAGGED {seg(r)[0]}/{seg(r)[1]}/{seg(r)[2]} TICKS {r["least8_tick96"]}/{r["least8_tick100"]}/{r["least8_tick210"]}')
    print(f'RANK8_LOWER CP_SAT_VERSION 9.15.6755 TEST <=21 STATUS INFEASIBLE MAXIMAL_TYPES {len(max8)} BRANCHES {s.num_branches} CONFLICTS {s.num_conflicts} WORKERS 1 SUBSOLVER max_lp RANDOM_SEED 0')
    print('RANK8_EXACT_MINIMUM 22 REASON exhaustive_infeasibility_at21 AND EXPLICIT_SIZE22_COVER')
    print('GRIDS 96:36482318832960 100:91205797082400 210:18241159416480; ticks are pair-specific numerators and are never cross-compared')
    print('SCOPE COMPLETE_FROZEN_RANK8_RANK9_RESPONSE_ATLAS_ONLY_NO_PHYSICAL_ENTRY_NO_LRC14')
    print('VERDICT PASS EXACT_MINIMA_MIXED10_RANK8_22')

if __name__=='__main__':main()
