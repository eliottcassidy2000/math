"""Independent literal controls for universal nonempty-context R4 dissolution."""
from pathlib import Path
from hashlib import sha256
import json
import arithmetic_dissolution as d

HERE = Path(__file__).resolve().parent
WORDS = [('A',2,3,0,1,'W'), ('W',3,0,'A',2,1),
         ('W',3,1,0,'A',2), ('A',2,0,3,1,'W')]


def first_scc_size(t):
    n = len(t)
    reach = [t[i] | (1 << i) for i in range(n)]
    for v in range(n):
        for u in range(n):
            if reach[u] >> v & 1:
                reach[u] |= reach[v]
    # All presentations have the standard HP, so vertex0 lies in first SCC.
    return sum(bool(reach[0] >> v & 1 and reach[v] & 1) for v in range(n))


def ledger(a,w,ap,wp):
    m,k = len(a),len(w)
    state = d.join(d.join(a,d.B),w)
    flags,steps,excesses = [],[],[]
    for word in WORDS:
        path = tuple(v for token in word for v in
                     (ap if token=='A' else tuple(m+4+x for x in wp)
                      if token=='W' else (m+token,)))
        d.need(d.is_path(state,path),'interior path has absent arc')
        flags.append(d.blue(state,path))
        excesses.append(state[path[0]].bit_count()+state[path[-1]].bit_count()
                        -(len(state)-1))
        steps.append(path)
        state=d.k_move(state,path)
    d.need(state==d.join(d.join(a,d.D),w),'interior labelled target differs')
    score=w[wp[0]].bit_count()
    expected=[m+2-score,score-m-1] if k>=2 else [m+1,-m]
    d.need(excesses[1:3]==expected,'middle endpoint formulas differ')
    d.need(not flags[2],'third interior move is blue')
    if flags[1]:
        # Relabel W by its current path before using the positional SCC helper.
        relabeled=tuple(sum(1<<j for j,v in enumerate(wp) if d.arc(w,u,v))
                        for u in wp)
        d.need(first_scc_size(relabeled)>=m+4,'blue P2 violates SCC lower bound')
    return flags,steps


def solve(a,w):
    target=d.join(d.join(a,d.D),w)
    target_active=any(not d.blue(target,p) for p in d.paths(target))
    standard_a,standard_w=tuple(range(len(a))),tuple(range(len(w)))
    initial_flags,_=ledger(a,w,standard_a,standard_w)
    route=None
    for ap in d.paths(a):
        for wp in d.paths(w):
            flags,steps=ledger(a,w,ap,wp)
            if not any(flags):
                route=steps
                break
        if route is not None:
            break
    d.need((route is not None)==target_active,'universal path-selection theorem fails')
    return dict(active=target_active,standard_blue=initial_flags,
                A=list(a),W=list(w),route=route)


def main():
    # Both formal ports distinct; then all singleton-port identifications.
    for a in ((0,),(2,0)):
        for w in ((0,),(2,0)):
            ledger(a,w,tuple(range(len(a))),tuple(range(len(w))))
    rows,witnesses=[],{}
    for m in range(1,4):
        for k in range(1,7):
            count=active=repairs=blue_p2=0
            for a in d.presentations(m):
                for w in d.presentations(k):
                    r=solve(a,w)
                    count+=1
                    active+=int(r['active'])
                    repairs+=int(r['active'] and any(r['standard_blue']))
                    blue_p2+=int(r['standard_blue'][1])
                    if r['active'] and r['standard_blue'][1]:
                        witnesses.setdefault('middle_path_selection',r)
                    if r['active'] and r['standard_blue'][0]:
                        witnesses.setdefault('source_endpoint_selection',r)
                    if r['active'] and r['standard_blue'][3]:
                        witnesses.setdefault('target_endpoint_selection',r)
            expected=2**((m-1)*(m-2)//2+(k-1)*(k-2)//2)
            d.need(count==expected,'complete presentation rectangle mismatch')
            rows.append(dict(A_order=m,W_order=k,presentations=count,
                             active_targets=active,standard_path_repairs=repairs,
                             standard_blue_second_moves=blue_p2))
    extra=[]
    for a,w in ((d.B,d.B), ((6,12,8,1),d.join((2,4,1),(6,12,8,1)))):
        r=solve(a,w)
        d.need(r['active'],'extra coupled-endpoint control is pure')
        extra.append(r)
        if r['standard_blue'][0]:witnesses.setdefault('source_endpoint_selection',r)
        if r['standard_blue'][3]:witnesses.setdefault('target_endpoint_selection',r)
    d.need('middle_path_selection' in witnesses,'no actual blue-middle hostile')
    # Exceptional endpoint normalization, independently replayed.
    source=d.join((0,),d.B)
    exceptional=[(0,2,3,4,1),(1,0,2,3,4)]
    state=source
    for p in exceptional:
        d.need(not d.blue(state,p),'n5 exceptional route is blue')
        state=d.k_move(state,p)
    d.need(state==(28,29,8,16,4),'n5 exceptional labelled target differs')
    report=dict(scope='all standard-HP presentations, A orders1..3 and W orders1..6',
                census=rows,branch_witnesses=witnesses,extra_controls=extra,
                n5_exception=dict(source=source,paths=exceptional,target=state),
                source_sha256_lf=sha256(Path(__file__).read_text().encode()).hexdigest(),
                utility_sha256_lf=sha256((HERE/'arithmetic_dissolution.py').read_text().encode()).hexdigest())
    output=json.dumps(report,indent=2)+'\n'
    (HERE/'arithmetic_interior.json').write_text(output,encoding='utf-8')
    print(output,end='')


if __name__=='__main__':
    main()
