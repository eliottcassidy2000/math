"""Audit the anti-involution colour law for arc-disjoint HP squares.

Universe: every standard-blue fixed-path presentation of orders 3..7,
and every directed Hamiltonian path arc-disjoint from the standard path.
All checks are explicit, and remain active under python -O.
"""
from collections import Counter
import json
from pathlib import Path
from detours_probe import G, paths

HERE=Path(__file__).resolve().parent


def opposite(a):
    full=(1<<len(a))-1
    return tuple(full^(1<<u)^row for u,row in enumerate(a))


def reverse(a,p):
    return opposite(G.complement(a,p))


def reflection(p):
    r=[None]*len(p)
    for u,v in zip(p,p[::-1]):r[u]=v
    return tuple(r)


def main():
    results=[];positive={};hostile={}
    for n in range(3,8):
        P=tuple(range(n));rP=reflection(P);ep=G.edge_set(P);stat=Counter()
        for a in G.blue_tilings(n):
            stat['blue_presentations']+=1
            b=reverse(a,P)
            for Q in paths(a):
                if ep&G.edge_set(Q):continue
                stat['arc_disjoint_pairs']+=1
                aq=reverse(a,Q)
                G.need(G.is_path(b,Q) and G.is_path(aq,P),'commuting directedness')
                G.need(reverse(b,Q)==reverse(aq,P),'commuting reversal identity')
                qb=G.blue(a,Q);pb_after=G.blue(aq,P);qb_after=G.blue(b,Q)
                same=reflection(Q)==rP
                G.need(pb_after==same,'exact P transport colour law')
                if qb:G.need(qb_after==same,'blue Q transport colour law')
                else:G.need(not same and not pb_after,'black Q forces middle black')
                pattern=f'{int(qb)}{int(pb_after)}{int(qb_after)}'
                stat['pattern_'+pattern]+=1
                if not qb and not qb_after:
                    words=[Q,P[::-1],Q[::-1]];state=a
                    for word in words:
                        G.need(G.is_path(state,word) and not G.blue(state,word),'black square detour')
                        state=G.complement(state,word)
                    G.need(state==G.complement(a,P),'exact labelled detour endpoint')
                    positive.setdefault('three_black',dict(rows=a,P=P,Q=Q))
                if qb and not same:
                    # The two partner states are connected through the double reversal.
                    G.need(not pb_after and not qb_after,'two-blue partner link')
                    positive.setdefault('two_black_partner_link',dict(rows=a,P=P,Q=Q))
                if qb and same:hostile.setdefault('all_blue_square',dict(rows=a,P=P,Q=Q))
        results.append(dict(n=n,counts=dict(stat)))
    # C5 proper-tail rotation is nonautomorphic but remains blue.
    a=tuple(sum(1<<v for v in range(5) if (v-u)%5 in (1,2)) for u in range(5))
    P=tuple(range(5));Q=(0,2,3,4,1)
    G.need(G.is_path(a,Q) and G.blue(a,P) and G.blue(a,Q),'C5 rotation hostile')
    G.need(G.induced(a,Q)!=a,'tail rotation must not be an automorphism')
    hostile['nonautomorphic_tail_rotation_can_be_blue']=dict(rows=a,P=P,Q=Q)
    out=dict(status='FINITE-EXACT audit of proved colour laws',results=results,
             positive_controls=positive,hostile_controls=hostile)
    (HERE/'detours_square.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
