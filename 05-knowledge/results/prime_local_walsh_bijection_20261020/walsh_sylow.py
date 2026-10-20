#!/usr/bin/env python3
"""Direct signed p-orbit and Sylow transporter controls, complete n<=6 atlas.
This bounded verifier exploits that every tournament Aut group here has order
1,3,5,9, so its p-Sylow is the whole group or the identity. It explicitly checks
that premise; no such shortcut is asserted for larger n.
"""
import json
from collections import Counter
from pathlib import Path
from walsh_support import need, permdata, transform, vp


def main():
    atlas=json.loads(Path('walsh_support.json').read_text()); result=[]
    for d in atlas:
        n=d['n']; T=d['tournaments']; G=d['even_graphs']; M=d['matrix']
        pairs,perms=permdata(n); sizes={t['aut'] for t in T}
        need(sizes<={1,3,5,9},'bounded Sylow shortcut premise')
        auts=[[rec for rec in perms if transform(t['mask'],rec[0])^rec[2]==t['mask']] for t in T]
        need(all(len(a)==t['aut'] for a,t in zip(auts,T)),'direct automorphism orders')
        signed_copies=[]
        for g in G:
            copies={}
            for images,inv,target_inv in perms:
                y=transform(g['mask'],images); s=(-1)**((g['mask']&inv).bit_count()%2)
                need(y not in copies or copies[y]==s,'consistent signed copies')
                copies[y]=s
            signed_copies.append(copies)
        checks=0; records=[]
        for p in [3,5,7]:
            matched={a:b for a,b,v in d['prime_data'][str(p)]['matching']}
            class_hist=Counter()
            for i,t in enumerate(T):
                power=p**vp(t['aut'],p)
                if power==1: continue
                need(power==t['aut'],'Aut itself is p-group in bounded scope')
                P=auts[i]
                # The sorted orbit sizes on vertices identify C3 actions and C3xC3 here.
                # Derive vertex action from pair maps? Reconstruct tuples separately.
                import itertools
                actual=[q for q,rec in zip(itertools.permutations(range(n)),perms) if rec in P]
                seen=set(); vertex_orbits=[]
                for v in range(n):
                    if v not in seen:
                        orb={q[v] for q in actual}; seen|=orb; vertex_orbits.append(len(orb))
                signature=(power,tuple(sorted(vertex_orbits)))
                class_hist[str(signature)]+=1
                for j,g in enumerate(G):
                    if M[i][j]%p==0: continue
                    fixed=[]; signed_fixed=0
                    for y,eps in signed_copies[j].items():
                        if all(transform(y,rec[0])==y for rec in P):
                            fixed.append(y); signed_fixed+=eps*(-1)**((t['mask']&y).bit_count()%2)
                    need(fixed,'p-unit fixed copy')
                    need((signed_fixed-M[i][j])%p==0,'localization congruence')
                    checks+=1
                    if matched[t['mask']]==g['mask']:
                        need(p**vp(g['aut'],p)==power,'same Sylow size')
                        records.append({'p':p,'tournament_mask':t['mask'],'canonical_graph_mask':g['mask'],
                            'fixed_graph_copy':min(fixed),'graph_transporter':list(next(q for q,rec in zip(itertools.permutations(range(n)),perms) if transform(g['mask'],rec[0])==min(fixed))),'sylow_order':power,'vertex_orbit_sizes':list(signature[1]),
                            'fixed_copy_count':len(fixed),'fixed_signed_sum':signed_fixed,'Phi':M[i][j]})
        result.append({'n':n,'nontrivial_p_unit_edges_checked':checks,'explicit_matched_transporters':records})
        print(f'n={n}: direct signed localization on {checks} nontrivial p-unit edges; {len(records)} matched Sylow transporters')
    Path('walsh_sylow.json').write_text(json.dumps(result,indent=2)+'\n')
    print('ALL SYLOW LOCALIZATION CHECKS PASSED')

if __name__=='__main__':main()
