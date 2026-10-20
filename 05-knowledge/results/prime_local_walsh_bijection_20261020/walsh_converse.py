#!/usr/bin/env python3
"""Exact converse/parity block refinement of walsh_support.py."""
from __future__ import annotations
import json
from collections import Counter
from pathlib import Path
import sympy as sp
from walsh_support import need, permdata, transform, cofactor_matching, vp


def main():
    source=json.loads(Path('walsh_support.json').read_text()); out=[]
    for d in source:
        n=d['n']; T=d['tournaments']; G=d['even_graphs']; M=d['matrix']; N=len(T)
        pairs,perms=permdata(n); top=(1<<len(pairs))-1
        codes={t['mask']:i for i,t in enumerate(T)}
        opposite=[]
        for t in T:
            mask=top^t['mask']; opposite.append(codes[min(transform(mask,a)^c for a,b,c in perms)])
        singles=[i for i in range(N) if opposite[i]==i]
        doubles=[(i,j) for i,j in enumerate(opposite) if i<j]
        erows=sorted(singles+[i for i,j in doubles]); orows=[i for i,j in doubles]
        ecols=[j for j,g in enumerate(G) if g['mask'].bit_count()%2==0]
        ocols=[j for j,g in enumerate(G) if g['mask'].bit_count()%2]
        need(len(erows)==len(ecols) and len(orows)==len(ocols),'block counts')
        E=[[M[i][j] for j in ecols] for i in erows]; O=[[M[i][j] for j in ocols] for i in orows]
        for i in range(N):
            for j,g in enumerate(G): need(M[opposite[i]][j]==(-1)**(g['mask'].bit_count()%2)*M[i][j],'converse sign')
        for name,A,rows,cols in [('even',E,erows,ecols),('odd',O,orows,ocols)]:
            mat=sp.Matrix(A) if A else sp.zeros(0,0)
            det=int(mat.det()); need(det and abs(det)&(abs(det)-1)==0,'block determinant')
            if len(pairs):
                weights=[T[i]['aut']*(2 if name=='even' and i in singles else 1) for i in rows]
                Dg=sp.diag(*[G[j]['aut'] for j in cols]); Dt=sp.diag(*weights)
                need(mat*Dg*mat.T==2**(len(pairs)-1)*Dt,'block weighted orthogonality')
        record={'n':n,'self_converse':[T[i]['mask'] for i in singles],
            'nonself_converse_pairs':[[T[i]['mask'],T[j]['mask']] for i,j in doubles],
            'even_block_rows':[T[i]['mask'] for i in erows],'even_block_columns':[G[j]['mask'] for j in ecols],
            'odd_block_rows':[T[i]['mask'] for i in orows],'odd_block_columns':[G[j]['mask'] for j in ocols],
            'even_matrix':E,'odd_matrix':O,'prime_data':{}}
        for p in [3,5,7]:
            em,_=cofactor_matching(E,p); om,_=cofactor_matching(O,p)
            matched={erows[i]:ecols[j] for i,j in em}
            matched.update({opposite[orows[i]]:ocols[j] for i,j in om})
            need(sorted(matched)==list(range(N)),'whole domain')
            need(sorted(matched.values())==list(range(N)),'whole codomain')
            need(all(M[i][j]%p for i,j in matched.items()),'support units')
            need(all(vp(T[i]['aut'],p)==vp(G[j]['aut'],p) for i,j in matched.items()),'valuation preserved')
            need(all(G[matched[i]]['mask'].bit_count()%2==0 for i in singles),'self converse parity')
            need(all((G[matched[i]]['mask'].bit_count()+G[matched[j]]['mask'].bit_count())%2 for i,j in doubles),'pair opposite parity')
            sh=Counter(vp(T[i]['aut'],p) for i in singles)
            gh=Counter()
            for g in G: gh[vp(g['aut'],p)]+=(-1)**(g['mask'].bit_count()%2)
            need(dict(sh)=={k:v for k,v in gh.items() if v},'signed valuation distribution')
            record['prime_data'][str(p)]={'matching':[[T[i]['mask'],G[j]['mask'],M[i][j]] for i,j in sorted(matched.items())],
                'self_converse_valuation_histogram':dict(sorted(sh.items()))}
        print(f'n={n}: self={len(singles)} nonself_pairs={len(doubles)} even-edge={len(ecols)} odd-edge={len(ocols)}; all parity+valuation+support controls pass')
        if n<=5:
            print('  p3 matching:',record['prime_data']['3']['matching'])
        out.append(record)
    Path('walsh_converse.json').write_text(json.dumps(out,indent=2)+'\n')
    print('ALL CONVERSE BLOCK CHECKS PASSED')

if __name__=='__main__':main()
