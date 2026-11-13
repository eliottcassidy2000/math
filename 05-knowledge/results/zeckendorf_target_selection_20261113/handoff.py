"""Literal arbitrary-module handoff between chirality and the size repair."""
from itertools import combinations, permutations
from pathlib import Path
import json
from chirality import require, tt, cyclic, profile, tournament

I = (46,8,58,48,35,2)


def inflate(frame, site, module):
    blocks = [module if i == site else (0,) for i in range(len(frame))]
    offsets = []; n = 0
    for block in blocks:
        offsets.append(n); n += len(block)
    rows = []
    for i, block in enumerate(blocks):
        external = sum(((1 << len(blocks[j]))-1) << offsets[j]
                       for j in range(len(frame)) if frame[i] >> j & 1)
        rows.extend((row << offsets[i]) | external for row in block)
    return tuple(rows), offsets


def birth(rows, root):
    n = len(rows); output = list(rows)
    output[root] |= 1 << n
    output.append(((1 << n)-1) ^ (1 << root))
    return tuple(output)


def relabel(rows, mapping):
    require(sorted(mapping) == list(range(len(rows))), 'bijection required')
    output = [0]*len(rows)
    for i,row in enumerate(rows):
        output[mapping[i]] = sum(1 << mapping[j] for j in range(len(rows)) if row >> j & 1)
    return tuple(output)


def run_coefficients(frame, site, kind):
    outside = [i for i in range(len(frame)) if i != site and not(kind == 'q' and i == 4)]
    answer = []
    for k in (1,2):
        total = 0
        for pos in combinations(range(k+len(outside)),len(outside)):
            for perm in permutations(outside):
                word = [site]*(k+len(outside))
                for p,label in zip(pos,perm): word[p] = label
                if kind == 's' and word[0] != 4: continue
                if kind == 'F':
                    cut = word.index(4)+1
                    paths = (word[:cut],word[cut:])
                else:
                    paths = (word,)
                if all(all((frame[x]>>y)&1 for x,y in zip(path,path[1:])) for path in paths):
                    total += 1
        answer.append(total)
    return answer


def main():
    coefficients = {str(slot):{kind:run_coefficients(I,slot,kind)
                    for kind in ('H','s','q','F')} for slot in (0,1)}
    require(coefficients['0'] == {'H':[17,9],'s':[3,0],'q':[3,0],'F':[21,9]}, 'old slot tensor')
    require(coefficients['1'] == {'H':[17,9],'s':[3,1],'q':[3,1],'F':[21,11]}, 'handoff slot tensor')
    bank = [tournament(n, code) for n in range(1,4)
            for code in range(1 << (n*(n-1)//2))]+[tt(4),(2,4,9,3)]
    controls = []
    for M in bank:
        m=len(M); cm,_,_,_ = profile(M);h=cm[1];p=2*h+(cm[2] if m>1 else 0)
        initial = cyclic((tt(2),M,tt(1)))
        nr = birth(birth(initial,0),m+3)
        slot1, off1 = inflate(I,1,M)
        # Original names r=0,a=1,M=2..m+1,b=m+2,x=m+3,y=m+4.
        mapping=[off1[4],off1[5]]+list(range(off1[1],off1[1]+m))+[off1[3],off1[0],off1[2]]
        require(relabel(nr,mapping)==slot1,'literal NR/frame handoff')
        slot0, off0 = inflate(I,0,M)
        # Re-express slot0 in the slot1 naming, retaining the SAME old M vertices.
        to_old=[None]*len(slot0)
        for j in range(m):to_old[off0[0]+j]=off1[1]+j
        to_old[off0[1]]=off1[0]
        for v in (2,3,4,5):to_old[off0[v]]=off1[v]
        migrated = relabel(slot0,to_old)
        changed={(i,j) for i in range(len(slot1)) for j in range(i+1,len(slot1))
                 if ((slot1[i]>>j)&1)!=((migrated[i]>>j)&1)}
        module=set(range(off1[1],off1[1]+m));u=off1[0]
        expected={tuple(sorted((i,j))) for i in module for j in (u,off1[2],off1[5])}
        expected.update(tuple(sorted((u,off1[j]))) for j in (2,5))
        require(changed==expected and len(changed)==3*m+2,'exact relocation cuts')
        root=off1[4]
        require(slot1[root]==migrated[root],'retained full old root star')
        a=profile(slot1,root)[2];b=profile(migrated,root)[2]
        require(a==(9*p-h,p+h,p+h,11*p-h),'incoming four-state')
        require(b==(9*p-h,3*h,3*h,9*p+3*h),'outgoing four-state')
        # Add the sole-exit helper: it beats only the now-singleton role1.
        final=list(migrated);new=len(final)
        for i in range(new):
            if i != u:final[i]|=1<<new
        final.append(1<<u);final=tuple(final)
        require(profile(final,root)[2]==(15*p-h,5*h,5*h,15*p+5*h),'composed repair target')
        # At the earliest stage the root keeps its arcs into every old M vertex.
        require(all(initial[0]>>j&1 for j in range(2,m+2)),'initial root/module port')
        controls.append({'module':M,'h_p':[h,p], 'counts':[3*p-h,9*p-h,15*p-h],
                         'incoming':a,'outgoing':b,'final':profile(final,root)[2],
                         'changed_old_edges':len(changed),'root':root,'slot1_NR_map':mapping})
    result={'status':'PASS','complete_run_coefficients':coefficients,'literal_modules':controls,
            'scope':'all labelled modules through3 plus two order4 controls; arbitrary-module proof via at-most-two runs',
            'reproduce':'python3 -B 05-knowledge/results/zeckendorf_target_selection_20261113/handoff.py'}
    Path(__file__).with_name('handoff.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS exact module relocation and composed chirality/neutral-return handoff:',len(controls),'literal modules')


if __name__ == '__main__':
    main()
