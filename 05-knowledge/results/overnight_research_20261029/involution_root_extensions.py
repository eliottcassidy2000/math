"""Complete marked index-two extension recurrence, with a retained involution."""
from itertools import product
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent
GATES=0


def need(ok,why):
    global GATES
    GATES+=1
    if not ok:raise ValueError(why)


def initial():
    # N=C2+C16; S(a,b)=(a+b,8a+5b). Character target C8 has S=3.
    n=32
    add=[[(a%2+b%2)%2+2*((a//2+b//2)%16) for b in range(n)] for a in range(n)]
    S=[(a%2+a//2)%2+2*((8*(a%2)+5*(a//2))%16) for a in range(n)]
    chars=[]
    for x,y in product((0,4),range(8)):
        f=[(x*(a%2)+y*(a//2))%8 for a in range(n)]
        if all(f[S[a]]==3*f[a]%8 for a in range(n)):chars.append(f)
    return {'add':add,'clock':S,'characters':chars,'chain':[]}


def extensions(A):
    add,S,chars=A['add'],A['clock'],A['characters']
    n=len(add)
    negatives=[next(b for b in range(n) if add[a][b]==0) for a in range(n)]
    doubles={add[a][a] for a in range(n)}
    covered,reps=set(),[]
    for v in range(n):
        if v not in covered:
            reps.append(v)
            covered.update(add[v][a] for a in doubles)
    for v in reps:
        difference=add[S[v]][negatives[v]]
        for w in range(n):
            if add[w][w]!=difference or add[w][S[w]]!=0:continue
            # 2x=v, Sx=x+w. These two equations are necessary and sufficient.
            newadd=[[add[add[a%n][b%n]][v] if a>=n and b>=n else add[a%n][b%n]
                     for b in range(2*n)] for a in range(2*n)]
            for a,b in product(range(2*n),repeat=2):
                newadd[a][b]+=n*((a//n+b//n)%2)
            newS=S+[add[S[a]][w]+n for a in range(n)]
            newchars=[]
            for f,z in product(chars,range(8)):
                if 2*z%8==f[v] and 2*z%8==f[w]:
                    newchars.append(f+[(y+z)%8 for y in f])
            yield {'add':newadd,'clock':newS,'characters':newchars,
                   'chain':A['chain']+[{'old_order':n,'root_v':v,'clock_w':w}]}


def verify_witness(A,f):
    add,S=A['add'],A['clock']
    n=len(add)
    need(sorted(set(f))==list(range(8)),'character is surjective, not just nonzero')
    need(all(S[S[a]]==a and f[S[a]]==3*f[a]%8 for a in range(n)),'joint involution and character clocks')
    for a,b in product(range(n),repeat=2):
        need(S[add[a][b]]==add[S[a]][S[b]],'literal clock is a group automorphism')
        need(f[add[a][b]]==(f[a]+f[b])%8,'literal full character homomorphism')
    old=initial()
    need(all(add[a][b]==old['add'][a][b] for a,b in product(range(32),repeat=2)),
         'all original group relations survive the successive attachments')
    need(S[:32]==old['clock'],'all original clock values retained')


def run():
    base=initial()
    need(not any(len(set(f))==8 for f in base['characters']),'original cubic pair has no primitive character')
    levels,frontier,witness=[],[base],None
    for depth in range(1,4):
        following=[]
        candidates=0
        for A in frontier:
            for B in extensions(A):
                candidates+=1
                f=next((f for f in B['characters'] if len(set(f))==8),None)
                if f is not None:
                    verify_witness(B,f)
                    witness={'depth':depth,'group_order':len(B['add']),'chain':B['chain'],
                             'add':B['add'],'clock':B['clock'],'surjective_character':f}
                    break
                following.append(B)
            if witness:break
        levels.append({'depth':depth,'tested_marked_presentations':candidates,
                       'complete_level':witness is None,'primitive_character_found':bool(witness)})
        if witness:break
        frontier=following
    need(witness is not None,'bounded explicit cubic completion witness')
    out={'status':'FINITE-EXACT complete smaller-depth exclusion and explicit attaining extension',
         'levels':levels,'witness':witness,'gates':GATES,
         'scope':'marked presentations, not isomorphism-class counts; C2 involution and C8 character'}
    (HERE/'involution_root_extensions.json').write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(json.dumps({'levels':levels,'witness_depth':witness['depth'],
                      'chain':witness['chain'],'gates':GATES}),flush=True)
    return out


if __name__=='__main__':run()
