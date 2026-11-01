"""Transition-table-only simulator; no numeric witness predicate is imported."""
import json
from pathlib import Path
from build_tm import delta

HERE=Path(__file__).resolve().parent

def encode(n,word):
    return '^'+'1'*n+'#'+'#'.join('1'*v for v in word)+'$'

def run(text,mode='path',limit=10000000,machine=None,trace=False):
    machine=machine or json.loads((HERE/(mode+'.tm.json')).read_text())
    if any(c not in machine['alphabets'][0] or c=='_' for c in text):
        raise ValueError('input must use the nonblank input alphabet')
    tapes=[dict(enumerate(text))]+[{} for _ in range(5)]
    heads=[0]*6; state=machine['start']; steps=0; visited={state}; lookup={}
    while state not in ['ACCEPT','REJECT']:
        if steps>=limit: raise RuntimeError('simulation step limit')
        letters=[t.get(p,'_') for t,p in zip(tapes,heads)]
        key=(state,tuple(letters))
        if key not in lookup: lookup[key]=delta(machine,state,letters)
        nxt,write,moves=lookup[key]
        for i in range(6):
            if write[i]=='_': tapes[i].pop(heads[i],None)
            else: tapes[i][heads[i]]=write[i]
            heads[i]+={'L':-1,'R':1,'S':0}[moves[i]]
        state=nxt; steps+=1; visited.add(state)
    out={'state':state,'steps':steps}
    if trace: out['states']=sorted(visited)
    return out

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(); p.add_argument('mode',choices=['path','cycle']); p.add_argument('n',type=int); p.add_argument('word',nargs='*',type=int)
    a=p.parse_args(); print(json.dumps(run(encode(a.n,a.word),a.mode),sort_keys=True))
