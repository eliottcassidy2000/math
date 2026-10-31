#!/usr/bin/env python3
"""Symbolic per-offset colors and independently readable child-member tails."""
import json
from pathlib import Path
HERE=Path(__file__).resolve().parent
TABLE=HERE.parent/'session124_recursive_square_prefix_20261031'/'paired'


def need(x,msg):
    if not x:raise RuntimeError(msg)


def profile(word,m):
    at=0;heads={};tails={};singletons={};readable={}
    for atom in word:
        kind,c=atom[:2]
        if kind=='single':singletons[c]=at%2;at+=1
        elif kind=='head':heads[c]=at%2;at+=1
        else:
            _,c,i,o=atom;length=m+i;start=int(kind=='tail')
            shift=(at-start if o==1 else at+length-1)%2
            need(c not in tails,'unique tail')
            tails[c]=shift;readable[c]=(at,length-start,i,o,start)
            if not start:heads[c]=shift
            at+=length-start
    return heads,tails,singletons,readable,at


states=channels=terminal_states=0
for parity in (0,1):
    for d in range(24,73):
        m=42+parity
        record=json.loads((TABLE/f'B_p{parity}_d{d}.json').read_text())
        for word in record['words']:
            need(word[:4]==[['single',x] for x in (1,3,6,19)],'fixed C head')
            need(word[4][0:2]==['block',4] and word[4][3]==1,'intact oriented prefix channel')
        a,b=[profile(w,m) for w in record['words']]
        need(a[:3]==b[:3],'per-offset colors agree across output members')
        need(a[4]==49*m+d and b[4]==49*m+d+1,'member sizes')
        for p in (a,b):
            need(set(p[0])==set(range(-24,25))==set(p[1]),'all heads and tails')
            need(set(p[2])==set(range(1,25)),'low interval')
            need(all(c in p[3] for c in range(25)),'every absolute channel is readable in each member')
            channels+=25
        states+=1
        terminal=json.loads((TABLE.parent/'cycles'/f'cycle_p{parity}_d{d}.json').read_text())
        word=terminal['word'];p=profile(word,m)
        need(word[:4]==[['single',x] for x in (1,3,6,19)] and word[4][0:2]==['block',4] and word[4][3]==1,'terminal C prefix channel')
        need(p[4]==49*m+d,'cycle length')
        need(word[:2]==[['single',1],['single',3]],'cycle orientation anchor')
        need(set(p[0])==set(range(-24,25))==set(p[1]),'cycle source ownership')
        need(set(p[2])==set(range(1,25)),'cycle low singletons')
        need(all(c in p[3] for c in range(25)),'all terminal child channels')
        terminal_states+=1
print(json.dumps({'status':'passed','complete_states':states,'readable_member_channels':channels,
                  'terminal_states':terminal_states,
                  'scope':'All98 Q49 weak-B pair states; symbolic parity and tail locations, without assuming equal colors in unrelated channels.'},indent=2,sort_keys=True))
