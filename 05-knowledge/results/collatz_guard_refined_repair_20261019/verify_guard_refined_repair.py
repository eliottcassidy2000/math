#!/usr/bin/env python3
"""Exact source-aware offset certificates: two cylinders and fifteen trees.
All arithmetic relevant to certification is integral/rational. No assertions,
network calls, external imports, or source-orbit census are used.
"""
from fractions import Fraction as F
import hashlib
import json

checks=0
def require(value,label):
    global checks
    checks+=1
    if not value:raise RuntimeError(label)

def branch(word):
    a,b,M=1,0,1
    for bit in word:
        if bit=='1':a,b=3*a,3*b+M
        elif bit!='0':raise ValueError(bit)
        M*=2
    return a,b,M

def cylinder(word):
    a,b,M=branch(word)
    return (-b*pow(a,-1,M))%M

def step(x):return (3*x+1)//2 if x%2 else x//2

def rational(x):return f'{x.numerator}/{x.denominator}'

K='110110110110110101011011101010'
V=K*4+'10'
blocks={'K':K,'V':V,'F':'10','H':'110101000'}
C=F(83829391903560839836419813590128166431,
    157489106280555890461927692448373459)
E=C+1
cap25,cap1000=int(25*E),int(1000*E)
require((cap25,cap1000)==(13332,533286),'carry budgets')

short=[(f'K^{k}{name}',K*k+blocks[name]) for k in (1,2,3) for name in ('F','H')]
short+=[('K^4H',K*4+blocks['H'])]
short+=[(f'VK^{k}{name}',V+K*k+blocks[name]) for k in (0,1,2) for name in ('F','H')]
short_receipt=[]
for name,word in short:
    a,b,M=branch(word);r=cylinder(word)
    require(a<M,'short first-crossing multiplier')
    require(r>0 and (M-a)*r>b,'short all-lift endpoint descent')
    short_receipt.append({'name':name,'steps':len(word),'residue':r,
                          'fixed_point':rational(F(b,M-a))})
require(len(short)==13,'complete displayed short-branch list')

# The two long source-prefix cylinders in the symbolic partition.
sources=[('A',K*5),('B',V+K*3)]
source_data={name:(word,cylinder(word)) for name,word in sources}
require(len(source_data['A'][0])==150 and len(source_data['B'][0])==212,
        'long prefix lengths')

# Direct integer iteration selects the prefix; the separate closed affine
# branch and exact congruence below authenticate it on every lawful lift.
def first_contraction(r,d,L):
    n=r+d;a,b,M=1,0,1;word=''
    for j in range(1,L+1):
        bit=n%2;word+=str(bit)
        if bit:a,b=3*a,3*b+M
        M*=2;n=step(n)
        if a<M:return j,a,b,M,word
    return None

def certify(source_word,r,d,found):
    j,a,b,M,word=found
    require(j<=len(source_word),'known source bits suffice')
    require(branch(word)==(a,b,M),'separate branch reconstruction')
    require((r+d-cylinder(word))%M==0,'lawful continuation guard')
    require((1<<len(source_word))%M==0,'all-lift modulus containment')
    require(r>0 and (M-a)*r>a*d+b,'exact original-source inequality')
    return F(a*d+b,M-a)

# Cache by first-contraction word/residue. The same parity class is checked
# once, while its worst offset checks the full interval of offsets in it.
# Each d still contributes to an exact complete-coverage count and digest.
classes={};exceptions=[];digest=hashlib.sha256();max25=0;maxbase=0
for d in range(1,cap1000+1):
    failed=False
    for name,(word,r) in source_data.items():
        found=first_contraction(r,d,len(word))
        if found is None:
            failed=True;digest.update(f'{name},{d},NONE\n'.encode());continue
        j,a,b,M,shiftword=found
        require((M-a)*r>a*d+b,'each offset lies above its exact threshold')
        key=(name,j,d%M)
        if key not in classes:
            threshold=certify(word,r,d,found)
            classes[key]=[d,d,1,found]
        else:
            classes[key][1]=d;classes[key][2]+=1
        digest.update(f'{name},{d},{j},{a},{b}\n'.encode())
        maxbase=max(maxbase,j)
        if d<=cap25:max25=max(max25,j)
    if failed:exceptions.append(d)
require(not any(d<=cap25 for d in exceptions),'all 25-episode offsets covered')
require(max25==183 and maxbase==183,'two-cylinder repair depths')
expected=[14196,15364,18436,49324,128308,215276,266092,288932,295660,
          310916,390388,394148,480032,522404,529652]
require(exceptions==expected,'complete exceptional offset list')
max_threshold=F(0);max_threshold25=F(0)
for (name,j,residue),(lo,hi,count,found) in classes.items():
    word,r=source_data[name];jj,a,b,M,w=found
    require(hi%M==residue and lo%M==residue,'class endpoint residues')
    require(count==(hi-lo)//M+1,'no holes in an offset class')
    require((M-a)*r>a*hi+b,'worst class offset descends for every source lift')
    max_threshold=max(max_threshold,F(a*hi+b,M-a))
    if lo<=cap25:
        last=lo+((min(hi,cap25)-lo)//M)*M
        max_threshold25=max(max_threshold25,F(a*last+b,M-a))
require(max_threshold25==F(810692053755261919,1738366812781),
        '25-episode maximum threshold')
require(max_threshold25<466400,'25-episode simple threshold')

# Finite complete prefix forests for the fifteen exceptional offsets.
# A node records a macroblock word. Both roots and all four children of each
# nonterminal node are included, making coverage checkable from the receipt.
def build_tree(d):
    stack=['V','K'];nodes={}
    while stack:
        name=stack.pop();word=''.join(blocks[s] for s in name)
        a,b,M=branch(word);r=cylinder(word)
        require(name not in nodes,'unique tree node')
        if a<M:
            require(r>0 and (M-a)*r>b,'early first crossing descends on all lifts')
            nodes[name]={'kind':'source_descent','steps':len(word)}
            continue
        found=first_contraction(r,d,len(word))
        if found is not None:
            threshold=certify(word,r,d,found)
            nodes[name]={'kind':'offset_repair','steps':found[0],
                         'source_bits':len(word),'source_residue':r,
                         'threshold':rational(threshold)}
            continue
        require(len(word)<600,'finite refinement depth budget')
        nodes[name]={'kind':'branch'}
        for letter in reversed(tuple(blocks)):stack.append(name+letter)
    require('K' in nodes and 'V' in nodes,'both source roots present')
    for name,row in nodes.items():
        if row['kind']=='branch':
            require(all(name+letter in nodes for letter in blocks),'all four children present')
        else:
            require(not any(name+letter in nodes for letter in blocks),'terminal leaf')
        if len(name)>1:require(nodes[name[:-1]]['kind']=='branch','connected complete forest')
    return nodes

trees={str(d):build_tree(d) for d in exceptions}
max_tree=max(row['steps'] for nodes in trees.values() for row in nodes.values()
             if row['kind']=='offset_repair')
require(max_tree==273,'refined worst continuation time')
max_source_bits=max(row['source_bits'] for nodes in trees.values() for row in nodes.values()
                    if row['kind']=='offset_repair')
print(json.dumps({'status':'PASS','explicit_checks':checks,
 'scope':'all-lift offset certificate for the specified first-crossing bank language; no arbitrary source coverage',
 'main_episode_limit':25,'main_offset_cap':cap25,'main_repair_steps':max25,
 'main_max_threshold':rational(max_threshold25),
 'computational_episode_limit':1000,'computational_offset_cap':cap1000,
 'computational_repair_steps':max(maxbase,max_tree),
 'base_class_count':len(classes),'base_row_digest':digest.hexdigest(),
 'base_max_threshold':rational(max_threshold),
 'long_source_cylinders':{name:{'length':len(w),'residue':r,'modulus':1<<len(w)}
                         for name,(w,r) in source_data.items()},
 'short_descending_branches':short_receipt,
 'exceptional_offsets':exceptions,
 'refinement_node_count':sum(map(len,trees.values())),
 'refinement_max_source_bits':max_source_bits,
 'refinement_trees':trees},indent=2,sort_keys=True))
