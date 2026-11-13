"""Literal eight-vertex escapes and an all-length primitive rooted return.

Only the named sources and constructions are counted; this is not a spectrum census.
"""
from pathlib import Path
from math import factorial, gcd
from fractions import Fraction
import argparse
import json

HERE = Path(__file__).resolve().parent
C = (56,117,129,22,100,76,141,59)
P = (1,0,4,5,6,2,7,3)

def require(ok, message):
    if not ok:
        raise RuntimeError(message)

def converse(t):
    full=(1<<len(t))-1
    return tuple(full^(1<<i)^row for i,row in enumerate(t))

def flip(t,i,j):
    r=list(t);r[i]^=1<<j;r[j]^=1<<i
    return tuple(r)

def native(t,path):
    require(sorted(path)==list(range(len(t))),'not a spanning path')
    require(all(t[i]>>j&1 for i,j in zip(path,path[1:])),'spine is not directed')
    r=list(converse(t))
    for i,j in zip(path,path[1:]):
        r[i]^=1<<j;r[j]^=1<<i
    return tuple(r)

def ends(t):
    n=len(t);full=(1<<n)-1
    f=[[0]*n for _ in range(full+1)]
    for i in range(n):
        f[1<<i][i]=1
    for mask in range(1,full+1):
        for i,count in enumerate(f[mask]):
            if not count:
                continue
            todo=t[i]&~mask
            while todo:
                bit=todo&-todo;todo-=bit
                f[mask|bit][bit.bit_length()-1]+=count
    return f

def covers(t):
    n=len(t);full=(1<<n)-1
    h=[sum(row) for row in ends(t)]
    f=[[0]*(full+1) for _ in range(n+1)];f[0][0]=1
    for k in range(1,n+1):
        for mask in range(1,full+1):
            sub=mask
            while sub:
                f[k][mask]+=h[sub]*f[k-1][mask^sub]
                sub=(sub-1)&mask
    return [f[k][full] for k in range(1,n+1)]

def slot_coefficients(profile):
    n=len(profile);out=[Fraction(0)]*(n+1);falling=[1]
    for k,count in enumerate(profile,1):
        new=[0]*(len(falling)+1)
        for j,v in enumerate(falling):
            new[j]-=(k-1)*v;new[j+1]+=v
        falling=new
        for j,v in enumerate(falling):
            out[j]+=Fraction(count*v,factorial(k))
    require(all(c.denominator==1 for c in out),'nonintegral slot coefficient')
    return [int(c) for c in out]

def summary(t):
    f=covers(t)
    return {'rows':t,'H':f[0],'ordered_covers':f,
            'raw_content':gcd(*f),
            'uniform_context_content':gcd(*(factorial(k)*x for k,x in enumerate(f))),
            'slot_coefficients':slot_coefficients(f)}

def root_state(t,root):
    full=(1<<len(t))-1;e=ends(t);s=ends(converse(t))
    h=[sum(row) for row in e];h[0]=1
    return (h[full],s[full][root],h[full^(1<<root)],
            sum(e[mask][root]*h[full^mask] for mask in range(1,full+1)))

def child(t,root,new_root=False):
    """Add x with root -> x -> all other old vertices."""
    n=len(t);r=list(t)+[((1<<n)-1)^(1<<root)]
    r[root]|=1<<n
    return tuple(r),n if new_root else root

def realize(k):
    require(k>=1,'k must be positive')
    t,root=converse(C),5
    for _ in range(k):
        t,root=child(t,root)
    t,root=child(t,root,True)
    t,root=child(t,root)
    return t,root

def transitive(n):
    return tuple(sum(1<<j for j in range(i+1,n)) for i in range(n))

def context(module,exterior,beats_module):
    n=len(module);m=len(exterior)
    r=list(module)+[row<<n for row in exterior]
    for j in range(m):
        if beats_module[j]:
            r[n+j]|=(1<<n)-1
        else:
            for i in range(n):
                r[i]|=1<<(n+j)
    return tuple(r)

def strong(t):
    reach=[row|(1<<i) for i,row in enumerate(t)]
    for k in range(len(t)):
        for i in range(len(t)):
            if reach[i]>>k&1:
                reach[i]|=reach[k]
    return all(row==(1<<len(t))-1 for row in reach)

def controls():
    sources={'C':C,'C2':flip(C,1,6),'C4':flip(C,2,3),
             'C24':flip(flip(C,1,6),2,3),'native':native(C,P)}
    out={name:summary(t) for name,t in sources.items()}
    require([out[n]['H'] for n in sources]==[513,515,517,495,513],'wrong source counts')
    require([out[n]['uniform_context_content'] for n in sources]==[27,5,1,9,3],
            'wrong context contents')
    require(all(strong(t) for t in sources.values()),'a source is not strong')
    for name,amplitude,degree in [('C2',2,4),('C4',4,2)]:
        delta=[a-b for a,b in zip(out[name]['slot_coefficients'],out['C']['slot_coefficients'])]
        expected=[0]*9;expected[degree]=amplitude
        require(delta==expected,'wrong full polynomial response')
    mixed=[out['C24']['slot_coefficients'][i]-out['C2']['slot_coefficients'][i]
           -out['C4']['slot_coefficients'][i]+out['C']['slot_coefficients'][i] for i in range(9)]
    require(mixed==[0,0,-16,0,-8,0,0,0,0],'wrong mixed response')
    # Named exterior root: all four coordinates obey the same module relation.
    response=[];selectors=[]
    for exterior,polarity in [((2,0),(False,True)),
                              (transitive(4),(False,True,False,True))]:
        old=context(C,exterior,polarity);before=root_state(old,8)
        for name,amplitude,degree in [('C2',2,4),('C4',4,2)]:
            after=root_state(context(sources[name],exterior,polarity),8)
            small=root_state(context(transitive(degree),exterior,polarity),degree)
            require(tuple(y-x for x,y in zip(before,after))==tuple(amplitude*x for x in small),
                    'external-root response relation failed')
            response.append({'exterior':exterior,'beats_module':polarity,'replacement':name,
                             'before':before,'after':after,'small_context':small,
                             'amplitude':amplitude})
        small4=root_state(context(transitive(4),exterior,polarity),4)
        small2=root_state(context(transitive(2),exterior,polarity),2)
        both=root_state(context(sources['C24'],exterior,polarity),8)
        require(both==tuple(h-6*a-12*b for h,a,b in zip(before,small4,small2)),
                'mixed external-root response failed')
        quartet=[before[0],before[0]+2*small4[0],before[0]+4*small2[0],both[0]]
        require(sorted(x%8 for x in quartet)==[1,3,5,7],'not a mod-eight selector')
        require(all(495*before[0]<=513*x<=517*before[0] for x in quartet),
                'selector escaped its universal relative band')
        selectors.append({'exterior':exterior,'beats_module':polarity,
                          'h_A_B':[before[0],small4[0],small2[0]],'four_counts':quartet,
                          'residues_mod8':[x%8 for x in quartet]})
    entry=root_state(converse(C),5)
    require(entry==(513,77,129,711),'wrong rooted entry')
    returned=[]
    for k in (1,2,3):
        t,root=realize(k);z=root_state(t,root)
        b=1746*(1<<k)+241;a=582*(1<<k)-17
        require(z==(b,a,a,b),'wrong returned state')
        require(gcd(b,a)==1 and strong(t),'nonprimitive or nonstrong return')
        returned.append({'k':k,'rows':t,'root':root,'state_h_s_q_F':z})
    powers=[pow(2,k,73) for k in range(9)]
    require(pow(2,9,73)==1 and len(set(powers))==9 and 28 not in powers,
            'prime-clock certificate failed')
    require(pow(28,9,73)==10 and pow(10,4,73)==72 and pow(28,8,73)==16,
            'eight-nine character certificate failed')
    receipt={'status':'Exact controls for the all-length and all-context proofs in PROOF.md',
             'sources':out,'native_path':P,'rooted_entry':entry,
             'external_root_controls':response,'selector_controls':selectors,
             'literal_returns':returned,
             'powers_of_two_mod73':powers,
             'eight_nine_character':{'modulus':73,'kernel_order':9,'image_order':8,
                                     'forbidden_residue':28,'character_value':10}}
    (HERE/'escape_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':'pass','source_counts':[out[n]['H'] for n in sources],
                      'context_contents':[out[n]['uniform_context_content'] for n in sources],
                      'return_counts':[r['state_h_s_q_F'][0] for r in returned]},indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--k',type=int)
    args=parser.parse_args()
    if args.k is None:
        controls()
    else:
        t,root=realize(args.k)
        print(json.dumps({'k':args.k,'rows':t,'root':root,
                          'H':1746*(1<<args.k)+241,
                          'state_h_s_q_F':[1746*(1<<args.k)+241,582*(1<<args.k)-17,
                                           582*(1<<args.k)-17,1746*(1<<args.k)+241]},indent=2))
