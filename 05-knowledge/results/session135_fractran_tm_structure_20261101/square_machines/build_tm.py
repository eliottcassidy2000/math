"""Build literal six-tape deterministic Turing machines, with no arithmetic oracle.

A wildcard row denotes finitely many ordinary transition rows. The emitted gzip
TSV expands every alphabet tuple, including explicit rejecting default rows.
"""
import gzip
import hashlib
import itertools
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
ALPHABETS=[['_','^','1','#','$'],['_','^','U','X','$']]+[['_','^','1'] for _ in range(4)]

class Builder:
    def __init__(self, mode):
        self.mode=mode
        self.rows={}
    def row(self,q,read,nxt,write=None,move=None):
        r=['*']*6; w=['=']*6; m=['S']*6
        for i,v in read.items(): r[i]=v
        for i,v in (write or {}).items(): w[i]=v
        for i,v in (move or {}).items(): m[i]=v
        self.rows.setdefault(q,[]).append([r,nxt,w,m])
    def rewind(self,name,tape,nxt,advance=False):
        for c in ALPHABETS[tape]:
            if c!='^': self.row(name,{tape:c},name,move={tape:'L'})
        self.row(name,{tape:'^'},nxt,move={tape:'R'} if advance else {})
    def transfer(self,name,src,dst,nxt):
        self.rewind(name,src,name+'_copy',True)
        self.row(name+'_copy',{src:'1',dst:'_'},name+'_copy',
                 {src:'_',dst:'1'},{src:'R',dst:'R'})
        self.row(name+'_copy',{src:'_'},name+'_reset')
        self.rewind(name+'_reset',src,nxt,True)
    def square(self,name,nxt):
        self.row(name,{4:'_',5:'_'},name+'_subtract',{5:'1'},{4:'L'})
        self.row(name+'_subtract',{4:'1',5:'1'},name+'_subtract',
                 {4:'_'},{4:'L',5:'R'})
        self.row(name+'_subtract',{4:'1',5:'_'},name+'_grow2',
                 {5:'1'},{5:'R'})
        self.row(name+'_grow2',{5:'_'},name+'_rewind',
                 {5:'1'},{5:'R'})
        self.rewind(name+'_rewind',5,name+'_subtract',True)
        self.row(name+'_subtract',{4:'^',5:'_'},name+'_clear',move={4:'R'})
        self.row(name+'_clear',{5:'_'},name+'_clear',move={5:'L'})
        self.row(name+'_clear',{5:'1'},name+'_clear',{5:'_'},{5:'L'})
        self.row(name+'_clear',{5:'^'},nxt,move={5:'R'})
    def build(self):
        self.row('start',{0:'^'},'n0',
                 {i:'^' for i in range(1,6)},{i:'R' for i in range(6)})
        for j in range(4):
            self.row('n'+str(j),{0:'1',1:'_'},'n'+str(min(j+1,3)),
                     {1:'U'},{0:'R',1:'R'})
            if j and (self.mode=='path' or j==3):
                self.row('n'+str(j),{0:'#',1:'_'},'initial_rewind',
                         {1:'$'},{0:'R'})
        self.rewind('initial_rewind',1,'first_begin')
        self.row('first_begin',{0:'1'},'first_read')
        self.row('first_read',{0:'1',2:'_',3:'_'},'first_bound',
                 {2:'1',3:'1'},{0:'R',1:'R',2:'R',3:'R'})
        for marker in ['U','X']:
            self.row('first_bound',{0:'1',1:marker},'first_read')
        for delim in ['#','$']:
            self.row('first_bound',{0:delim,1:'U'},'boundary',{1:'X'})
        self.row('boundary',{0:'#'},'next_begin',move={0:'R'})
        self.row('next_begin',{0:'1'},'visited_rewind')
        self.rewind('visited_rewind',1,'previous_transfer')
        self.transfer('previous_transfer',2,4,'current_read')
        self.row('current_read',{0:'1',2:'_',4:'_'},'current_bound',
                 {2:'1',4:'1'},{0:'R',1:'R',2:'R',4:'R'})
        for marker in ['U','X']:
            self.row('current_bound',{0:'1',1:marker},'current_read')
        for delim in ['#','$']:
            self.row('current_bound',{0:delim,1:'U'},'edge_square',{1:'X'})
        self.square('edge_square','boundary')
        self.row('boundary',{0:'$'},'end_check',move={0:'R'})
        if self.mode=='cycle':
            self.row('end_check',{0:'_'},'close_previous')
            self.transfer('close_previous',2,4,'close_first')
            self.transfer('close_first',3,4,'closing_square')
            self.square('closing_square','final_rewind')
        else:
            self.row('end_check',{0:'_'},'final_rewind')
        self.rewind('final_rewind',1,'final_scan',True)
        self.row('final_scan',{1:'X'},'final_scan',move={1:'R'})
        self.row('final_scan',{1:'$'},'ACCEPT')
        machine={'kind':'deterministic-six-tape-Turing-machine',
                 'mode':self.mode,'start':'start','accept':'ACCEPT','reject':'REJECT',
                 'blank':'_','alphabets':ALPHABETS,
                 'default':'REJECT with unchanged symbols and stationary heads',
                 'tapes':['input','named_vertex_marks','previous_vertex',
                          'first_vertex','edge_sum','odd_subtrahend'],
                 'rules':self.rows}
        validate(machine)
        return machine

def delta(machine,state,letters):
    matches=[]
    for read,nxt,write,move in machine['rules'][state]:
        if all(a=='*' or a==b for a,b in zip(read,letters)):
            matches.append((nxt,tuple(b if a=='=' else a for a,b in zip(write,letters)),tuple(move)))
    if len(matches)>1: raise ValueError('overlapping transitions')
    return matches[0] if matches else ('REJECT',tuple(letters),tuple('S'*6))

def validate(machine):
    states=set(machine['rules'])|{'ACCEPT','REJECT'}
    for state,rows in machine['rules'].items():
        for read,nxt,write,move in rows:
            if nxt not in states: raise ValueError('unknown target '+nxt)
            for i,(r,w,m) in enumerate(zip(read,write,move)):
                if r not in machine['alphabets'][i]+['*'] or w not in machine['alphabets'][i]+['='] or m not in 'LRS':
                    raise ValueError('invalid symbol')
        for a,b in itertools.combinations(rows,2):
            if all(x==y or x=='*' or y=='*' for x,y in zip(a[0],b[0])):
                raise ValueError('overlap at '+state)
    return True

def expand(machine,destination):
    """Materialize one ordinary delta table from the published finite rules."""
    concrete=0; h=hashlib.sha256()
    with destination.open('wb') as raw:
        with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0) as f:
            for state in sorted(machine['rules']):
                for letters in itertools.product(*machine['alphabets']):
                    nxt,write,move=delta(machine,state,letters)
                    line=(state+'\t'+''.join(letters)+'\t'+nxt+'\t'+''.join(write)+'\t'+''.join(move)+'\n').encode()
                    f.write(line); h.update(line); concrete+=1
    return {'concrete_delta_rows':concrete,'expanded_sha256':h.hexdigest(),
            'gzip_bytes':destination.stat().st_size}

def emit(mode):
    machine=Builder(mode).build()
    path=HERE/(mode+'.tm.json')
    path.write_text(json.dumps(machine,sort_keys=True,indent=2)+'\n')
    # Expanded delta is an ordinary finite transition table, not a macro VM.
    data=expand(machine,HERE/(mode+'.tm.tsv.gz'))
    return {'mode':mode,'nonhalting_states':len(machine['rules']),
            'symbolic_rows':sum(map(len,machine['rules'].values())),
            'json_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),**data}

if __name__=='__main__':
    results=[emit(mode) for mode in ['path','cycle']]
    (HERE/'TM_BUILD.json').write_text(json.dumps(results,indent=2,sort_keys=True)+'\n')
    print(json.dumps(results,sort_keys=True))
