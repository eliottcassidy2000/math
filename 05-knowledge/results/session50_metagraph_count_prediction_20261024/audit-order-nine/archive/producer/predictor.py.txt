#!/usr/bin/env python3
"""Order-nine prediction from complete parent actions and a small regular sector.
Inputs: authenticated ordinary parents through eight. The bulk pass only
extends and refines; the exceptional regular sector uses parent-eight reroot
identification. No target class count, target representative atlas, or target
canonicalizer is read or implemented.
"""
import argparse, hashlib, itertools, json, sys
from collections import Counter
from pathlib import Path
CHECKS=Counter();OPS=Counter()
ATLAS_HASH='1cd97666c0ee6294bf379b6d77f95bea5e24f8f5b27749ef804c46a45892ac7b'
PREFIX_HASH='df16881db1060e0423ed73658bdcb5923266f839681b74ac5b1bddd5a232f3ad'
def require(c,m):
    CHECKS[m]+=1
    if not c:raise RuntimeError(m)
def partition(rows):
    n=len(rows);cells={}
    for v,r in enumerate(rows):cells.setdefault(r.bit_count(),[]).append(v)
    blocks=[tuple(cells[d]) for d in sorted(cells)]
    while True:
        masks=[sum(1<<v for v in block) for block in blocks];new=[]
        for block in blocks:
            groups={}
            for v in block:
                signature=tuple((rows[v]&mask).bit_count() for mask in masks)
                groups.setdefault(signature,[]).append(v)
            new.extend(tuple(groups[s]) for s in sorted(groups))
        if len(new)==len(blocks):return tuple(new)
        blocks=new

def root_cell(rows):
    # This independent bulk evaluator refines explicit vertex blocks and may
    # stop when the distinguished root is alone. C++ uses full color vectors.
    OPS['child_refinement_calls']+=1
    cells={}
    for v,r in enumerate(rows):cells.setdefault(r.bit_count(),[]).append(v)
    blocks=list(cells.values());root=len(rows)-1
    if len(cells[rows[root].bit_count()])==1:return 1
    while True:
        OPS['refinement_rounds']+=1
        masks=[sum(1<<v for v in block) for block in blocks];new=[];size=None
        for block in blocks:
            if len(block)==1:new.append(block);continue
            groups={}
            for v in block:
                signature=tuple((rows[v]&mask).bit_count() for mask in masks)
                groups.setdefault(signature,[]).append(v)
            for group in groups.values():
                if root in group:size=len(group)
                new.append(group)
        if size==1:return 1
        if len(new)==len(blocks):
            require(size is not None,'stable distinguished root block found')
            return size
        blocks=new

def canonical_parent(rows):
    require(len(rows)==8,'only order-eight parent canonicalization');OPS['parent_canonicalizations']+=1
    best=None
    for pieces in itertools.product(*(itertools.permutations(block) for block in partition(rows))):
        order=sum(pieces,());candidate=tuple(sum(((rows[v]>>w)&1)<<j for j,w in enumerate(order)) for v in order)
        OPS['parent_canonical_label_candidates']+=1
        if best is None or candidate<best:best=candidate
    return best

def regular_sector(records):
    parents=[tuple(record['rows']) for record in records if set(v.bit_count() for v in record['rows'])=={3,4}]
    parentset=set(parents);trace=0
    for rows in parents:
        child=list(rows)+[0]
        for v in range(8):
            if rows[v].bit_count()==3:child[v]|=1<<8
            else:child[8]|=1<<v
        require(all(row.bit_count()==4 for row in child),'regular extension forced by degree');OPS['regular_child_extensions']+=1
        for root in range(9):
            order=[v for v in range(9) if v!=root]
            deleted=tuple(sum(((child[v]>>w)&1)<<j for j,w in enumerate(order)) for v in order);OPS['regular_child_deletions']+=1
            key=canonical_parent(deleted);require(key in parentset,'regular reroot stays in parent domain');trace+=key==rows
    require(trace%9==0,'regular trace multiple of nine')
    classes=trace//9
    return {'pointed_states':len(parents),'trace':trace,'classes':classes,'error':9*(len(parents)-classes)}

def known_error(record):
    rows=record['rows'];actions=record['actions'];blocks=partition(rows);cellsize={v:len(b) for b in blocks for v in b};remaining=set(range(len(rows)));error=0
    while remaining:
        v=min(remaining);orbit={g[v] for g in actions};require(orbit<=remaining,'known parent action orbit partition');remaining-=orbit
        error+=cellsize[v]-len(orbit)
    return error

def predict(atlas,prefix):
    require(atlas['maximum_parent_order']==8,'parent atlas order eight');require(prefix['maximum_parent_order']==7,'known prefix order seven')
    records=atlas['parents'];levels=prefix['levels'];counts=[len(level['parents']) for level in levels]
    lower_errors=[sum(known_error(r) for r in level['parents']) for level in levels]
    require(all(e==0 for e in lower_errors[:7]),'computed lower-order exactness')
    seven_error=lower_errors[7];regular=regular_sector(records)
    # Structural theorem: the only nonregular strong-nine contributions are
    # one C3 substitution for each regular-seven class and the rigid part of
    # the three-triangle transitive-majority matching family.
    triangle_cells=3;edges=triangle_cells*(triangle_cells-1)//2
    rigid_triangle_classes=(2**edges-2**(triangle_cells-1))*3**(edges-triangle_cells)
    triangle_error=2*triangle_cells*triangle_cells*rigid_triangle_classes
    padding_coefficient=sum(counts[i]*counts[2-i] for i in range(3))
    correction=regular['error']+seven_error+triangle_error+padding_coefficient*seven_error
    states=upper=0;hist=Counter();perparent=[]
    for record in records:
        rows=tuple(record['rows']);actions=record['actions'];seen=set();localhist=Counter();localstates=localupper=0
        for mask in range(1<<len(rows)):
            if mask in seen:continue
            orbit={sum(1<<g[v] for v in range(len(rows)) if mask>>v&1) for g in actions}
            require(mask==min(orbit) and not orbit&seen,'attachment orbit unique representative');seen.update(orbit)
            child=tuple(row|(((mask>>v)&1)<<len(rows)) for v,row in enumerate(rows))+(((1<<len(rows))-1)^mask,)
            OPS['bulk_child_extensions']+=1
            cell=root_cell(child);require(cell%2==1,'odd stable root cell');hist[cell]+=1;localhist[cell]+=1
            states+=1;upper+=cell;localstates+=1;localupper+=cell
        require(len(seen)==1<<len(rows),'all attachment subsets covered')
        perparent.append([localstates,localupper]+[localhist[k] for k in range(1,10,2)])
    numerator=upper-correction;require(numerator>0 and numerator%9==0,'corrected trace positive multiple of nine')
    return {'status':'PASS','target_order':9,'input_parent_classes':len(records),'pointed_extension_states':states,'stable_root_cell_histogram':sorted(hist.items()),'upper_refinement_trace':upper,'known_prefix_counts':counts,'known_prefix_errors':lower_errors,'regular_sector':regular,'nonregular_strong_nine_seven_cell_error':seven_error,'nonregular_strong_nine_triangle_family_classes':rigid_triangle_classes,'nonregular_strong_nine_triangle_error':triangle_error,'lower_scc_padding_error':padding_coefficient*seven_error,'total_error_correction':correction,'exact_unpointing_trace':numerator,'predicted_class_count':numerator//9,'child_canonicalizations':0,'general_child_reroots':0,'operations':dict(OPS),'checks':dict(CHECKS),'per_parent':perparent}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',default='normal.json');args=parser.parse_args();here=Path(__file__).resolve().parent
    source=Path(__file__).resolve();input1=here/'parents-eight.json';input2=here/'archive/parents-through-seven.json.txt';out=(here/args.output).resolve()
    allowed={str(source),str(input1),str(input2)}
    def audit(event,arguments):
        if event=='open':
            target,mode,flags=arguments;require(isinstance(target,(str,bytes)),'no descriptor-based input')
            path=str(Path(target).resolve());writing=mode is not None and any(c in mode for c in 'wax+')
            require(path==str(out) if writing else path in allowed,'prediction I/O boundary')
        require(not event.startswith(('subprocess.','socket.','ctypes.')),'no external prediction operations')
    sys.addaudithook(audit)
    raw1=input1.read_bytes();raw2=input2.read_bytes()
    require(hashlib.sha256(raw1).hexdigest()==ATLAS_HASH,'pinned complete order-eight atlas');require(hashlib.sha256(raw2).hexdigest()==PREFIX_HASH,'pinned complete lower prefix')
    result=predict(json.loads(raw1),json.loads(raw2));result['source_sha256']=hashlib.sha256(source.read_bytes()).hexdigest()
    out.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n');print(json.dumps({k:result[k] for k in ('status','pointed_extension_states','upper_refinement_trace','total_error_correction','predicted_class_count','regular_sector')},sort_keys=True))
if __name__=='__main__':main()
