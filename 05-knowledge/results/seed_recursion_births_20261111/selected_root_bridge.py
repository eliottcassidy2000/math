"""Independent full-coefficient check of two incoming selected-root ladders."""
from pathlib import Path
import hashlib,json,sys
import cones

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'endpoint_coverage_compiler_20261111'))
import verify as iv


def main():
    original=[2,24,3,5,13];base=[24,1,3,6,14];p=[1,0,2,3,4]
    iv.require(all(bool(original[i]>>j&1)==bool(base[p[i]]>>p[j]&1) for i in range(5) for j in range(i)),'literal marked isomorphism')
    records=[]
    for mask,selected,port,slope in [(5,[1,2],1,[89,90,30]),(24,[3,4],3,[147,150,50])]:
        child,_,labels=cones.expand(base,[],[port],[2])
        marks=[labels.index(v) for v in selected]
        fresh=[i for i,v in enumerate(labels) if v==port][1]
        delta_values=[];u_values=[]
        for a in range(1,len(child)+1):
            expanded,es,_=cones.expand(child,marks,[fresh],[a])
            h,w,u=iv.endpoint_state(expanded,es)
            delta_values.append(h-w);u_values.append(u)
        d=iv.invert(delta_values);u=iv.invert(u_values)
        iv.require(d==[1]+[0]*5 and u==slope+[0]*3,'full seed and slope certificate')
        h,w,first_u=iv.endpoint_state(base,selected)
        iv.require(h==11 and h-w==1 and first_u==(29 if mask==5 else 47),'initial deletion root')
        records.append(dict(original_rows=original,original_selected_mask=mask,isomorphism=p,
                            incoming_rows=base,incoming_selected=selected,selected_port=port,
                            initial_U=first_u,child_rows=child,child_selected=marks,fresh_stable_port=fresh,
                            defect_coefficients=d,slope_coefficients=u,
                            unit_seeds_physically_recounted=len(child)+1))
    digest=lambda path:hashlib.sha256(path.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
    report=dict(credit='Session260 checkpoint 4315a74fe0, constructors/RESEARCH.md; independently recovered coefficients here',
                source_sha256_lf=digest(Path(__file__)),shared_input_manifest_sha256_lf=digest(HERE/'SOURCE_PINS.json'),records=records)
    target=HERE/'selected_root_bridge.json'
    if target.exists():iv.require(json.loads(target.read_text())==report,'selected-root replay unchanged')
    target.write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    print(json.dumps(dict(roots=[r['initial_U'] for r in records],
                         stable_children=[r['slope_coefficients'] for r in records],physical_states=14),indent=2))


if __name__=='__main__':main()
