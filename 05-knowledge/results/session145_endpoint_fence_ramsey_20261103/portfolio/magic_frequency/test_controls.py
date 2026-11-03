#!/usr/bin/env python3
"""Rejection tests and mathematical hostiles; all checks survive python -O."""
from fractions import Fraction as Q
import argparse, copy, json
from pathlib import Path
import verify as v

def run():
    rejected=[]
    good=v.verify()
    changes={
      'one_grid_entry_changed':lambda d:d['target'].__setitem__(0,109),
      'frequency_center_mistated':lambda d:d.__setitem__('frequency_mean',528),
      'trace_and_order_maps_confused':lambda d:d.__setitem__('trace_map',d['order_map']),
      'raw_labels_declared_magic':lambda d:d.__setitem__('frequency_magic_arrangements',8),
      'carry_defect_erased':lambda d:d['near_line_sums'].__setitem__(7,1665),
      'magic_function_dimension_wrong':lambda d:d.__setitem__('magic_function_dimension',4),
      'degree_seven_mode_missed':lambda d:d['polynomial_constraint_ranks_degree_0_to_8'].__setitem__(7,6),
      'projection_error_wrong':lambda d:d.__setitem__('projection_squared_error',0),
      'rotation_declared_magic':lambda d:d['order_map_rotation_line_sums'].__setitem__(7,378),
    }
    for name,mutate in changes.items():
        d=copy.deepcopy(good);mutate(d)
        try:v.validate_bundle(d)
        except ValueError:rejected.append(name)
        else:raise ValueError('hostile accepted: '+name)
    checks={
      'raw_frequency_rank_grid_not_magic':not v.magic(tuple(v.FREQUENCIES[a-1] for a in v.RANKS)),
      'same_placement_traces_not_magic':not v.magic(v.T),
      'squared_target_not_magic':not v.magic(tuple(x*x for x in v.TARGET)),
      'literal_cubic_trace_not_magic':not v.magic(tuple(t**3 for t in v.T)),
      'symmetric_cube_trace_not_magic':not v.magic(tuple(t**3-10*t for t in v.T)),
      'U_only_loses_injectivity':len(set(v.U))==3 and v.magic(v.U),
      'equal_weights_collide':len({u+w for u,w in zip(v.U,v.V)})<9,
      'order_boundary_collides':len({u+2*w for u,w in zip(v.U,v.V)})<9,
      'fiber_preserving_18_mutation_breaks_magic':not v.magic((v.TARGET[0]+18,)+v.TARGET[1:]),
      'center_shift_preserves_magic':v.magic(tuple(x+1 for x in v.TARGET)),
      'nonconstant_affine_raw_not_magic':not v.magic(tuple(7*v.FREQUENCIES[a-1]-13 for a in v.RANKS)),
    }
    for name,ok in checks.items():v.require(ok,name)
    # Exhaustively test the injectivity and monotonicity criteria in a bounded integer box.
    for p in range(-12,13):
      for q in range(-12,13):
        vals=[v.order_map(a,c=0,p=p,q=q) for a in range(1,10)]
        criterion=p!=0 and q!=0 and Q(p,q) not in (Q(-2),Q(-1),Q(-1,2),Q(1,2),Q(1),Q(2))
        v.require((len(set(vals))==9)==criterion,'injectivity criterion')
        v.require(all(x<y for x,y in zip(vals,vals[1:]))==(p>0 and q>2*p),'monotonicity criterion')
    return {'status':'PASS','rejected_mutations':rejected,'mathematical_hostiles':checks,'integer_parameter_pairs_tested':625}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output');args=p.parse_args()
    out=json.dumps(run(),indent=2,sort_keys=True)+'\n'
    if args.output:Path(args.output).write_text(out)
    print(out,end='')
