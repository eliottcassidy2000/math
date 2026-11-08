"""Replay the six published exact physical-oracle smoke tests."""
import json
from verify import find_state, graph, state

cases = [((1, 1, 1, 1), True), ((1, 0, 2, 1), True),
         ((2, 2, 2, 2), False), ((1, 1, 0, 1), False),
         ((3, 1, 3, 1), True), ((1, 0, 1, 1), False)]
rows = []
for target, expected in cases:
    witness = find_state(target)
    if (witness is not None) != expected:
        raise RuntimeError(('oracle answer', target, witness))
    if witness is not None and state(graph(witness['n'], witness['code'])) != target:
        raise RuntimeError(('literal witness', target, witness))
    rows.append({'state': list(target), 'expected_physical': expected, 'witness': witness})
print(json.dumps({'status': 'FINITE-EXACT PASS', 'cases': rows}, indent=2))
