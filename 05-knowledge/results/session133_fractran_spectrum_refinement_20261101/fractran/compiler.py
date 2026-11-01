"""Session133 refinement of the audited Session132 counter compiler.

Arithmetic shortcuts preserve N-valued semantics. Jump threading preserves register
traces up to silent JUMP steps, but inserts guards whenever threading would form
an illegal same-control transition. All emitted fractions are static and reduced.
"""
import ast
import json
from pathlib import Path
from baseline_compiler import (Compiler as BaselineCompiler, require, fractions,
    run_register, run_fractran_literal, run_fractran_exponents)


def thread_jumps(program):
    """Remove unreachable/silent states, never cancel a control-prime guard."""
    code = program['instructions']
    def resolve(pc):
        seen = set()
        while code[pc][0] == 'JUMP':
            require(pc not in seen, 'pure JUMP cycle unsupported by optimizer')
            seen.add(pc)
            pc = code[pc][1]
        return pc
    graph = {}
    # Iterative discovery avoids Python recursion limits on large fixed programs.
    entry = resolve(0)
    pending = [entry]
    while pending:
        pc = pending.pop()
        if pc in graph:
            continue
        ins = code[pc]
        targets = ins[2:] if ins[0] in ('INC', 'DECJZ') else ()
        new = []
        for target in targets:
            dest = resolve(target)
            if dest == pc:
                bridge = ('bridge', pc)
                graph[bridge] = ('JUMP', pc)
                new.append(bridge)
            else:
                new.append(dest)
                pending.append(dest)
        graph[pc] = (ins[0], ins[1], *new) if targets else ins
    keys = [entry] + sorted((k for k in graph if isinstance(k, int) and k != entry))
    keys += sorted(k for k in graph if isinstance(k, tuple))
    # Put the HALT last for a transparent artifact dictionary.
    keys = [k for k in keys if graph[k][0] != 'HALT'] + [k for k in keys if graph[k][0] == 'HALT']
    mapping = {k: i for i, k in enumerate(keys)}
    out = []
    for k in keys:
        ins = graph[k]
        if ins[0] in ('INC', 'DECJZ'):
            out.append((ins[0], ins[1], *(mapping[t] for t in ins[2:])))
        elif ins[0] == 'JUMP':
            out.append(('JUMP', mapping[ins[1]]))
        else:
            out.append(ins)
    refined = {**program, 'instructions': out}
    fractions(refined)  # Reduced-guard contract is checked, not merely asserted.
    return refined


class Compiler(BaselineCompiler):
    def statements(self, body):
        for stmt in body:
            if (isinstance(stmt, ast.Assign) and len(stmt.targets) == 1
                    and isinstance(stmt.targets[0], ast.Name)
                    and isinstance(stmt.value, ast.Constant)
                    and type(stmt.value.value) is int and stmt.value.value >= 0):
                self.set(self.reg(stmt.targets[0].id), stmt.value.value)
            elif isinstance(stmt, ast.AugAssign) and isinstance(stmt.op, ast.Add):
                require(isinstance(stmt.target, ast.Name), 'name augmented assignment only')
                target = self.reg(stmt.target.id)
                if isinstance(stmt.value, ast.Constant):
                    require(type(stmt.value.value) is int and stmt.value.value >= 0,
                            'nonnegative integer only')
                    for _ in range(stmt.value.value):
                        self.inc(target)
                else:
                    temporary = self.expr(stmt.value)
                    self.move_add(temporary, target)
                    self.release(temporary)
            else:
                super().statements([stmt])

    def compile(self, source):
        return thread_jumps(super().compile(source))


def write_artifacts(program, stem):
    """Write fixed ordered fraction shards and IR shards, each strictly < 100KB."""
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    f = fractions(program)
    metadata = {k: v for k, v in program.items() if k != 'instructions'}
    metadata['data_primes'] = f['data_primes']
    metadata['instruction_parts'] = []
    # Shards align at instruction boundaries; order is specified in metadata.
    chunks = []
    chunk = []
    size = 2
    for ins in program['instructions']:
        line = json.dumps(ins, separators=(',', ':'))
        if size + len(line) + 1 > 70000:
            chunks.append(chunk)
            chunk = []
            size = 2
        chunk.append(ins)
        size += len(line) + 1
    if chunk:
        chunks.append(chunk)
    for j, chunk in enumerate(chunks):
        path = stem.parent / f'{stem.name}.instructions.{j:02d}.json'
        path.write_text(json.dumps(chunk, separators=(',', ':')) + '\n')
        metadata['instruction_parts'].append(path.name)
    metadata['state_prime_parts'] = []
    for j, start in enumerate(range(0, len(f['state_primes']), 8000)):
        path = stem.parent / f'{stem.name}.states.{j:02d}.json'
        path.write_text(json.dumps(f['state_primes'][start:start+8000], separators=(',', ':')) + '\n')
        metadata['state_prime_parts'].append(path.name)
    chunks = []
    chunk = ''
    for a, b, _, _ in f['fractions']:
        line = f'{a}/{b}\n'
        if len(chunk) + len(line) > 80000:
            chunks.append(chunk)
            chunk = ''
        chunk += line
    if chunk:
        chunks.append(chunk)
    metadata['fraction_parts'] = []
    for j, chunk in enumerate(chunks):
        path = stem.parent / f'{stem.name}.fractions.{j:02d}.txt'
        path.write_text(chunk)
        metadata['fraction_parts'].append(path.name)
    (stem.parent / f'{stem.name}.meta.json').write_text(json.dumps(metadata, indent=2, sort_keys=True) + '\n')
    return dict(registers=program['register_count'], states=len(program['instructions']),
                fractions=len(f['fractions']), fraction_parts=len(chunks),
                fraction_bytes=sum(map(len, chunks)), input_prime=f['data_primes'][program['input']],
                output_prime=f['data_primes'][program['output']], entry_prime=f['state_primes'][0],
                halt_prime=f['state_primes'][-1])
