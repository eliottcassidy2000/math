#!/usr/bin/env python3
"""Exact valuation-matrix analysis of the frozen literal 22-fraction list."""
from fractions import Fraction
from pathlib import Path
import json
from verify import load_program, PRIMES, STATES, TRANSFERS, MAIN, AUX, CHECKPOINT, require, formula, SEED

HERE = Path(__file__).resolve().parent


def factor(n):
    row = []
    for p in PRIMES:
        e = 0
        while n % p == 0:
            e += 1
            n //= p
        row.append(e)
    require(n == 1, 'unknown prime')
    return row


def rank(matrix):
    a = [list(map(Fraction, row)) for row in matrix]
    pivots = []
    for col in range(len(a[0])):
        pivot = next((i for i in range(len(pivots), len(a)) if a[i][col]), None)
        if pivot is None:
            continue
        row = len(pivots)
        a[row], a[pivot] = a[pivot], a[row]
        divisor = a[row][col]
        a[row] = [x / divisor for x in a[row]]
        for i in range(len(a)):
            if i != row and a[i][col]:
                scale = a[i][col]
                a[i] = [x - scale*y for x, y in zip(a[i], a[row])]
        pivots.append(col)
    return len(pivots)


def determinant(matrix):
    a = [list(map(Fraction, row)) for row in matrix]
    result = Fraction(1)
    for col in range(len(a)):
        pivot = next((i for i in range(col, len(a)) if a[i][col]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != col:
            a[pivot], a[col] = a[col], a[pivot]
            result *= -1
        value = a[col][col]
        result *= value
        for i in range(col+1, len(a)):
            scale = a[i][col] / value
            a[i] = [x-scale*y for x,y in zip(a[i],a[col])]
    return result


def strong(vertices, edges):
    for root in vertices:
        visited, todo = {root}, [root]
        while todo:
            v = todo.pop()
            for a, b in edges:
                if a == v and b not in visited:
                    visited.add(b)
                    todo.append(b)
        if visited != set(vertices):
            return False
    return True


def main():
    columns = [[a-b for a,b in zip(factor(num),factor(den))] for num,den in load_program()]
    matrix = [list(row) for row in zip(*columns)]
    r = rank(matrix)
    require(r == 22 and len(matrix) == 23)
    conservation = [int(p in STATES) for p in PRIMES]
    require(all(sum(a*b for a,b in zip(conservation,column)) == 0 for column in columns))
    weights = {2:12, 3:10, 5:8, 7:6, 11:18, 13:16, 17:14, 19:8,
               23:0, 29:-7, 31:-8, 37:-6, 41:-7, 43:-5, 47:-6,
               53:-4, 59:-5, 61:-3, 67:-4, 71:-2, 73:-3, 79:-1, 83:-2}
    grading = [weights[p] for p in PRIMES]
    increments = [sum(w*d for w,d in zip(grading,column)) for column in columns]
    require(increments == [1]*22, 'microstep grading is not exact')
    seed_grade = sum(w*e for w,e in zip(grading,factor(SEED)))
    require(seed_grade == 108)
    checkpoints = json.loads((HERE/'verification.json').read_text())['checkpoints']
    totals = []
    for row in checkpoints:
        n = row['n']
        exact = sum(weights[p]*formula(n+j) for j,p in enumerate((2,3,5,7,11,13)))+8*n-108
        require(exact == row['rational_steps_from_seed'], 'closed-form step count fails')
        totals.append(exact)
    trace = json.loads((HERE/'first_transition_trace.json').read_text())
    for row in trace:
        grade = sum(w*e for w,e in zip(grading,factor(int(row['integer']))))
        require(grade == seed_grade+row['step'], 'literal trace grading fails')
    minor = [row for p,row in zip(PRIMES,matrix) if p != CHECKPOINT]
    det = determinant(minor)
    require(abs(det) == 1)
    control_edges = [(CHECKPOINT, MAIN[0])]
    for i in range(7):
        control_edges += [(MAIN[i], AUX[i]), (AUX[i], MAIN[i]),
                          (MAIN[i], MAIN[i+1] if i<6 else CHECKPOINT)]
    transfer_edges = [(source,dest) for source,dests in TRANSFERS for dest in dests]
    require(strong(STATES,control_edges))
    require(strong((2,3,5,7,11,13,17),transfer_edges))
    result = {'matrix_row_primes': list(PRIMES), 'matrix': matrix,
              'shape': [23,22], 'rank_over_Q': r, 'right_kernel_dimension': 0,
              'left_kernel_dimension': 1, 'left_kernel_generator': conservation,
              'microstep_grading_weights': weights, 'column_grading_increments': increments,
              'seed_grade': seed_grade, 'checkpoint_steps_from_grading': totals,
              'literal_first_transition_grading_checked': True,
              'minor_deleted_row_prime': CHECKPOINT, 'minor_determinant': int(det),
              'control_graph_edges': control_edges, 'control_graph_SCC_count': 1,
              'transfer_graph_edges': transfer_edges, 'transfer_graph_SCC_count': 1,
              'scope': 'No nonempty product of fractions equals 1, hence no finite periodic orbit from any positive integer. This does not classify termination or output behavior of arbitrary multi-control seeds.'}
    (HERE/'geometry.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS: 23x22 valuation matrix has rank22; right kernel0; left kernel is total control exponent.')
    print('Deleted-23 row minor determinant:', det)
    print('Integer grading increments exactly1 on all22 fractions; seed grade108; n24 grade-derived steps:', totals[-1])
    print('Both semantic control graph and seven-register transfer graph are strongly connected.')


if __name__ == '__main__':
    main()
