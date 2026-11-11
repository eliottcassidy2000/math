#!/usr/bin/env python3
"""Independent iterative-g audit of the compiler's odd-multiplier matrices."""
import ast
from hashlib import sha256
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
PRODUCER = HERE.parent / 'compiler/check.py'
MULTIPLIERS = (1, 3, 5, 7, 9, 15, 17, 21, 25, 35)
ADDRESSES = (0, 1, 2, 7, 18, 123, 1024, 170807)


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def iterative_g(n):
    """Read the recurrence downwards; no producer g, state or reader used."""
    need(n >= 0, 'negative address')
    total = 1
    while n > 1:
        half, bit = divmod(n, 2)
        if bit == 0:
            total += half + (half % 2)
        n = half
    return total


def expected_state(multiplier, address):
    return tuple(iterative_g(multiplier * address + r)
                 for r in range(multiplier)) + (
                     address, int(address % 2 == 0), address % 2)


def multiply(matrix, vector):
    need(len(matrix) == len(vector), 'matrix/vector dimension')
    return tuple(sum(matrix[i][j] * vector[j]
                     for j in range(len(vector)))
                 for i in range(len(matrix)))


def read_address(multiplier, address, matrices):
    state = expected_state(multiplier, 0)
    digits = []
    rest = address
    while rest:
        rest, digit = divmod(rest, 2)
        digits.append(digit)
    for digit in reversed(digits):
        state = multiply(matrices[digit], state)
    return state


def producer_matrix_constructor():
    # Extract only the operation under audit. Do not import the producer's
    # module, its recursive g, its binary reader, or its graph helper imports.
    tree = ast.parse(PRODUCER.read_text())
    definitions = [node for node in tree.body
                   if isinstance(node, ast.FunctionDef)
                   and node.name == 'digit_matrix']
    need(len(definitions) == 1, 'unique producer matrix constructor')
    module = ast.Module(body=definitions, type_ignores=[])
    namespace = {'need': need}
    exec(compile(module, str(PRODUCER), 'exec'), namespace)
    return namespace['digit_matrix']


def main():
    digit_matrix = producer_matrix_constructor()
    state_checks = transition_checks = channel_checks = 0
    for multiplier in MULTIPLIERS:
        matrices = tuple(digit_matrix(multiplier, digit) for digit in (0, 1))
        initial = expected_state(multiplier, 0)
        need(multiply(matrices[0], initial) == initial, 'leading-zero invariant')
        for address in ADDRESSES:
            state = read_address(multiplier, address, matrices)
            need(state == expected_state(multiplier, address),
                 f'iterative direct state m={multiplier}, n={address}')
            state_checks += 1
            for digit in (0, 1):
                destination = 2 * address + digit
                need(multiply(matrices[digit], state)
                     == expected_state(multiplier, destination),
                     f'whole matrix transition m={multiplier}, n={address}, a={digit}')
                transition_checks += 1
                for residue in range(multiplier):
                    quotient, remainder = divmod(multiplier * digit + residue, 2)
                    predicted = state[quotient]
                    if remainder == 0:
                        predicted += (multiplier * address + quotient
                                      + ((address + quotient) % 2))
                    need(predicted == iterative_g(multiplier * destination + residue),
                         f'carry channel m={multiplier}, n={address}, a={digit}, r={residue}')
                    channel_checks += 1
    need((state_checks, transition_checks, channel_checks) == (80, 160, 2208),
         'declared audit universe')

    # Oddness is essential to replacing (m*n+j)%2 with (n+j)%2.
    even_input = expected_state(2, 1)
    even_generated = multiply(digit_matrix(2, 0), even_input)
    even_actual = expected_state(2, 2)
    need(even_generated == (6, 3, 2, 1, 0), 'specified even hostile generated state')
    need(even_actual == (5, 3, 2, 1, 0), 'specified even hostile actual state')
    need(even_generated != even_actual, 'oddness hostile')

    print(json.dumps({
        'status': 'PASS',
        'scope': 'finite independent carry audit; no global target coverage assertion',
        'multipliers': MULTIPLIERS,
        'addresses': ADDRESSES,
        'checks': {'direct_states': state_checks,
                   'matrix_transitions': transition_checks,
                   'carry_channels': channel_checks,
                   'leading_zero_states': len(MULTIPLIERS)},
        'even_multiplier_hostile': {
            'multiplier': 2, 'address': 1, 'digit': 0,
            'generated_state': even_generated, 'actual_state': even_actual,
            'first_failed_identity': '(m*n+j)%2 = (n+j)%2'},
        'producer_sha256': sha256(PRODUCER.read_bytes()).hexdigest(),
        'script_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
    }, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
