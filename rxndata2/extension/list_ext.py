from itertools import permutations, product


def generate_patterns(lst):
    """Generate all order permutations of the given list and return them as a list of strings joined by '.'"""
    return ['.'.join(p) for p in permutations(lst)]


def generate_combinations(lst1, lst2):
    """Generate all combinations of ordered patterns from two lists"""
    patterns1 = generate_patterns(lst1)
    patterns2 = generate_patterns(lst2)
    return list(product(patterns1, patterns2))
