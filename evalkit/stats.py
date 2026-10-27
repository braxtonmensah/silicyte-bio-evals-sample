"""Pure-Python statistics used by builders and graders.

Each function here is cross-checked against SciPy inside the builders (see tasks/*/build.py), so a
bug in either implementation fails CI instead of silently changing a key.
"""
import math


def average_ranks(values):
    """1-based ranks; tied values share the mean of their positions."""
    n = len(values)
    order = sorted(range(n), key=lambda i: values[i])
    ranks = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and values[order[j + 1]] == values[order[i]]:
            j += 1
        rank = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = rank
        i = j + 1
    return ranks


def pearson(x, y):
    """Pearson correlation, or None when either input has zero variance."""
    n = len(x)
    if n != len(y) or n < 2:
        return None
    mx = math.fsum(x) / n
    my = math.fsum(y) / n
    sxx = math.fsum((a - mx) ** 2 for a in x)
    syy = math.fsum((b - my) ** 2 for b in y)
    if sxx == 0.0 or syy == 0.0:
        return None
    sxy = math.fsum((a - mx) * (b - my) for a, b in zip(x, y))
    return sxy / math.sqrt(sxx * syy)


def spearman(x, y):
    """Spearman rank correlation with average ranks for ties; None if undefined."""
    return pearson(average_ranks(x), average_ranks(y))


def log_comb(n, k):
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def hypergeom_sf_ge(k, N, K, n):
    """P(X >= k) for X ~ Hypergeometric(population N, K successes, n draws)."""
    lo = max(0, n - (N - K))
    hi = min(K, n)
    if k <= lo:
        return 1.0
    if k > hi:
        return 0.0
    denom = log_comb(N, n)
    terms = [math.exp(log_comb(K, i) + log_comb(N - K, n - i) - denom) for i in range(k, hi + 1)]
    return min(1.0, math.fsum(terms))


def benjamini_hochberg(pvalues):
    """Benjamini-Hochberg step-up q-values with monotonicity enforced, capped at 1."""
    m = len(pvalues)
    order = sorted(range(m), key=lambda i: (pvalues[i], i))
    q = [0.0] * m
    running = 1.0
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        running = min(running, pvalues[i] * m / rank)
        q[i] = min(running, 1.0)
    return q
