"""First longitudinal instability for exp(-r**p). Python 3.9+ and NumPy.

Run: python verify_first_mode.py
q is lattice wave number. q=0 below denotes the acoustic limit after
division by 1-cos(q), not a translating chain.
"""

from decimal import Decimal, localcontext
from math import ceil, fsum, pi, sin, sqrt

import numpy as np


def weights(a, p, cutoff):
    last = ceil(cutoff**(1 / p) / a)
    m = np.arange(1, last + 1)
    r = a * m
    z = r**p
    coeff = p * r**(p - 2) * (p * z - (p - 1)) * np.exp(-z)
    return m, np.where(m % 2, coeff, -coeff)


def modulus(a, p, q, cutoff=56):
    m, coeff = weights(a, p, cutoff)
    ratio = m if q == 0 else np.sin(m * q / 2) / sin(q / 2)
    return fsum(coeff * ratio**2)


def all_modes(a, p, n):
    m, coeff = weights(a, p, 56)
    by_residue = np.bincount(m % n, weights=coeff, minlength=n)
    cosine_sum = np.fft.rfft(by_residue).real
    j = np.arange(1, n // 2 + 1)
    result = (fsum(coeff) - cosine_sum[1:]) / (2 * np.sin(pi * j / n)**2)
    # Avoid subtracting two nearly equal sums for the longest waves.
    for index in range(min(8, len(result))):
        q = 2 * pi * (index + 1) / n
        result[index] = modulus(a, p, q)
    return result


def bisect_zero(function, lo, hi):
    if not function(lo) < 0 < function(hi):
        raise ValueError("bracket does not straddle an instability")
    for _ in range(42):
        midpoint = (lo + hi) / 2
        if function(midpoint) < 0:
            lo = midpoint
        else:
            hi = midpoint
    return (hi + lo) / 2


def decimal_acoustic(a, p):
    with localcontext() as ctx:
        ctx.prec = 38
        gap, power = Decimal(str(a)), Decimal(str(p))
        total = Decimal(0)
        for m in range(1, ceil(65**(1 / p) / a) + 1):
            r = m * gap
            z = (power * r.ln()).exp()
            term = power * z * (power * z - (power - 1)) * (-z).exp() / gap**2
            total += term if m % 2 else -term
        return float(total)


def main():
    print("delta    acoustic gap    gap/sqrt(2 delta)   finite gaps (N=4,16,128)")
    for delta in (0.1, 0.03, 0.01, 0.003, 0.001, 0.0001):
        p = 1 + delta
        scale = sqrt(2 * delta)
        lo, hi = 0.45 * scale, 1.7 * scale
        acoustic = bisect_zero(lambda a: modulus(a, p, 0), lo, hi)
        finite = []
        for n in (4, 16, 128):
            def lowest(a):
                return float(np.min(all_modes(a, p, n)))
            onset = bisect_zero(lowest, lo, hi)
            assert np.argmin(all_modes(onset, p, n)) == 0
            finite.append(onset)
        assert abs(modulus(acoustic, p, 0, 70)) < 5e-8
        assert abs(decimal_acoustic(acoustic, p)) < 5e-8
        print(f"{delta:<7g} {acoustic:13.8f} {acoustic/scale:17.7f}"
              + "".join(f" {a:.8f}" for a in finite))

    # Independent all-mode sign scan, at successively smaller gaps.
    for delta in (0.01, 0.001, 0.0001):
        p = 1 + delta
        n = 512
        for a in np.geomspace(7.0, 0.6 * sqrt(2 * delta), 180):
            b = all_modes(a, p, n)
            if np.min(b) < 0:
                assert np.argmin(b) == 0
                q = 2 * pi / n
                assert abs(b[8] - modulus(a, p, 9 * q)) < 1e-8
                print(f"N={n}, delta={delta:g}, first negative sampled gap={a:.7f}, j=1")
                break
        else:
            raise AssertionError("scan did not reach instability")


if __name__ == "__main__":
    main()
