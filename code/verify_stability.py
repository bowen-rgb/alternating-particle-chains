"""Rational certificates and independent numerical checks. Python 3.9+."""
from fractions import Fraction as F
from functools import lru_cache
from math import cos, exp, fsum, pi, sin

SCALE = 10**42


def ceil(q):
    return -(-q.numerator // q.denominator)


@lru_cache(maxsize=None)
def exponential_bounds(z):
    """Enclose exp(-z), z >= 0, by alternating Taylor and outward squaring."""
    z = F(z)
    assert z >= 0
    halves = 0
    while z > 1:
        z /= 2
        halves += 1
    term = total = F(1)
    for j in range(1, 37):
        term *= -z / j
        total += term
    # The even partial sum is an upper bound and the preceding odd one a lower.
    lower = max(0, int((total - term) * SCALE))
    upper = min(SCALE, ceil(total * SCALE))
    for _ in range(halves):
        lower = lower * lower // SCALE
        upper = (upper * upper + SCALE - 1) // SCALE
    return F(lower, SCALE), F(upper, SCALE)


def add_scaled(bounds, coefficient, item):
    a, b = item
    if coefficient < 0:
        a, b = b, a
    return bounds[0] + coefficient*a, bounds[1] + coefficient*b


def threshold_bounds(b, kind):
    b = F(b)
    assert 1 <= b <= 3 and kind in ("D", "insertion")
    constant = 0 if kind == "D" else -1
    bounds = (F(constant), F(constant))
    for m in range(1, 11):
        if kind == "D":
            coeff = m*m*(2*b*m*m-1)
        else:
            coeff = 4*b*m*m+2
        bounds = add_scaled(bounds, (-1)**(m+1)*coeff, exponential_bounds(b*m*m))
    # |coefficient| <= 16 m^4; successive envelope ratios < 1/2 for m >= 11.
    tail = F(32 * 11**4, 2**121)
    return bounds[0]-tail, bounds[1]+tail


def enclose_energy(sites, b):
    # Positions are in mean-spacing units. Include the self term, then remove it.
    b = F(b)
    sites = tuple(map(F, sites))
    n = len(sites)
    assert n % 2 == 0 and all(0 <= x < n for x in sites)
    bounds = (F(1), F(1))
    for i, x in enumerate(sites):
        for j, y in enumerate(sites):
            coeff = F(-(-1)**(i+j), n)
            for k in range(-3, 4):
                bounds = add_scaled(bounds, coeff, exponential_bounds(b*(x-y+k*n)**2))
    assert b >= 1 and n >= 2
    # For |k| >= 4 each distance >= (|k|-1)n. Use exp(-z) <= 2^-z.
    tail = F(4*n, 2**(9*n*n))
    return bounds[0]-tail, bounds[1]+tail


def signed_root_brackets():
    for kind, left, right, sign in [
        ("D", "1.08182608", "1.08182609", 1),
        ("insertion", "2.47416700", "2.47416702", -1),
    ]:
        a = threshold_bounds(left, kind)
        b = threshold_bounds(right, kind)
        if sign == 1:
            assert a[1] < 0 < b[0]
        else:
            assert b[1] < 0 < a[0]
        print(f"{kind}: ({left}, {right}), certified signs")


def proof_constants():
    r = F(4, 3)**6 * exponential_bounds(7)[1]
    assert r < F(6, 1000)
    assert F(729, 3) * exponential_bounds(5)[1] / (1-r) < 2
    assert F(11, 5) * 729 * exponential_bounds(9)[1] / (1-r) < F(21, 100)
    assert F(4, 5) * exponential_bounds(F(11, 10))[0] > F(26, 100)
    b = F(11, 10)
    assert (2*b-1) > 4*(8*b-1)*exponential_bounds(3*b)[1]
    print("Analytic tail and derivative estimates verified by rational arithmetic")


def symbol(b, theta):
    return fsum((-1)**(m+1)*(2*b*m*m-1)*exp(-b*m*m)*2*sin(m*theta/2)**2
                for m in range(1, 90))


def image_energy(b, sites):
    n = len(sites)
    return 1 - fsum((-1)**(i+j)*exp(-b*(x-y+k*n)**2)
                    for i, x in enumerate(sites) for j, y in enumerate(sites)
                    for k in range(-3, 4))/n


def exponential_identity():
    q = tuple(map(F, ("1/2", "3/5", "4/7", "1/4")))
    n = len(q)
    response = []
    for i in range(n):
        weight = F(1)
        value = F(0)
        for k in range(n):
            current = q[(i+k) % n]
            value += weight * (1-current)
            weight *= -current
        response.append(value / (1-weight))
    t = [(1-z)/(1+z) for z in q]
    remainder = F(0)
    for i,z in enumerate(q):
        nxt = response[(i+1) % n]
        assert response[i]+z*nxt == 1-z
        remainder += t[i]*(response[i]-nxt)**2 / 4
    assert sum(t)-sum(response) == remainder and remainder > 0
    print("Exact gap identity verified on a nonuniform rational cell")


def independent_checks():
    for n, b, mode in [(16, 1.05, 1), (12, 1.5, 3), (6, 1.5, 3)]:
        theta = 2*pi*mode/n
        u = [cos(theta*j) for j in range(n)]
        step = 0.0002
        initial = image_energy(b, list(range(n)))
        displaced = [image_energy(b, [j+s*step*v for j,v in enumerate(u)])
                     for s in (-1, 1)]
        observed = (fsum(displaced)-2*initial)/step**2
        predicted = 8*b*symbol(b, theta)*fsum(v*v for v in u)/n
        assert abs(predicted-observed) < 0.000003, (predicted, observed)
        print(f"N={n}, b={b}, mode={mode}: Hessian {predicted:.8f}, difference {observed:.8f}")

    a = enclose_energy(range(6), "1.5")
    c = enclose_energy(["0", ".74", ".76", "1.5", "3", "4.5"], "1.5")
    assert F(".44130555") < a[0] <= a[1] < F(".44130557")
    assert F(".39161235") < c[0] <= c[1] < F(".39161238")
    assert a[0]-c[1] > F(".04969")
    print(f"Uniform {float(a[0]):.12f}, competitor {float(c[0]):.12f}, certified gain > 0.04969")


if __name__ == "__main__":
    exponential_identity()
    signed_root_brackets()
    proof_constants()
    independent_checks()
