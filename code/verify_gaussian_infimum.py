"""Independent numerical checks for en/gaussian-infimum.md (standard library)."""

import cmath
import math
import random


def parameters(length, alpha):
    return math.sqrt(math.pi / alpha) / length, math.exp(-math.pi**2 / (alpha * length**2))


def image_energy(sites, length, alpha):
    """Sum the original directed interactions, skipping only a site's zero image."""
    radius = 1 + math.ceil(math.sqrt(50 / alpha) / length)
    return math.fsum(
        (-1 if (i + j) % 2 == 0 else 1)
        * math.exp(-alpha * (x - y + k * length) ** 2)
        for i, x in enumerate(sites)
        for j, y in enumerate(sites)
        for k in range(-radius, radius + 1)
        if i != j or k != 0
    ) / len(sites)


def modes(q):
    return range(1, 2 + math.ceil(math.sqrt(55 / -math.log(q))))


def dipole_energy(count, c, q):
    return 1 - 8 * c * math.fsum(q ** (n * n) for n in modes(q) if n % 2) / count


def fourier_energy(sites, length, c, q):
    angles = [math.tau * x / length for x in sites]
    def structure(n):
        return sum((-1) ** j * cmath.exp(1j * n * t) for j, t in enumerate(angles))
    return 1 - 2 * c * math.fsum(q ** (n*n) * abs(structure(n))**2 for n in modes(q)) / len(sites)


def semicircle_distance(sites, length):
    """Overlap is piecewise affine in the arc start, with breakpoints at sites and antipodes."""
    angles = [math.tau * x / length for x in sites]
    arcs = list(zip(angles[::2], angles[1::2]))
    mass = math.fsum(b - a for a, b in arcs)
    starts = [t for x in angles for t in (x, (x - math.pi) % math.tau)]
    def overlap(start):
        return math.fsum(max(0, min(b, start + math.pi + shift) - max(a, start + shift))
                         for a, b in arcs for shift in (0, -math.tau))
    return max(0, mass + math.pi - 2 * max(map(overlap, starts)))


def spacing_trial(count, length, h):
    m = count // 2
    return [0.] + [length/2 + j*h for j in range(1-m, m)]


def trial_coefficient(count, length, c, q):
    """Coefficient of physical h^2 for the one-site plus (N-1)-site trial."""
    m = count // 2
    curvature = 4*c*math.fsum(n*n*q**(n*n) for n in modes(q) if n % 2)
    return m*(m-1)/count * curvature * (math.tau / length)**2


def close(a, b, tol=2e-12):
    if abs(a - b) > tol:
        raise AssertionError((a, b, tol))


def balanced_trial(count, length, h):
    m = count // 2
    k = (m - 1) // 2
    l = m - 1 - k
    # Translate the positive group so that its first site is at zero.
    first = [j*h for j in range(2*k+1)]
    second = [length/2 + (j+k)*h for j in range(-l, l+1)]
    return first + second


def grid_minimum(count, mesh, length, alpha):
    """Exhaust all grid cells up to translation (a positive site is at zero)."""
    from itertools import combinations
    h = length/mesh
    kernel = [math.fsum(math.exp(-alpha*(i*h+k*length)**2) for k in range(-9, 10))
              for i in range(mesh)]
    pairs = [(i, j, (-1)**(i+j)) for i in range(count) for j in range(i)]
    lowest = math.inf
    for tail in combinations(range(1, mesh), count-1):
        sites = (0,) + tail
        cross = sum(sign*kernel[sites[i]-sites[j]] for i, j, sign in pairs)
        lowest = min(lowest, 1-kernel[0]-2*cross/count)
    return lowest


def check_balanced_clusters(length, alpha, c, q):
    for count in (4, 6, 8, 10, 12):
        m = count//2
        expected = trial_coefficient(count, length, c, q) * (m*m//2) / (m*(m-1))
        h = .001
        e = image_energy(balanced_trial(count, length, h), length, alpha)
        ratio = (e - dipole_energy(count, c, q))/h**2
        close(ratio/expected, 1, 4e-5)
        print(f"N={count:2d}  balanced excess/h^2={ratio:.8f}  sharp limit={expected:.8f}")
    for count, meshes in ((4, (16, 32, 64, 128)), (6, (16, 24, 32)), (8, (16,))):
        for mesh in meshes:
            h = length/mesh
            minimum = grid_minimum(count, mesh, length, alpha)
            trial = image_energy(balanced_trial(count, length, h), length, alpha)
            close(minimum, trial)
            print(f"Grid N={count} mesh={mesh:3d} grid min={minimum:.12f} balanced trial={trial:.12f}")


def check():
    length, alpha = 4., .25
    c, q = parameters(length, alpha)
    for sites in ([0, 1, 2, 3], [0, 1.5, 2, 2.5], [0, 2, 2, 2]):
        direct = image_energy(sites, length, alpha)
        close(direct, fourier_energy(sites, length, c, q))
        print(f"{sites}    E = {direct:.15f}")
    bound = dipole_energy(4, c, q)
    close(image_energy([0, 2, 2, 2], length, alpha), bound)
    close(image_energy([0, .75, .75, 2], length, alpha), bound)
    close(image_energy([0, 2], length, alpha), dipole_energy(2, c, q))
    rng = random.Random(20261126)
    count_checked = 0
    for target_q in (q, .20, 1/3):
        a = math.pi**2 / (-math.log(target_q)*length**2)
        c1, q1 = parameters(length, a)
        for count in (2, 4, 6, 8, 12):
            infimum = dipole_energy(count, c1, q1)
            for _ in range(45):
                sites = sorted(length*rng.random() for _ in range(count))
                e = image_energy(sites, length, a)
                close(e, fourier_energy(sites, length, c1, q1))
                d = semicircle_distance(sites, length)
                lower = 32*c1*(q1/3)/count * math.sin(d/8)**2
                if e - infimum < lower - 2e-12:
                    raise AssertionError(("stability", sites, lower, e-infimum))
                if count >= 4:
                    gaps = [b-x for x, b in zip(sites, sites[1:] + [sites[0]+length])]
                    if d + 2e-12 < math.tau*min(gaps)/length:
                        raise AssertionError(("gap geometry", d, gaps))
                count_checked += 1
    print(f"{count_checked} deterministic random cells passed both energy and stability checks.")
    for count in (4, 6, 8):
        coeff = trial_coefficient(count, length, c, q)
        for h in (.04, .02, .01, .005, .0025):
            energy = image_energy(spacing_trial(count, length, h), length, alpha)
            ratio = (energy - dipole_energy(count, c, q)) / h**2
            print(f"N={count:2d}  h={h:.4f}  excess/h^2={ratio:.8f}  limit={coeff:.8f}")
        close(ratio/coeff, 1, 2e-4)
    check_balanced_clusters(length, alpha, c, q)
    narrow_alpha = 16.
    narrow_c, narrow_q = parameters(length, narrow_alpha)
    narrow_uniform = image_energy([0, 1, 2, 3], length, narrow_alpha)
    if narrow_uniform >= dipole_energy(4, narrow_c, narrow_q):
        raise AssertionError("The narrow Gaussian comparison failed.")
    print("Narrow Gaussian comparison passed.")


if __name__ == "__main__":
    check()
