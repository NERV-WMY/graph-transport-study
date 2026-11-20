"""Reproduce the exclusion-current benchmarks from the configuration generator."""

from itertools import combinations

import numpy as np


def build_generator(channels, rail_edges, mixing_rate):
    """Return the row generator and a mask for counted left injections.

    The first half of the bits encode entrance sites. The second half
    encodes exit sites in the same channel order.
    """
    states = 1 << (2 * channels)
    generator = np.zeros((states, states))
    injections = np.zeros((states, states), dtype=bool)
    bonds = [(i, channels + i, 1.0) for i in range(channels)]
    bonds.extend((i, j, mixing_rate) for i, j in rail_edges)

    for state in range(states):
        for i in range(channels):
            entrance = 1 << i
            exit_site = 1 << (channels + i)
            if not state & entrance:
                target = state | entrance
                generator[state, target] += 1
                injections[state, target] = True
            if state & exit_site:
                generator[state, state ^ exit_site] += 1

        for i, j, rate in bonds:
            if (state >> i & 1) != (state >> j & 1):
                target = state ^ (1 << i) ^ (1 << j)
                generator[state, target] += rate

    np.fill_diagonal(generator, -generator.sum(axis=1))
    return generator, injections


def current_statistics(generator, injections):
    """Solve the stationary and Poisson equations for the variance rate."""
    stationary_system = generator.T.copy()
    stationary_system[-1] = 1
    normalization = np.zeros(len(generator))
    normalization[-1] = 1
    stationary = np.linalg.solve(stationary_system, normalization)

    drift = (generator * injections).sum(axis=1)
    current = stationary @ drift
    # Subtracting the stationary projection fixes the corrector's mean.
    corrector = np.linalg.solve(generator - stationary[None, :], current - drift)
    increments = injections + corrector[None, :] - corrector[:, None]
    variance = stationary @ (generator * increments**2).sum(axis=1)
    return current, variance


def main():
    tolerance = {"rtol": 1e-11, "atol": 1e-12}
    for rate in (0, 1, 3):
        generator, injections = build_generator(2, [(0, 1)], rate)
        current, variance = current_statistics(generator, injections)
        expected_noise = 8 / 27 - 2 * rate / (27 * (6 + 5 * rate))
        np.testing.assert_allclose(
            [current, variance], [2 / 3, expected_noise], **tolerance
        )
        print(f"Four sites, g={rate}: I={current:.9f}, variance rate={variance:.9f}")
        if rate == 1:
            # Suppress counted transitions while retaining the escape rates.
            killed = np.where(injections, 0, generator)
            eigenvalue = np.linalg.eigvals(killed).real.max()
            np.testing.assert_allclose(eigenvalue, -3 + np.sqrt(3), **tolerance)
            print(f"Four sites, g=1: zero-count decay rate={-eigenvalue:.9f}")

    graphs = [
        (f"C{m}", m, [(i, (i + 1) % m) for i in range(m)], expected)
        for m, expected in [(3, 14 / 99), (4, 964 / 6831), (5, 11 / 78)]
    ]
    graphs.append(("K4", 4, list(combinations(range(4), 2)), 53 / 378))
    for name, channels, edges, expected_noise in graphs:
        rails = edges + [(i + channels, j + channels) for i, j in edges]
        generator, injections = build_generator(channels, rails, 1)
        current, variance = current_statistics(generator, injections)
        np.testing.assert_allclose(
            [current, variance / channels], [channels / 3, expected_noise], **tolerance
        )
        print(f"{name} on both rails, g=1: variance rate/channel={variance / channels:.9f}")


if __name__ == "__main__":
    main()
