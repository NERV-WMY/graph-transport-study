"""Verify finite quantum absorption examples from the note with NumPy."""

import numpy as np


def mean_lifetime(hamiltonian, detectors, loss_rate, dephasing_rate=0.0):
    """Solve the integrated equation for an initially maximally mixed state."""
    n = hamiltonian.shape[0]
    projector = np.zeros((n, n))
    projector[detectors, detectors] = 1
    operator = np.empty((n * n, n * n), dtype=complex)

    # Each column is the integrated operator applied to one matrix unit.
    for column in range(n * n):
        unit = np.zeros((n, n))
        unit.flat[column] = 1
        commutator = hamiltonian @ unit - unit @ hamiltonian
        loss = projector @ unit + unit @ projector
        off_diagonal = unit - np.diag(unit.diagonal())
        operator[:, column] = (
            1j * commutator + loss_rate / 2 * loss + dephasing_rate * off_diagonal
        ).ravel()

    integrated = np.linalg.solve(operator, (np.eye(n) / n).ravel())
    return np.trace(integrated.reshape(n, n))


def check(label, computed, expected):
    np.testing.assert_allclose(computed, expected, rtol=1e-10, atol=1e-12, err_msg=label)
    print(f"{label}: {computed.real:.10f} (expected {expected:.10f})")


def main():
    # The two four-site benchmarks use unit squared detector hopping.
    n = 4
    sigma = 1.0
    loss_rate = 8 / 3
    path_weights = sigma * np.sqrt(np.arange(n - 1, 0, -1) / (n - 1))
    path = np.diag(path_weights, 1) + np.diag(path_weights, -1)
    check("Coherent path", mean_lifetime(path, [0], loss_rate), 3)

    roots = np.polynomial.hermite_e.hermeroots([0, 0, 0, 1])
    leaf_hopping = sigma / np.sqrt(n - 1)
    hermite_graph = np.zeros((n, n), dtype=complex)
    hermite_graph[0, 1:] = hermite_graph[1:, 0] = leaf_hopping
    for j in range(n - 1):
        for k in range(j + 1, n - 1):
            hermite_graph[j + 1, k + 1] = 2j * leaf_hopping / (roots[k] - roots[j])
            hermite_graph[k + 1, j + 1] = hermite_graph[j + 1, k + 1].conjugate()
    check(
        "Dephased Hermite graph",
        mean_lifetime(hermite_graph, [0], loss_rate, dephasing_rate=1),
        33 / 8,
    )

    # Check both branches of the fixed-rate optimum for two absorbing sites.
    for loss_rate in (1.0, 2 * np.sqrt(5) * sigma):
        minimum_y = loss_rate**2 + 2 * sigma**2
        if loss_rate < np.sqrt(6) * sigma:
            optimal_y = sigma * np.sqrt(8 * minimum_y)
            sink_hopping = np.sqrt((optimal_y - minimum_y) / 4)
            expected = (10 + 4 * np.sqrt(2) * np.sqrt(minimum_y) / sigma) / (12 * loss_rate)
        else:
            sink_hopping = 0.0
            expected = np.sqrt(5) / (3 * sigma)

        weights = [sigma, sink_hopping]
        two_sink_path = np.diag(weights, 1) + np.diag(weights, -1)
        check(
            f"Three-site, two-sink system (loss rate {loss_rate:.6g})",
            mean_lifetime(two_sink_path, [1, 2], loss_rate),
            expected,
        )


if __name__ == "__main__":
    main()
