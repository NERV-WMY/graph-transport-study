# Exact limits in quantum absorption and exclusion transport

A research note on quantum absorption and current fluctuations in the symmetric simple exclusion process on graphs.

[Read the note (PDF)](graph_transport_study.pdf)

## Quantum absorption

The initial state is maximally mixed on $n$ sites. Write $\kappa$ for the loss rate per detector and $\sigma^2$ for the sum of squared hoppings across the detector boundary.

**One detector.**

$$
\begin{aligned}
\tau &\ge \frac{n}{\kappa}+\frac{(\kappa+2\gamma)(n-1)^2}{4n\sigma^2},\\
\tau_{\ast} &=\frac{n-1}{\sigma}+\frac{\gamma(n-1)^2}{2n\sigma^2}.
\end{aligned}
$$

$\tau$ is the mean absorption time and $\gamma$ is the site dephasing rate. A weighted path reaches the bound at $\gamma=0$, and a Hermite-root graph reaches it with dephasing. The minimum  $\tau_{\ast}$ uses $\kappa_{\ast}=2n\sigma/(n-1)$.

**Several coherent detectors.**

$$
\tau \ge \frac{n}{r\kappa}+\frac{\kappa(n-r)^2}{4n\sigma^2}.
$$

For $r$ detectors and zero dephasing, equality is attainable exactly when $r$ divides $n$.

## Exclusion current at full bias

**Mean and noise.**

$$
\begin{aligned}
I &= \sum_{\mathrm{bonds}} w_{ij}(h_i-h_j)^2,\\
\mathcal V &= \frac I3+\frac23\sum_{\mathrm{bonds}}w_{ij}(h_i-h_j)^4
              +2s^{\mathsf{T}}A_2^{-1}s .
\end{aligned}
$$

$h_i$ is the stationary occupation and $w_{ij}$ is the bond rate. The sums include reservoir contacts. $\mathcal V$ is the long-time variance per unit time of the left injection count.

The covariance term uses $s_{ij}=w_{ij}(h_i-h_j)^2$ on internal pairs and the two-particle Laplacian $A_2$ with killing at the reservoirs.

**Small voltage drops.**

$$
\frac13\le \frac{\mathcal V}{I}\le \frac{1}{3(1-\epsilon)},
\qquad \epsilon=\max_{\mathrm{bonds}} |h_i-h_j|<1.
$$

For $I>0$, the Fano factor approaches $1/3$ as the largest positive-rate bond drop $\epsilon$ tends to zero.

Equal-occupation links preserve $I$ and can lower $\mathcal V$. For any internal exchange rate $g$, the long-time cumulant generating function obeys $\partial_g\Lambda(\lambda)\le0$ at negative tilt $\lambda$.

## Mixing parallel two-site channels

**Finite budget.** For $m\ge2$ channels with unit contact and internal rates and total added mixing rate $W$, set $a=2W/[m(m-1)]$.

$$
\frac{\mathcal V_{\min}}{m}
   =\frac{11}{81}+\frac{2(3+a)}{81[6+(3m-1)a]}.
$$

The minimum is attained exactly when each pair of channels has combined mixing rate $a$ across both rails.

**Sparse mixing.** For the same vertex- and edge-transitive mixing graph on both rails, let $L$ be its weighted Laplacian.

$$
\begin{aligned}
\frac{\mathcal V}{m}&=\frac{11}{81}+\frac{2}{81(2+\beta)},\\
\beta&=\frac{2b_2b_6}{b_2+b_6},\\
b_\mu&=\left[\frac1m\mathrm{Tr}\big((\mu I_m+2L)^{-1}\big)\right]^{-1}-\mu.
\end{aligned}
$$

Here $I_m$ is the identity matrix and $\beta=0$ in the absence of mixing. At fixed degree and edge rate, transitive mixing graphs with growing girth approach the regular-tree noise bound.

## Validate the benchmarks

The scripts solve the finite absorption and exclusion models directly and compare the results with the formulas in the note.

    python -m pip install numpy
    python validate_absorption.py
    python validate_exclusion.py
