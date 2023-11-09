# Clock-like signature attribution

See Methods in _Disentangling sources of clock-like mutations in germline and soma_, [`doi.org/10.1101/2023.09.07.556720`](https://doi.org/10.1101/2023.09.07.556720).

We assume that independent distributions govern the signature compositions of the age-dependent and age-independent mutations and denote them $P_a(s)$ and $P_b(s)$ respectively. 

Given linear dependence of the number of mutations on age, $\mathbb{E}[y|x]=ax+b$, the signature attribution at age $x$ is given by

$$
P(s|x) = \frac{ax}{ax+b} P_a(s) + \frac{b}{ax+b} P_b(s),
$$

and for Poisson regression with a log link function (not clock-like), $\log\mathbb{E}[y|x]=ax+b$,

$$
P(s|x) = (1-e^{-ax}) P_a(s) + e^{-ax} P_b(s),
$$

The code in `inference.py` allows to find the two distributions, $P_{a,b}(s)$, using expectation-maximization.
