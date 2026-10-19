# Tensor Attack on HFE

In this repository we present versions of a new attack that can be used against the following signature schemes:

- https://eprint.iacr.org/2021/084 (Key recovery)

- https://eprint.iacr.org/2024/1706 (Key distinguish in smaller parameters, which I think can be also used for key recovery, and perhaps can be extended for full parameters)

- https://eprint.iacr.org/2026/1650 (Key distinguish, which I think can be also used for key recovery)

- https://eprint.iacr.org/2026/404 (Key recovery)

This is a project I've been working on since August.

The idea behind the attack was obtained without any AI assistance.

## Notes

- Even though in the article it is specified the rank bound of $n(\binomial{q}{p}-\binomial{q-d}{p-d})$ for the flattening of a tensor obtained from the learning low degree polynomials, it actually works against any family of linear polynomials provided they are sparse (that is, if we pick a set $S$ of size at most $d$, and sample functions on the form $\sum_{i \in S} a_i x^{q^i}$).

## References

- https://eprint.iacr.org/2026/290

- https://arxiv.org/abs/2411.14344

- Multivariate Public Key Cryptosystems, Jintai Ding and Albrecht Petzoldt and Dieter S. Schmidt


