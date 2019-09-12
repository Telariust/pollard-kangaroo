# pollard-kangaroo

Pollard, kangaroo method, solving discrete logarithm problem (DLP) using pseudorandom walks.
Its the kangaroo method (of [11, p. 13]) using distinguished points.
Runtime expected of 2w^(1/2) group operations.

[11] P. C. van Oorschot and M. J. Wiener, Parallel collision search with cryptanalytic applications, J. Cryptology, #12 (1999)

# How pollard-kangaroo works, the Tame and Wild kangaroos, is a simple explanation.

Suppose there is pubkeyX, unknow privkeyX, but privkeyX is in range w=[L..U]
The keys have a property - if we increase pubkey by S, then its privkey will also increase by S.
We start step-by-step to increase pubkeyX by S(i), keeping sumX(S(i)). This is a Wild kangaroo.
We select a random privkeyT from range [L..U], compute pubkeyT.
We start step-by-step to increment pubkeyT by S(i) while maintaining sumT(S(i)). This is a Tame kangaroo.
The size of the jump S(i) is determined by the x coordinate of the current point, so if a Wild or Tame kangaroo lands on one point, their paths will merge.
(we are concerned with pseudo random walks whose next step is determined by the current position)
Thanks to the Birthday Paradox (Kruskal's card trick), their paths will one day meet.
Knowing each traveled path, privkeyX is calculated.
The number of jumps is approximately 2*(w^1/2) group operations, which is much less than a full search w.

# articles

base
https://andrea.corbellini.name/2015/05/17/elliptic-curve-cryptography-a-gentle-introduction/
hand translate to ru/ua (recommend instead translate.google)
https://habr.com/ru/post/335906/

0) best of the best, all in 1, epic,  2012
Chapter 14. Factoring and Discrete Logarithms using Pseudorandom Walks
https://www.math.auckland.ac.nz/~sgal018/crypto-book/ch14.pdf

1) with reference to old
J. M. Pollard, “Kangaroos, monopoly and discrete logarithms,” Journal of Cryptology, #13 (2000).
https://web.northeastern.edu/seigen/11Magic/KruskalsCount/PollardKangarooMonopoly.pdf
(good dir web.northeastern.edu/seigen/11Magic/KruskalsCount/)

2) about parallelism problems
P. C. van Oorschot and M. J. Wiener, Parallel collision search with cryptanalytic applications, J. Cryptology, #12 (1999)
https://people.scs.carleton.ca/~paulv/papers/JoC97.pdf
