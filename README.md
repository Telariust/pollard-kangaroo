# Pollard-Kangaroo

Pollard, kangaroo method, solving discrete logarithm problem (DLP) using pseudorandom walks.
Its the kangaroo method (of [11, p. 13]) using distinguished points.
Runtime expected of 2w<sup>1/2</sup> group operations.

[11] P. C. van Oorschot and M. J. Wiener, Parallel collision search with cryptanalytic applications, J. Cryptology, #12 (1999)

# Feature:

 - singlecore, multicore
 - python2/3 compatibility (print/time/input/xrange/IntDiv)
 - auto adaptation under environment
 - raw python/coincurve+cffi/gmpy2
 - some checks, debug, more info, minor accelerate
 - (un)compress/custom/random pubkey
 - format time: 0y 0m 0d 00:00:00s 000ms
 - format prefix SI: Kilo/Mega/Giga/Tera/Peta/Exa
 - heuristic prognosis (avg)time and (avg)jump per bits
 - repair kangaroos, that falls into the trap of same sequences
 - location privkey in keyspace on the strip
 - percent status progress, lost/left time
 - support arbitrary range (keyspace) start:end
 - profiles settings/algorithms

Expected in the future
 - named argv
 - precision S(i) set without problem of odd of pow2
 - <s>neural network for user entertainment speaking in the voice of Jarvis</s> or anything else

Multicore options:
 1) naive parallelism, split equally the search range
 2) parallelism by Oorschot&Wiener
 - has problem of collisions between members of the same herd
 3) parallelism by Pollard (best choice)
 - coprime/odd numbers U=V=2p+/-1, U+V<=Ncores
 - without problem of collisions between members of the same herd
 
*Acceleration on several devices (i.e. with shared RAM) for 2) and 3) (w<sup>1/2</sup>/Ncores) can only be achieved using the pool.
Filtered distinguished points should be sent to the pool in a single hashtable.

*removed support coincurve lib, reason: if u can install coincurve - that u can install gmpy2

# Compilation

Its python implementation, need install python 2 or 3.

Raw python is more x10 slower, recommended install module gmpy2.

how to install gmpy2 on windows
 1) download file .whl from https://www.lfd.uci.edu/~gohlke/pythonlibs/
 2) put to root dir of python
 3) install as:
 
Python27>python -m pip install gmpy2-2.0.8-cp27-cp27m-win_amd64.whl

Python37>python -m pip install gmpy2-2.0.8-cp37-cp37m-win_amd64.whl

# Benchmark libs
Algo: 1 Tame + 1 Wild with distinguished points,  expected of 2w<sup>1/2</sup> group operations

1core i5-2540, win7x64, python2x64
```
[lib#raw]		 7.1K j/s
[lib#coincurve]		24.5K j/s
[lib#coincurve+cffi]	35.4K j/s
[lib#gmpy2]		89.4K j/s
```

1core i5-2540, win7x64, python3x64
```
[lib#raw]		 8.9K j/s
[lib#coincurve]		31.2K j/s
[lib#coincurve+cffi]	43.2K j/s
[lib#gmpy2]		97.8K j/s
```

# Heuristic prognose
Algo: 1 Tame + 1 Wild with distinguished points,  expected of 2w<sup>1/2</sup> group operations

avg stat after 100tests

1core i5-2540, python37x64 + gmpy2, win7x64

[i] 97.8K j/s;..

|   W   |jump avg/2w^(1/2)| time                         avg/2w^(1/2)                         |
|:-----:|:---------------:|:-----------------------------------------------------------------:|
| 2^020 |    1.8K/   2.0K |              0d 00:00:00s 018ms /              0d 00:00:00s 020ms |
| 2^030 |   57.3K/  65.5K |              0d 00:00:00s 586ms /              0d 00:00:00s 670ms |
|>2^031 |   81.1K/  92.7K |              0d 00:00:00s 829ms /              0d 00:00:00s 948ms |
| 2^032 |  114.6K/ 131.1K |              0d 00:00:01s 173ms /              0d 00:00:01s 341ms |
| 2^040 |    1.8M/   2.1M |              0d 00:00:18s 777ms /              0d 00:00:21s 471ms |
| 2^050 |   58.7M/  67.1M |              0d 00:10:00s 886ms /              0d 00:11:27s 086ms |
| 2^060 |    1.9G/   2.1G |              0d 05:20:28s 355ms /              0d 06:06:26s 759ms |
| 2^070 |   60.1G/  68.7G |              7d 02:55:07s 378ms /              8d 03:26:16s 292ms |
| 2^080 |    1.9T/   2.2T |          7m 17d 21:23:56s 096ms /          8m 20d 14:00:41s 355ms |
| 2^090 |   61.5T/  70.4T |     20y  3m  2d 12:45:55s 095ms /     23y  1m 28d 16:22:03s 390ms |
| 2^100 |    2.0P/   2.3P |    648y  2m 21d 00:29:23s 067ms /    741y  2m 17d 19:45:48s 481ms |
| 2^105 |   11.1P/  12.7P |   3.7Ky 10m 29d 06:43:07s 581ms /   4.2Ky 11m 12d 16:12:57s 514ms |
| 2^110 |   63.0P/  72.1P |  20.7Ky  2m 12d 15:40:18s 166ms /  23.7Ky 11m  0d 08:25:51s 416ms |
| 2^120 |    2.0E/   2.3E | 663.8Ky  5m 14d 21:29:41s 312ms / 759.0Ky  4m 11d 05:47:25s 328ms |


# Discussion
more info

https://bitcointalk.org/index.php?topic=5173445.msg52473992#msg52473992

https://bitcointalk.org/index.php?topic=5166284.msg52318676#msg52318676

# How pollard-kangaroo works, the Tame and Wild kangaroos, is a simple explanation.

Suppose there is pubkeyX, unknow privkeyX, but privkeyX is in range w=[L..U]. 
The keys have a property - if we increase pubkey by S, then its privkey will also increase by S. 
We start step-by-step to increase pubkeyX by S(i), keeping sumX(S(i)). This is a Wild kangaroo. 
We select a random privkeyT from range [L..U], compute pubkeyT. 
We start step-by-step to increment pubkeyT by S(i) while maintaining sumT(S(i)). This is a Tame kangaroo. 
The size of the jump S(i) is determined by the x coordinate of the current point, so if a Wild or Tame kangaroo lands on one point, their paths will merge. 
(we are concerned with pseudo random walks whose next step is determined by the current position) 
Thanks to the Birthday Paradox (Kruskal's card trick), their paths will one day meet. 
Knowing each traveled path (sumX and sumT), privkeyX is calculated. 
The number of jumps is approximately 2w<sup>1/2</sup> group operations, which is much less than a full search w. 

# Articles

base
https://andrea.corbellini.name/2015/05/17/elliptic-curve-cryptography-a-gentle-introduction/

hand translate to ru/ua (recommend instead translate.google)
https://habr.com/ru/post/335906/


- [1] Best of the best, all in 1, epic,  2012

Chapter 14. Factoring and Discrete Logarithms using Pseudorandom Walks 
https://www.math.auckland.ac.nz/~sgal018/crypto-book/ch14.pdf

- [2] with reference to old

J. M. Pollard, “Kangaroos, monopoly and discrete logarithms,” Journal of Cryptology, #13 (2000) 
https://web.northeastern.edu/seigen/11Magic/KruskalsCount/PollardKangarooMonopoly.pdf

(good dir web.northeastern.edu/seigen/11Magic/KruskalsCount/)

- [3] About parallelism problems

P. C. van Oorschot and M. J. Wiener, Parallel collision search with cryptanalytic applications, J. Cryptology, #12 (1999) 
https://people.scs.carleton.ca/~paulv/papers/JoC97.pdf
