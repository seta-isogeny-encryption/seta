# Séta

This code implements the isogeny-based scheme Séta.

(C) 2021, The Séta team. MIT license.

## Dependencies

The code depends on the latest stable version of the [PARI/GP
library](http://pari.math.u-bordeaux.fr/), 2.12.0, and on
[GMP](https://gmplib.org/) 6.1.2, which is also an optional dependency
of PARI/GP and is typically installed along with it.

## Supported platforms

The code compiles and runs on 64-bit Linux and MacOS.

## Compile

To compile and test the code, run

```
make
make check
```

The tests typically take 5-10 minutes.

## Benchmark

First build the benchmark

```
make benchs
```

then run it

```
./build/bench_seta
```

WARNING: the benchmark takes more than 10 hours!
