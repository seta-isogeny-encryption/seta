CC=clang
CFLAGS=-Wall -Wextra -pedantic -std=gnu99 -I./include -I/usr/local/include
CFLAGSLINK=-lpari -lm -L/usr/local/lib/ 
DEBUG_FLAGS=-g
BENCH_FLAGS=-DNDEBUG -O3 -Os -march=native -mtune=native

default: all

# The implemented primes
PRIME_LIST=p6983 p6334792777 p8426067021
PRIME?=p8426067021

# "Library" object files depended upon by all executables
LIBFP  = fp2.o constants.o precomputed.o steps_tunecycles.o fp.o curve.o
LIBCOM = toolbox.o ideal.o klpt.o idiso.o rng.o poly.o mont.o tedwards.o steps.o uint.o\
	isogenies.o isomorphism.o two_walks.o mitm.o seta.o qfsolve_factor.o #verif.o
LIB = $(LIBFP) $(LIBCOM)

build/obj/$(PRIME)/precomputed.o: src/$(PRIME)/precomputed_generated.h

$(LIBFP:%=build/obj/$(PRIME)/%): build/obj/%.o: src/%.c
	@mkdir -p $(@D)
	$(CC) $< $(CFLAGS) $(DEBUG_FLAGS) -c -o $@

$(LIBCOM:%=build/obj/%): build/obj/%.o: src/%.c
	@mkdir -p $(@D)
	$(CC) $< $(CFLAGS) $(DEBUG_FLAGS) -c -o $@

# Search path for building executables
vpath %.c src:src/$(PRIME)
vpath %.o build/obj:build/obj/$(PRIME)

# Benchmarks
BENCHS=$(patsubst bench/%.c,build/bench_%,$(wildcard bench/*.c))

$(BENCHS): build/bench_%: bench/%.c $(LIB:%.o=%.c)
	@mkdir -p $(@D)
	$(CC) $^ $(CFLAGS) $(BENCH_FLAGS) $(CFLAGSLINK) -lgmp -o $@

# Tests
TESTS=$(patsubst test/%.c,build/test_%,$(wildcard test/*.c))

$(TESTS): build/test_%: test/%.c $(LIB)
	@mkdir -p $(@D)
	$(CC) $^ $(CFLAGS) $(DEBUG_FLAGS) $(CFLAGSLINK) -lgmp -o $@

# Precomputed values
EXES=build/precomp

build/precomp: precomp.c fp2.o constants.o fp.o curve.o poly.o mont.o tedwards.o uint.o steps.o rng.o steps_tunecycles.o
	@mkdir -p $(@D)
	$(CC) $^ $(CFLAGS) $(DEBUG_FLAGS) $(CFLAGSLINK) -lgmp -o $@

src/$(PRIME)/precomputed_generated.h: | build/precomp
	$(MAKE) precompute

precompute:
	./build/precomp > src/$(PRIME)/precomputed_generated.h


# Velusqrt Tuning

build/tunecycles_$(PRIME): tunecycles.c isogenies.c mont.c fp2.c uint.c fp.c constants.c rng.c poly.c steps.c steps_default.c
	@mkdir -p $(@D)
	$(CC) $^ $(CFLAGS) $(BENCH_FLAGS) $(CFLAGSLINK) -o $@

src/$(PRIME)/tunecycles.out: build/tunecycles_$(PRIME)
	# 8 minutes on 1.9GHz Kaby Lake
	time ./$< > $@

tune: src/tune2c src/$(PRIME)/tunecycles.out
	./src/tune2c < src/$(PRIME)/tunecycles.out > src/$(PRIME)/steps_tunecycles.c



# Run tests
$(TESTS:%=%_run): %_run: %
	@echo
	./$^

check: $(TESTS:%=%_run)

# Run benchmarks
$(BENCHS:build/%=%.tsv): %.tsv: build/%
	./$^ >> $@

benchmark: $(BENCHS:build/%=%.tsv)

# Phony targets
tests: $(TESTS)
benchs: $(BENCHS)
lib: tests benchs
all: $(EXES) lib

distclean:
	rm -r build

.PHONY: distclean all lib tests benchs benchmark check tune default precompute
