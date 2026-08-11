# SQIsign C implementation using Qlapoty


This library is a C implementation of SQIsign, modified to use the new Qlapoty algorithm for the ideal-to-isogeny translation step. 

It was adapted from [the Qlapoti repo](https://www.github.com/KULeuven-COSIC/Qlapoti) accessed on Mai 18, 2026.

---------

## Differences to the the C-implementation of Qlapoti

Differences compared to the state of [the Qlapoti repo](https://www.github.com/KULeuven-COSIC/Qlapoti) accessed on Mai 18, 2026.

### Ideal to isogeny 

- Adapted the benchmarks file `test/qlapoti_normeq_benchmarks.c` to benchmark Qlapoty's norm equation solver by default and Qlapoti's if the argument `--qlapoti` is provided. Also add some cpu warmup iterations to the benchmark function.
- Make the function `dim2id2iso_ideal_to_isogeny_clapotis` in file `dim2id2so.c` use Qlapoty's norm equation solver instead of Qlapoti's.
- Remove the use of Weil pairings, and replace where still needed by Tate pairings. This affects `test/dim2id2iso_benchmarks.c`, `id2iso.h` and `dim2id2so.c`. 

### Quaternions

- Move the `ibz_xgcd` function from a dedicated file in the `hnf` subfolder to the `intbig.c` file and adapt headers and cmake files accordingly.
- Add the Qlapoty norm equation solver and its subfunctions in the file `qlapoti.c`, and tests for them in `test/qlapoti.c`.
- Add the Qlapoty norm equation solver to the header file `quaternion.h` and some of its subfunctions to `internal_quaternion_headers/internal.h`.

### Others

- This README is entirely new, and the README of the Qlapoti is copied to `Qlapoti_README.md`


## Difference to the SQIsign NIST round 3 submission

The difference to the SQIsign round 2 NIST submission is the above difference from our code to Qlapoti added to the difference between Qlapoti and the SQisign round 2 NIST submission. The latter difference is described in the readme of the C-implementation of Qlapoti which was renamed to `Qlapoti_README.md` here.

---------

## Replicating our experimental results

### First setup and compilation

- Use a machine meeting the requirements of SQIsign's Round 2 NIST submission as stated in the requirements section of the `SQIsign_README.md` file.
- Create a folder `build/` inside C-implementation
- Inside `build/`, run `cmake -DSQISIGN_BUILD_TYPE=ref -DCMAKE_BUILD_TYPE=Release ..` (optionally choose your C compiler by the flag `-DCMAKE_C_COMPILER`), then run `make`. No errors nor warnings should show.
- For a full test of all SQIsign subfunctions, run `make test`. It is expected that the `KAT`tests fail, as qlapoty changes the way randomness is used in the scheme compared to the NIST round 2 submission (these tests also already fail in the `KULeuven-COSIC/Qlapoti` repository).
- For rerunning the precomputations, run `make precomp` in the build folder (after `cmake` and `make`). This requires a very recent version of SageMath (10.5 for example). To use the new precomp files, re-run `make` afterwards.
- To use or compile the code in other ways, please use the SQIsign NIST Round 2 README, provided in the `SQIsign_README.md` file.

### Replicating measures

- For benchmarks of the full SQIsign signature, go into `build/apps` and run `./benchmark_lvl1 --iterations=<number of iterations>`. Change lvl1 to lvl3 or lvl5 for the other levels.
- For benchmarking the equation solving part of Qlapoty and Qlapoti separately from the isogeny computations, run `./sqisign_id2iso_benchmark_qlapoti_normeq_lvl1 --iterations=<number of iterations>` in `build/src/id2iso/ref/lvl1/test`. Change lvl1 to lvl3 or lvl5 for the other levels. Add the argument `--qlapoti`to benchmark Qlapoti's norm equation solver instead of Qlapoty's.

### Other benchmarks available

- For benchmarks of `idiso` only, go into `build/src/id2iso/ref/lvl1/test` and run `./sqisign_id2iso_benchmark_dim2id2iso_lvl1 --iterations=<number of iterations>`. Change lvl1 to lvl3 or lvl5 for the other levels.

### Comparing to Qlapoti and SQIsign

For comparison to Qlapoti norm equation on the norm equation benchmarks alone, it is sufficient to pass the argument `--qlapoti` to `sqisign_id2iso_benchmark_qlapoti_normeq_lvl1`. For comparison to Qlapoti's ideal-to-isogeny (without our removal of pairings) and its impact on SQIsign, clone the `KULeuven-COSIC/Qlapoti` repository in github and run the exact same benchmark files there and in this repository. 

For comparison full scheme or signing benchmraks to SQIsign's round 2 NIST submission, clone the `SQIsign/the-sqisign` repo's version from April 2025 and run the full signature benchmarks there by following the same steps as in this repo. Benchmark files for the subroutines do not exist in that repository.
