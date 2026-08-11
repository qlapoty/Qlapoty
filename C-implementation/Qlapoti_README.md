# SQIsign C implementation using Qlapoti


This library is a C implementation of SQIsign, modified to use the new Qlapoti algorithm for the ideal-to-isogeny translation step. 

---------

## Differences to the SQIsign NIST round 2 submission

Differences compared to the state of [the SQIsign team's repo](https://www.github.com/SQIsign/the-sqisign) accessed on April 30, 2025.

### Ideal to isogeny 

- Replaced the implementation of the function `dim2id2iso_ideal_to_isogeny_clapotis` in file `src/id2iso/ref/lvlx/dim2id2so.c` by a new Qlapoti-based version and removing the dependencies of the old version. This results in an almost full rewrite off the file.
- In folder `src/id2iso/ref/lvlx/test` adaptation of the tests in `dim2id2iso_test.c`, `test_dim2id2iso.c` and `test_id2iso_test.c`, addition of benchmarking files `dim2id2iso_benchmarks.c` and `qlapoti_normeq_benchmarks.c`
- Adaptation of `src/id2iso/ref/lvlx/lvlx_test.cmake` to compile the new benchmarking files.

### Precomputations

- Addition of a new constant for `QUAT_cornacchia_extended_params` in `src/id2iso/ref/include/id2iso.h`.
- In the `scripts/precomp` folder, adapt the files `parameters.py`, `precompute_quaternion_constants.sage` and `precompute_quaternion_data.sage` to remove the computation of additional orders and add the computation of parameters for the extended cornacchia algoritm.
- In the `src/precomp folder`, overwrite all files with the output of the scripts in scripts/precomp by running `make precomp`

### Quaternions

- In the `src/quaternion/ref/generic` folder, add a file `qlapoti.c` containing the norm equation solving Qlapoti algorithm. Add tests for all functions in this file in a file in `test/qlapoti.c` in the same folder, and make sure these tests are called from `test_quaternions.c`. Adapt `CMakeLists.txt` in this folder and the `test/` folder it contains to compile these files. Add the public-facing `quat_qlapoti` function to `quaternion.h`, as well as the type definition for `ibz_cornacchia_extended_params_t` it requires and the headers of a few functions useful for testing in `id2iso` (`quat_alg_elem_set`, `quat_lideal_equals`, `quat_lideal_mul`, `quat_alg_elem_equal`).
- In the files `intbig.h` and `intbig.c`, the argument to `ibz_two_adic` was made constant.
- Qlapoti needs some additional lower-level functions, which are added in the files where they fit best within the `src/quaternion/ref/generic` folder. For each of them, tests were added in the corresponding file in the `src/quaternion/ref/generic/test/` folder, and headers in one of the header files in  `src/quaternion/ref/internal_quaternion_headers/`, if not in `quaternion.h`.
  - `dim2.c ibz_2x2_mul`, `ibz_mat_2x2_inv_with_det_as_denom`
  - `intbig.c` `ibz_cornacchia_extended` and its dependencies `ibz_cornacchia_extended_prime_loop`, `ibz_complex_mul_by_complex_power`, `ibz_complex_mul`
  - The file `rationals.c` was moved from the hnf subfolder to the `quaternion/ref/generic` base folder, and all headers, inclusions and CMakeLists adapted accordingly.

### Others

- For compiling, one line in `.cmake/target` was changed to "if (UNIX AND NOT APPLE)" from "if (UNIX)".
- This README is entirely new, and the README of the SQIsign NIST round 2 implementation is copied to `SQIsign_README.md`


---------

## Replicating our experimental results

### First setup and compilation

- Use a machine meeting the requirements of SQIsign's Round 2 NIST submission as stated in the requirements section of the `SQIsign_README.md` file.
- Create a folder `build/` inside C-implementation
- Inside `build/`, run `cmake -DSQISIGN_BUILD_TYPE=ref -DCMAKE_BUILD_TYPE=Release ..` (optionally choose your C compiler by the flag `-DCMAKE_C_COMPILER`), then run `make`. No errors nor warnings should show.
- For rerunning the precomputations, run `make precomp` in the build folder (after `cmake` and `make`). This requires a very recent version of SageMath (10.5 for example). To use your new precomp files, re-run `make` afterwards.
- To use or compile the code in other ways, please use the SQIsign NIST Round 2 README, provided in the `SQIsign_README.md` file.

### Replicating our measures of Qlapoti

- For benchmarks of the full SQIsign signature (as in Table 7 of the [paper](https://eprint.iacr.org/2025/1604.pdf)), go into `build/apps` and run `./benchmark_lvl1 --iterations=<number of iterations>`. Change lvl1 to lvl3 or lvl5 for the other levels.
- For heap memory usage (as in Table 8 of the [paper](https://eprint.iacr.org/2025/1604.pdf)), go into `build/test` and then run `valgrind --tool=massif ./sqisign_test_scheme_lvl1`. Visualize the result by calling `ms_print` on the output file `massif.out.<process-id>`. For averaging, it is recommended to script generating and parsing the outputs. Change lvl1 to lvl3 or lvl5 for the other levels.
- For the executable sizes reported in Table 9 of the [paper](https://eprint.iacr.org/2025/1604.pdf), measure the disk space used by the files `sqisign_test_scheme_lvl1`, `sqisign_test_scheme_lvl3` and `sqisign_test_scheme_lvl5` in `build/test`.

### Other benchmarks available

- For benchmarks of `idiso` only, go into `build/src/id2iso/ref/lvl1/test` and run `./sqisign_id2iso_benchmark_dim2id2iso_lvl1 --iterations=<number of iterations>`. Change lvl1 to lvl3 or lvl5 for the other levels.
- We also provide a tool for benchmarking the equation solving part of Qlapoti separately from the isogeny computations. This can be done by running `./sqisign_id2iso_benchmark_qlapoti_normeq_lvl1 --iterations=<number of iterations>` in `build/src/id2iso/ref/lvl1/test`. Change lvl1 to lvl3 or lvl5 for the other levels.

### Comparing to SQIsign

For comparison, follow the exact same procedure in the `SQIsign/the-sqisign` repo's version from April 2025. The files `sqisign_id2iso_benchmark_qlapoti_normeq` and `sqisign_id2iso_benchmark_dim2id2iso_lvl1` do not exist in that version. While the former is meaningless in that context, the latter can be copied into that codebase and made to compile with only minimal changes to some id2iso and quaternion files (essentially modifying `lvlx_test.cmake` in the `src/id2iso/ref/lvlx/` folder to compile the new file, and make the function it depends on available to it, either by copying there headers into public header files if they are not yet there, or by copying them from the Qlapoti code).
