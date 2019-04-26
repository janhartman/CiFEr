/*
 * Copyright (c) 2018 XLAB d.o.o.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <stdio.h>
#include <stdlib.h>
#include "test.h"
#include "internal/common.h"

int main(int argc, char *argv[]) {
    if (cfe_init()) {
        perror("Insufficient entropy available for random generation\n");
        return CFE_ERR_INIT;
    }

    MunitSuite ddh_suites[] = {
            //ddh_suite,
            //ddh_multi_suite,
            lwe_suite,
            //damgard_suite,
            //damgard_multi_suite,
            //lwe_fully_secure_suite,
            //paillier_suite,
            {NULL, NULL, NULL, 0, MUNIT_SUITE_OPTION_NONE}
    };

    MunitSuite abe_suites[] = {
            gpsw_suite,
            fame_suite,
            {NULL, NULL, NULL, 0, MUNIT_SUITE_OPTION_NONE}
    };

    if (argc < 2) {
        printf("no arguments passed\n");
        return 1;
    }

    int err = 0;

    if (argv[1][0] == 'd') {
        if (argc < 5) {
            printf("not enough arguments\n");
            return 1;
        }

        n_bench = strtol(argv[2], NULL, 10);
        l_bench = strtol(argv[3], NULL, 10);
        b_bench = strtol(argv[4], NULL, 10);
        mpz_init_set_ui(bound_bench, b_bench);

        MunitSuite main_suite = {(char *) "", NULL, ddh_suites, 1, MUNIT_SUITE_OPTION_NONE};
        err = munit_suite_main(&main_suite, NULL, 0, argv);
        mpz_clear(bound_bench);
    } else if (argv[1][0] == 'a') {
        if (argc < 3) {
            printf("not enough arguments\n");
            return 1;
        }

        l_bench = strtol(argv[2], NULL, 10);

        MunitSuite main_suite = {(char *) "", NULL, abe_suites, 1, MUNIT_SUITE_OPTION_NONE};
        err = munit_suite_main(&main_suite, NULL, 0, argv);
    } else {
        printf("invalid argument: %s\n", argv[1]);
    }

    return err;
}
