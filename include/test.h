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

#ifndef CIFER_TEST_H
#define CIFER_TEST_H

#include <munit.h>
#include <gmp.h>

MunitSuite ddh_suite;
MunitSuite damgard_suite;
MunitSuite ddh_multi_suite;
MunitSuite damgard_multi_suite;
MunitSuite lwe_suite;
MunitSuite lwe_fully_secure_suite;
MunitSuite paillier_suite;
MunitSuite policy_suite;
MunitSuite gpsw_suite;
MunitSuite fame_suite;

size_t n_bench;
size_t l_bench;
size_t b_bench;
mpz_t bound_bench;

#endif
