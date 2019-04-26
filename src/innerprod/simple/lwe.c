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

#include "innerprod/simple/lwe.h"

#include "internal/prime.h"
#include "sample/normal_double.h"
#include "sample/uniform.h"

// Calculates the center function t(x) = floor(x*q/p) % q for a vector x.
void center(cfe_lwe *s, cfe_vec *t, cfe_vec *x) {
    mpz_t t_i, x_i;
    mpz_inits(t_i, x_i, NULL);

    for (size_t i = 0; i < t->size; i++) {
        cfe_vec_get(x_i, x, i);
        mpz_mul(t_i, x_i, s->q);
        mpz_fdiv_q(t_i, t_i, s->p);
        mpz_mod(t_i, t_i, s->q);
        cfe_vec_set(t, t_i, i);
    }

    mpz_clears(t_i, x_i, NULL);
}

// Initializes scheme struct with the desired configuration
// and configures public parameters for the scheme.
cfe_error cfe_lwe_init(cfe_lwe *s, size_t l, mpz_t bound_x, mpz_t bound_y, size_t n) {
    cfe_error err = CFE_ERR_NONE;
    clock_t start = clock();

    size_t sec = 128;
    size_t i = 10;

    double b = sec / 0.265;
    double delta = pow(pow(3.14159 * b, 1/b) * b / (2 * 3.14159 * 2.71828), 1 / (2 * b - 2));

    s->l = l;
    s->n = n;
    s->A.mat = NULL;
    mpz_init_set(s->bound_x, bound_x);
    mpz_init_set(s->bound_y, bound_y);

    mpf_t p_f, q_f, bound_x_f, bound_y_f, sigma_prime, tmp, sigma, one;
    mpf_inits(s->sigma_q, p_f, q_f, bound_x_f, bound_y_f, sigma_prime, tmp, sigma, one, NULL);

    mpz_t l_z;
    mpz_inits(s->p, s->q, l_z, NULL);
    mpz_set_ui(l_z, l);

    bool safe = false;
    while (!safe) {
        if (cfe_get_prime(s->q, i, false)) {
            cfe_lwe_free(s);
            err = CFE_ERR_PRIME_GEN_FAILED;
            goto cleanup;
        }

        s->m = (n + l + 1) * i + 2 * sec;
        mpz_mul(s->p, bound_x, bound_y);
        mpz_mul(s->p, s->p, l_z);

        mpf_set_z(p_f, s->p);
        mpf_set_z(q_f, s->q);
        mpf_set_z(bound_x_f, bound_x);
        mpf_set_z(bound_y_f, bound_y);
        mpf_div(sigma, q_f, p_f);
        mpf_div(sigma, sigma, bound_y_f);
        mpf_set_d(tmp, 4 * sqrt(l * s->m * sec));
        mpf_div(sigma, sigma, tmp);
        mpf_set_d(tmp, sqrt(l));
        mpf_div(sigma_prime, sigma, tmp);
        mpf_div(sigma_prime, sigma_prime, bound_x_f);

        double sigma_prime_f = mpf_get_d(sigma_prime);
        double q_f_f = mpf_get_d(q_f);
        safe = true;

        for (size_t m_test = n; m_test < s->m; m_test++) {
            double d = (double) n + (double) m_test;
            double left = sigma_prime_f * sqrt(b);
            double right = pow(delta, 2 * b - d - 1) * pow(q_f_f, ((double) m_test) / d );

            if (left < right) {
                safe = false;
                break;
            }
        }
        
        i++;
    }

    mpf_trunc(sigma, sigma);
    mpf_set_ui(one, 1);
    mpf_add(sigma, sigma, one);
    mpf_set(s->sigma_q, sigma);

    // Create a random m*n matrix A
    // The matrix is a public parameter of the scheme.
    cfe_mat_init(&s->A, s->m, s->n);
    cfe_uniform_sample_mat(&s->A, s->q);

    cleanup:
    mpz_clears(l_z, NULL);
    mpf_clears(p_f, q_f, bound_x_f, bound_y_f, sigma_prime, tmp, sigma, one, NULL);
    clock_t end = clock();
    double time = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("\nLWE_Init %.6f\n", time);

    return err;
}

void cfe_lwe_sec_key_init(cfe_mat *SK, cfe_lwe *s) {
    cfe_mat_init(SK, s->n, s->l);
}

void cfe_lwe_pub_key_init(cfe_mat *PK, cfe_lwe *s) {
    cfe_mat_init(PK, s->m, s->l);
}

// Generates a secret key for the scheme.
// The key is represented by a matrix with dimensions n*l whose
// elements are random values from the interval [0, q).
void cfe_lwe_generate_sec_key(cfe_mat *SK, cfe_lwe *s) {
    clock_t start = clock();
    cfe_uniform_sample_mat(SK, s->q);
    clock_t end = clock();
    double time = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("\nLWE_KeyGen1 %.6f\n", time);

}

// Generates a public key for the scheme.
// Public key is a matrix of m*l elements.
cfe_error cfe_lwe_generate_pub_key(cfe_mat *PK, cfe_lwe *s, cfe_mat *SK) {
    clock_t start = clock();

    if (SK->rows != s->n || SK->cols != s->l) {
        return CFE_ERR_MALFORMED_SEC_KEY;
    }

    cfe_error err = CFE_ERR_NONE;

    // Initialize and fill noise matrix E with m*l samples
    cfe_mat E;
    cfe_mat_init(&E, s->m, s->l);

    cfe_normal_double sampler;
    mpf_t one;
    mpf_init_set_ui(one, 1);

    if (cfe_normal_double_init(&sampler, s->sigma_q, s->n, one)) {
        err = CFE_ERR_PUB_KEY_GEN_FAILED;
        goto cleanup;
    }
    cfe_normal_double_sample_mat(&E, &sampler);

    // Calculate public key as PK = (A * SK + E) % q
    cfe_mat_mul(PK, &s->A, SK);
    cfe_mat_mod(PK, PK, s->q);
    cfe_mat_add(PK, PK, &E);
    cfe_mat_mod(PK, PK, s->q);

    cleanup:
    cfe_mat_free(&E);
    mpf_clear(one);
    cfe_normal_double_free(&sampler);
    clock_t end = clock();
    double time = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("\nLWE_KeyGen2 %.6f\n", time);

    return err;
}

void cfe_lwe_fe_key_init(cfe_vec *sk_y, cfe_lwe *s) {
    cfe_vec_init(sk_y, s->n);
}

// Derives a secret key sk_y for decryption of inner product of y and
// a secret operand.
// Secret key is a linear combination of input vector y and master secret key.
cfe_error cfe_lwe_derive_key(cfe_vec *sk_y, cfe_lwe *s, cfe_mat *SK, cfe_vec *y) {
    clock_t start = clock();
    if (!cfe_vec_check_bound(y, s->bound_y)) {
        return CFE_ERR_BOUND_CHECK_FAILED;
    }
    if (SK->rows != s->n || SK->cols != s->l) {
        return CFE_ERR_MALFORMED_SEC_KEY;
    }

    cfe_mat_mul_vec(sk_y, SK, y);
    cfe_vec_mod(sk_y, sk_y, s->q);

    clock_t end = clock();
    double time = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("\nLWE_KeyDerive %.6f\n", time);

    return CFE_ERR_NONE;
}

void cfe_lwe_ciphertext_init(cfe_vec *ct, cfe_lwe *s) {
    cfe_vec_init(ct, s->n + s->l);
}

// Encrypts vector x using public key PK.
// The resulting ciphertext is stored in vector ct.
cfe_error cfe_lwe_encrypt(cfe_vec *ct, cfe_lwe *s, cfe_vec *x, cfe_mat *PK) {
    clock_t start = clock();

    if (!cfe_vec_check_bound(x, s->bound_x)) {
        return CFE_ERR_BOUND_CHECK_FAILED;
    }
    if (PK->rows != s->m || PK->cols != s->l) {
        return CFE_ERR_MALFORMED_PUB_KEY;
    }
    if (x->size != s->l) {
        return CFE_ERR_MALFORMED_INPUT;
    }

    // A vector comprising the first n elements of the cipher
    cfe_vec ct_0;
    cfe_vec_init(&ct_0, s->n);

    // A vector comprising the last l elements of the cipher
    cfe_vec ct_last;
    cfe_vec_init(&ct_last, s->l);

    // Create a random vector comprised of m 0s and 1s
    mpz_t two;
    mpz_init_set_ui(two, 2);
    cfe_vec r;
    cfe_vec_init(&r, s->m);
    cfe_uniform_sample_vec(&r, two);

    // Obtain first n elements of the cipher as
    // ct_0 = A_transposed * r
    cfe_mat A_t;
    cfe_mat_init(&A_t, s->n, s->m);
    cfe_mat_transpose(&A_t, &s->A);
    cfe_mat_mul_vec(&ct_0, &A_t, &r);

    // Calculate coordinates ct_last_i = <pk_i, r> + t(x_i) mod q
    // We can obtain the vector of dot products <pk_i, r> as PK
    // transposed * r
    // Function t(x) is denoted as the center function
    cfe_mat PK_t;
    cfe_mat_init(&PK_t, s->l, s->m);
    cfe_mat_transpose(&PK_t, PK);
    cfe_mat_mul_vec(&ct_last, &PK_t, &r);

    cfe_vec t;
    cfe_vec_init(&t, s->l);
    center(s, &t, x);
    cfe_vec_add(&ct_last, &ct_last, &t);
    cfe_vec_mod(&ct_last, &ct_last, s->q);

    // Construct the final ciphertext vector by joining both parts
    cfe_vec_join(ct, &ct_0, &ct_last);

    // Cleanup
    mpz_clear(two);
    cfe_vec_frees(&t, &ct_last, &ct_0, &r, NULL);
    cfe_mat_frees(&A_t, &PK_t, NULL);

    clock_t end = clock();
    double time = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("\nLWE_Encrypt %.6f\n", time);

    return CFE_ERR_NONE;
}

// Decrypts the ciphertext ct.
// res will hold the decrypted inner product <x,y>
// sk_y is the derived secret key for decryption of <x,y>
// y is plaintext input vector
cfe_error cfe_lwe_decrypt(mpz_t res, cfe_lwe *s, cfe_vec *ct, cfe_vec *sk_y, cfe_vec *y) {
    clock_t start = clock();
    if (!cfe_vec_check_bound(y, s->bound_y)) {
        return CFE_ERR_BOUND_CHECK_FAILED;
    }
    if (sk_y->size != s->n) {
        return CFE_ERR_MALFORMED_FE_KEY;
    }
    if (y->size != s->l) {
        return CFE_ERR_MALFORMED_INPUT;
    }
    if (ct->size != s->n + s->l) {
        return CFE_ERR_MALFORMED_CIPHER;
    }

    // Break down the ciphertext vector into
    // ct_0     which holds first n elements of the cipher, and
    // ct_last  which holds last n elements of the cipher
    cfe_vec ct_0, ct_last;
    cfe_vec_init(&ct_0, s->n);
    cfe_vec_init(&ct_last, s->l);
    cfe_vec_extract(&ct_0, ct, 0, s->n);
    cfe_vec_extract(&ct_last, ct, s->n, s->l);

    // Calculate d = <y, ct_last> - <ct_0, sk_y>
    mpz_t d;    // will hold the decrypted message
    mpz_t prod; // temporary variable for holding dot products
    mpz_inits(d, prod, NULL);

    cfe_vec_dot(prod, y, &ct_last);
    mpz_mod(d, prod, s->q);

    cfe_vec_dot(prod, &ct_0, sk_y);
    mpz_mod(prod, prod, s->q);

    mpz_sub(d, d, prod);
    mpz_mod(d, d, s->q);

    // Return the plaintext res, where res is such that
    // d - center(m) % q is closest to 0.
    // half_q = floor(q/2)
    mpz_t half_q;
    mpz_init(half_q);
    mpz_fdiv_q_ui(half_q, s->q, 2);

    if (mpz_cmp(d, half_q) == 1) {
        mpz_sub(d, d, s->q);
    }

    mpz_mul(d, d, s->p);
    mpz_add(d, d, half_q);
    mpz_fdiv_q(d, d, s->q);

    mpz_set(res, d); // set the value of decrypted message as the result

    // Cleanup
    mpz_clears(d, prod, half_q, NULL);
    cfe_vec_frees(&ct_0, &ct_last, NULL);
    clock_t end = clock();
    double time = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("\nLWE_Decrypt %.6f\n", time);

    return CFE_ERR_NONE;
}

// Frees the memory allocated for configuration of the scheme.
void cfe_lwe_free(cfe_lwe *s) {
    mpz_clears(s->p, s->q, s->bound_x, s->bound_y, NULL);
    mpf_clear(s->sigma_q);

    if (s->A.mat != NULL) {
        cfe_mat_free(&s->A);
    }
}
