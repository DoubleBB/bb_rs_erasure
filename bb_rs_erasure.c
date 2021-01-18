/*
 Reed-Solomon code implementation for erasure coding  v.1.0.3

 Copyright (c) 2019-2021 Bela Bodecs   (bodecsb#vivanet.hu)


   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
 */

#include <bb_rs_erasure.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>


// in GF(2**m) adding is the XORing the two component and adding is the same as substraction
static inline uint8_t gf_add(uint8_t a, uint8_t b) {

  return a^b;
}


// multiply two numbers in GF aritmetic
static inline uint8_t gf_mul(rs_ctx const * const rs, uint8_t a, uint8_t b) {

  return rs->mul_table[a * GF_ORDER + b];

  // this is an exp table based slower version
  //
  //  if (a == 0 || b == 0) return 0;
  //  //exp table is oversized so no need for modulo by  GF_MUL_ORDER on the sum
  //  return rs->exp_table[ rs->log_table[a] + rs->log_table[b] ];
}


// calculating a - b in GF
static inline uint8_t gf_div(rs_ctx const * const rs, uint8_t a, uint8_t b) {


  return rs->div_table[a * GF_ORDER + b];

  // this is an exp table based slower version
  //
  //  if (a == 0 || b == 0)  // b == 0 is an error
  //    return 0;
  //  // we add GF_MUL_ORDER to the sum to avoid negativ result, alpha**GF_MUL_ORDER == alpha,
  //  // by definition because th eorder of alpha is GF_MUL_ORER
  //  // a - b == GF_MUL_ORDER - b + a
  //  int16_t exp = GF_MUL_ORDER - rs->log_table[b] + rs->log_table[a];
  //  //exp table is oversized so no need for modulo by  GF_MUL_ORDER
  //  return rs->exp_table[ exp ];
}



// free up all GF aritmetic helper table allocated earlier by init_gf_tables
static void free_gf_tables(rs_ctx * rs) {
  if (!rs)
    return;

  if (!rs->external_tables) {
    free(rs->tables); // If a null pointer is passed as argument, no action occurs.
  } //if

  rs->tables = NULL;
  rs->exp_table = NULL;
  rs->log_table = NULL;
  rs->div_table = NULL;
  rs->mul_table = NULL;
}


uint8_t set_external_tables(rs_ctx * rs, clone_data_t * external_tables) {

  size_t size_of_gf_tables = (3 * GF_ORDER + 2 * GF_ORDER * GF_ORDER) * sizeof(uint8_t);
  size_t size_of_g = (rs->n - rs->k + 1) * sizeof(uint8_t);
  size_t size_of_generator_matrix = rs->n * rs->k * sizeof(uint8_t);
  size_t size_of_generator_matrix_t = rs->k * rs->n * sizeof(uint8_t);
  size_t size_of_parity_matrix = rs->n * (rs->n - rs->k) * sizeof(uint8_t);

  size_t size_of_data = size_of_gf_tables + size_of_g + size_of_generator_matrix + size_of_generator_matrix_t + size_of_parity_matrix;

  if (external_tables->k != rs->k || external_tables->n != rs->n ||
      external_tables->size_of_data != size_of_data || !external_tables->data)
    return 0;


  uint8_t * p = external_tables->data;
  rs->tables = p;
  p += size_of_gf_tables;
  rs->g = p;
  p += size_of_g;
  rs->generator_matrix = p;
  p += size_of_generator_matrix;
  rs->generator_matrix_t = p;
  p += size_of_generator_matrix_t;
  rs->parity_matrix = p;


  return 1;
}


clone_data_t * clone_internal_tables(rs_ctx const * restrict const rs) {


  if (!rs || !rs->tables ||
      !rs->g || !rs->parity_matrix ||
      !rs->generator_matrix || !rs->generator_matrix_t)
    return NULL;

  clone_data_t * clone_data = NULL;

  clone_data = (clone_data_t *) calloc(1, sizeof(clone_data_t));
  if (!clone_data)
    return NULL;

  clone_data->n = rs->n;
  clone_data->k = rs->k;

  size_t size_of_gf_tables = (3 * GF_ORDER + 2 * GF_ORDER * GF_ORDER) * sizeof(uint8_t);
  size_t size_of_g = (rs->n - rs->k + 1) * sizeof(uint8_t);
  size_t size_of_generator_matrix = rs->n * rs->k * sizeof(uint8_t);
  size_t size_of_generator_matrix_t = rs->k * rs->n * sizeof(uint8_t);
  size_t size_of_parity_matrix = rs->n * (rs->n - rs->k) * sizeof(uint8_t);


  clone_data->data = (uint8_t *) malloc(size_of_gf_tables +
                                        size_of_g +
                                        size_of_generator_matrix +
                                        size_of_generator_matrix_t +
                                        size_of_parity_matrix);
  if (!clone_data->data) {
    free(clone_data);
    return NULL;
  } // if
  clone_data->size_of_data = size_of_gf_tables + size_of_g + size_of_generator_matrix + size_of_generator_matrix_t + size_of_parity_matrix;

  uint8_t * p = clone_data->data;
  memcpy(p, rs->tables, size_of_gf_tables);
  p += size_of_gf_tables;
  memcpy(p, rs->g, size_of_g);
  p += size_of_g;
  memcpy(p, rs->generator_matrix, size_of_generator_matrix);
  p += size_of_generator_matrix;
  memcpy(p, rs->generator_matrix_t, size_of_generator_matrix_t);
  p += size_of_generator_matrix_t;
  memcpy(p, rs->parity_matrix, size_of_parity_matrix);

  return clone_data;
}




// allocate and fill-up all GF aritmetic helper table
static uint8_t init_gf_tables(rs_ctx * rs) {

  size_t exp_table_size = 2 * GF_ORDER * sizeof(uint8_t);
  size_t log_table_size = GF_ORDER * sizeof(uint8_t);
  size_t div_table_size = GF_ORDER * GF_ORDER * sizeof(uint8_t);
  size_t mul_table_size = GF_ORDER * GF_ORDER * sizeof(uint8_t) ;
  size_t size_of_sum_gf_table_memories = exp_table_size + log_table_size + div_table_size + mul_table_size;

  if (!rs->external_tables) {
    rs->tables = (uint8_t *) calloc(1, size_of_sum_gf_table_memories);
  } //  if

  if (!rs->tables)
    return 0;

  // max power value of adding two alphas as primitive element is GF_MUL_ORDER + GF_MUL_ORDER
  // exp_table contains alpha powers
  rs->exp_table = rs->tables; // (uint8_t *) calloc(2 * GF_ORDER, sizeof(uint8_t));
  // exp values range between 0 and GF_MUL_ORDER-1  alpha**GF_MUL_ORDER = 1
  // rs->exp_table[i] == alpha**i mod GF_PRIMITIVE_POLYNOMIAL (i. power of alpha modulo by GF_PRIMITIVE_POLYNOMIAL)
  // we use the alpha == 2 as primitive multiplicative generator element
  // rs->exp_table[0] == 1, rs->exp_table[2] == 2 etc.
  rs->log_table = &rs->tables[exp_table_size]; // (uint8_t *) calloc(GF_ORDER, sizeof(uint8_t));
  rs->div_table = &rs->tables[exp_table_size + log_table_size]; // (uint8_t *) calloc(GF_ORDER * GF_ORDER, sizeof(uint8_t));
  rs->mul_table = &rs->tables[exp_table_size + log_table_size + div_table_size]; // (uint8_t *) calloc(GF_ORDER * GF_ORDER, sizeof(uint8_t));

  if (rs->external_tables)
    return 1; // no need any table calculation


  uint16_t i,j,x;  // helper variables

  // exp_table[0] = 1 because it is alpha**0 (b1), exp_table[1] = b10 because it is alpha**1,
  // exp_table[2] = b100 because it is alpha**2 and so on until exp_table[7],
  // but exp_table[8] == lower bits of primitive polinomial coeffs because using the substitution that prim.poly(alpha) == 0,
  // so alpha**8 == prim.poly(alpha) - alpha**8 (in binary -b == +b)
  x = 1; // we use the alpha == 2 as primitive multiplicative generator element  (in GF(2**F_POWER) alpha**GF_MUL_ORDER == 1)
         // so multiplying by alpha is the same as multiplying by 2 and then modulo by PRIMITIVE_POLYNOMIAL
         // multiplying by 2 is very easy/efficient (left shift by 1 position)
  for(i=0; i<GF_ORDER; i++) {

    // to easily get the answer for question:  "what is the number equals with the i-th power of alpha in GF(2**8)?"
    rs->exp_table[i] = x;


    // to easily get the answer for question:  "which power of alpha is x  in GF(2**8)?"
    rs->log_table[x] = i;

    // pre-calculating x for next table item
    x <<= 1; // multiplyig by 2 == shift 1 bit to left
    if (x & GF_ORDER)
      x ^= GF_PRIMITIVE_POLYNOMIAL; // substituting the 2**8 element with primitive polynom lower elements.
                                    // it also clears the highest bit
                                    // this substitution comes from the fact that
                                    // alpha is root of primitive polynomial by definition
                                    // substituting occurs by adding each bit positions. XOR is the add method in GF(2**8)

    // to easily calculate multiplication of two numbers of GF(2**8) we oversize exp table by 2
    // the mutiplication is the same as adding the power values:  a^x * a^y = a^(x+y)
    // in our case x = 0...7 and y = 0...7
    // but the "upper half" of the matrix is the same as the lower half,
    // because multplicative order means that after this many consecutive powers the generator values will be the same again (cyclic ring)
    // raising generator element to power of MULt.order is the genertor element again. We never can get zero, so MUL.ORDEr is GF_ORDER-1
    // a^MUL_ORDER == a
    // so initialization of the whole matrix is based on calculating only the lower half
    rs->exp_table[i + 1 + GF_MUL_ORDER] = x; // do not forget that x is power of i+1 now
//    printf("(%u:) "BYTE_TO_BINARY_PATTERN" => (%u:) "BYTE_TO_BINARY_PATTERN"\n", i, BYTE_TO_BINARY(rs->exp_table[i]), i+1, BYTE_TO_BINARY(x));
  } // for
  rs->log_table[0] = GF_MUL_ORDER; // this would remain uninitialized because log(0) does not exist
                                   // (we may setting it to any value)


  // to easily calculate multiplication of two numbers of GF(2**8) we oversize exp table by 2
  // but the matrix is simmetrical
//  for(t=GF_ORDER; t<2*GF_ORDER; t++)
//    if (rs->exp_table[t] != rs->exp_table[t - GF_MUL_ORDER]) {
//      printf("oversize error at position %u", t);
//      return 0;
//    } //if

  // a full multiplication table
//  for(i=0; i<=GF_MUL_ORDER; i++)
//    for(j=0; j<=GF_MUL_ORDER; j++)
//      if (i != 0 && j != 0)

  // a full multiplication and division table
  // J == 0 and  i == 0 values remain zero, because  zero times anything is zero
  for(i=1; i<=GF_MUL_ORDER; i++)
    for(j=1; j<=GF_MUL_ORDER; j++) {

      // theoretically need to do mod by GF_MUL_ORDER, but this  is unnecessary here because size of exp_table is 2 * GF_ORDER
      // intentionally to avoid/save this relative expensive mod operation here
      rs->mul_table[i * GF_ORDER + j] = rs->exp_table[ rs->log_table[i] + rs->log_table[j] ];

      // we add GF_MUL_ORDER to avoid negativ result, alpha**GF_MUL_ORDER == alpha,
      // by definition because the order of alpha is GF_MUL_ORDER
      // a - b == GF_MUL_ORDER - b + a
      rs->div_table[i * GF_ORDER + j] = rs->exp_table[ (int16_t)GF_MUL_ORDER - rs->log_table[j] + rs->log_table[i] ];
    } // for

  return 1;
}


//  divides c[n]*x**n + c[n-1]*x**(n-1) + ... c[0]  by (x-a)
//  result placed into c[] coefficient array,  c[n] becames always zero
//  returns the remainder
//  basic idea comes from following examples:
//   Ax + B  / (x-a) == A and remainder B+aA
//   Ax2 + Bx + C / (x-a) == Ax + (B+aA) and remainder:  C+a(B+aA)
//   Ax3 + Bx2 + Cx + D / (x-a) == Ax2 + (B+aA)x + (C+a(B+aA)) and remainder D+a(C+a(B+aA))
//   and so on ... so the result and the remainder may be calculated paralel in n steps
//  parameter n is the highest power in polynom, NOT the size of array!!!!
uint8_t gf_simple_polynom_divison_in_place(rs_ctx * rs, uint8_t * c, uint8_t n, uint8_t a) {
  int16_t i;
  uint8_t prev_remainder;
  uint8_t remainder = c[n];
  c[n] = 0;
  for(i=n-1; i>=0; i--) {
    prev_remainder = c[i];
    c[i] = remainder;
    remainder = (a && remainder) ? gf_add(prev_remainder, gf_mul(rs, remainder, a)) : prev_remainder;
  } // for
  return remainder;
}



// long division method as by pencil and paper in school but with GF aritmetic
// given the n+1 coefficients of a polynomial of degree n in u[0 .. n], and the m+1 coefficients
// of another polynomial of degree m in v[0 .. m], divide the polynomial u by the polynomial v
// resulting a quotient polynomial whose coefficients are returned in q[0 .. n-m], and a
// remainder polynomial whose coefficients are returned in r[0 .. n].
// The elements in r[0 .. m .. n] are always returned as zero because highest non-zero power coefficient
// of r is m-1 or less, but we use r[] for calculations, so its size must be n+1
// m is NOT greater than n
static void gf_polynom_div(rs_ctx const * const rs, uint8_t * u, uint8_t n, uint8_t * v, uint8_t m, uint8_t * q, uint8_t * r) {
  int16_t i,j; // i and j are intentionally signed and bigger than uint8_t

  memcpy(r,u,n+1); // init remainder as full u polynomial

  // the result created in n-m+1 steps
  for(i=n-m; i>=0; i--)
    if (r[m+i]) { // skip this step if r[i+m] is zero

      // calculate quotient of current power coefficients in remainder by dividing the actual coefficients
      q[i] = gf_div(rs, r[m+i], v[m]); // r[m+i] / v[m]

      // substract the q[i]*v from the remainder
      if (q[i]) // skip whole step if q[i] is zero
        for (j=m+i-1; j>=i; j--)
          if (v[j-i]) // skip if v[j-i] is zero
            GF_ADD_INTO(r[j], gf_mul(rs, q[i], v[j-i]));  // r[j] - q[k]*v[j-k]  (addition is the same as substraction)

      r[m+i] = 0; // set this remainder coefficient to zero

    } else
      q[i] = 0;

}


// create a systematic k*n sized R-S Generator Matrix (GM) and its transponse
// basic idea: let i-th row in a k*n sized systematic generator matrix represents something like x**(i+k) + s(x) where highest power of s(x) is n-k
// So calculate s(x) for each k rows as x**i mod g(x) (for i: n-k .... n-1), where g(x) is the R-S generator polynom
// this way each row will be divisible by g(x) and when calculate a c(x) codeword to an u(x) message -
// as c(x) = u(x)*GM() will be codewords because c(x) mod g(x) == 0 and
// c(x) values will be systematic because the Generator Matrix format is (I(k) + R(n-k,k))
//
// these will be the rows of k x n generator matrix: coefficients of x**j + (x**j mod g(x)) (k items and (n-k) items)
//  1. row: x**(n-1) + (x**(n-1) mod g(x))               1 0 0 ... 0  R1(n-k-1) R1(n-k-2) ... R1(0)
//  2. row: x**(n-2) + (x**(n-2) mod g(x))               0 1 0 ... 0  R2(n-k-1) R2(n-k-2) ... R2(0)
//    ...
//  k. row: x**(n-k) + (x**(n-k) mod g(x))               0 0 0 ... 1  Rk(n-k-1) Rk(n-k-2) ... Rk(0)
//
static uint8_t * create_rs_generator_matrix(rs_ctx * rs) {
  uint16_t i,j;
  uint8_t * p = (uint8_t *) calloc(GF_ORDER, sizeof(uint8_t)); // max possible size
  uint8_t * q = (uint8_t *) calloc(GF_ORDER, sizeof(uint8_t)); // max possible size
  uint8_t * r = (uint8_t *) calloc(GF_ORDER, sizeof(uint8_t)); // max possible size

  rs->generator_matrix = (uint8_t *) calloc( rs->k * rs->n, sizeof(uint8_t));
  rs->generator_matrix_t = (uint8_t *) calloc( rs->n * rs->k, sizeof(uint8_t));

  if (!p || !q || !r || !rs->generator_matrix || !rs->generator_matrix_t) {
    free(p);
    free(q);
    free(r);
    free(rs->generator_matrix);
    free(rs->generator_matrix_t);
    rs->generator_matrix = NULL;
    rs->generator_matrix_t = NULL;
    return NULL;
  } // if

  uint8_t row = rs->k;
  for (i=rs->n - rs->k; i<rs->n; i++) {
    memset(p, 0, rs->n * sizeof(uint8_t));
    p[i]=1;
    gf_polynom_div(rs, p, i, rs->g, rs->n - rs->k, q, r);
    row--;
    // copy the result remainder coefficients into the last n-k position of the actual matrix row
    for(j=0; j<(rs->n - rs->k); j++)
      rs->generator_matrix[ row * rs->n + rs->k + j ] = r[j];
    rs->generator_matrix[ row * rs->n + row ] = 1;  // this will be the systematix part
  } // for

  free(p);
  free(q);
  free(r);

  // also create the transponse of generator matrix
  for(j=0; j<rs->k; j++)
    for(i=0; i<rs->n; i++)
      rs->generator_matrix_t[i * rs->k + j] = rs->generator_matrix[j * rs->n + i];

  return rs->generator_matrix;
}


// H * Gt == 0 by definition,
// return
//    0 on success,
//    1 on error
uint8_t check_parity_and_generator_matrix(rs_ctx * rs) {

  uint16_t i,j,m;
  uint8_t error_flag = 0;

  uint8_t * s = calloc(rs->k * (rs->n - rs->k), sizeof(uint8_t));
  if (!s)
    return 1;

  // calculate HxGt matrix
  for(j=0; j<(rs->n - rs->k); j++)
    for(m=0; m<rs->k; m++)
      for(i=0;i<rs->n;i++)
        GF_ADD_INTO(s[j * rs->k + m], gf_mul(rs, rs->parity_matrix[ j * rs->n + i ], rs->generator_matrix[ m * rs->n + i] ));


  // check HxGt result matrix for all zero values
  for(i=0; i<(rs->n - rs->k) && !error_flag; i++)
    for(j=0; j<rs->k && !error_flag; j++)
      if (s[i*rs->k + j])
        error_flag = 1;

  free(s);
  return error_flag;
}


//  calculate  the H (n-k,k) sized parity matrix of systematic generator matrix.
//  Basic idea: if Gt is the transponse of G generator matrix, then H x Gt == 0 by definition for any linear code
//  So if G is systematic then looks like this: G == (I(k),B(k,n-k) and thus H == (A(n-k,k),I(n-k))
//  From this quaity comes that A == - Bt   where Bt is the transponse of B matrix and I is an identity matrix
static uint8_t * create_rs_parity_matrix(rs_ctx * rs) {

  uint16_t i,j;

  rs->parity_matrix = (uint8_t *) calloc( rs->n * (rs->n - rs->k), sizeof(uint8_t));

  if (!rs->parity_matrix)
    return NULL;


  for(i=0; i<(rs->n - rs->k); i++) {
    // fill up the identitiy sub-matrix part
    rs->parity_matrix[i * rs->n + rs->k + i] = 1;

    // fill up the transponse of generator matrix sub-matrix part
    for(j=0; j<rs->k; j++)
      rs->parity_matrix[ i* rs->n + j ] = rs->generator_matrix[ j * rs->n + rs->k + i ];
  } // for

  return rs->parity_matrix;
}


// Given the n+1 coefficients of a polynomial of degree n in u[0..n], and the m+1 coefficients
// of another polynomial of degree m in v[0..m], multiple the two polynomials as u * v over GF into r[0 ... n*m]
static void gf_polynom_mul(rs_ctx * rs, uint8_t * u, uint8_t n, uint8_t * v, uint8_t m, uint8_t * r) {

  uint16_t i,j;

  // init the target space
  memset(r,0,(m*n+1)*sizeof(uint8_t));

  for (i=0; i<=m; i++)
    for (j=0; j<=n && v[i]; j++) // skip if v[i] == 0
      if (u[j])
        GF_ADD_INTO(r[i+j], gf_mul(rs, u[j], v[i])); // r[i+j] = gf_add(r[i+j], gf_mul(rs, u[j], v[i]));

}


// calculate Reed-Solomon simplest code generator polynom as (x-alpha**i)(x-alpha**i+1)...(x-alpha**(i+n-k))
// we choose i=1 case, so first tag is (x-alpha), the 2nd tag is (x-alpha**2) and so on ...
// because  we used the alpha == 2 as primitive multiplicative generator element to create exp_table
// any n-k consecutive powers of alpha is usable
static uint8_t * create_rs_generator_polynom(rs_ctx * rs) {
  uint16_t i, j;
  uint8_t s[2];
  uint8_t * g = (uint8_t *) calloc(rs->n - rs->k + 1, sizeof(uint8_t));
  uint8_t * t = (uint8_t *) calloc(rs->n - rs->k + 1, sizeof(uint8_t));

  if (!g || !t) {
    free(g);
    free(t);
    return NULL;
  } // if

  // starting by (x-a)
  g[0] = rs->exp_table[1]; // this is alpha itself
  g[1] = 1;

  for(i=2; i<=(rs->n - rs->k); i++) {
    for(j=0; j<i; j++)
      t[j] = g[j]; // copy current g coeffitients into t

    // put (x-a**î) into s
    s[0] = rs->exp_table[i];
    s[1] = 1;
    gf_polynom_mul(rs,t,i-1,s,1,g); // calculate s * t into g
  } // for
  rs->g = g;
  free(t);

  return rs->g;
}



// Gauss elimination and then substitution in an augmented (A(kxk)|I(k)) matrix
// assume that input matrix is not singular
static uint8_t gf_inverse_matrix_in_place(rs_ctx * rs, uint8_t * source_matrix, uint8_t nb_rows, uint8_t nb_columns) {

  uint16_t current_row = 0; // Initialization of the current row
  uint16_t current_col = 0; // Initialization of the current column
  int16_t i,j,t;  // intentionally signed numbers
  uint8_t * result_matrix = (uint8_t *) calloc( nb_rows * nb_columns, sizeof(uint8_t));
  uint8_t * tmp_swap_row = (uint8_t *) calloc(nb_columns, sizeof(uint8_t));
  if (!result_matrix || !tmp_swap_row) {
    free(result_matrix);
    free(tmp_swap_row);
    return 0;
  } // if

  // init result matrix as an identity matrix
  for(i=0; i<nb_rows; i++)
    result_matrix[i * nb_columns + i] = 1;

  // go trough the source matrix
  while (current_row < nb_rows && current_col < nb_columns) {

    // Find the pivot in source matrix at/below current row in curent col
    // in general Gauss elimination process pivot item is the abs_max to achive best numerical stability but
    // here we can choose pivot as the first non-null item because we are in GF => no need for improved numerical stability
    uint16_t pivot_row = current_row;
    uint8_t pivot_val = source_matrix[current_row * nb_columns + current_col];

    // if current row holds the pivot we do not need the for loop block at all
    if (pivot_val==0)
      for(i=current_row+1; i<nb_rows; i++)
        if (source_matrix[nb_columns * i + current_col] != 0) {
          pivot_row = i;
          pivot_val = source_matrix[nb_columns * i + current_col];
        } //if


    if (pivot_val==0) {
      // No pivot found in this column, get to next column
      current_col++;

    } else {
      // if pivot_row != current_row then swap rows item by item: current_row <=> idx_pivot_row in both matrix
      if (current_row!=pivot_row) {
        memcpy(tmp_swap_row, &source_matrix[current_row*nb_columns], nb_columns*sizeof(uint8_t));
        memcpy(&source_matrix[current_row*nb_columns], &source_matrix[pivot_row*nb_columns], nb_columns*sizeof(uint8_t));
        memcpy(&source_matrix[pivot_row*nb_columns], tmp_swap_row, nb_columns*sizeof(uint8_t));

        memcpy(tmp_swap_row, &result_matrix[current_row*nb_columns], nb_columns*sizeof(uint8_t));
        memcpy(&result_matrix[current_row*nb_columns], &result_matrix[pivot_row*nb_columns], nb_columns*sizeof(uint8_t));
        memcpy(&result_matrix[pivot_row*nb_columns], tmp_swap_row, nb_columns*sizeof(uint8_t));
      } // if

      // a slower version for row swap
      //      uint8_t tmp;
      //      for(i=0; i<nb_columns && current_row!=pivot_row; i++) {
      //        tmp = source_matrix[current_row*nb_columns + i];
      //        source_matrix[current_row*nb_columns + i] = source_matrix[pivot_row*nb_columns + i];
      //        source_matrix[pivot_row*nb_columns + i] = tmp;
      //
      //        tmp = result_matrix[current_row*nb_columns + i];
      //        result_matrix[current_row*nb_columns + i] = result_matrix[pivot_row*nb_columns + i];
      //        result_matrix[pivot_row*nb_columns + i] = tmp;
      //      } // for

      // Do for all rows below pivot in both matrix
      for(i = current_row + 1; i < nb_rows; i++) {
        uint8_t f = gf_div(rs, source_matrix[i * nb_columns + current_col] , source_matrix[current_row * nb_columns + current_col]);

        // Do for all elements in current row in result matrix
        // (skip if f==0 : no need to add zero (each gf_mul result would be zero))
        if (f!=0)
          for(j = 0; j < nb_columns; j++)
            if (result_matrix[current_row * nb_columns + j] != 0) // no need to add zero (gf_mul result is zero)
              GF_ADD_INTO(result_matrix[i * nb_columns + j], gf_mul(rs, result_matrix[current_row * nb_columns + j], f));
              // subs is the same as add in GF with base prime 2

        // Fill with zeros the lower part of pivot column in source matrix
        source_matrix[i * nb_columns + current_col] = 0;

        // Do for all remaining elements in current row in source matrix
        // (skip if f==0 : no need to add zero (each gf_mul result would be zero))
        if (f!=0)
          for(j = current_col + 1; j < nb_columns; j++)
            if (source_matrix[current_row * nb_columns + j] != 0)  // no need to add zero (gf_mul result is zero)
              GF_ADD_INTO(source_matrix[i * nb_columns + j], gf_mul(rs, source_matrix[current_row * nb_columns + j], f));
              // subs is the same as add in GF with base prime 2

      } // for

      // Increase current row and column
      current_row++;
      current_col++;
    } // if
  } // if

  // now we have the matrix in row echelon form and the corresponding inverse may be solved by back substitution

  // and now do the back substitution
  for(i=nb_rows-1;i>=0;i--) {
    // div all items in the result matrix at i-th row by source_matrix[i * nb_columns + i] value
    if (source_matrix[i * nb_columns + i] != 1)  // division by 1 is unneeded because alpha**0 = 1
      for(t=0;t<nb_columns;t++)
        if (result_matrix[i * nb_columns + t] != 0)
          result_matrix[i * nb_columns + t] = gf_div(rs, result_matrix[i * nb_columns + t], source_matrix[i * nb_columns + i]);

    // in each column add the multiplied value of i-th row to each upper row
    for(j=i-1;j>=0;j--) {
      if (source_matrix[j * nb_columns + i] != 0) // no need to add zero to anything
        for(t=0;t<nb_columns;t++)
          if (result_matrix[i * nb_columns + t] != 0)  // multiply by zero is zero, so no need to add zero to any item
            GF_ADD_INTO(result_matrix[j * nb_columns + t], gf_mul(rs, result_matrix[i * nb_columns + t], source_matrix[j * nb_columns + i]));
            // subs is the same as add in GF with base prime 2


      source_matrix[j * nb_columns + i] = 0; // just for clarity
    } // for
    source_matrix[i * nb_columns + i] = 1;  // just for clarity

  } // for i

  // copy the result matrix into the source
  memcpy(source_matrix, result_matrix, nb_rows*nb_columns*sizeof(uint8_t));
  //  A slower version for copiing
  //  for(i=0;i<nb_rows;i++)
  //    for(j=0;j<nb_columns;j++)
  //      source_matrix[i * nb_columns + j] = result_matrix[i * nb_columns + j];
  //

  free(tmp_swap_row);
  free(result_matrix);
  return 1;
}


// copy k different columns of n columns of generator matrix according to provided column indexes
// they must be all different ones but may be in any order
// first rs->k values of col_indexes array will be used
// returns a pointer to the new allocated and filled quadratic matrix with size rs->k or NULL in case of error
// caller must free up the allocated matrix,
// returns NULL if provided indexes is not in ascending order or there are same indexes multiple times
static uint8_t * create_rs_reduced_gen_matrix(rs_ctx * rs, uint8_t * col_indexes) {

  uint16_t i,j;

  uint8_t * used_indexes = (uint8_t *) calloc( rs->k, sizeof(uint8_t));  // local temp table for filtering duplicates
  uint8_t * reduced_generator_matrix = (uint8_t *) calloc( rs->k * rs->k, sizeof(uint8_t));
  if (!reduced_generator_matrix || !used_indexes) {
    free(used_indexes);
    free(reduced_generator_matrix);
    return NULL;
  } // if

  for(i=0; i<rs->k; i++) {
    if (col_indexes[i]>=rs->n || used_indexes[i]) {
      // invalid index value or this index value is a duplicate of a previous one
      free(used_indexes);
      free(reduced_generator_matrix);
      return NULL;
    } // if
    used_indexes[i] = 1; // marking this index value as used

    for(j=0; j<rs->k; j++) // copy the whole generator matrix column
      reduced_generator_matrix[j*rs->k + i] = rs->generator_matrix[j * rs->n + col_indexes[i]];

  } // for

  free(used_indexes);

  return reduced_generator_matrix;
}



// creates a reduced generator matrix based on provided column indexes and invert that matrix and finally transponse it
// returns NULL on error or pointer to the result matrix
// caller must free up the result matrix
// number of col_indexes must be exactly rs->k and must be all different numbers in any order
// to calculate all the original data words based on any k items of c(x) use: u(x) = ck(x)D(k)
// where D is the kxk sized decode matrix and ck(x) is the choosen k items of n in c(x)
uint8_t * create_rs_decode_matrix(rs_ctx * rs, uint8_t * col_indexes) {

  uint16_t i,j;

  // result quadratic matrix with size rs->k
  // choosing any k columns of n columns of R-S generator matrix the resulting sub-matrix is still non-singular
  // it is the same effect as calculating only the survived k code words of n ck(x)=u(k)G(k*k)
  uint8_t * reduced_generator_matrix = create_rs_reduced_gen_matrix(rs, col_indexes);
  uint8_t * inv_reduced_generator_matrix_t = NULL;

  if (!reduced_generator_matrix)
    return NULL;


  if (!gf_inverse_matrix_in_place(rs, reduced_generator_matrix, rs->k, rs->k)) {
    free(reduced_generator_matrix);
    return NULL;
  } // if


  inv_reduced_generator_matrix_t = (uint8_t *) calloc(rs->k * rs->k, sizeof(uint8_t));
  if (!inv_reduced_generator_matrix_t) {
    free(reduced_generator_matrix);
    return NULL;
  } // if

  // we need the transponse of inverted reduced_generator_matrix for easier calculations
  for(j=0; j<rs->k; j++)
    for(i=0; i<rs->k; i++)
      inv_reduced_generator_matrix_t[i * rs->k + j] = reduced_generator_matrix[j * rs->k + i];

  free(reduced_generator_matrix);
  return inv_reduced_generator_matrix_t;
}


// free-up rs and its all members
void rs_free(rs_ctx * rs) {
  if (!rs)
    return;

  free_gf_tables(rs);
  if (!rs->external_tables) {
    free(rs->g);
    free(rs->generator_matrix);
    free(rs->generator_matrix_t);
    free(rs->parity_matrix);
  } // if
  rs->g = NULL;
  rs->generator_matrix = NULL;
  rs->generator_matrix_t = NULL;
  rs->parity_matrix = NULL;

  free(rs);
}


static rs_ctx * internal_rs_init(uint8_t n, uint8_t k, clone_data_t * external_tables) {

  rs_ctx * rs = NULL;
  if (external_tables) {
    k = external_tables->k;
    n = external_tables->n;
  } //if

  if (n <= k)
    return NULL;

  rs = (rs_ctx *) calloc(1, sizeof(rs_ctx));
  if (!rs)
    return NULL;

  // initilaize/reset rs structure
  rs->k = k;
  rs->n = n;
  rs->g = NULL;
  rs->exp_table = NULL;
  rs->log_table = NULL;
  rs->div_table = NULL;
  rs->mul_table = NULL;
  rs->generator_matrix = NULL;
  rs->generator_matrix_t = NULL;
  rs->parity_matrix = NULL;
  rs->tables = NULL;
//  rs->size_of_sum_gf_table_memories = 0;
  rs->external_tables = 0;

  if (external_tables) {
    rs->external_tables = 1;
    if (!set_external_tables(rs, external_tables)) {
      rs_free(rs);
      return NULL;
    } // if
  } // if
  if (!init_gf_tables(rs)) {
    rs_free(rs);
    return NULL;
  } // if

  if (external_tables)
    return rs;

  if (!create_rs_generator_polynom(rs)) { // alloc and init
    rs_free(rs);
    return NULL;
  } // if

  if (!create_rs_generator_matrix(rs)) { // alloc and init
    rs_free(rs);
    return NULL;
  } // if

  if (!create_rs_parity_matrix(rs)) { // alloc and init
    rs_free(rs);
    return NULL;
  } // if

  return rs;

}


rs_ctx * rs_init_with_external_tables(clone_data_t * external_tables) {
  return internal_rs_init(0, 0, external_tables);
}

rs_ctx * rs_init(uint8_t n, uint8_t k) {
  return internal_rs_init(n, k, NULL);
}



// calculate original data words at specified erasure index position (i.e. original code block
// position) based on k different code words
// require the decode matrix (transponse of the inverted reduced generator matrix)
// c contains k different code words of n in ascending order leaving out the erasure data items
// according to inv_reduced_generator_matrix_t creation
uint8_t rs_decode(rs_ctx const * restrict const rs, uint8_t const * restrict const inv_reduced_generator_matrix_t,
                  uint8_t const * restrict const c, uint8_t const * restrict const req_indexes, uint8_t nb_req_indexes,
                  uint8_t * restrict const v) {

  int16_t i,j;
  uint8_t calculated_index_count = 0;
  uint16_t register row_start;
  register uint8_t a;

  for(j=0; j<nb_req_indexes; j++)
    if (req_indexes[j]<rs->k) {
      row_start = req_indexes[j] * rs->k; // this was the reason for transponsing the inverted matrix
                                           // to optimize matrix item access by pre computing row start index
      a = 0;
      // multiply the reduced c received message polynom with proper column of the inverted reduced generator matrix
      for(i=0; i<rs->k; i++)
        GF_ADD_INTO(a, gf_mul(rs, c[i], inv_reduced_generator_matrix_t[ row_start + i ] )); // j
      v[j] = a;
      calculated_index_count ++;

    } else {
      v[j] = 0; // to avoid uninitialized state
    } // if

  return calculated_index_count;
}

// H x c == 0 by definition, so check it!
// c: full length (n) code word
// e: resulted syndrome (n-k sized vector)
// return the count of non-zero syndrome positions (0...n-k)
uint8_t rs_calculate_syndrome(rs_ctx const * restrict const rs, uint8_t const * restrict const c, uint8_t * restrict const e) {

  uint16_t j,i;
  uint8_t register a;
  uint16_t register row_start;

  // no need for init e vector because all of its members will be set

  uint8_t syndrome_error_count = 0;

  for(j=0; j<(rs->n - rs->k); j++) {
    a = 0; // must be zeroing each loop
    row_start = j * rs->n;
    for(i=0;i<rs->n;i++)
      GF_ADD_INTO(a, gf_mul(rs, rs->parity_matrix[ row_start + i ], c[i] ));
    e[j] = a;
    if (a) // should be zero
      syndrome_error_count ++;
  } // for

  return syndrome_error_count;
}


// because we use systematic generator matrix calculate only the last n-k items of the code word but only required ones
// required index values may be in any order but values must be between n-k ... n-1  (because on the first k position there are the systematic input code words)
// u contains the message block items in normal neutral ascending order
uint8_t rs_encode(rs_ctx const * restrict const rs, uint8_t const * restrict const u,
                  uint8_t const * restrict const req_indexes, const uint8_t nb_req_indexes, uint8_t * restrict const r) {

  int16_t  i, j;
  uint16_t register row_start;

  uint8_t calculated_index_count = 0;
  uint8_t register a;

  // go through all required indexes
  for(j=0; j<nb_req_indexes; j++)
    if (req_indexes[j]<rs->n && req_indexes[j]>=rs->k) {
     // multiply the u message polynom with the n-k+index column of the generator matrix
      row_start = req_indexes[j]*rs->k;
      a = 0; // must be initilized az zero

      for(i=0; i<rs->k; i++)
        GF_ADD_INTO(a, gf_mul(rs, u[i], rs->generator_matrix_t[row_start + i]));
        // using transponse of generator matrix to save a multiplication by row_start

      r[j] = a;
      calculated_index_count ++;

    } else {
      r[j] = 0; // to avoid uninitialized state
    } // if

  return calculated_index_count;
}



uint8_t rs_decode_block(rs_ctx const * restrict const rs, uint8_t const * restrict const inv_reduced_generator_matrix_t,
                        uint8_t * restrict * c, const uint32_t block_length,
                        uint8_t const * restrict const req_indexes, uint8_t nb_req_indexes, uint8_t * * v) {

  int16_t i,j;
  uint32_t m;

  uint16_t gm_row_start;
  uint8_t * restrict mul_table_row;
  uint8_t * restrict ci, * restrict vj;

  uint16_t mod_block_length = (block_length >> 4) << 4 ; // block_length - block_length % 16;
  uint8_t calculated_index_count = 0;

  // go through all required indexes
  for(j=0; j<nb_req_indexes; j++) {

    // initialize all result  values to zero
    memset(v[j],0,block_length*sizeof(uint8_t));

    if (req_indexes[j]<rs->k) {

      gm_row_start = req_indexes[j] * rs->k; // this was the reason for transponing the inverted matrix
                                             // to optimize matrix item access by pre computing row start index

      // multiply the reduced c received message polynom with proper column of the inverted reduced generator matrix
      for(i=0; i<rs->k; i++) {
        mul_table_row = &rs->mul_table[ inv_reduced_generator_matrix_t[gm_row_start + i] * GF_ORDER ];
        for(m=0, ci = c[i], vj = v[j]; m < mod_block_length; m+=16) {
          GF_ADD_INTO(vj[ 0], mul_table_row[ci[ 0]]);
          GF_ADD_INTO(vj[ 1], mul_table_row[ci[ 1]]);
          GF_ADD_INTO(vj[ 2], mul_table_row[ci[ 2]]);
          GF_ADD_INTO(vj[ 3], mul_table_row[ci[ 3]]);
          GF_ADD_INTO(vj[ 4], mul_table_row[ci[ 4]]);
          GF_ADD_INTO(vj[ 5], mul_table_row[ci[ 5]]);
          GF_ADD_INTO(vj[ 6], mul_table_row[ci[ 6]]);
          GF_ADD_INTO(vj[ 7], mul_table_row[ci[ 7]]);
          GF_ADD_INTO(vj[ 8], mul_table_row[ci[ 8]]);
          GF_ADD_INTO(vj[ 9], mul_table_row[ci[ 9]]);
          GF_ADD_INTO(vj[10], mul_table_row[ci[10]]);
          GF_ADD_INTO(vj[11], mul_table_row[ci[11]]);
          GF_ADD_INTO(vj[12], mul_table_row[ci[12]]);
          GF_ADD_INTO(vj[13], mul_table_row[ci[13]]);
          GF_ADD_INTO(vj[14], mul_table_row[ci[14]]);
          GF_ADD_INTO(vj[15], mul_table_row[ci[15]]);
          vj+=16;
          ci+=16;
        } // for m

        // remaining items
        for(m=mod_block_length; m<block_length; m++) {
          GF_ADD_INTO(*vj, mul_table_row[*ci]);
          vj++;
          ci++;
        } // for

        // un-optimized version
        //  GF_ADD_INTO(v[j], gf_mul_inline(rs, c[i], inv_reduced_generator_matrix_t[ row_start + i ] )); // j

      } // for i
      calculated_index_count ++;

    } else {
      // this required index is out of scope, rsult vector values remain zero
    } // if
  } // for
  return calculated_index_count;
}


// this is the most speed critical part
// u contains the input data and the result R-S code words will placed into r
// u contains k different pointers to block_length sized data vectors
// r contains n-k different pointers to block_length sized data vectors
// see rs_encode for other restrictions
uint8_t rs_encode_block(rs_ctx const * restrict const rs, uint8_t *  *  u, const uint32_t block_length,
                        uint8_t const * restrict const req_indexes, const uint8_t nb_req_indexes, uint8_t *  * r) {

  int16_t  i, j;
  uint32_t m;
  uint16_t gm_row_start;
  register uint8_t * restrict mul_table_row;
  uint8_t  *  restrict ui, * restrict rj;

  uint16_t mod_block_length = (block_length >> 4) << 4 ; // fast way to clculate  block_length - block_length % 16;
  uint8_t calculated_index_count = 0;


  // go through all required indexes
  for(j=0; j<nb_req_indexes; j++) {
    // initialize all result  values to zero (r is a ponter array)
    memset(r[j],0,block_length*sizeof(uint8_t));

    // calculate only fec indexes (k ... n-1)
    if (req_indexes[j] >= rs->n  ||  req_indexes[j] < rs->k)
      continue; // nothing to do, this required index value is out of scope


    // hint: take the values of u as coefficients of a polinom and
    //       multiply the u message polynom with the n-k+index column of the generator matrix

    // to speed up things, we exploit the fact that u values and r values are on consecutive memory places
    gm_row_start = req_indexes[j]*rs->k; // pre-calculate here the used generator matrix row
    for(i=0; i<rs->k; i++) { // go through all matrix row coeficients

      // because the gen.matrix coefficients tipically are non-zero values,
      // do not try to optimize by checking for zero values because it would waste of time

      // one of the two items in the multiplication will be in this multiplcation table row
      // this was the reason to create transposed gen.matrix, to go through on row values
      mul_table_row = &rs->mul_table[ rs->generator_matrix_t[gm_row_start + i] * GF_ORDER];

      // go through all input/output code words
      // calculate multiple code words "paralell" to optimize memory caching
      // unroling 16 times is a good compromise (32 gave no benefit)
      for(m=0, ui = u[i], rj = r[j]; m < mod_block_length; m+=16) {
        GF_ADD_INTO(rj[ 0], mul_table_row[ui[ 0]]);
        GF_ADD_INTO(rj[ 1], mul_table_row[ui[ 1]]);
        GF_ADD_INTO(rj[ 2], mul_table_row[ui[ 2]]);
        GF_ADD_INTO(rj[ 3], mul_table_row[ui[ 3]]);
        GF_ADD_INTO(rj[ 4], mul_table_row[ui[ 4]]);
        GF_ADD_INTO(rj[ 5], mul_table_row[ui[ 5]]);
        GF_ADD_INTO(rj[ 6], mul_table_row[ui[ 6]]);
        GF_ADD_INTO(rj[ 7], mul_table_row[ui[ 7]]);
        GF_ADD_INTO(rj[ 8], mul_table_row[ui[ 8]]);
        GF_ADD_INTO(rj[ 9], mul_table_row[ui[ 9]]);
        GF_ADD_INTO(rj[10], mul_table_row[ui[10]]);
        GF_ADD_INTO(rj[11], mul_table_row[ui[11]]);
        GF_ADD_INTO(rj[12], mul_table_row[ui[12]]);
        GF_ADD_INTO(rj[13], mul_table_row[ui[13]]);
        GF_ADD_INTO(rj[14], mul_table_row[ui[14]]);
        GF_ADD_INTO(rj[15], mul_table_row[ui[15]]);
        rj+=16;
        ui+=16;
      } // for m

      // remaining items (if block length is not multiple of 16)
      // rj and ui has been initialized in previous for even when mod_block_length == 0
      for(m = mod_block_length; m < block_length; m++) {
        GF_ADD_INTO(*rj, mul_table_row[*ui]);
        rj++;
        ui++;
      } // for

      // un-optimized version:
      // GF_ADD_INTO(r[j][m], gf_mul_inline(rs, u[i][m], rs->generator_matrix_t[row_start + i]));

    } // for i
    calculated_index_count ++;

  } // for j
  return calculated_index_count;
}


