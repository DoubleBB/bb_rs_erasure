/*
 Reed-Solomon code implementation for erasure coding  v.1.0

 Copyright (c) 2019 Bela Bodecs   (bodecsb#vivanet.hu)


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


#include <stdint.h>
#ifndef _bb_rs_erasure_h_
#define _bb_rs_erasure_h_

// Galoise Field (GF) of prime 2 with power of 8  a.k.a. GF(2**8)
#define GF_ORDER 256 // 2**8 - number of elements in this GF
#define GF_MUL_ORDER 255 // multiplicative order of primitive element is GF_ORDER - 1


// in GF(2**m) adding is the XORing the two component ans adding is the same as substraction
// for efficiency reasons, it is better to use as macro than as function
#define GF_ADD_INTO(a, b) a ^= b


#define GF_PRIMITIVE_POLYNOMIAL 0b100011101 // choosen irreducible primitive polynomial to generate elements of GF

/*
 * all 16 possible primitive polynomials for GF(2**8) (as binary and decimal number)
 *
 * https://link.springer.com/content/pdf/bbm%3A978-1-4615-1509-8%2F1.pdf
 * http://users.ece.cmu.edu/~koopman/lfsr/8.txt
 *
 *
 *   100011101   (285)
 *   100101011   (299)
 *   100101101   (301)
 *   101001101   (333)
 *   101011111   (351)   inefficient
 *   101100011   (355)
 *   101100101   (357)
 *   101101001   (361)
 *   101110001   (369)
 *   110000111   (391)
 *   110001101   (397)
 *   110101001   (425)
 *   111000011   (451)
 *   111001111   (463)   inefficient
 *   111100111   (487)   inefficient
 *   111110101   (501)   inefficient
 *
 *
 */



/*
 * Some words about Reed-Solomon codes
 *
 *   I know you know more about Reed-Solomon codes than these short introduction,
 *   but I wrote it to clarify my terminology, used throughout the code comments
 *   and function descriptions. Because there is a finite chance that I use
 *   different words for that you are familiar with.
 *
 *   A Reed-Solomon code word over GF(2**8) consists of n bytes,
 *   and the original message consist of k bytes. Encoding means,
 *   you create the output n bytes based on the input k bytes.
 *   This implementation use systematic generator matrices, so the first k bytes
 *   of the output n bytes are always the same as the input k bytes.
 *   So encoding practically means to produce n-k bytes to the input message.
 *   somethins these n-k bytes are called forward error correction info.
 *
 *   Due to properties of Reed-Solomon code, any k bytes of the codeword
 *   may be used to re-create all the n bytes of the whole codeword,
 *   and thus recover k bytes of the original message.
 *   The only requirement for you when you lost bytes of codewrod,
 *   is to know what the k bytes original positions were
 *   in the original codeword. The lost data byte positions are called erasure errors.
 *   So if you lost up to n-k bytes of the codeword,
 *   no problem, using the forward error correction part, it is always possible
 *   to recover the original code word data.
 *   Reed-Solomon code is even capable to correct if some bytes in your codeword altered,
 *   but this implementation only able to signal the altering, not correc them.
 *   This implementation use GF(2**8) for Reed-Solomon code aritmetic.
 *   So the code word items are bytes and you may choose value of n between 1 and 255 freely.
 *   Value of k may be betwwen 1 and n-1. How to choose values for n and k?
 *   It depends on your use case. But rule of thumb, the bigger the chance you lost
 *   data, choose more difference between n and k. Bigger the n-k value,
 *   You will have more error correction info appended to the original data.
 *
 *
 * Typical usage steps
 *
 *   Choose you Reed-Solomon code parameters as 1<n<256 and 0<k<n and call init function.
 *   Let n=15 and k=12.
 *   This means that fec values are 3 bytes appended to the 12 bytes long input messages
 *   A full codewrod length will be 15 bytes long. This implentation use systematic
 *   generator matrix, so in each codeword on positions 0...11
 *   the original message bytes will be, and index positions 12..14 hold fec values.
 *
 *   uint8_t rs_n=15, rs_k=12;
 *   rs_ctx * rs = rs_init(rs_n, rs_k);
 *
 *   To calculate a R-S code word vector last forward error corrector part:
 *   (code word bytes on the first k position are the same as input message bytes,
 *    it means no need to calculate them, but only the last n-k bytes)
 *
 *   uint8_t r[3], u[12] = { 'h', 'e', 'l', 'l', 'o', ' ', 'w', 'o', 'r', 'l', 'd', '!'};
 *   uint8_t w[3] = {12,13,14};
 *
 *   uint8_t created_rs_code_word_count = rs_encode(rs, u, w, rs->n - rs->k, r);
 *
 *   There maybe cases when you do not want to create all the fec bytes. rs_encode function
 *   makes it possible to choos which fec index positons to calculate (see w vector)
 *
 *   Put the whole code word into c array:
 *   uint8_t c[15];
 *
 *   // put the original message on the first k positions
 *   memcpy(c,u,12);
 *
 *   append the r0, r1, r2 values to the message
 *   memcpy(&c[12],r,3)
 *
 *   (To save these two memory copies we could declare u as 15 bytes long array and set r as a pointer to the 13th bytes)
 *
 *
 *   Now c holds the whole codeword.
 *   Now, Let's check the whole code word wheter it is error free:
 *
 *   uint8_t e[3];
 *   uint8_t syndrome_is_ok = rs_check(rs, c, e) {
 *
 *
 *   Recalculate some missing message bytes based on any 12 byte of the original message:
 *   Let' assume we lost only 2 bytes of codeword according to the followings:
 *
 *   'h', ?, 'l', ?, 'o', ' ', w, 'o', 'r', 'l', 'd', '!', r0, r1, r2
 *
 *   Put the really available first 12 bytes of codeword  into v
 *   uint8_t v[12] = { 'h', 'l', 'o', ' ', 'w', 'o', 'r', 'l', 'd', '!', r0, r1};
 *   (in fact you may choose any 12 of available 13 data bytes)
 *
 *   We lost data on index positions 1 and 3 and require both of them to recalculate:
 *   uint8_t z[2] = {1,3}}
 *
 *   First we need to calculate a new decode matrix based on available data indexes:
 *   (choose any 12 of available 13 data index values but same ones as you filled up v vector)
 *   Please notice that this decode matrix depends only on available indexes and not the actual data,
 *   so you don not need to recreate it as long as you decode the same erasure pattern
 *
 *   uint8_t avail_col_indexes[12] = {0,2,4,5,6,7,8,9,10,11,12,13,14}
 *   uint8_t * inv_reduced_generator_matrix_t = create_rs_decode_matrix(rs, avail_col_indexes);
 *
 *   uint8_t recalculated_count rs_decode(rs, v, inv_reduced_generator_matrix_t, v, z, 2, m)
 *
 *   m will hold the required missing values as follows:
 *   m[0] == 'e'
 *   m[1] == 'l'
 *
 *   Sometimes you do not want to re-create all the missing bytes. rs_decode function
 *   makes it possible to choose which index positons to recalculate (see z vector)
 *
 *   Do mot forget to free the decode matrix:
 *   free(inv_reduced_generator matrix);
 *
 *   To free all the rs data:
 *   rs_free(rs)
 *
 *
 *   If you have several messages to encode them at once use rs_encode_block and
 *   rs_decode_block to decode code vectors vectors  accordingly.
 */


// Holds parameters and lookup tables
typedef struct rs_ctx_struct {
  uint8_t n;
  uint8_t k;
  uint8_t * g; // coefficients of R-S generator polinom x**0 ... x**(n-k)
  uint8_t * generator_matrix; // will hold k x n systematic generator matrix for Reed-Solomon code
                              // with parameters (n,k) over GF(2^8)
  uint8_t * generator_matrix_t; // transponse of generator matrix
                                // (working on rows is more optimized than column operations on generator matrix)
  uint8_t * parity_matrix; // n x k systematic "parity check" matrix for Reed-Solomon code
                           // with parameters (n,k) over GF(2^8)


  // to speed up GF calculations we use table lookups
  uint8_t * exp_table;  // size of GF_MUL_ORDER
  uint8_t * log_table;  // size of 2 * GF_MIL_ORDER
  uint8_t * div_table;  // size of GF_ORDER * GF_ORDER
  uint8_t * mul_table;  // size of GF_ORDER * GF_ORDER
} rs_ctx;


// H * Gt == 0 by definition, so calculate it and check,
// return 0 on failure and 1 if all result elements are zero
uint8_t check_parity_and_generator_matrix(rs_ctx * rs);



// creates a reduced generator matrix based on provided column indexes and
// then invert that matrix and finally transponse it
// number of provided column indexes is k,
// their order is indifferent and may be any k of n possible index values
// (but you must provide data blocks to decode function in this very same order)
// returns NULL on error or pointer to the result matrix
// caller must free up the result matrix
// to calculate all the original data words based on any k items of c(x) useing u(x) = ck(x)D(k)
// where D is the kxk sized decode matrix and ck(x) is the choosen k items of n in c(x)
uint8_t * create_rs_decode_matrix(rs_ctx * rs, uint8_t * col_indexes);


// deallocates all data allocated by rs_init
void rs_free(rs_ctx * rs);

// allocate and init all rs internal data tables for encoding
rs_ctx * rs_init(uint8_t n, uint8_t k);




// calculate original data words at specified erasure req_index position
// (i.e. original code block position) based on k different code words into v
// require the decode matrix (transponse of the inverted reduced generator matrix) created by create_rs_decode_matrix()
// c contains k different code words of n leaving out the erasure data items,
// their order according to inv_reduced_generator_matrix_t creation
// return number of calculated code words, on success returns nb_req_indexes
uint8_t rs_decode(rs_ctx const * restrict const rs, uint8_t const * restrict const inv_reduced_generator_matrix_t,
                  uint8_t const * restrict const c, uint8_t const * restrict const req_indexes, uint8_t nb_req_indexes,
                  uint8_t * restrict const v);


// same as rs_decode but works on data arrays
uint8_t rs_decode_block(rs_ctx const * restrict const rs, uint8_t const * restrict const inv_reduced_generator_matrix_t,
                  uint8_t * restrict * c, const uint32_t block_length,
                  uint8_t const * restrict const req_indexes, uint8_t nb_req_indexes, uint8_t * * v);




// check code word values by calculate n-k sized syndrome vector,
// Each syndrome will be zero in case of error free code word
// in input c holds the n code word values,
// on return the result vector e will hold the n-k syndrome values
// return 0 if any syndrome value is non-zero
// return 1 if all syndrome values are zero
uint8_t rs_check(rs_ctx const * restrict const rs, uint8_t const * restrict const c, uint8_t * restrict const e);



// calculate the required code word values based on message values into r
// (codeword values count is n, so its indexes go through 0 ... n-1)
// because we use systematic generator matrix,  only the last n-k items of the code word values may be required to calculate
// (because on the first k position there are the systematic input message word values,
// no need to calculate anything on the 0 .. k elements)
// req_indexes holds the required code word index values, they may be in any order but values must be between n-k ... n-1
// u contains the message block values in normal neutral ascending order
// nb_req_indexes must be between 1 ... n-k, and size of r vector must be equal or greater than nb_req_indexes
// returns number of  calculated R-S code words, on success return value is nb_req_indexes
uint8_t rs_encode(rs_ctx const * restrict const rs, uint8_t const * restrict const u,
                  uint8_t const * restrict const req_indexes, const uint8_t nb_req_indexes, uint8_t * restrict const r);


// same as rs_encode but works on data arrays
uint8_t rs_encode_block(rs_ctx const * restrict const rs, uint8_t  *  * u, const uint32_t block_length,
                        uint8_t const * restrict const req_indexes, const uint8_t nb_req_indexes, uint8_t * * r);

#endif
