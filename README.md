# Reed-Solomon code implementation for erasure coding  
  
  
  I know you know more about Reed-Solomon codes than these short introduction,
  but I wrote it to clarify my terminology, used throughout the code comments
  and function descriptions. Because there is a finite chance that I use
  different words for that you are familiar with.  
   
  A Reed-Solomon code word over GF(2^8) consists of n bytes,
  and the original message consist of k bytes. Encoding means,
  you create the output n bytes based on the input k bytes.
  This implementation use systematic generator matrices, so the first k bytes
  of the output n bytes are always the same as the input k bytes.
  So encoding practically means to produce n-k bytes to the input message.
  somethins these n-k bytes are called forward error correction info.  
   
  Due to properties of Reed-Solomon code, any k bytes of the codeword
  may be used to re-create all the n bytes of the whole codeword,
  and thus recover k bytes of the original message.
  The only requirement for you when you lost bytes of codewrod,
  is to know what the k bytes original positions were
  in the original codeword. The lost data byte positions are called erasure errors.
  So if you lost up to n-k bytes of the codeword,
  no problem, using the forward error correction part, it is always possible
  to recover the original code word data.
  Reed-Solomon code is even capable to correct if some bytes in your codeword altered,
  but this implementation only able to signal the altering, not correc them.
  This implementation use GF(2^8) for Reed-Solomon code aritmetic.
  So the code word items are bytes and you may choose value of n between 1 and 255 freely.
  Value of k may be betwwen 1 and n-1. How to choose values for n and k?
  It depends on your use case. But rule of thumb, the bigger the chance you lost
  data, choose more difference between n and k. Bigger the n-k value,
  You will have more error correction info appended to the original data.  
   
   
 <h>Typical usage steps</h>  
   
  Choose you Reed-Solomon code parameters as 1<n<256 and 0<k<n and call init function.  
  Let n=15 and k=12.  
  This means that fec values are 3 bytes appended to the 12 bytes long input messages.
  A full codeword length will be 15 bytes long. This implentation use systematic
  generator matrix, so in each codeword on positions 0 .. 11  
  the original message bytes will be, and index positions 12 .. 14 hold fec values.  
   
>    uint8_t rs_n=15, rs_k=12;  
>    rs_ctx rs = rs_init(rs_n, rs_k);  
   
  To calculate a R-S code word vector last forward error corrector part:
  (code word bytes on the first k position are the same as input message bytes,
   it means no need to calculate them, but only the last n-k bytes)  
   
>    uint8_t r[3], u[12] = { 'h', 'e', 'l', 'l', 'o', ' ', 'w', 'o', 'r', 'l', 'd', '!'};  
>    uint8_t w[3] = {12, 13, 14};  
>    
>    uint8_t created_rs_code_word_count = rs_encode(rs, u, w, rs->n - rs->k, r);  
   
  There maybe cases when you do not want to create all the fec bytes. rs_encode function
  makes it possible to choos which fec index positons to calculate (see w vector)  
   
  Put the whole code word into c array:  
>    uint8_t c[15];  
   
  Put the original message on the first k positions:    
>    memcpy(c, u, 12);  
   
  Append the r0, r1, r2 values to the message  
>    memcpy(&c[12], r, 3);  
 
  Now c holds the whole codeword.  
  Let's check the whole code word wheter it is error free:  
   
>    uint8_t e[3];  
>    uint8_t syndrome_is_ok = rs_check(rs, c, e);  
   
   
    Let' assume we lost only 2 bytes of codeword according to the followings:  
   
  'h', ?, 'l', ?, 'o', ' ', w, 'o', 'r', 'l', 'd', '!', r0, r1, r2  
    
  Recalculate some missing message bytes based on any 12 byte of the original message:  
  Put the really available first 12 bytes of codeword  into v:    
>    uint8_t v[12] = { 'h', 'l', 'o', ' ', 'w', 'o', 'r', 'l', 'd', '!', r0, r1};  
  (in fact you may choose any 12 of available 13 data bytes)  
   
  We lost data on index positions 1 and 3 and require both of them to recalculate:  
>    uint8_t z[2] = {1,3}}  
 
  First we need to calculate a new decode matrix based on available data indexes:  
  (choose any 12 of available 13 data index values but same ones as you filled up v vector)  
  Please notice that this decode matrix depends only on available indexes and not the actual data,
  so you don not need to recreate it as long as you decode the same erasure pattern.  
   
>    uint8_t avail_col_indexes[12] = {0,2,4,5,6,7,8,9,10,11,12,13,14}  
>    uint8_t inv_reduced_generator_matrix_t = create_rs_reduced_gen_matrix(rs, avail_col_indexes);
> 
>    uint8_t recalculated_count rs_decode(rs, v, inv_reduced_generator_matrix_t, v, z, 2, m)
 
  m will hold the required missing values as follows:  
  m[0] == 'e'  
  m[1] == 'l'  
 
  Sometimes you do not want to re-create all the missing bytes. rs_decode function
  makes it possible to choose which index positons to recalculate (see z vector)  
   
  Do not forget to free the decode matrix:  
>    free(inv_reduced_generator matrix);;  
   
  To free all the rs data:  
>    rs_free(rs);  
   
   
  If you have several messages to encode them at once use <em>rs_encode_block</em> and
  <em>rs_decode_block</em> to decode code vectors vectors  accordingly.  
