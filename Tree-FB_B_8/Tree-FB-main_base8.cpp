#include <stdio.h>
#include <unistd.h>
#include "./functional_bootstrap.cpp"

/*
  g++ Tree-FB-main_base8.cpp -O3 -std=c++11 -funroll-all-loops -march=native -I"../tfhe/src/include/" -L"../tfhe/build/libtfhe/" -ltfhe-spqlios-fma -lm -g -o TreeFBbase8
*/

void test_LUT_error(fbt_integer * output, fbt_integer * input, fbt_context * context){

 tfhe_functional_bootstrap_woKS_FFT(&context->output[0], context->bootstrap_key, &context->tlwe_temp[0], &input->lwe_samples[0], context->torus_base);
 lweCopy(&output->lwe_samples[0], &context->output[0], &context->params->tgsw_params->tlwe_params->extracted_lweparams);

}

void test_LUT_6_bit_to_6_bit(fbt_integer * output, fbt_integer * input, fbt_context * context){
  const TFheGateBootstrappingParameterSet * params = context->params;
  int f_l3a[3] = {fbt_op_load_hardcoded_LUT, 0, 0};       // this LUT is hardcoded at line 28 of functional_bootstrap.cpp and implements the identity
  int f_l3b[4] = {fbt_op_generate_mv_b_LUT, 0, 0, 128};   // load the LUTs in polynomials
  int f_l3c[1] = {fbt_op_mv_f_bootstrap_init};            // Carpov multi value bootstrap first phase
  int f_l3d[4] = {fbt_op_mv_f_bootstrap, 0, 0, f_l3b[3]/context->torus_base};       // Carpov multi value bootstrap second phase
  int f_l3e[1] = {fbt_op_next_tree_level};                // next level of the tree evaluation
  
  int f_l2a[4] = {fbt_op_tlwe_keyswitch, 0, 0, f_l3d[3]/context->torus_base};       // tlwe to trlwe packing.
  int f_l2b[4] = {fbt_op_f_bootstrap, 0, 0, f_l2a[3]};          // 12 functional bootstraps
  
  int * function[7] = {f_l3a, f_l3b, f_l3c, f_l3d, f_l3e, f_l2a, f_l2b};

  fbt_main_loop(function, 7, input->lwe_samples, context);            // Execute everything defined in "function"
    
  lweCopy(&output->lwe_samples[0], &context->output[0], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&output->lwe_samples[1], &context->output[1], &params->tgsw_params->tlwe_params->extracted_lweparams);
}

void test_LUT_8_bit_to_8_bit(fbt_integer * output, fbt_integer * input, fbt_context * context){
  const TFheGateBootstrappingParameterSet * params = context->params;
  int f_l3a[3] = {fbt_op_load_hardcoded_LUT, 0, 1};       // this LUT is hardcoded at line 28 of functional_bootstrap.cpp and implements the identity
  int f_l3b[4] = {fbt_op_generate_mv_b_LUT, 0, 0, 768};   // load the LUTs in polynomials
  int f_l3c[1] = {fbt_op_mv_f_bootstrap_init};            // Carpov multi value bootstrap first phase
  int f_l3d[4] = {fbt_op_mv_f_bootstrap, 0, 0, f_l3b[3]/context->torus_base};       // Carpov multi value bootstrap second phase
  int f_l3e[1] = {fbt_op_next_tree_level};                // next level of the tree evaluation
  
  int f_l2a[4] = {fbt_op_tlwe_keyswitch, 0, 0, f_l3d[3]/context->torus_base};       // tlwe to trlwe packing.
  int f_l2b[4] = {fbt_op_f_bootstrap, 0, 0, f_l2a[3]};          // functional bootstraps
  int f_l2c[1] = {fbt_op_next_tree_level};                // next level of the tree evaluation

  //To perform the computation correctly, it is necessary to copy the ciphertext to meet the packaging quantity requirements.
  int f_copy[4] = {100, 12, 2, 4};

  int f_l0a[4] = {fbt_op_tlwe_keyswitch, 0, 0, 3};        // tlwe to trlwe packing
  int f_l0b[4] = {fbt_op_f_bootstrap, 0, 0, f_l0a[3]};           // functional bootstraps

  int * function[11] = {f_l3a, f_l3b, f_l3c, f_l3d, f_l3e, f_l2a, f_l2b, f_l2c, f_copy, f_l0a, f_l0b};

  fbt_main_loop(function, 11, input->lwe_samples, context);            // Execute everything defined in "function"
    
  lweCopy(&output->lwe_samples[0], &context->output[0], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&output->lwe_samples[1], &context->output[1], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&output->lwe_samples[2], &context->output[2], &params->tgsw_params->tlwe_params->extracted_lweparams);
}

void test_LUT_10_bit_to_10_bit(fbt_integer * output, fbt_integer * input, fbt_context * context){
  const TFheGateBootstrappingParameterSet * params = context->params;
  int f_l4a[3] = {fbt_op_load_hardcoded_LUT, 0, 2};       // this LUT is hardcoded at line 28 of functional_bootstrap.cpp and implements the identity
  int f_l4b[4] = {fbt_op_generate_mv_b_LUT, 0, 0, 4096};   // load the LUTs in polynomials
  int f_l4c[1] = {fbt_op_mv_f_bootstrap_init};            // Carpov multi value bootstrap first phase
  int f_l4d[4] = {fbt_op_mv_f_bootstrap, 0, 0, f_l4b[3]/context->torus_base};       // Carpov multi value bootstrap second phase
  int f_l4e[1] = {fbt_op_next_tree_level};                // next level of the tree evaluation
  
  int f_l3a[4] = {fbt_op_tlwe_keyswitch, 0, 0, f_l4d[3]/context->torus_base};       // tlwe to trlwe packing.
  int f_l3b[4] = {fbt_op_f_bootstrap, 0, 0, f_l3a[3]};          // 64 functional bootstraps
  int f_l3c[1] = {fbt_op_next_tree_level};                // next level of the tree evaluation

  int f_l2a[4] = {fbt_op_tlwe_keyswitch, 0, 0, f_l3a[3]/context->torus_base};        // tlwe to trlwe packing
  int f_l2b[4] = {fbt_op_f_bootstrap, 0, 0, f_l2a[3]};           // 16 functional bootstraps
  int f_l2c[1] = {fbt_op_next_tree_level};                // next level of the tree evaluation

  int f_copy[4] = {100, 8, 3, 2};

  int f_l1a[4] = {fbt_op_tlwe_keyswitch, 0, 0, 4};        // tlwe to trlwe packing
  int f_l1b[4] = {fbt_op_f_bootstrap, 0, 0, 4};           // 4 functional bootstraps

  int * function[14] = {f_l4a, f_l4b, f_l4c, f_l4d, f_l4e, f_l3a, f_l3b, f_l3c, f_l2a, f_l2b, f_l2c, f_copy, f_l1a, f_l1b};

  fbt_main_loop(function, 14, input->lwe_samples, context);            // Execute everything defined in "function"
    
  lweCopy(&output->lwe_samples[0], &context->output[0], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&output->lwe_samples[1], &context->output[1], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&output->lwe_samples[2], &context->output[2], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&output->lwe_samples[3], &context->output[3], &params->tgsw_params->tlwe_params->extracted_lweparams);
}

void test_LUT_12_bit_to_12_bit(fbt_integer * output, fbt_integer * input, fbt_context * context){
  const TFheGateBootstrappingParameterSet * params = context->params;
  int f_l5a[3] = {fbt_op_load_hardcoded_LUT, 0, 3};       // this LUT is hardcoded at line 28 of functional_bootstrap.cpp and implements the identity
  int f_l5b[4] = {fbt_op_generate_mv_b_LUT, 0, 0, 16384};   // load the LUTs in polynomials
  int f_l5c[1] = {fbt_op_mv_f_bootstrap_init};            // Carpov multi value bootstrap first phase
  int f_l5d[4] = {fbt_op_mv_f_bootstrap, 0, 0, f_l5b[3]/context->torus_base};       // Carpov multi value bootstrap second phase
  int f_l5e[1] = {fbt_op_next_tree_level};                // next level of the tree evaluation
  
  int f_l4a[4] = {fbt_op_tlwe_keyswitch, 0, 0, f_l5d[3]/context->torus_base};       // tlwe to trlwe packing.
  int f_l4b[4] = {fbt_op_f_bootstrap, 0, 0, f_l4a[3]};          // 64 functional bootstraps
  int f_l4c[1] = {fbt_op_next_tree_level};                // next level of the tree evaluation

  int f_l3a[4] = {fbt_op_tlwe_keyswitch, 0, 0, f_l4a[3]/context->torus_base};       // tlwe to trlwe packing.
  int f_l3b[4] = {fbt_op_f_bootstrap, 0, 0, f_l3a[3]};          // 64 functional bootstraps
  int f_l3c[1] = {fbt_op_next_tree_level};                // next level of the tree evaluation

  int f_l2a[4] = {fbt_op_tlwe_keyswitch, 0, 0, f_l3a[3]/context->torus_base};        // tlwe to trlwe packing
  int f_l2b[4] = {fbt_op_f_bootstrap, 0, 0, f_l2a[3]};           // 16 functional bootstraps

  int * function[13] = {f_l5a, f_l5b, f_l5c, f_l5d, f_l5e, f_l4a, f_l4b, f_l4c,f_l3a, f_l3b, f_l3c, f_l2a, f_l2b};

  fbt_main_loop(function, 13, input->lwe_samples, context);            // Execute everything defined in "function"
    
  lweCopy(&output->lwe_samples[0], &context->output[0], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&output->lwe_samples[1], &context->output[1], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&output->lwe_samples[2], &context->output[2], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&output->lwe_samples[3], &context->output[3], &params->tgsw_params->tlwe_params->extracted_lweparams);
}

#ifndef NUM_EXE
#define NUM_EXE 10
#endif

int main(int argc, char const *argv[])
{
  printf("Generating keys... This may take a few minutes... (~1m in our machine)\n");fflush(stdout);
#ifdef P6_4_6_3
  TFheGateBootstrappingParameterSet* params = new_TFHE_parameters(6, 4, 6, 3);
#else
  TFheGateBootstrappingParameterSet* params = new_TFHE_parameters(5, 5, 6, 2);
#endif

  unsigned int seed;
  generate_random_bytes(8, (uint8_t *) &seed);
  srand(seed);


  TFheFunctionalBootstrappingSecretKeySet* key = new_random_functional_bootstrapping_secret_keyset(params, 3);
  fbt_context * context = fbt_init(key->cloud_key->params, key->cloud_key->bkFFT, key->tlweKS, key->log_torus_base, 6144);
  int64_t size = sizeof(Torus32) * key->tlweKS->N * key->tlweKS->k * key->tlweKS->N * key->tlweKS->t * ((1 << (key->tlweKS->basebit)) - 1) * context->torus_base; 
  printf("switch key size is %ld.\n", size);

 

 /* 6-bit LUT */
  int8_t lut_input_clear = rand()&0x3F;
  fbt_integer * lut_input = fbt_new_encrypted_integer(lut_input_clear, 3, (LweKey *) key->tfhe_keys->lwe_key, context);
  fbt_integer * lut_output = fbt_new_encrypted_integer(lut_input_clear, 3, key->extracted_key, context);
  printf("\nRunning 6-bit-to-6-bit LUT evaluation %d times. The LUT encodes the identity function.\n", NUM_EXE);
  MEASURE_TIME_AND_PRINT(
  for (size_t i = 0; i < NUM_EXE; i++)
    {
      test_LUT_6_bit_to_6_bit(lut_output, lut_input, context);
      
    }
  );
  printf("Result LUT(%d): %d\n", lut_input_clear, fbt_integer_decrypt(lut_output, key->extracted_key, context));


 /*8-bit LUT*/
  int16_t lut_input_clear_8 = rand()&0xFF;
  fbt_integer * lut_input_8 = fbt_new_encrypted_integer(lut_input_clear_8, 3, (LweKey *) key->tfhe_keys->lwe_key, context);
  fbt_integer * lut_output_8 = fbt_new_encrypted_integer(lut_input_clear_8, 3, key->extracted_key, context);
  printf("\nRunning 8-bit-to-8-bit LUT evaluation %d times. The LUT encodes the identity function.\n", NUM_EXE);
  MEASURE_TIME_AND_PRINT(
  for (size_t i = 0; i < NUM_EXE; i++)
    {
      test_LUT_8_bit_to_8_bit(lut_output_8, lut_input_8, context);
    }
  );
  printf("Result LUT(%d): %d\n", lut_input_clear_8, fbt_integer_decrypt(lut_output_8, key->extracted_key, context));

 /*10-bit LUT*/
  int16_t lut_input_clear_10 = rand()&0x3FF;
  fbt_integer * lut_input_10 = fbt_new_encrypted_integer(lut_input_clear_10, 4, (LweKey *) key->tfhe_keys->lwe_key, context);
  fbt_integer * lut_output_10 = fbt_new_encrypted_integer(lut_input_clear_10, 4, key->extracted_key, context);
  printf("\nRunning 10-bit-to-10-bit LUT evaluation %d times. The LUT encodes the identity function.\n", NUM_EXE);
  MEASURE_TIME_AND_PRINT(
  for (size_t i = 0; i < NUM_EXE; i++)
    {
      test_LUT_10_bit_to_10_bit(lut_output_10, lut_input_10, context);
    }
  );
  printf("Result LUT(%d): %d\n", lut_input_clear_10, fbt_integer_decrypt(lut_output_10, key->extracted_key, context));

 /*12-bit LUT*/
  int16_t lut_input_clear_12 = rand()&0xFFF;
  fbt_integer * lut_input_12 = fbt_new_encrypted_integer(lut_input_clear_12, 5, (LweKey *) key->tfhe_keys->lwe_key, context);
  fbt_integer * lut_output_12 = fbt_new_encrypted_integer(lut_input_clear_12, 5, key->extracted_key, context);
  printf("\nRunning 12-bit-to-12-bit LUT evaluation %d times. The LUT encodes the identity function.\n", NUM_EXE);
  MEASURE_TIME_AND_PRINT(
  for (size_t i = 0; i < NUM_EXE; i++)
    {
      test_LUT_12_bit_to_12_bit(lut_output_12, lut_input_12, context);
    }
  );
  printf("Result LUT(%d): %d\n", lut_input_clear_12, fbt_integer_decrypt(lut_output_12, key->extracted_key, context));

  return 0;
}
