#include <stdio.h>
#include <unistd.h>
#include "BootManyLUT.cpp"

/* 
export LD_LIBRARY_PATH="/home/liu/cuFHE/cufhe/include/"  export LD_LIBRARY_PATH="/home/liu/cuFHE/cufhe/lib" g++ test_api_gpu.cu -I"/home/liu/cuFHE/cufhe/include/" -L"/home/liu/cuFHE/cufhe/lib" -o main
export LD_LIBRARY_PATH="/home/liu/cuFHE/cufhe/lib"

step1: 
  g++ Tree-BML-main_base8.cpp -O3 -std=c++11 -funroll-all-loops -march=native -I"../tfhe/src/include/" -L"../tfhe/build/libtfhe/" -ltfhe-spqlios-fma -lm -g -o TreeBMLmainbase8
step2: 
  ./TreeBMLmain 2
*/

void test_LUT_6_bit_to_6_bit_BML(fbt_integer * output, fbt_integer * input, fbt_context_BML * context, const int theta){
  const TFheGateBootstrappingParameterSet * params = context->params;
  // theta = 1
  int f_l3a[3] = {fbt_op_load_hardcoded_LUT, 0, 0}; //Let LUTs into the lut_input integer array in order
  int f_l3b[4] = {fbt_op_generate_mv_b_LUT, 0, 0, 128}; //Let lut_input into a polynomial according to rules
  int f_l3c[1] = {fbt_op_mv_f_bootstrap_init};  //Changing the integer coefficient polynomial to a noiseless TLWE ciphertext
  int f_l3d[4] = {fbt_op_mv_f_bootstrap, 0, 0, 16}; //MVB
  int f_l3e[1] = {fbt_op_next_tree_level};

  int f_l2a[4] = {fbt_op_tlwe_keyswitch_BML_conclusion, 0, 0, f_l3d[3]};       // tlwe to trlwe packing key switching.
  
  int f_l2b[6] = {fbt_op_f_bootstrap_BML, 0, 0, f_l2a[3]/((1<<theta)*context->torus_base), 1, 0};

   
  if (theta == 2)
  {       
    f_l2b[3] = 1; f_l2b[4] = 1; f_l2b[5] = 1;
  }
  else if(theta == 3)
  {
    f_l2b[3] = 1; f_l2b[4] = 1; f_l2b[5] = 1;
  }
  int * function[7] = {f_l3a, f_l3b, f_l3c, f_l3d, f_l3e, f_l2a, f_l2b};
  fbt_main_loop_BML(function, 7, input->lwe_samples, context, theta);

  lweCopy(&output->lwe_samples[0], &context->output[0], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&output->lwe_samples[1], &context->output[1], &params->tgsw_params->tlwe_params->extracted_lweparams);
  //lweCopy(&output->lwe_samples[2], &context->output[2], &params->tgsw_params->tlwe_params->extracted_lweparams);

}

void test_LUT_8_bit_to_8_bit_BML(fbt_integer * output, fbt_integer * input, fbt_context_BML * context, const int theta){
  const TFheGateBootstrappingParameterSet * params = context->params;
  int f_l3a[3] = {fbt_op_load_hardcoded_LUT, 0, 1};
  int f_l3b[4] = {fbt_op_generate_mv_b_LUT, 0, 0, 768};
  int f_l3c[1] = {fbt_op_mv_f_bootstrap_init};
  int f_l3d[4] = {fbt_op_mv_f_bootstrap, 0, 0, f_l3b[3]/context->torus_base};
  int f_l3e[1] = {fbt_op_next_tree_level};

  int f_l2a[4] = {fbt_op_tlwe_keyswitch_BML_conclusion, 0, 0, f_l3d[3]}; 
  int f_l2b[6] = {fbt_op_f_bootstrap_BML, 0, 0, f_l2a[3]/((1<<theta)*context->torus_base), theta, 0};  
  int f_l2c[1] = {fbt_op_next_tree_level}; 
  
  int f_copy[4] = {100, 12, 2, 4};

  int f_l1a[4] = {fbt_op_tlwe_keyswitch_BML_conclusion, 0, 0, 24};
  int f_l1b[6] = {fbt_op_f_bootstrap_BML, 0, 0, 2, 0, 1};

  if (theta == 2)
  {       
    f_l1b[3] = 1; f_l1b[4] = theta; f_l1b[5] = 0;
  }
  else if(theta == 3)
  {
    f_l2b[3] = 2; f_l2b[4] = 2; f_l2b[5] = 1;
    f_l1b[3] = 1; f_l1b[4] = 2; f_l1b[5] = 1;
  }
  int * function[11] = {f_l3a, f_l3b, f_l3c, f_l3d, f_l3e, f_l2a,f_l2b, f_l2c, f_copy, f_l1a, f_l1b};


  fbt_main_loop_BML(function, 11, input->lwe_samples, context, theta);

  lweCopy(&output->lwe_samples[0], &context->output[0], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&output->lwe_samples[1], &context->output[1], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&output->lwe_samples[2], &context->output[2], &params->tgsw_params->tlwe_params->extracted_lweparams);
}


void test_LUT_10_bit_to_10_bit_BML(fbt_integer * output, fbt_integer * input, fbt_context_BML * context, const int theta){
  const TFheGateBootstrappingParameterSet * params = context->params;
  int f_l4a[3] = {fbt_op_load_hardcoded_LUT, 0, 2};
  int f_l4b[4] = {fbt_op_generate_mv_b_LUT, 0, 0, 4096};
  int f_l4c[1] = {fbt_op_mv_f_bootstrap_init};
  int f_l4d[4] = {fbt_op_mv_f_bootstrap, 0, 0, f_l4b[3]/context->torus_base};
  int f_l4e[1] = {fbt_op_next_tree_level};

  int f_l3a[4] = {fbt_op_tlwe_keyswitch_BML_conclusion, 0, 0, f_l4d[3]}; 
  int f_l3b[6] = {fbt_op_f_bootstrap_BML, 0, 0, f_l3a[3]/((1<<theta)*context->torus_base), theta, 0};         
  int f_l3c[1] = {fbt_op_next_tree_level}; 

  int f_l2a[4] = {fbt_op_tlwe_keyswitch_BML_conclusion, 0, 0, f_l3a[3]/context->torus_base}; 
  int f_l2b[6] = {fbt_op_f_bootstrap_BML, 0, 0, f_l2a[3]/((1<<theta)*context->torus_base), theta, 0};
  int f_l2c[1] = {fbt_op_next_tree_level};  

  int f_copy[4] = {100, 8, 3, 2};

  int f_l1a[4] = {fbt_op_tlwe_keyswitch_BML_conclusion, 0, 0, 32};
  int f_l1b[6] = {fbt_op_f_bootstrap_BML, 0, 0, f_l1a[3]/((1<<theta)*context->torus_base), theta, 0};

  if (theta == 2)
  {       
    f_l1b[3] = 1; f_l1b[4] = theta; f_l1b[5] = 0;
  }
  else if(theta == 3)
  {
    f_l1b[3] = 1; f_l1b[4] = 2; f_l1b[5] = 1;
  }

  int * function[14] = {f_l4a, f_l4b, f_l4c, f_l4d, f_l4e, f_l3a, f_l3b, f_l3c, f_l2a, f_l2b, f_l2c, f_copy, f_l1a, f_l1b};


  fbt_main_loop_BML(function, 14, input->lwe_samples, context, theta);

  lweCopy(&output->lwe_samples[0], &context->output[0], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&output->lwe_samples[1], &context->output[1], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&output->lwe_samples[2], &context->output[2], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&output->lwe_samples[3], &context->output[3], &params->tgsw_params->tlwe_params->extracted_lweparams);
}

void test_LUT_12_bit_to_12_bit_BML(fbt_integer * output, fbt_integer * input, fbt_context_BML * context, const int theta){
  const TFheGateBootstrappingParameterSet * params = context->params;
  int f_l5a[3] = {fbt_op_load_hardcoded_LUT, 0, 3};
  int f_l5b[4] = {fbt_op_generate_mv_b_LUT, 0, 0, 16384};
  int f_l5c[1] = {fbt_op_mv_f_bootstrap_init};
  int f_l5d[4] = {fbt_op_mv_f_bootstrap, 0, 0, f_l5b[3]/context->torus_base};
  int f_l5e[1] = {fbt_op_next_tree_level};

  int f_l4a[4] = {fbt_op_tlwe_keyswitch_BML_conclusion, 0, 0, f_l5d[3]};       
  int f_l4b[6] = {fbt_op_f_bootstrap_BML, 0, 0, f_l4a[3]/((1<<theta)*context->torus_base), theta, 0};       
  int f_l4c[1] = {fbt_op_next_tree_level};               

  int f_l3a[4] = {fbt_op_tlwe_keyswitch_BML_conclusion, 0, 0, f_l4a[3]/context->torus_base};        
  int f_l3b[6] = {fbt_op_f_bootstrap_BML, 0, 0, f_l3a[3]/((1<<theta)*context->torus_base), theta, 0};
  int f_l3c[1] = {fbt_op_next_tree_level};              

  int f_l2a[4] = {fbt_op_tlwe_keyswitch_BML_conclusion, 0, 0, f_l3a[3]/context->torus_base};    
  int f_l2b[6] = {fbt_op_f_bootstrap_BML, 0, 0, f_l2a[3]/((1<<theta)*context->torus_base), theta, 0};
  
  if (theta == 2)
  {       
    f_l2b[3] = 1; f_l2b[4] = theta; f_l2b[5] = 0;
  }
  else if(theta == 3)
  {
    f_l2b[3] = 1; f_l2b[4] = 2; f_l2b[5] = 1;
  }

  int * function[13] = {f_l5a, f_l5b, f_l5c, f_l5d, f_l5e, f_l4a, f_l4b, f_l4c, f_l3a, f_l3b, f_l3c, f_l2a, f_l2b};


  fbt_main_loop_BML(function, 13, input->lwe_samples, context, theta);

  lweCopy(&output->lwe_samples[0], &context->output[0], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&output->lwe_samples[1], &context->output[1], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&output->lwe_samples[2], &context->output[2], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&output->lwe_samples[3], &context->output[3], &params->tgsw_params->tlwe_params->extracted_lweparams);
}

/**/

#ifndef NUM_EXE
#define NUM_EXE 10
#endif

int main(int argc, char const *argv[])
{
  const int theta = atoi(argv[1]);;
  printf("Generating keys... This may take a few minutes...\nNow theta = %d\n", theta);fflush(stdout);

#ifdef P6_4_6_3
  TFheGateBootstrappingParameterSet* params = new_TFHE_parameters(6, 4, 6, 3);
#else
  TFheGateBootstrappingParameterSet* params = new_TFHE_parameters(5, 5, 6, 2);
#endif


  unsigned int seed;
  generate_random_bytes(8, (uint8_t *) &seed);
  srand(seed);

 
  TFheFunctionalBootstrappingSecretKeySet_BML* key_BML = new_random_functional_bootstrapping_secret_keyset_BML(params, 3, theta);
  fbt_context_BML * context_BML = fbt_init_BML(key_BML->cloud_key->params, key_BML->cloud_key->bkFFT, key_BML->tlweKS, key_BML->log_torus_base, 6144);
  
  int64_t size = sizeof(Torus32) * key_BML->tlweKS->N * key_BML->tlweKS->k * key_BML->tlweKS->N * key_BML->tlweKS->t * ((1 << (key_BML->tlweKS->basebit)) - 1); 
  printf("switch key size is %ld bit.\n", size);

 

  /*6-bit*/
  int lut_input_clear_new = rand()&0x3F;
  fbt_integer * lut_input_new = fbt_new_encrypted_integer_BML(lut_input_clear_new, 3, (LweKey *) key_BML->tfhe_keys->lwe_key, context_BML);
  fbt_integer * lut_output_new = fbt_new_encrypted_integer_BML(lut_input_clear_new, 3, key_BML->extracted_key, context_BML);
  printf("\nRunning 6-bit-to-6-bit LUT evaluation %d times. The LUT encodes the identity function.------------------\n", NUM_EXE);
  MEASURE_TIME_AND_PRINT(
  for (size_t i = 0; i < NUM_EXE; i++)
    {
      test_LUT_6_bit_to_6_bit_BML(lut_output_new, lut_input_new, context_BML, theta);
    }
  );
  printf("Result LUT(%d): %d\n", lut_input_clear_new, fbt_integer_decrypt_BML(lut_output_new, key_BML->extracted_key, context_BML)); 



  int16_t lut_input_clear_new_8 = rand()&0xFF;
  fbt_integer * lut_input_new_8 = fbt_new_encrypted_integer_BML(lut_input_clear_new_8, 3, (LweKey *) key_BML->tfhe_keys->lwe_key, context_BML);
  fbt_integer * lut_output_new_8 = fbt_new_encrypted_integer_BML(lut_input_clear_new_8, 3, key_BML->extracted_key, context_BML);
  printf("\nRunning 8-bit-to-8-bit LUT evaluation %d times. The LUT encodes the identity function.------------------\n", NUM_EXE);
  MEASURE_TIME_AND_PRINT(
  for (size_t i = 0; i < NUM_EXE; i++)
    {
      test_LUT_8_bit_to_8_bit_BML(lut_output_new_8, lut_input_new_8, context_BML, theta);
    }
  );
  printf("Result LUT(%d): %d\n", lut_input_clear_new_8, fbt_integer_decrypt_BML(lut_output_new_8, key_BML->extracted_key, context_BML)); 


  int16_t lut_input_clear_new_10 = rand()&0x3FF;
  fbt_integer * lut_input_new_10 = fbt_new_encrypted_integer_BML(lut_input_clear_new_10, 4, (LweKey *) key_BML->tfhe_keys->lwe_key, context_BML);
  fbt_integer * lut_output_new_10 = fbt_new_encrypted_integer_BML(lut_input_clear_new_10, 4, key_BML->extracted_key, context_BML);
  printf("\nRunning 10-bit-to-10-bit LUT evaluation %d times. The LUT encodes the identity function.------------------\n", NUM_EXE);
  MEASURE_TIME_AND_PRINT(
  for (size_t i = 0; i < NUM_EXE; i++)
    {
      test_LUT_10_bit_to_10_bit_BML(lut_output_new_10, lut_input_new_10, context_BML, theta);
    }
  );
  printf("Result LUT(%d): %d\n", lut_input_clear_new_10, fbt_integer_decrypt_BML(lut_output_new_10, key_BML->extracted_key, context_BML)); 



  int16_t lut_input_clear_new_12 = rand()&0xFFF;
  fbt_integer * lut_input_new_12 = fbt_new_encrypted_integer_BML(lut_input_clear_new_12, 5, (LweKey *) key_BML->tfhe_keys->lwe_key, context_BML);
  fbt_integer * lut_output_new_12 = fbt_new_encrypted_integer_BML(lut_input_clear_new_12, 5, key_BML->extracted_key, context_BML);
  printf("\nRunning 12-bit-to-12-bit LUT evaluation %d times. The LUT encodes the identity function.------------------\n", NUM_EXE);
  MEASURE_TIME_AND_PRINT(
  for (size_t i = 0; i < NUM_EXE; i++)
    {
      test_LUT_12_bit_to_12_bit_BML(lut_output_new_12, lut_input_new_12, context_BML, theta);
    }
  );
  printf("Result LUT(%d): %d\n", lut_input_clear_new_12, fbt_integer_decrypt_BML(lut_output_new_12, key_BML->extracted_key, context_BML)); 



  
  return 0;
}
