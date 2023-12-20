#include <stdio.h>
#include <unistd.h>
#include "BootManyLUT.cpp"

/* 
export LD_LIBRARY_PATH="/home/liu/cuFHE/cufhe/include/"  export LD_LIBRARY_PATH="/home/liu/cuFHE/cufhe/lib" g++ test_api_gpu.cu -I"/home/liu/cuFHE/cufhe/include/" -L"/home/liu/cuFHE/cufhe/lib" -o main
export LD_LIBRARY_PATH="/home/liu/cuFHE/cufhe/lib"

step1: 
  g++ Tree-BML-mainMulthread.cpp -O3 -std=c++11 -funroll-all-loops -march=native -I"../tfhe/src/include/" -L"../tfhe/build/libtfhe/" -ltfhe-spqlios-fma -lm -g -o TreeBMLmainMultithread -pthread
step2: 
  ./TreeBMLmain 2
*/

void test_LUT_6_bit_to_6_bit_BML(fbt_integer * output, fbt_integer * input, fbt_context_BML * context, const int thet){
  const TFheGateBootstrappingParameterSet * params = context->params;
  int f_l3a[3] = {fbt_op_load_hardcoded_LUT, 0, 0}; //Let LUTs into the lut_input integer array in order
  int f_l3b[4] = {fbt_op_generate_mv_b_LUT, 0, 0, 192}; //Let lut_input into a polynomial according to rules
  int f_l3c[1] = {fbt_op_mv_f_bootstrap_init};  //Changing the integer coefficient polynomial to a noiseless TLWE ciphertext
  int f_l3d[4] = {fbt_op_mv_f_bootstrap, 0, 0, 48}; //MVB
  int f_l3e[1] = {fbt_op_next_tree_level};

  int f_l2a[4] = {fbt_op_tlwe_keyswitch_BML_conclusion, 0, 0, f_l3d[3]};       // tlwe to trlwe packing key switching.
  
  int f_l2b[6] = {fbt_op_f_bootstrap_BML, 0, 0, f_l2a[3]/((1<<thet)*context->torus_base), thet, 0};
  int f_l2c[1] = {fbt_op_next_tree_level};                // next level of the tree evaluation
  
  int f_l1a[4] = {fbt_op_tlwe_keyswitch_BML_conclusion, 0, 0, f_l2a[3]/context->torus_base}; 
  int f_l1b[6] = {fbt_op_f_bootstrap_BML, 0, 0, 1, 2, 1}; 
  /*
    if f_l1b[5] == 0, 
      carry out f_l1b[3] times bootstrappings with thet
    else{
      carry out f_l1b[3] - 1 times bootstrappings with thet
      carry out 1 times bootstrapping with thet = f_lb[4]
    } 
      
  */      
  if (thet == 1)
  {       
    f_l1b[3] = 2; f_l1b[4] = 0; f_l1b[5] = 1;
  }
  else if(thet == 3)
  {
    f_l2b[3] = 2; f_l2b[4] = 2; f_l2b[5] = 1;
    f_l1b[3] = 1; f_l1b[4] = 2; f_l1b[5] = 1;
  }
  int * function[10] = {f_l3a, f_l3b, f_l3c, f_l3d, f_l3e, f_l2a,f_l2b, f_l2c, f_l1a, f_l1b};
  fbt_main_loop_BML(function, 10, input->lwe_samples, context, thet);

  lweCopy(&output->lwe_samples[0], &context->output[0], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&output->lwe_samples[1], &context->output[1], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&output->lwe_samples[2], &context->output[2], &params->tgsw_params->tlwe_params->extracted_lweparams);

}

void test_LUT_8_bit_to_8_bit_BML(fbt_integer * output, fbt_integer * input, fbt_context_BML * context, const int thet){
  const TFheGateBootstrappingParameterSet * params = context->params;
  int f_l3a[3] = {fbt_op_load_hardcoded_LUT, 0, 1};
  int f_l3b[4] = {fbt_op_generate_mv_b_LUT, 0, 0, 1024};
  int f_l3c[1] = {fbt_op_mv_f_bootstrap_init};
  int f_l3d[4] = {fbt_op_mv_f_bootstrap, 0, 0, f_l3b[3]/context->torus_base};
  int f_l3e[1] = {fbt_op_next_tree_level};

  int f_l2a[4] = {fbt_op_tlwe_keyswitch_BML_conclusion, 0, 0, f_l3d[3]}; 
  int f_l2b[6] = {fbt_op_f_bootstrap_BML, 0, 0, f_l2a[3]/((1<<thet)*context->torus_base), thet, 0};  
  int f_l2c[1] = {fbt_op_next_tree_level}; 

  int f_l1a[4] = {fbt_op_tlwe_keyswitch_BML_conclusion, 0, 0, f_l2a[3]/context->torus_base};
  int f_l1b[6] = {fbt_op_f_bootstrap_BML, 0, 0, f_l1a[3]/((1<<thet)*context->torus_base), thet, 0};
  int f_l1c[1] = {fbt_op_next_tree_level}; 

  int f_l0a[4] = {fbt_op_tlwe_keyswitch_BML_conclusion, 0, 0, f_l1a[3]/context->torus_base};  
  int f_l0b[6] = {fbt_op_f_bootstrap_BML, 0, 0, 1, thet, 0};


  if (thet == 1)
  {       
    f_l0b[3] = 2; f_l0b[4] = thet; f_l0b[5] = 0;
  }
  else if(thet == 3)
  {
    f_l0b[3] = 1; f_l0b[4] = 2; f_l0b[5] = 1;
  }
  int * function[13] = {f_l3a, f_l3b, f_l3c, f_l3d, f_l3e, f_l2a,f_l2b, f_l2c, f_l1a, f_l1b, f_l1c, f_l0a, f_l0b};


  fbt_main_loop_BML(function, 13, input->lwe_samples, context, thet);

  lweCopy(&output->lwe_samples[0], &context->output[0], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&output->lwe_samples[1], &context->output[1], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&output->lwe_samples[2], &context->output[2], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&output->lwe_samples[3], &context->output[3], &params->tgsw_params->tlwe_params->extracted_lweparams);
}


void test_LUT_10_bit_to_10_bit_BML(fbt_integer * output, fbt_integer * input, fbt_context_BML * context, const int thet){
  const TFheGateBootstrappingParameterSet * params = context->params;
  int f_l4a[3] = {fbt_op_load_hardcoded_LUT, 0, 2};
  int f_l4b[4] = {fbt_op_generate_mv_b_LUT, 0, 0, 5120};
  int f_l4c[1] = {fbt_op_mv_f_bootstrap_init};
  int f_l4d[4] = {fbt_op_mv_f_bootstrap, 0, 0, f_l4b[3]/context->torus_base};
  int f_l4e[1] = {fbt_op_next_tree_level};

  int f_l3a[4] = {fbt_op_tlwe_keyswitch_BML_conclusion, 0, 0, f_l4d[3]}; 
  int f_l3b[6] = {fbt_op_f_bootstrap_BML, 0, 0, f_l3a[3]/((1<<thet)*context->torus_base), thet, 0};         
  int f_l3c[1] = {fbt_op_next_tree_level}; 

  int f_l2a[4] = {fbt_op_tlwe_keyswitch_BML_conclusion, 0, 0, f_l3a[3]/context->torus_base}; 
  int f_l2b[6] = {fbt_op_f_bootstrap_BML, 0, 0, f_l2a[3]/((1<<thet)*context->torus_base), thet, 0};
  int f_l2c[1] = {fbt_op_next_tree_level};  

  int f_l1a[4] = {fbt_op_tlwe_keyswitch_BML_conclusion, 0, 0, f_l2a[3]/context->torus_base};
  int f_l1b[6] = {fbt_op_f_bootstrap_BML, 0, 0, f_l1a[3]/((1<<thet)*context->torus_base), thet, 0};
  int f_l1c[1] = {fbt_op_next_tree_level};  

  int f_l0a[4] = {fbt_op_tlwe_keyswitch_BML_conclusion, 0, 0, f_l1a[3]/context->torus_base};
  int f_l0b[6] = {fbt_op_f_bootstrap_BML, 0, 0, 2, 0, 1};

  if (thet == 1)
  {       
    f_l0b[3] = 3; f_l0b[4] = 0; f_l0b[5] = 1;
  }
  else if(thet == 3)
  {
    f_l1b[3] = 3; f_l1b[4] = 2; f_l1b[5] = 1;
    f_l0b[3] = 1; f_l0b[4] = thet; f_l0b[5] = 1;
  }

  int * function[16] = {f_l4a, f_l4b, f_l4c, f_l4d, f_l4e, f_l3a, f_l3b, f_l3c, f_l2a, f_l2b, f_l2c, f_l1a, f_l1b, f_l1c, f_l0a, f_l0b};


  fbt_main_loop_BML(function, 16, input->lwe_samples, context, thet);

  lweCopy(&output->lwe_samples[0], &context->output[0], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&output->lwe_samples[1], &context->output[1], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&output->lwe_samples[2], &context->output[2], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&output->lwe_samples[3], &context->output[3], &params->tgsw_params->tlwe_params->extracted_lweparams);
}

void test_LUT_12_bit_to_12_bit_BML(fbt_integer * output, fbt_integer * input, fbt_context_BML * context, const int thet){
  const TFheGateBootstrappingParameterSet * params = context->params;
  int f_l5a[3] = {fbt_op_load_hardcoded_LUT, 0, 3};
  int f_l5b[4] = {fbt_op_generate_mv_b_LUT, 0, 0, 24576};
  int f_l5c[1] = {fbt_op_mv_f_bootstrap_init};
  int f_l5d[4] = {fbt_op_mv_f_bootstrap, 0, 0, f_l5b[3]/context->torus_base};
  int f_l5e[1] = {fbt_op_next_tree_level};

  int f_l4a[4] = {fbt_op_tlwe_keyswitch_BML_conclusion, 0, 0, f_l5d[3]};       
  int f_l4b[6] = {fbt_op_f_bootstrap_BML, 0, 0, f_l4a[3]/((1<<thet)*context->torus_base), thet, 0};       
  int f_l4c[1] = {fbt_op_next_tree_level};               

  int f_l3a[4] = {fbt_op_tlwe_keyswitch_BML_conclusion, 0, 0, f_l4a[3]/context->torus_base};        
  int f_l3b[6] = {fbt_op_f_bootstrap_BML, 0, 0, f_l3a[3]/((1<<thet)*context->torus_base), thet, 0};
  int f_l3c[1] = {fbt_op_next_tree_level};              

  int f_l2a[4] = {fbt_op_tlwe_keyswitch_BML_conclusion, 0, 0, f_l3a[3]/context->torus_base};    
  int f_l2b[6] = {fbt_op_f_bootstrap_BML, 0, 0, f_l2a[3]/((1<<thet)*context->torus_base), thet, 0};
  int f_l2c[1] = {fbt_op_next_tree_level};            

  int f_l1a[4] = {fbt_op_tlwe_keyswitch_BML_conclusion, 0, 0, f_l2a[3]/context->torus_base};     
  int f_l1b[6] = {fbt_op_f_bootstrap_BML, 0, 0, f_l1a[3]/((1<<thet)*context->torus_base), thet, 0};
  int f_l1c[1] = {fbt_op_next_tree_level};   

  int f_l0a[4] = {fbt_op_tlwe_keyswitch_BML_conclusion, 0, 0, f_l1a[3]/context->torus_base};  
  int f_l0b[6] = {fbt_op_f_bootstrap_BML, 0, 0, 2, 1, 1};

  if (thet == 1)
  {       
    f_l0b[3] = 3; f_l0b[4] = thet; f_l0b[5] = 0;
  }
  else if(thet == 3)
  {
    f_l0b[3] = 1; f_l0b[4] = thet; f_l0b[5] = 1;
  }

  int * function[19] = {f_l5a, f_l5b, f_l5c, f_l5d, f_l5e, f_l4a, f_l4b, f_l4c, f_l3a, f_l3b, f_l3c, f_l2a, f_l2b, f_l2c, f_l1a, f_l1b, 
                        f_l1c, f_l0a, f_l0b};


  fbt_main_loop_BML(function, 19, input->lwe_samples, context, thet);

  lweCopy(&output->lwe_samples[0], &context->output[0], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&output->lwe_samples[1], &context->output[1], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&output->lwe_samples[2], &context->output[2], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&output->lwe_samples[3], &context->output[3], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&output->lwe_samples[4], &context->output[4], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&output->lwe_samples[5], &context->output[5], &params->tgsw_params->tlwe_params->extracted_lweparams);
}

/**/
void test_32_bit_integer_comparison_new_ks(LweSample * output, fbt_integer * input1, fbt_integer * input2, fbt_context_BML * context, const int thet){
  const TFheGateBootstrappingParameterSet * params = context->params;
  int f_l0a[3] = {fbt_op_load_hardcoded_LUT, 0, 4};
  int f_l0b[4] = {fbt_op_generate_b_LUT, 0, 0, 4};
  int f_l0c[4] = {fbt_op_f_bootstrap, 0, 0, 1};
  int f_l0d[4] = {fbt_op_load_literal, 1, 1, 3};
  int f_l0e[1] = {fbt_op_next_tree_level};

  int f_l1a[4] = {fbt_op_tlwe_keyswitch_BML_conclusion, 0, 0, context->torus_base};
  int f_l1b[4] = {fbt_op_f_bootstrap, 0, 0, 1};
  int f_l1c[4] = {fbt_op_load_literal, 1, 1, 3};
  int f_l1d[1] = {fbt_op_next_tree_level};

  int * function[64] = {f_l0a, f_l0b, f_l0c, f_l0d, f_l0e, f_l1a, f_l1b, f_l1c, f_l1d, f_l1a, f_l1b, f_l1c, f_l1d,
                        f_l1a, f_l1b, f_l1c, f_l1d, f_l1a, f_l1b, f_l1c, f_l1d, f_l1a, f_l1b, f_l1c, f_l1d,
                        f_l1a, f_l1b, f_l1c, f_l1d, f_l1a, f_l1b, f_l1c, f_l1d, f_l1a, f_l1b, f_l1c, f_l1d,
                        f_l1a, f_l1b, f_l1c, f_l1d, f_l1a, f_l1b, f_l1c, f_l1d, f_l1a, f_l1b, f_l1c, f_l1d,
                        f_l1a, f_l1b, f_l1c, f_l1d, f_l1a, f_l1b, f_l1c, f_l1d, f_l1a, f_l1b, f_l1c, f_l1d,
                        f_l1a, f_l1b, f_l1c};


  for (size_t i = 0; i < 16; i++){
    lweCopy(&context->lwe_temp_in[i], &input1->lwe_samples[i], params->in_out_params);
    lweSubTo(&context->lwe_temp_in[i], &input2->lwe_samples[i], params->in_out_params);
  }

  fbt_main_loop_BML(function, 64, context->lwe_temp_in, context, thet);

  lweCopy(output, &context->output[0], &params->tgsw_params->tlwe_params->extracted_lweparams);
}

// relu
void relu_8bit_BML(fbt_integer * output, fbt_integer * input, fbt_context_BML * context, int thet){
  const TFheGateBootstrappingParameterSet * params = context->params;
  int f_l3a[4] = {fbt_op_tlwe_keyswitch_BML_conclusion, 0, 0, 16};
  int f_l3b[6] = {fbt_op_f_bootstrap_BML, 0, 0, 1, 2, 0};

  int * function[2] = {f_l3a, f_l3b};

  lweKeySwitch(context->lwe_temp_in, context->bootstrap_key->ks, &input->lwe_samples[3]);

  lweNoiselessTrivial(&context->input[2], 0, &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweNoiselessTrivial(&context->input[3], 0, &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweNoiselessTrivial(&context->input[6], 0, &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweNoiselessTrivial(&context->input[7], 0, &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweNoiselessTrivial(&context->input[10], 0, &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweNoiselessTrivial(&context->input[11], 0, &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweNoiselessTrivial(&context->input[12], 0, &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweNoiselessTrivial(&context->input[13], modSwitchToTorus32(1,2*context->torus_base), &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweNoiselessTrivial(&context->input[14], 0, &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweNoiselessTrivial(&context->input[15], 0, &params->tgsw_params->tlwe_params->extracted_lweparams);

  lweCopy(&context->input[0], &input->lwe_samples[0], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&context->input[1], &input->lwe_samples[0], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&context->input[4], &input->lwe_samples[1], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&context->input[5], &input->lwe_samples[1], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&context->input[8], &input->lwe_samples[2], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&context->input[9], &input->lwe_samples[2], &params->tgsw_params->tlwe_params->extracted_lweparams);

  fbt_main_loop_BML(function, 2, context->lwe_temp_in, context, thet);

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
  const int thet = atoi(argv[1]);;
  printf("Generating keys... This may take a few minutes...\nNow theta = %d\n", thet);fflush(stdout);

#ifdef P6_4_6_3
  TFheGateBootstrappingParameterSet* params = new_TFHE_parameters(6, 4, 6, 3);
#else
  TFheGateBootstrappingParameterSet* params = new_TFHE_parameters(5, 5, 6, 2);
#endif


  unsigned int seed;
  generate_random_bytes(8, (uint8_t *) &seed);
  srand(seed);

 
  TFheFunctionalBootstrappingSecretKeySet_BML* key_BML = new_random_functional_bootstrapping_secret_keyset_BML(params, 2, thet);
  fbt_context_BML * context_BML = fbt_init_BML(key_BML->cloud_key->params, key_BML->cloud_key->bkFFT, key_BML->tlweKS, key_BML->log_torus_base, 6144);
  
  int64_t size = sizeof(Torus32) * key_BML->tlweKS->N * key_BML->tlweKS->k * key_BML->tlweKS->N * key_BML->tlweKS->t * ((1 << (key_BML->tlweKS->basebit)) - 1); 
  printf("switch key size is %ld bit.\n", size);

 /*Relu*/
  fbt_integer * number1_8bit, * output_8bit;
  int8_t number1_8bit_clear = rand()&0xFF;
  number1_8bit = fbt_new_encrypted_integer_BML(number1_8bit_clear, 4, key_BML->extracted_key, context_BML);
  output_8bit = fbt_new_encrypted_integer_BML(0, 4, key_BML->extracted_key, context_BML);
  printf("\nRunning 8-bit ReLU %d times\n", NUM_EXE);
  MEASURE_TIME_AND_PRINT(
  for (size_t i = 0; i < NUM_EXE; i++)
    {
      relu_8bit_BML(output_8bit, number1_8bit, context_BML, thet);
    }
  );
  printf("Result ReLU(%d): %d\n", number1_8bit_clear, fbt_integer_decrypt_BML(output_8bit, key_BML->extracted_key, context_BML));


/* 32-bit integer comparison */
  fbt_integer * number1_32bit, * number2_32bit;
  uint32_t number1_32bit_clear = rand();
  uint32_t number2_32bit_clear = rand();
  number1_32bit = fbt_new_encrypted_integer_BML(number1_32bit_clear, 16, (LweKey *) key_BML->tfhe_keys->lwe_key, context_BML);
  number2_32bit = fbt_new_encrypted_integer_BML(number2_32bit_clear, 16, (LweKey *) key_BML->tfhe_keys->lwe_key, context_BML);

  LweSample * output = new_LweSample(&params->tgsw_params->tlwe_params->extracted_lweparams);
  printf("\nRunning 32-bit integer comparison %d times\n", NUM_EXE);
  MEASURE_TIME_AND_PRINT(
  for (size_t i = 0; i < NUM_EXE; i++)
    {
      test_32_bit_integer_comparison_new_ks(output, number1_32bit, number2_32bit, context_BML, thet);
    }
  );
  int intComp_result = modSwitchFromTorus32(lwePhase(output, key_BML->extracted_key), 2*context_BML->torus_base);
  printf("Result Compare(%u, %u): ", number1_32bit_clear, number2_32bit_clear);
  if(intComp_result == 1) printf("%u > %u\n", number1_32bit_clear, number2_32bit_clear);
  if(intComp_result == 7) printf("%u < %u\n", number1_32bit_clear, number2_32bit_clear);
  if(intComp_result == 0) printf("%u == %u\n", number1_32bit_clear, number2_32bit_clear);


  /*6-bit*/
  int lut_input_clear_new = rand()&0x3F;
  fbt_integer * lut_input_new = fbt_new_encrypted_integer_BML(lut_input_clear_new, 4, (LweKey *) key_BML->tfhe_keys->lwe_key, context_BML);
  fbt_integer * lut_output_new = fbt_new_encrypted_integer_BML(lut_input_clear_new, 4, key_BML->extracted_key, context_BML);
  printf("\nRunning 6-bit-to-6-bit LUT evaluation %d times. The LUT encodes the identity function.------------------\n", NUM_EXE);
  MEASURE_TIME_AND_PRINT(
  for (size_t i = 0; i < NUM_EXE; i++)
    {
      test_LUT_6_bit_to_6_bit_BML(lut_output_new, lut_input_new, context_BML, thet);
    }
  );
  printf("Result LUT(%d): %d\n", lut_input_clear_new, fbt_integer_decrypt_BML(lut_output_new, key_BML->extracted_key, context_BML)); 



  int8_t lut_input_clear_new_8 = rand()&0xFF;
  fbt_integer * lut_input_new_8 = fbt_new_encrypted_integer_BML(lut_input_clear_new_8, 4, (LweKey *) key_BML->tfhe_keys->lwe_key, context_BML);
  fbt_integer * lut_output_new_8 = fbt_new_encrypted_integer_BML(lut_input_clear_new_8, 4, key_BML->extracted_key, context_BML);
  printf("\nRunning 8-bit-to-8-bit LUT evaluation %d times. The LUT encodes the identity function.------------------\n", NUM_EXE);
  MEASURE_TIME_AND_PRINT(
  for (size_t i = 0; i < NUM_EXE; i++)
    {
      test_LUT_8_bit_to_8_bit_BML(lut_output_new_8, lut_input_new_8, context_BML, thet);
    }
  );
  printf("Result LUT(%d): %d\n", lut_input_clear_new_8, fbt_integer_decrypt_BML(lut_output_new_8, key_BML->extracted_key, context_BML)); 


  int16_t lut_input_clear_new_10 = rand()&0x3FF;
  fbt_integer * lut_input_new_10 = fbt_new_encrypted_integer_BML(lut_input_clear_new_10, 6, (LweKey *) key_BML->tfhe_keys->lwe_key, context_BML);
  fbt_integer * lut_output_new_10 = fbt_new_encrypted_integer_BML(lut_input_clear_new_10, 6, key_BML->extracted_key, context_BML);
  printf("\nRunning 10-bit-to-10-bit LUT evaluation %d times. The LUT encodes the identity function.------------------\n", NUM_EXE);
  MEASURE_TIME_AND_PRINT(
  for (size_t i = 0; i < NUM_EXE; i++)
    {
      test_LUT_10_bit_to_10_bit_BML(lut_output_new_10, lut_input_new_10, context_BML, thet);
    }
  );
  printf("Result LUT(%d): %d\n", lut_input_clear_new_10, fbt_integer_decrypt_BML(lut_output_new_10, key_BML->extracted_key, context_BML)); 



  int16_t lut_input_clear_new_12 = rand()&0xFFF;
  fbt_integer * lut_input_new_12 = fbt_new_encrypted_integer_BML(lut_input_clear_new_12, 8, (LweKey *) key_BML->tfhe_keys->lwe_key, context_BML);
  fbt_integer * lut_output_new_12 = fbt_new_encrypted_integer_BML(lut_input_clear_new_12, 8, key_BML->extracted_key, context_BML);
  printf("\nRunning 12-bit-to-12-bit LUT evaluation %d times. The LUT encodes the identity function.------------------\n", NUM_EXE);
  MEASURE_TIME_AND_PRINT(
  for (size_t i = 0; i < NUM_EXE; i++)
    {
      test_LUT_12_bit_to_12_bit_BML(lut_output_new_12, lut_input_new_12, context_BML, thet);
    }
  );
  printf("Result LUT(%d): %d\n", lut_input_clear_new_12, fbt_integer_decrypt_BML(lut_output_new_12, key_BML->extracted_key, context_BML)); 



  
  return 0;
}
