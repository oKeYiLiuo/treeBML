#include <stdio.h>
#include <unistd.h>
#include "BootManyLUT.cpp"

/* 
export LD_LIBRARY_PATH="/home/liu/cuFHE/cufhe/include/"  export LD_LIBRARY_PATH="/home/liu/cuFHE/cufhe/lib" g++ test_api_gpu.cu -I"/home/liu/cuFHE/cufhe/include/" -L"/home/liu/cuFHE/cufhe/lib" -o main
export LD_LIBRARY_PATH="/home/liu/cuFHE/cufhe/lib"

step1: 
  g++ main.cpp -O3 -std=c++11 -funroll-all-loops -march=native -I"../tfhe/src/include/" -L"../tfhe/build/libtfhe/" -ltfhe-spqlios-fma -lm -g -o main -pthread
step2: 
  ./main 2
*/

void MVB_timeSurvey(fbt_integer * output, fbt_integer * input, fbt_context_BML * context, const int thet){
  const TFheGateBootstrappingParameterSet * params = context->params;
  int f_l3a[3] = {fbt_op_load_hardcoded_LUT, 0, 0}; //Let LUTs into the lut_input integer array in order
  int f_l3b[4] = {fbt_op_generate_mv_b_LUT, 0, 0, 192}; //Let lut_input into a polynomial according to rules
  int f_l3c[1] = {fbt_op_mv_f_bootstrap_init};  //Changing the integer coefficient polynomial to a noiseless TLWE ciphertext
  int f_l3d[4] = {fbt_op_mv_f_bootstrap, 0, 0, 48}; //MVB
  /*
    if f_l1b[5] == 0, 
      carry out f_l1b[3] times bootstrappings with thet
    else{
      carry out f_l1b[3] - 1 times bootstrappings with thet
      carry out 1 times bootstrapping with thet = f_lb[4]
    } 
      
  */      
  int * function[4] = {f_l3a, f_l3b, f_l3c, f_l3d};
  fbt_main_loop_BML(function, 4, input->lwe_samples, context, thet);
}

void parallel_MVB_timeSurvey(fbt_integer * output, fbt_integer * input, fbt_context_BML * context, const int thet){
  const TFheGateBootstrappingParameterSet * params = context->params;
  int f_l3a[3] = {fbt_op_load_hardcoded_LUT, 0, 0}; //Let LUTs into the lut_input integer array in order
  int f_l3b[4] = {fbt_op_generate_mv_b_LUT, 0, 0, 192}; //Let lut_input into a polynomial according to rules
  int f_l3c[1] = {fbt_op_mv_f_bootstrap_init};  //Changing the integer coefficient polynomial to a noiseless TLWE ciphertext
  int f_l3d[4] = {fbt_op_mv_f_bootstrap_parallel, 0, 0, 48}; //MVB
  /*
    if f_l1b[5] == 0, 
      carry out f_l1b[3] times bootstrappings with thet
    else{
      carry out f_l1b[3] - 1 times bootstrappings with thet
      carry out 1 times bootstrapping with thet = f_lb[4]
    } 
      
  */      
  int * function[4] = {f_l3a, f_l3b, f_l3c, f_l3d};
  fbt_main_loop_BML(function, 4, input->lwe_samples, context, thet);

}

void BML_timeSurvey(fbt_integer * output, fbt_integer * input, fbt_context_BML * context, const int thet){
  const TFheGateBootstrappingParameterSet * params = context->params;
  
  int f_l2b[6] = {fbt_op_f_bootstrap_BML, 0, 0, 1, thet, 0};  

  
  int * function[1] = {f_l2b};


  fbt_main_loop_BML(function, 1, input->lwe_samples, context, thet);

}


void KS_timeSurvey(fbt_integer * output, fbt_integer * input, fbt_context_BML * context, const int thet){
  const TFheGateBootstrappingParameterSet * params = context->params;


  int f_l3a[4] = {fbt_op_tlwe_keyswitch_BML_conclusion, 0, 0, (1<<thet)*context->torus_base}; 
  

  int * function[1] = {f_l3a};


  fbt_main_loop_BML(function, 1, input->lwe_samples, context, thet);
}

void BlindRotate_timeSurveyFunction(fbt_integer * output, fbt_integer * input, fbt_context_BML * context, const int thet){
  const TFheGateBootstrappingParameterSet * params = context->params;


  int f_l3a[6] = {BlindRotate_timeSurvey, 0, 0, 1, thet, 0}; 
  

  int * function[1] = {f_l3a};


  fbt_main_loop_BML(function, 1, input->lwe_samples, context, thet);
}

#ifndef NUM_EXE
#define NUM_EXE 1
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
  
  /*MVB_timeSurvey*/
  int lut_input_clear_new_MVB = rand()&0x3F;
  fbt_integer * lut_input_new_MVB = fbt_new_encrypted_integer_BML(lut_input_clear_new_MVB, 4, (LweKey *) key_BML->tfhe_keys->lwe_key, context_BML);
  fbt_integer * lut_output_new_MVB = fbt_new_encrypted_integer_BML(lut_input_clear_new_MVB, 4, key_BML->extracted_key, context_BML);
  printf("\nRunning MVB_timeSurvey %d times. \n", NUM_EXE);
  MEASURE_TIME_AND_PRINT(
  for (size_t i = 0; i < NUM_EXE; i++)
    {
      MVB_timeSurvey(lut_output_new_MVB, lut_input_new_MVB, context_BML, thet);
    }
  );

  /*parallel_MVB_timeSurvey*/
  int lut_input_clear_new_parallel_MVB = rand()&0x3F;
  fbt_integer * lut_input_new_parallel_MVB = fbt_new_encrypted_integer_BML(lut_input_clear_new_parallel_MVB, 4, (LweKey *) key_BML->tfhe_keys->lwe_key, context_BML);
  fbt_integer * lut_output_new_parallel_MVB = fbt_new_encrypted_integer_BML(lut_input_clear_new_parallel_MVB, 4, key_BML->extracted_key, context_BML);
  printf("\nRunning parallel_MVB_timeSurvey %d times. \n", NUM_EXE);
  MEASURE_TIME_AND_PRINT(
  for (size_t i = 0; i < NUM_EXE; i++)
    {
      parallel_MVB_timeSurvey(lut_output_new_parallel_MVB, lut_input_new_parallel_MVB, context_BML, thet);
    }
  );

  /*BlindRotate_timeSurvey*/
  int lut_input_clear_new_BlindRotate = rand()&0x3F;
  fbt_integer * lut_input_new_BlindRotate = fbt_new_encrypted_integer_BML(lut_input_clear_new_BlindRotate, 4, (LweKey *) key_BML->tfhe_keys->lwe_key, context_BML);
  fbt_integer * lut_output_new_BlindRotate = fbt_new_encrypted_integer_BML(lut_input_clear_new_BlindRotate, 4, key_BML->extracted_key, context_BML);
  printf("\nRunning BlindRotate_timeSurvey %d times. \n", NUM_EXE);
  MEASURE_TIME_AND_PRINT(
  for (size_t i = 0; i < NUM_EXE; i++)
    {
      BlindRotate_timeSurveyFunction(lut_output_new_BlindRotate, lut_input_new_BlindRotate, context_BML, thet);
    }
  );

  /*BML_timeSurvey*/
  int lut_input_clear_new_BML = rand()&0x3F;
  fbt_integer * lut_input_new_BML = fbt_new_encrypted_integer_BML(lut_input_clear_new_BML, 4, (LweKey *) key_BML->tfhe_keys->lwe_key, context_BML);
  fbt_integer * lut_output_new_BML = fbt_new_encrypted_integer_BML(lut_input_clear_new_BML, 4, key_BML->extracted_key, context_BML);
  printf("\nRunning BML_timeSurvey %d times. \n", NUM_EXE);
  MEASURE_TIME_AND_PRINT(
  for (size_t i = 0; i < NUM_EXE; i++)
    {
      BML_timeSurvey(lut_output_new_BML, lut_input_new_BML, context_BML, thet);
    }
  );


  

  /*KS_timeSurvey*/
  int lut_input_clear_new_KS = rand()&0x3F;
  fbt_integer * lut_input_new_KS = fbt_new_encrypted_integer_BML(lut_input_clear_new_KS, 4, (LweKey *) key_BML->tfhe_keys->lwe_key, context_BML);
  fbt_integer * lut_output_new_KS = fbt_new_encrypted_integer_BML(lut_input_clear_new_KS, 4, key_BML->extracted_key, context_BML);
  printf("\nRunning KS_timeSurvey %d times. \n", NUM_EXE);
  MEASURE_TIME_AND_PRINT(
  for (size_t i = 0; i < NUM_EXE; i++)
    {
      KS_timeSurvey(lut_output_new_KS, lut_input_new_KS, context_BML, thet);
    }
  );



  
  return 0;
}
