#include <stdio.h>
#include <unistd.h>
#include "BootManyLUT.cpp"

/* 
export LD_LIBRARY_PATH="/home/liu/cuFHE/cufhe/include/"  export LD_LIBRARY_PATH="/home/liu/cuFHE/cufhe/lib" g++ test_api_gpu.cu -I"/home/liu/cuFHE/cufhe/include/" -L"/home/liu/cuFHE/cufhe/lib" -o main
export LD_LIBRARY_PATH="/home/liu/cuFHE/cufhe/lib"

step1: 
  g++ Tree-BML-mainPipeline.cpp -O3 -std=c++11 -funroll-all-loops -march=native -I"../tfhe/src/include/" -L"../tfhe/build/libtfhe/" -ltfhe-spqlios-fma -lm -g -o TreeBMLmainPipeline -pthread
step2: 
  ./TreeBMLmainPipeline 2
*/

void test_LUT_bit_to_bit_BML_pipeline(vector<fbt_integer *> &outputNums, vector<fbt_integer *> &inputNums, vector<fbt_context_BML *> &contextNums, 
int thet, const int product, int * functionValue, const int layers, vector<int> function_KS_i_3, vector<int> function_BML_i_3, vector<int> function_BML_i_4, vector<int> function_BML_i_5){
  const TFheGateBootstrappingParameterSet * params = contextNums[0]->params;
  
  const TLweParams * tlwe_params = contextNums[0]->params->tgsw_params->tlwe_params;
  const LweParams * lwe_params = contextNums[0]->params->in_out_params;
  const LweParams * extracted_params = &tlwe_params->extracted_lweparams;
  vector<int> current_input(product, 0);
  const int log_carpov_base = 1;
  vector<TLweSample *> TV_0_nums;
  vector<LweSample *> temp_nums;
  vector<IntPolynomial *> lut_clear_nums(product);
  int * lut_input, number_of_luts;

  int function_fbt_op_load_hardcoded_LUT_2 = functionValue[0];//第几个查找表
  int function_fbt_op_generate_mv_b_LUT_3 = functionValue[1];//查找表大小
  int function_fbt_op_mv_f_bootstrap_1 = functionValue[2];
  int function_fbt_op_mv_f_bootstrap_3 = functionValue[3];

  lut_input = BML_LUT_list[function_fbt_op_load_hardcoded_LUT_2];

  //公用参数
  const LweBootstrappingKeyFFT * bootstrap_key = contextNums[0]->bootstrap_key;
  const int torus_base = contextNums[0]->torus_base;
  //fbt_op_generate_mv_b_LUT
  for(int i = 0; i < product; ++i){
      IntPolynomial * luts;
      number_of_luts = generate_polynomial_LUT(&luts, lut_input, function_fbt_op_generate_mv_b_LUT_3, torus_base, tlwe_params->N);
      lut_clear_nums[i] = new_IntPolynomial_array((contextNums[i]->log_torus_base/log_carpov_base)*(function_fbt_op_generate_mv_b_LUT_3/torus_base), tlwe_params->N);
      for (size_t j = 0; j < function_fbt_op_generate_mv_b_LUT_3/torus_base; j++){
        carpov_factorization(&lut_clear_nums[i][(contextNums[i]->log_torus_base/log_carpov_base)*j], &luts[j], contextNums[i]->log_torus_base, log_carpov_base);
      }
      delete_IntPolynomial_array(number_of_luts, luts);
  }

  
  //fbt_op_mv_f_bootstrap_init
  for(int i = 0; i < product; ++i){
    TV_0_nums.push_back(new_TLweSample(tlwe_params));
  }
  thread Mythread[product];
  for(int i = 0; i < product; ++i){
    Mythread[i] = thread(tfhe_carpov_multivalue_bootstrap_init,
    TV_0_nums[i], bootstrap_key, &(inputNums[i]->lwe_samples)[current_input[i]], torus_base);
  }
  for(int i = 0; i < product; ++i){
    Mythread[i].join();
  }

  //fbt_op_mv_f_bootstrap
  thread Mythread_mv[product];
  for(int i = 0; i < product; ++i){
    //Mythread_mv[i] = thread(carpov_second_phase_with_base_composition, &(contextNums[i])->output[function_fbt_op_mv_f_bootstrap_1], lut_clear_nums[i],
    //TV_0_nums[i], function_fbt_op_mv_f_bootstrap_3, tlwe_params, contextNums[0]->log_torus_base, log_carpov_base);
    carpov_second_phase_with_base_composition(&(contextNums[i])->output[function_fbt_op_mv_f_bootstrap_1], lut_clear_nums[i],
    TV_0_nums[i], function_fbt_op_mv_f_bootstrap_3, tlwe_params, contextNums[0]->log_torus_base, log_carpov_base);
  }
  for(int i = 0; i < product; ++i){
    //Mythread_mv[i].join();
  }
  for(int i = 0; i < product; ++i){
    delete_IntPolynomial_array(number_of_luts, lut_clear_nums[i]);
  }


  //PipelineCaiculationCore
  for (size_t i = 0; i < product + layers; i++)
  {
    thread Mythread_core[layers];
    for(int j = 0; j < layers; ++j){
      if(i - j >= 0 && i - j < product){
        Mythread_core[j] = thread(PipelineCaiculationCore, inputNums[i - j]->lwe_samples, contextNums[i - j], ++current_input[i - j],
        thet, function_KS_i_3[j], function_BML_i_3[j], function_BML_i_4[j], function_BML_i_5[j]);
      }
    }
    for(int j = 0; j < layers; ++j){
      if(i - j >= 0 && i - j < product)
        Mythread_core[j].join();
    }
  }
  

  for(int i = 0; i < product; ++i){
    lweCopy(&outputNums[i]->lwe_samples[0], &contextNums[i]->output[0], &params->tgsw_params->tlwe_params->extracted_lweparams);
    lweCopy(&outputNums[i]->lwe_samples[1], &contextNums[i]->output[1], &params->tgsw_params->tlwe_params->extracted_lweparams);
    lweCopy(&outputNums[i]->lwe_samples[2], &contextNums[i]->output[2], &params->tgsw_params->tlwe_params->extracted_lweparams);

  }
}



#ifndef NUM_EXE
#define NUM_EXE 10
#endif

int NUM_PRODUCT =  10;

int main(int argc, char const *argv[])
{
  const int thet = atoi(argv[1]);
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
  


  for(int sss = 0; sss < 4; sss++){
    if(sss == 0){
      /*6-bit*/
  vector<fbt_integer *> lut_input_new_nums;
  vector<fbt_integer *> lut_out_new_nums;
  vector<int> lut_input_clear_new_nums;
  vector<fbt_context_BML *> context_BML_nums;
  int functionValue[4] = {0, 192, 0, 48};
  const int layers = 2;
  const int torus_base = 4;
  vector<int> function_KS_i_3 = {functionValue[3], functionValue[3]/torus_base};
  vector<int> function_BML_i_3 = {functionValue[3]/((1<<thet)*torus_base), 1};
  vector<int> function_BML_i_4 = {thet, 2};
  vector<int> function_BML_i_5 = {0, 1};

  if(thet == 1){
    function_BML_i_3[1] = 2;
    function_BML_i_4[1] = 0;
    function_BML_i_5[1] = 1;
  } else if(thet == 3){
    function_BML_i_3[0] = 2;
    function_BML_i_4[0] = 2;
    function_BML_i_5[0] = 1;

    function_BML_i_3[1] = 1;
    function_BML_i_4[1] = 2;
    function_BML_i_5[1] = 1;
  } 

  for(int i = 0; i < NUM_PRODUCT; ++i){
     fbt_context_BML * context_BML = fbt_init_BML(key_BML->cloud_key->params, key_BML->cloud_key->bkFFT, key_BML->tlweKS, key_BML->log_torus_base, 6144);
    context_BML_nums.push_back(context_BML);

    int lut_input_clear_new = rand()&0x3F;
    lut_input_clear_new_nums.push_back(lut_input_clear_new);

     fbt_integer * lut_input_new = fbt_new_encrypted_integer_BML(lut_input_clear_new, 4, (LweKey *) key_BML->tfhe_keys->lwe_key, context_BML);
    lut_input_new_nums.push_back(lut_input_new);

     fbt_integer * lut_output_new = fbt_new_encrypted_integer_BML(lut_input_clear_new, 4, key_BML->extracted_key, context_BML);
    lut_out_new_nums.push_back(lut_output_new);
  }
  printf("\nRunning 6-bit-to-6-bit LUT evaluation %d product. The LUT encodes the identity function.------------------\n", NUM_PRODUCT);
  MEASURE_TIME_AND_PRINT(
    test_LUT_bit_to_bit_BML_pipeline(lut_out_new_nums, lut_input_new_nums, context_BML_nums, thet, NUM_PRODUCT, functionValue, layers, function_KS_i_3, 
    function_BML_i_3, function_BML_i_4, function_BML_i_5);
  );
  for(int i = 0; i < NUM_PRODUCT; ++i){
    printf("Result LUT(%d): %d\n", lut_input_clear_new_nums[i], fbt_integer_decrypt_BML(lut_out_new_nums[i], key_BML->extracted_key, context_BML_nums[i])); 
  }
    }
    else if(sss == 1){
      /*8-bit*/
  vector<fbt_integer *> lut_input_new_nums;
  vector<fbt_integer *> lut_out_new_nums;
  vector<int> lut_input_clear_new_nums;
  vector<fbt_context_BML *> context_BML_nums;
  int functionValue[4] = {1, 1024, 0, 256};
  const int layers = 3;
  const int torus_base = 4;
  vector<int> function_KS_i_3 = {functionValue[3], functionValue[3]/torus_base, functionValue[3]/(torus_base*torus_base)};
  vector<int> function_BML_i_3 = {functionValue[3]/((1<<thet)*torus_base), functionValue[3]/(torus_base*(1<<thet)*torus_base), 1};
  vector<int> function_BML_i_4 = {thet, thet, thet};
  vector<int> function_BML_i_5 = {0, 0, 0};


  if(thet == 1){
    function_BML_i_3[2] = 2;
    function_BML_i_4[2] = 1;
    function_BML_i_5[2] = 0;
  } else if(thet == 3){
    function_BML_i_3[2] = 1;
    function_BML_i_4[2] = 2;
    function_BML_i_5[2] = 1;
  } 

  for(int i = 0; i < NUM_PRODUCT; ++i){
     fbt_context_BML * context_BML = fbt_init_BML(key_BML->cloud_key->params, key_BML->cloud_key->bkFFT, key_BML->tlweKS, key_BML->log_torus_base, 6144);
    context_BML_nums.push_back(context_BML);

    int lut_input_clear_new = rand()&0xFF;
    lut_input_clear_new_nums.push_back(lut_input_clear_new);

     fbt_integer * lut_input_new = fbt_new_encrypted_integer_BML(lut_input_clear_new, 4, (LweKey *) key_BML->tfhe_keys->lwe_key, context_BML);
    lut_input_new_nums.push_back(lut_input_new);

     fbt_integer * lut_output_new = fbt_new_encrypted_integer_BML(lut_input_clear_new, 4, key_BML->extracted_key, context_BML);
    lut_out_new_nums.push_back(lut_output_new);
  }
  printf("\nRunning 8-bit-to-8-bit LUT evaluation %d product. The LUT encodes the identity function.------------------\n", NUM_PRODUCT);
  MEASURE_TIME_AND_PRINT(
    test_LUT_bit_to_bit_BML_pipeline(lut_out_new_nums, lut_input_new_nums, context_BML_nums, thet, NUM_PRODUCT, functionValue, layers, function_KS_i_3, 
    function_BML_i_3, function_BML_i_4, function_BML_i_5);
  );
  for(int i = 0; i < NUM_PRODUCT; ++i){
    printf("Result LUT(%d): %d\n", lut_input_clear_new_nums[i], fbt_integer_decrypt_BML(lut_out_new_nums[i], key_BML->extracted_key, context_BML_nums[i])); 
  }

    }
    else if(sss == 2){
      /*10-bit*/
  vector<fbt_integer *> lut_input_new_nums;
  vector<fbt_integer *> lut_out_new_nums;
  vector<int> lut_input_clear_new_nums;
  vector<fbt_context_BML *> context_BML_nums;
  int functionValue[4] = {2, 5120, 0, 1280};
  const int layers = 4;
  const int torus_base = 4;
  vector<int> function_KS_i_3 = {functionValue[3], functionValue[3]/torus_base, functionValue[3]/(torus_base*torus_base), functionValue[3]/(torus_base*torus_base*torus_base)};
  vector<int> function_BML_i_3 = {functionValue[3]/((1<<thet)*torus_base), functionValue[3]/(torus_base*(1<<thet)*torus_base), functionValue[3]/(torus_base*torus_base*(1<<thet)*torus_base), 2};
  vector<int> function_BML_i_4 = {thet, thet, thet, 0};
  vector<int> function_BML_i_5 = {0, 0, 0, 1};

  if(thet == 1){
    function_BML_i_3[3] = 3;
    function_BML_i_4[3] = 0;
    function_BML_i_5[3] = 1;
  } else if(thet == 3){
    function_BML_i_3[2] = 3;
    function_BML_i_4[2] = 2;
    function_BML_i_5[2] = 1;

    function_BML_i_3[3] = 1;
    function_BML_i_4[3] = 3;
    function_BML_i_5[3] = 1;
  } 

  for(int i = 0; i < NUM_PRODUCT; ++i){
     fbt_context_BML * context_BML = fbt_init_BML(key_BML->cloud_key->params, key_BML->cloud_key->bkFFT, key_BML->tlweKS, key_BML->log_torus_base, 6144);
    context_BML_nums.push_back(context_BML);

    int16_t lut_input_clear_new = rand()&0x3FF;
    lut_input_clear_new_nums.push_back(lut_input_clear_new);

     fbt_integer * lut_input_new = fbt_new_encrypted_integer_BML(lut_input_clear_new, 6, (LweKey *) key_BML->tfhe_keys->lwe_key, context_BML);
    lut_input_new_nums.push_back(lut_input_new);

     fbt_integer * lut_output_new = fbt_new_encrypted_integer_BML(lut_input_clear_new, 6, key_BML->extracted_key, context_BML);
    lut_out_new_nums.push_back(lut_output_new);
  }
  printf("\nRunning 10-bit-to-10-bit LUT evaluation %d product. The LUT encodes the identity function.------------------\n", NUM_PRODUCT);
  MEASURE_TIME_AND_PRINT(
    test_LUT_bit_to_bit_BML_pipeline(lut_out_new_nums, lut_input_new_nums, context_BML_nums, thet, NUM_PRODUCT, functionValue, layers, function_KS_i_3, 
    function_BML_i_3, function_BML_i_4, function_BML_i_5);
  );
  for(int i = 0; i < NUM_PRODUCT; ++i){
    printf("Result LUT(%d): %d\n", lut_input_clear_new_nums[i], fbt_integer_decrypt_BML(lut_out_new_nums[i], key_BML->extracted_key, context_BML_nums[i])); 
  }

    }
    else if(sss == 3){
      /*12-bit*/
  vector<fbt_integer *> lut_input_new_nums;
  vector<fbt_integer *> lut_out_new_nums;
  vector<int> lut_input_clear_new_nums;
  vector<fbt_context_BML *> context_BML_nums;
  int functionValue[4] = {3, 24576, 0, 6144};
  const int layers = 5;
  const int torus_base = 4;
  vector<int> function_KS_i_3 = {functionValue[3], functionValue[3]/torus_base, functionValue[3]/(torus_base*torus_base), functionValue[3]/(torus_base*torus_base*torus_base), functionValue[3]/(torus_base*torus_base*torus_base*torus_base), functionValue[3]/(torus_base*torus_base*torus_base*torus_base*torus_base)};
  vector<int> function_BML_i_3 = {functionValue[3]/((1<<thet)*torus_base), functionValue[3]/(torus_base*(1<<thet)*torus_base), functionValue[3]/(torus_base*torus_base*(1<<thet)*torus_base), functionValue[3]/(torus_base*torus_base*torus_base*(1<<thet)*torus_base), 2};
  vector<int> function_BML_i_4 = {thet, thet, thet, thet, 1};
  vector<int> function_BML_i_5 = {0, 0, 0, 0, 1};

  if(thet == 1){
    function_BML_i_3[4] = 3;
    function_BML_i_4[4] = 1;
    function_BML_i_5[4] = 0;
  } else if(thet == 3){
    function_BML_i_3[4] = 1;
    function_BML_i_4[4] = 3;
    function_BML_i_5[4] = 1;
  }

  for(int i = 0; i < NUM_PRODUCT; ++i){
     fbt_context_BML * context_BML = fbt_init_BML(key_BML->cloud_key->params, key_BML->cloud_key->bkFFT, key_BML->tlweKS, key_BML->log_torus_base, 6144);
    context_BML_nums.push_back(context_BML);

    int lut_input_clear_new = rand()&0xFFF;
    lut_input_clear_new_nums.push_back(lut_input_clear_new);

     fbt_integer * lut_input_new = fbt_new_encrypted_integer_BML(lut_input_clear_new, 8, (LweKey *) key_BML->tfhe_keys->lwe_key, context_BML);
    lut_input_new_nums.push_back(lut_input_new);

     fbt_integer * lut_output_new = fbt_new_encrypted_integer_BML(lut_input_clear_new, 8, key_BML->extracted_key, context_BML);
    lut_out_new_nums.push_back(lut_output_new);
  }
  printf("\nRunning 12-bit-to-12-bit LUT evaluation %d product. The LUT encodes the identity function.------------------\n", NUM_PRODUCT);
  MEASURE_TIME_AND_PRINT(
    test_LUT_bit_to_bit_BML_pipeline(lut_out_new_nums, lut_input_new_nums, context_BML_nums, thet, NUM_PRODUCT, functionValue, layers, function_KS_i_3, 
    function_BML_i_3, function_BML_i_4, function_BML_i_5);
  );
  for(int i = 0; i < NUM_PRODUCT; ++i){
    printf("Result LUT(%d): %d\n", lut_input_clear_new_nums[i], fbt_integer_decrypt_BML(lut_out_new_nums[i], key_BML->extracted_key, context_BML_nums[i])); 
  }

    }
  }
  

  /*
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
  */


  
  return 0;
}
