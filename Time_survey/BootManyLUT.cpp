#include "../Tree-BML/libBML.cpp"
#include<thread>

enum TreeBML_operation {
  fbt_op_tlwe_keyswitch_BML_conclusion = 13,
  fbt_op_f_bootstrap_BML = 14,
  BlindRotate_timeSurvey = 15,
  fbt_op_mv_f_bootstrap_parallel = 16
  } TreeBML_operation;

typedef struct fbt_context_BML {
  const TFheGateBootstrappingParameterSet * params;
  const LweBootstrappingKeyFFT * bootstrap_key;
  const TLweKeySwitchKey_BML * tlweKS;
  LweSample * input, * output;
  TLweSample * tlwe_temp;
  LweSample * lwe_temp_in;
  int log_torus_base, torus_base;
} fbt_context_BML;

fbt_context_BML * fbt_init_BML(const TFheGateBootstrappingParameterSet * params,
                       const LweBootstrappingKeyFFT * bootstrap_key,
                       const TLweKeySwitchKey_BML * tlweKS, int log_torus_base,
                       const int size_of_temps){
  fbt_context_BML * context = (fbt_context_BML *) malloc(sizeof(fbt_context_BML));
  context->params = params;
  context->bootstrap_key = bootstrap_key;
  context->tlweKS = tlweKS;

  context->log_torus_base = log_torus_base;
  context->torus_base = 1 << log_torus_base;

  const TLweParams * tlwe_params = context->params->tgsw_params->tlwe_params;
  const LweParams * extracted_params = &tlwe_params->extracted_lweparams;
  context->input = new_LweSample_array(size_of_temps, extracted_params);
  context->output = new_LweSample_array(size_of_temps, extracted_params);
  context->tlwe_temp = new_TLweSample_array(size_of_temps, tlwe_params);
  context->lwe_temp_in = new_LweSample_array(size_of_temps, params->in_out_params);
  return context;
}

void carpov_second_phase_with_base_composition_Parallel(LweSample * results, IntPolynomial * tv_f, TLweSample * rotated_tv_0, int q, const TLweParams * params, int log_torus_base, const int log_carpov_base);

int lut_test_6bit[192] = {0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,
                      0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3};

int * BML_LUT_list[2] = {lut_test_6bit, lut_int_comp_base4};

void fbt_main_loop_BML(int ** function, const int size, LweSample * input, fbt_context_BML * context, const int thet){
  const TLweParams * tlwe_params = context->params->tgsw_params->tlwe_params;
  const LweParams * lwe_params = context->params->in_out_params;
  const LweParams * extracted_params = &tlwe_params->extracted_lweparams;
  int current_input = 0;
  const int log_carpov_base = 1;
  
  TLweSample * TV_0 = new_TLweSample(tlwe_params);
  LweSample * temp;
  IntPolynomial * lut_clear;
  int * lut_input, number_of_luts;

  for (size_t i = 0; i < size; i++)
  {
    switch (function[i][0])
    {

    case fbt_op_mv_f_bootstrap:
      carpov_second_phase_with_base_composition(&context->output[function[i][1]], lut_clear, TV_0, function[i][3], tlwe_params, context->log_torus_base, log_carpov_base);
      delete_IntPolynomial_array(number_of_luts, lut_clear);
      break;
    
    case fbt_op_mv_f_bootstrap_init:
      tfhe_carpov_multivalue_bootstrap_init(TV_0, context->bootstrap_key, &input[current_input], context->torus_base);
      break;

    case fbt_op_generate_mv_b_LUT:{
      IntPolynomial * luts;
      number_of_luts = generate_polynomial_LUT(&luts, lut_input, function[i][3], context->torus_base, tlwe_params->N);
      lut_clear = new_IntPolynomial_array((context->log_torus_base/log_carpov_base)*(function[i][3]/context->torus_base), tlwe_params->N);
      for (size_t j = 0; j < function[i][3]/context->torus_base; j++){
        carpov_factorization(&lut_clear[(context->log_torus_base/log_carpov_base)*j], &luts[j], context->log_torus_base, log_carpov_base);
      }
      delete_IntPolynomial_array(number_of_luts, luts);
    }
      break;

    case fbt_op_load_hardcoded_LUT:
      lut_input = BML_LUT_list[function[i][2]];
      break;

      
    case fbt_op_next_tree_level:{
      current_input+=1;
      temp = context->input;
      context->input = context->output;
      context->input->current_variance = context->output->current_variance;
      context->output = temp;

    }
      break;


    case fbt_op_f_bootstrap_BML:{
      thread mythread[function[i][3]];
      if(function[i][5] == 0){//no padding or divisi
        for (size_t j = 0; j < function[i][3]; j++)
        {
          mythread[j] = thread(tfhe_functional_bootstrap_woKS_FFT_BML, 
          &context->output[function[i][1]+ j * (1<<thet)], 
          context->bootstrap_key, &context->tlwe_temp[function[i][2] + j], 
          &input[current_input], context->torus_base, thet);
        }
        for (size_t j = 0; j < function[i][3]; j++)
        {
          mythread[j].join();
        }
      }
      else{
        for (size_t j = 0; j < function[i][3] - 1; j++)
        {
          mythread[j] = thread(tfhe_functional_bootstrap_woKS_FFT_BML, 
          &context->output[function[i][1]+ j * (1<<thet)], 
          context->bootstrap_key, &context->tlwe_temp[function[i][2] + j], 
          &input[current_input], context->torus_base, thet);
        }
        mythread[function[i][3] - 1] = thread(tfhe_functional_bootstrap_woKS_FFT_BML,
        &context->output[function[i][1]+ (function[i][3] - 1) * (1<<thet)], 
        context->bootstrap_key, &context->tlwe_temp[function[i][2] + (function[i][3] - 1)], 
        &input[current_input], context->torus_base, function[i][4]);
        for (size_t j = 0; j < function[i][3]; j++)
        {
          mythread[j].join();
        }
      }
    }
      break;
  
    case fbt_op_tlwe_keyswitch_BML_conclusion:{
      const int full_times = function[i][3]/((1<<thet)*context->torus_base);
      int residue = function[i][3] - full_times*((1<<thet)*context->torus_base);
      thread mythreadKS[full_times + (residue == 0 ? 0 : 1)];

      for (size_t j = 0; j < full_times; j++){
        mythreadKS[j] = thread(public_functional_tlweKeySwitch_BML_FULLThetB, 
        &context->tlwe_temp[(function[i][1] + j)], 
        context->tlweKS, &context->input[(function[i][2] + j*(1<<thet)*context->torus_base)], context->torus_base, thet);
        //public_functional_tlweKeySwitch_BML_FULLThetB(&context->tlwe_temp[(function[i][1] + j)], context->tlweKS, &context->input[(function[i][2] + j*(1<<thet)*context->torus_base)], context->torus_base, thet);
      }
      if(residue != 0){
        if((1<<thet)*context->torus_base % residue == 0){
          mythreadKS[full_times] = thread(public_functional_tlweKeySwitch_BML_divisible, &context->tlwe_temp[function[i][1] + full_times], context->tlweKS, &context->input[function[i][2] + full_times*((1<<thet)*context->torus_base)], context->torus_base, thet, residue);
          //public_functional_tlweKeySwitch_BML_divisible(&context->tlwe_temp[function[i][1] + full_times], context->tlweKS, &context->input[function[i][2] + full_times*((1<<thet)*context->torus_base)], context->torus_base, thet, residue);
        }
        else{
          mythreadKS[full_times] = thread(public_functional_tlweKeySwitch_BML_padding_then_divis, &context->tlwe_temp[function[i][1] + full_times], context->tlweKS, &context->input[function[i][2] + full_times*((1<<thet)*context->torus_base)], context->torus_base, thet, residue);
          //public_functional_tlweKeySwitch_BML_padding_then_divis(&context->tlwe_temp[function[i][1] + full_times], context->tlweKS, &context->input[function[i][2] + full_times*((1<<thet)*context->torus_base)], context->torus_base, thet, residue);
        }
      }
      for(auto k = 0; k < (full_times + (residue == 0 ? 0 : 1)); ++k){
        mythreadKS[k].join();
      }
    }
    break;

    case fbt_op_load_literal:
      for (size_t j = 0; j < function[i][3]; j++){
        lweNoiselessTrivial(&context->output[function[i][1] + j], modSwitchToTorus32(function[i][2], 2*context->torus_base), extracted_params);
      }
      break;

    case fbt_op_f_bootstrap:
    {
      for (size_t j = 0; j < function[i][3]; j++){
        tfhe_functional_bootstrap_woKS_FFT(&context->output[function[i][1] + j], context->bootstrap_key, &context->tlwe_temp[function[i][2] + j], &input[current_input], context->torus_base);
      }
    }
      break;
      
    case BlindRotate_timeSurvey:
    {
      //tfhe_functional_bootstrap_woKS_FFT(&context->output[function[i][1]], context->bootstrap_key, &context->tlwe_temp[function[i][2]], &input[current_input], context->torus_base);
      const TGswSampleFFT *bk = context->bootstrap_key->bkFFT;
      const TLweParams *accum_params = context->bootstrap_key->bk_params->tlwe_params;
      int32_t _2N = 2 * context->bootstrap_key->bk_params->tlwe_params->N;
      
      int32_t barb = (&input[current_input])->b%_2N;
      int32_t n = context->bootstrap_key->in_out_params->n;
      TLweSample *acc = new_TLweSample(accum_params);
      int32_t *bara = new int32_t[context->bootstrap_key->bk_params->tlwe_params->N];
      TLweSample *v = &context->tlwe_temp[function[i][2]];
      for (int32_t i = 0; i < n; i++) {
        bara[i] = rand()%(_2N);
      }


      // testvector = X^{2N-barb}*v
      if (barb != 0){
        for (size_t i = 0; i <= accum_params->k; i++){
          torusPolynomialMulByXai(&acc->a[i], _2N - barb, &v->a[i]);
        }
      }else{
        for (size_t i = 0; i <= accum_params->k; i++){
          torusPolynomialCopy(&acc->a[i], &v->a[i]);
        }
      } 
      
      // Blind rotation
      tfhe_blindRotate_FFT(acc, bk, bara, n, context->bootstrap_key->bk_params);

      //delete[] bara;
      //delete_TLweSample(acc);
      break;
    }

    case fbt_op_mv_f_bootstrap_parallel:{
      carpov_second_phase_with_base_composition_Parallel(&context->output[function[i][1]], lut_clear, TV_0, function[i][3], tlwe_params, context->log_torus_base, log_carpov_base);
      delete_IntPolynomial_array(number_of_luts, lut_clear);
      break;
    }

    default:
      break;
    }
  }
}



fbt_integer * fbt_new_encrypted_integer_BML(int value, int digits, LweKey * key, fbt_context_BML * context){
  fbt_integer * result = (fbt_integer*)malloc(sizeof(fbt_integer));
  result->lwe_samples = encryptInteger(value, digits, context->log_torus_base, key);
  result->digits = digits;
  result->lwe_params = key->params;
  result->log_torus_base = context->log_torus_base;
  return result;
}

int fbt_integer_decrypt_BML(fbt_integer * input, LweKey * key, fbt_context_BML * context){
  return decryptInteger(input->lwe_samples, input->digits, context->log_torus_base, key);
}


void parallel_CMulCipher(LweSample * results, IntPolynomial * tv_f, TLweSample * rotated_tv_0, int q, const TLweParams * params, int log_torus_base, const int log_carpov_base, const int number_of_elements,
TLweSample * result_tlwe, LweSample * res_temp, int i){
  lweNoiselessTrivial(&results[i], 0, &params->extracted_lweparams);
  for (size_t j = 0; j < number_of_elements; j++){
    tLweNoiselessTrivialT(result_tlwe, (const Torus32) 0, params);
    //printf("*********************************************\n");
    //MEASURE_TIME_AND_PRINT(tLweAddMulRTo(result_tlwe, &tv_f[i*number_of_elements + j], rotated_tv_0, params););
    //printf("*********************************************\n");
    tLweAddMulRTo(result_tlwe, &tv_f[i*number_of_elements + j], rotated_tv_0, params);
    for (size_t s = 0; s < (1<<(j*log_carpov_base)); s++){
      tLweExtractLweSampleIndex(res_temp, result_tlwe, s, &params->extracted_lweparams, params);
      lweAddTo(&results[i], res_temp, &params->extracted_lweparams);
    }
  }    
}

void carpov_second_phase_with_base_composition_Parallel(LweSample * results, IntPolynomial * tv_f, TLweSample * rotated_tv_0, int q, const TLweParams * params, int log_torus_base, const int log_carpov_base){
  TLweSample * result_tlwe = new_TLweSample(params);
  LweSample * res_temp = new_LweSample(&params->extracted_lweparams);

  const int number_of_elements = log_torus_base/log_carpov_base;
  thread Mythread[q];

  for (size_t i = 0; i < q; i++){
    Mythread[i] = thread(parallel_CMulCipher, results, tv_f, rotated_tv_0, q, params, log_torus_base, log_carpov_base, number_of_elements,
    result_tlwe, res_temp, i); 
  }
  for (size_t i = 0; i < q; i++){
    Mythread[i].join(); 
  }
  delete_TLweSample(result_tlwe);
  delete_LweSample(res_temp);
}