#include <stdio.h>
#include <unistd.h>
#include "./functional_bootstrap.cpp"

/*
If you want to calculate 8-bit, 10-bit or 12-bit LUT, 
  g++ Tree-FB-main.cpp -O3 -std=c++11 -funroll-all-loops -march=native -I"../tfhe/src/include/" -L"../tfhe/build/libtfhe/" -ltfhe-spqlios-fma -lm -g -o TreeFB
*/

void test_LUT_error(fbt_integer * output, fbt_integer * input, fbt_context * context){

 tfhe_functional_bootstrap_woKS_FFT(&context->output[0], context->bootstrap_key, &context->tlwe_temp[0], &input->lwe_samples[0], context->torus_base);
 lweCopy(&output->lwe_samples[0], &context->output[0], &context->params->tgsw_params->tlwe_params->extracted_lweparams);

}

void test_LUT_6_bit_to_6_bit(fbt_integer * output, fbt_integer * input, fbt_context * context){
  const TFheGateBootstrappingParameterSet * params = context->params;
  int f_l3a[3] = {fbt_op_load_hardcoded_LUT, 0, 0};       // this LUT is hardcoded at line 28 of functional_bootstrap.cpp and implements the identity
  int f_l3b[4] = {fbt_op_generate_mv_b_LUT, 0, 0, 192};   // load the LUTs in polynomials
  int f_l3c[1] = {fbt_op_mv_f_bootstrap_init};            // Carpov multi value bootstrap first phase
  int f_l3d[4] = {fbt_op_mv_f_bootstrap, 0, 0, 48};       // Carpov multi value bootstrap second phase
  int f_l3e[1] = {fbt_op_next_tree_level};                // next level of the tree evaluation
  
  int f_l2a[4] = {fbt_op_tlwe_keyswitch, 0, 0, 12};       // tlwe to trlwe packing.
  int f_l2b[4] = {fbt_op_f_bootstrap, 0, 0, 12};          // 12 functional bootstraps
  int f_l2c[1] = {fbt_op_next_tree_level};                // next level of the tree evaluation

  int f_l1a[4] = {fbt_op_tlwe_keyswitch, 0, 0, 3};        // tlwe to trlwe packing
  int f_l1b[4] = {fbt_op_f_bootstrap, 0, 0, 3};           // 3 functional bootstraps
  int * function[10] = {f_l3a, f_l3b, f_l3c, f_l3d, f_l3e, f_l2a, f_l2b, f_l2c, f_l1a, f_l1b};

  fbt_main_loop(function, 10, input->lwe_samples, context);            // Execute everything defined in "function"
    
  lweCopy(&output->lwe_samples[0], &context->output[0], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&output->lwe_samples[1], &context->output[1], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&output->lwe_samples[2], &context->output[2], &params->tgsw_params->tlwe_params->extracted_lweparams);
}

void test_LUT_8_bit_to_8_bit(fbt_integer * output, fbt_integer * input, fbt_context * context){
  const TFheGateBootstrappingParameterSet * params = context->params;
  int f_l3a[3] = {fbt_op_load_hardcoded_LUT, 0, 6};       // this LUT is hardcoded at line 28 of functional_bootstrap.cpp and implements the identity
  int f_l3b[4] = {fbt_op_generate_mv_b_LUT, 0, 0, 1024};   // load the LUTs in polynomials
  int f_l3c[1] = {fbt_op_mv_f_bootstrap_init};            // Carpov multi value bootstrap first phase
  int f_l3d[4] = {fbt_op_mv_f_bootstrap, 0, 0, f_l3b[3]/context->torus_base};       // Carpov multi value bootstrap second phase
  int f_l3e[1] = {fbt_op_next_tree_level};                // next level of the tree evaluation
  
  int f_l2a[4] = {fbt_op_tlwe_keyswitch, 0, 0, f_l3d[3]/context->torus_base};       // tlwe to trlwe packing.
  int f_l2b[4] = {fbt_op_f_bootstrap, 0, 0, f_l2a[3]};          // 64 functional bootstraps
  int f_l2c[1] = {fbt_op_next_tree_level};                // next level of the tree evaluation

  int f_l1a[4] = {fbt_op_tlwe_keyswitch, 0, 0, f_l2a[3]/context->torus_base};        // tlwe to trlwe packing
  int f_l1b[4] = {fbt_op_f_bootstrap, 0, 0, f_l1a[3]};           // 16 functional bootstraps
  int f_l1c[1] = {fbt_op_next_tree_level};                // next level of the tree evaluation

  int f_l0a[4] = {fbt_op_tlwe_keyswitch, 0, 0, f_l1a[3]/context->torus_base};        // tlwe to trlwe packing
  int f_l0b[4] = {fbt_op_f_bootstrap, 0, 0, f_l0a[3]};           // 4 functional bootstraps

  int * function[13] = {f_l3a, f_l3b, f_l3c, f_l3d, f_l3e, f_l2a, f_l2b, f_l2c, f_l1a, f_l1b, f_l1c, f_l0a, f_l0b};

  fbt_main_loop(function, 13, input->lwe_samples, context);            // Execute everything defined in "function"
    
  lweCopy(&output->lwe_samples[0], &context->output[0], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&output->lwe_samples[1], &context->output[1], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&output->lwe_samples[2], &context->output[2], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&output->lwe_samples[3], &context->output[3], &params->tgsw_params->tlwe_params->extracted_lweparams);
}

void test_LUT_10_bit_to_10_bit(fbt_integer * output, fbt_integer * input, fbt_context * context){
  const TFheGateBootstrappingParameterSet * params = context->params;
  int f_l4a[3] = {fbt_op_load_hardcoded_LUT, 0, 7};       // this LUT is hardcoded at line 28 of functional_bootstrap.cpp and implements the identity
  int f_l4b[4] = {fbt_op_generate_mv_b_LUT, 0, 0, 5120};   // load the LUTs in polynomials
  int f_l4c[1] = {fbt_op_mv_f_bootstrap_init};            // Carpov multi value bootstrap first phase
  int f_l4d[4] = {fbt_op_mv_f_bootstrap, 0, 0, f_l4b[3]/context->torus_base};       // Carpov multi value bootstrap second phase
  int f_l4e[1] = {fbt_op_next_tree_level};                // next level of the tree evaluation
  
  int f_l3a[4] = {fbt_op_tlwe_keyswitch, 0, 0, f_l4d[3]/context->torus_base};       // tlwe to trlwe packing.
  int f_l3b[4] = {fbt_op_f_bootstrap, 0, 0, f_l3a[3]};          // 64 functional bootstraps
  int f_l3c[1] = {fbt_op_next_tree_level};                // next level of the tree evaluation

  int f_l2a[4] = {fbt_op_tlwe_keyswitch, 0, 0, f_l3a[3]/context->torus_base};        // tlwe to trlwe packing
  int f_l2b[4] = {fbt_op_f_bootstrap, 0, 0, f_l2a[3]};           // 16 functional bootstraps
  int f_l2c[1] = {fbt_op_next_tree_level};                // next level of the tree evaluation

  int f_l1a[4] = {fbt_op_tlwe_keyswitch, 0, 0, f_l2a[3]/context->torus_base};        // tlwe to trlwe packing
  int f_l1b[4] = {fbt_op_f_bootstrap, 0, 0, f_l1a[3]};           // 4 functional bootstraps
  int f_l1c[1] = {fbt_op_next_tree_level};                // next level of the tree evaluation

  int f_l0a[4] = {fbt_op_tlwe_keyswitch, 0, 0, f_l1a[3]/context->torus_base};        // tlwe to trlwe packing
  int f_l0b[4] = {fbt_op_f_bootstrap, 0, 0, f_l0a[3]};           // 4 functional bootstraps


  int * function[16] = {f_l4a, f_l4b, f_l4c, f_l4d, f_l4e, f_l3a, f_l3b, f_l3c, f_l2a, f_l2b, f_l2c, f_l1a, f_l1b, f_l1c, f_l0a, f_l0b};

  fbt_main_loop(function, 16, input->lwe_samples, context);            // Execute everything defined in "function"
    
  lweCopy(&output->lwe_samples[0], &context->output[0], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&output->lwe_samples[1], &context->output[1], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&output->lwe_samples[2], &context->output[2], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&output->lwe_samples[3], &context->output[3], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&output->lwe_samples[4], &context->output[4], &params->tgsw_params->tlwe_params->extracted_lweparams);

}

void test_LUT_12_bit_to_12_bit(fbt_integer * output, fbt_integer * input, fbt_context * context){
  const TFheGateBootstrappingParameterSet * params = context->params;
  int f_l5a[3] = {fbt_op_load_hardcoded_LUT, 0, 8};       // this LUT is hardcoded at line 28 of functional_bootstrap.cpp and implements the identity
  int f_l5b[4] = {fbt_op_generate_mv_b_LUT, 0, 0, 24576};   // load the LUTs in polynomials
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
  int f_l2c[1] = {fbt_op_next_tree_level};                // next level of the tree evaluation

  int f_l1a[4] = {fbt_op_tlwe_keyswitch, 0, 0, f_l2a[3]/context->torus_base};        // tlwe to trlwe packing
  int f_l1b[4] = {fbt_op_f_bootstrap, 0, 0, f_l1a[3]};           // 4 functional bootstraps
  int f_l1c[1] = {fbt_op_next_tree_level};                // next level of the tree evaluation

  int f_l0a[4] = {fbt_op_tlwe_keyswitch, 0, 0, f_l1a[3]/context->torus_base};        // tlwe to trlwe packing
  int f_l0b[4] = {fbt_op_f_bootstrap, 0, 0, f_l0a[3]};           // 4 functional bootstraps


  int * function[19] = {f_l5a, f_l5b, f_l5c, f_l5d, f_l5e, f_l4a, f_l4b, f_l4c,f_l3a, f_l3b, f_l3c, f_l2a, f_l2b, f_l2c, f_l1a, f_l1b, f_l1c, f_l0a, f_l0b};

  fbt_main_loop(function, 19, input->lwe_samples, context);            // Execute everything defined in "function"
    
  lweCopy(&output->lwe_samples[0], &context->output[0], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&output->lwe_samples[1], &context->output[1], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&output->lwe_samples[2], &context->output[2], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&output->lwe_samples[3], &context->output[3], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&output->lwe_samples[4], &context->output[4], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&output->lwe_samples[5], &context->output[5], &params->tgsw_params->tlwe_params->extracted_lweparams);

}


void test_32_bit_integer_comparison(LweSample * output, fbt_integer * input1, fbt_integer * input2, fbt_context * context){
  const TFheGateBootstrappingParameterSet * params = context->params;
  int f_l0a[3] = {fbt_op_load_hardcoded_LUT, 0, 2};
  int f_l0b[4] = {fbt_op_generate_b_LUT, 0, 0, 4};
  int f_l0c[4] = {fbt_op_f_bootstrap, 0, 0, 1};
  int f_l0d[4] = {fbt_op_load_literal, 1, 1, 3};
  int f_l0e[1] = {fbt_op_next_tree_level};

  int f_l1a[4] = {fbt_op_tlwe_keyswitch, 0, 0, 1};
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

  fbt_main_loop(function, 64, context->lwe_temp_in, context);

  lweCopy(output, &context->output[0], &params->tgsw_params->tlwe_params->extracted_lweparams);
}

// max
void max_8bit(fbt_integer * output, fbt_integer * input1, fbt_integer * input2, fbt_context * context){
  const TFheGateBootstrappingParameterSet * params = context->params;
  int f_l0a[3] = {fbt_op_load_hardcoded_LUT, 0, 2};
  int f_l0b[4] = {fbt_op_generate_b_LUT, 0, 0, 4};
  int f_l0c[4] = {fbt_op_f_bootstrap, 0, 0, 1};
  int f_l0d[4] = {fbt_op_load_literal, 1, 1, 3};
  int f_l0e[1] = {fbt_op_next_tree_level};

  int f_l1a[4] = {fbt_op_tlwe_keyswitch, 0, 0, 1};
  int f_l1b[4] = {fbt_op_f_bootstrap, 0, 0, 1};
  int f_l1c[4] = {fbt_op_load_literal, 1, 1, 3};
  int f_l1d[1] = {fbt_op_next_tree_level};


  int * function[16] = {f_l0a, f_l0b, f_l0c, f_l0d, f_l0e, f_l1a, f_l1b, f_l1c, f_l1d, f_l1a, f_l1b, f_l1c, f_l1d,
                        f_l1a, f_l1b, f_l1d};

  int * bt[3] = {f_l1a, f_l1b, f_l1d};

  int f_l2a[4] = {fbt_op_tlwe_keyswitch, 0, 0, 4};
  int f_l2b[4] = {fbt_op_f_bootstrap, 0, 0, 4};

  int * select[2] = {f_l2a, f_l2b};

  for (size_t i = 0; i < 4; i++){
    lweCopy(&context->input[i], &input1->lwe_samples[i], &params->tgsw_params->tlwe_params->extracted_lweparams);
    lweSubTo(&context->input[i], &input2->lwe_samples[i], &params->tgsw_params->tlwe_params->extracted_lweparams);
    lweKeySwitch(&context->lwe_temp_in[i], context->bootstrap_key->ks, &context->input[i]);
  }
  

  fbt_main_loop(function, 16, context->lwe_temp_in, context);

  lweKeySwitch(context->lwe_temp_in, context->bootstrap_key->ks, &input1->lwe_samples[3]);
  lweCopy(&context->input[1], &context->input[0], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweNegate(&context->input[2], &context->input[0], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&context->input[3], &context->input[2], &params->tgsw_params->tlwe_params->extracted_lweparams);
  fbt_main_loop(bt, 3, context->lwe_temp_in, context);

  lweKeySwitch(context->lwe_temp_in, context->bootstrap_key->ks, &input2->lwe_samples[3]);
  lweCopy(&context->input[1], &context->input[0], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweNegate(&context->input[2], &context->input[0], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&context->input[3], &context->input[2], &params->tgsw_params->tlwe_params->extracted_lweparams);
  fbt_main_loop(bt, 2, context->lwe_temp_in, context);

  context->output[0].b += modSwitchToTorus32(1, 8);
  lweKeySwitch(context->lwe_temp_in, context->bootstrap_key->ks, &context->output[0]);

  for (size_t j = 0; j < 4; j++){
    lweCopy(&context->input[0 + 4*j], &input2->lwe_samples[j], &params->tgsw_params->tlwe_params->extracted_lweparams);
    lweCopy(&context->input[1 + 4*j], &input1->lwe_samples[j], &params->tgsw_params->tlwe_params->extracted_lweparams);
    lweCopy(&context->input[2 + 4*j], &input1->lwe_samples[j], &params->tgsw_params->tlwe_params->extracted_lweparams);
    lweNoiselessTrivial(&context->input[3 + 4*j], 0, &params->tgsw_params->tlwe_params->extracted_lweparams);
  }
  
  fbt_main_loop(select, 2, context->lwe_temp_in, context);
  
  lweCopy(&output->lwe_samples[0], &context->output[0], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&output->lwe_samples[1], &context->output[1], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&output->lwe_samples[2], &context->output[2], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&output->lwe_samples[3], &context->output[3], &params->tgsw_params->tlwe_params->extracted_lweparams);
}


// relu
void relu_8bit(fbt_integer * output, fbt_integer * input, fbt_context * context){
  const TFheGateBootstrappingParameterSet * params = context->params;
  int f_l3a[3] = {fbt_op_load_hardcoded_LUT, 0, 5};
  int f_l3b[4] = {fbt_op_generate_b_LUT, 3, 0, 4};
  int f_l3c[4] = {fbt_op_tlwe_keyswitch, 0, 0, 3};
  int f_l3d[4] = {fbt_op_f_bootstrap, 0, 0, 4};

  int * function[4] = {f_l3a, f_l3b, f_l3c, f_l3d};

  lweKeySwitch(context->lwe_temp_in, context->bootstrap_key->ks, &input->lwe_samples[3]);

  lweNoiselessTrivial(&context->input[2], 0, &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweNoiselessTrivial(&context->input[3], 0, &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweNoiselessTrivial(&context->input[6], 0, &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweNoiselessTrivial(&context->input[7], 0, &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweNoiselessTrivial(&context->input[10], 0, &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweNoiselessTrivial(&context->input[11], 0, &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&context->input[0], &input->lwe_samples[0], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&context->input[1], &input->lwe_samples[0], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&context->input[4], &input->lwe_samples[1], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&context->input[5], &input->lwe_samples[1], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&context->input[8], &input->lwe_samples[2], &params->tgsw_params->tlwe_params->extracted_lweparams);
  lweCopy(&context->input[9], &input->lwe_samples[2], &params->tgsw_params->tlwe_params->extracted_lweparams);

  fbt_main_loop(function, 4, context->lwe_temp_in, context);

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


  TFheFunctionalBootstrappingSecretKeySet* key = new_random_functional_bootstrapping_secret_keyset(params, 2);
  fbt_context * context = fbt_init(key->cloud_key->params, key->cloud_key->bkFFT, key->tlweKS, key->log_torus_base, 6144);
  int64_t size = sizeof(Torus32) * key->tlweKS->N * key->tlweKS->k * key->tlweKS->N * key->tlweKS->t * ((1 << (key->tlweKS->basebit)) - 1) * context->torus_base; 
  printf("switch key size is %ld.\n", size);

 /*Relu*/
  fbt_integer * number1_8bit, * output_8bit;
  int8_t number1_8bit_clear = rand()&0xFF;
  number1_8bit = fbt_new_encrypted_integer(number1_8bit_clear, 4, key->extracted_key, context);
  output_8bit = fbt_new_encrypted_integer(0, 4, key->extracted_key, context);

  printf("\nRunning 8-bit ReLU %d times\n", NUM_EXE);
  MEASURE_TIME_AND_PRINT(
  for (size_t i = 0; i < NUM_EXE; i++)
    {
      relu_8bit(output_8bit, number1_8bit, context);
    }
  );
  printf("Result ReLU(%d): %d\n", number1_8bit_clear, fbt_integer_decrypt(output_8bit, key->extracted_key, context));


/* 32-bit integer comparison */
  fbt_integer * number1_32bit, * number2_32bit; // 32-bit (16-digit) numbers
  uint32_t number1_32bit_clear = rand();
  uint32_t number2_32bit_clear = rand();
  number1_32bit = fbt_new_encrypted_integer(number1_32bit_clear, 16, (LweKey *) key->tfhe_keys->lwe_key, context);
  number2_32bit = fbt_new_encrypted_integer(number2_32bit_clear, 16, (LweKey *) key->tfhe_keys->lwe_key, context);

  LweSample * output = new_LweSample(&params->tgsw_params->tlwe_params->extracted_lweparams);
  printf("\nRunning 32-bit integer comparison %d times\n", NUM_EXE);
  MEASURE_TIME_AND_PRINT(
  for (size_t i = 0; i < NUM_EXE; i++)
    {
      test_32_bit_integer_comparison(output, number1_32bit, number2_32bit, context);
    }
  );
  int intComp_result = modSwitchFromTorus32(lwePhase(output, key->extracted_key), 2*context->torus_base);
  printf("Result Compare(%u, %u): ", number1_32bit_clear, number2_32bit_clear);
  if(intComp_result == 1) printf("%u > %u\n", number1_32bit_clear, number2_32bit_clear);
  if(intComp_result == 7) printf("%u < %u\n", number1_32bit_clear, number2_32bit_clear);
  if(intComp_result == 0) printf("%u == %u\n", number1_32bit_clear, number2_32bit_clear);
  

 /* 6-bit LUT */
  int8_t lut_input_clear = rand()&0x3F;
  fbt_integer * lut_input = fbt_new_encrypted_integer(lut_input_clear, 4, (LweKey *) key->tfhe_keys->lwe_key, context);
  fbt_integer * lut_output = fbt_new_encrypted_integer(lut_input_clear, 4, key->extracted_key, context);
  printf("\nRunning 6-bit-to-6-bit LUT evaluation %d times. The LUT encodes the identity function.\n", NUM_EXE);
  MEASURE_TIME_AND_PRINT(
  for (size_t i = 0; i < NUM_EXE; i++)
    {
      test_LUT_6_bit_to_6_bit(lut_output, lut_input, context);
      
    }
  );
  printf("Result LUT(%d): %d\n", lut_input_clear, fbt_integer_decrypt(lut_output, key->extracted_key, context));


 /*8-bit LUT*/
  int8_t lut_input_clear_8 = rand()&0xFF;
  fbt_integer * lut_input_8 = fbt_new_encrypted_integer(lut_input_clear_8, 4, (LweKey *) key->tfhe_keys->lwe_key, context);
  fbt_integer * lut_output_8 = fbt_new_encrypted_integer(lut_input_clear_8, 4, key->extracted_key, context);
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
  fbt_integer * lut_input_10 = fbt_new_encrypted_integer(lut_input_clear_10, 6, (LweKey *) key->tfhe_keys->lwe_key, context);
  fbt_integer * lut_output_10 = fbt_new_encrypted_integer(lut_input_clear_10, 6, key->extracted_key, context);
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
  fbt_integer * lut_input_12 = fbt_new_encrypted_integer(lut_input_clear_12, 7, (LweKey *) key->tfhe_keys->lwe_key, context);
  fbt_integer * lut_output_12 = fbt_new_encrypted_integer(lut_input_clear_12, 7, key->extracted_key, context);
  printf("\nRunning 12-bit-to-12-bit LUT evaluation %d times. The LUT encodes the identity function.\n", NUM_EXE);
  MEASURE_TIME_AND_PRINT(
  for (size_t i = 0; i < NUM_EXE; i++)
    {
      test_LUT_12_bit_to_12_bit(lut_output_12, lut_input_12, context);
    }
  );
  printf("Result LUT(%d): %d\n", lut_input_clear_12, fbt_integer_decrypt(lut_output_12, key->extracted_key, context));


  
  /* Neural network functions */
  /* fbt_integer * number1_8bit, * number2_8bit, * output_8bit; // 8-bit (4-digit) numbers
  int8_t number1_8bit_clear = rand()&0xFF, number2_8bit_clear = rand()&0xFF;
  number1_8bit = fbt_new_encrypted_integer(number1_8bit_clear, 4, key->extracted_key, context);
  number2_8bit = fbt_new_encrypted_integer(number2_8bit_clear, 4, key->extracted_key, context);
  output_8bit = fbt_new_encrypted_integer(0, 4, key->extracted_key, context);

  
  printf("\nRunning 8-bit addition %d times\n", NUM_EXE);
  MEASURE_TIME_AND_PRINT(
    for (size_t i = 0; i < 10; i++)
    {
      fbt_integer_add(output_8bit, number1_8bit, number2_8bit, context);
    }
  );
  printf("Result (%d + %d): %d\n", number1_8bit_clear, number2_8bit_clear, fbt_integer_decrypt(output_8bit, key->extracted_key, context));

  printf("\nRunning 8-bit max %d times\n", NUM_EXE);
  MEASURE_TIME_AND_PRINT(
  for (size_t i = 0; i < NUM_EXE; i++)
    {
      max_8bit(output_8bit, number1_8bit, number2_8bit, context);
    }
  );
  printf("Result max(%d, %d): %d\n", number1_8bit_clear, number2_8bit_clear,  fbt_integer_decrypt(output_8bit, key->extracted_key, context));

  printf("\nRunning 8-bit ReLU %d times\n", NUM_EXE);
  MEASURE_TIME_AND_PRINT(
  for (size_t i = 0; i < NUM_EXE; i++)
    {
      relu_8bit(output_8bit, number1_8bit, context);
    }
  );
  printf("Result ReLU(%d): %d\n", number1_8bit_clear, fbt_integer_decrypt(output_8bit, key->extracted_key, context)); */

  return 0;
}
