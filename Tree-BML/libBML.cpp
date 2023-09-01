#include<iostream>
#include<fstream>
using namespace std;
#include "../FBT-TFHE-main/functional_bootstrap.cpp"
/*
    new defined begin
*/
typedef struct TLweKeySwitchKey_BML
{
  TLweSample *  ks0_raw;
  TLweSample ** ks1_raw;
  TLweSample *** ks;
  int32_t N, k, basebit, t;
  const TLweParams* out_params;
} TLweKeySwitchKey_BML;
/*new defines end*/

int decryptInteger_output_error_value(LweSample * sample, int digits, const int log_torus_base, LweKey * key){
  int result = 0;
  const int torus_base = 1 << log_torus_base;

  ofstream outfile; // output error
  outfile.open("output-erro-value.txt", ios_base::app|ios_base::out); // output error
  for (int i = digits - 1; i >= 0; i--){
    result <<= log_torus_base;
    const int tmp = modSwitchFromTorus32(lweSymDecrypt(&sample[i], key, torus_base*2), torus_base*2);

    const int32_t n = key->params->n;
    Torus32 axs = 0;
    const Torus32 *__restrict a = (&sample[i])->a;
    const int32_t * __restrict k = key->key;
    for (int32_t i = 0; i < n; ++i) 
	   axs += a[i]*k[i]; 
    Torus32 erro =  (&sample[i])->b - axs - lweSymDecrypt(&sample[i], key, torus_base*2); // output error
    outfile << erro << endl;    // output error


    result += tmp;
  }
  
	outfile.close(); // output error
  return (result << (32 - log_torus_base*digits)) >> (32 - log_torus_base*digits);
}

int32_t redefined_modSwitchFromTorus32(Torus32 phase, int32_t Msize, int thet){
    uint64_t interv = ((UINT64_C(1)<<63)/Msize)*2; // width of each intervall
    uint64_t half_interval = interv/2; // begin of the first intervall
    uint64_t phase64 = (uint64_t(phase)<<32) + half_interval;

    return ((phase64/interv)>>thet)<<thet;
}

/*
  compress（2^thet）*B 
*/
int generate_polynomial_LUT_BML(IntPolynomial ** result, int * LUT, int size, int torus_base, int N,int thet){
  *result = new_IntPolynomial_array(size/((1<<thet)*torus_base), N);
  int slice_size = N/torus_base;
  int redundance_times = slice_size/(1<<thet);
  for (int i = 0; i < size/((1<<thet)*torus_base); i++)
  {
    
    for (size_t j = 0; j < N/(1<<thet); j++){
      for(size_t k = 0; k < (1<<thet); k++){
        (*result)[i].coefs[j * (1<<thet) + k] = LUT[i*((1<<thet)*torus_base) + (j/redundance_times) + k*torus_base];
      }
    }

  }
  
  return size/((1<<thet)*torus_base);
}


void tfhe_functional_blindRotateAnd_Many_Extract_FFT(LweSample *result,
                                           const TLweSample *v,
                                           const TGswSampleFFT *bk,
                                           const int32_t barb,
                                           const int32_t *bara,
                                           const int32_t n,
                                           const TGswParams *bk_params,
                                           int thet){

    const TLweParams *accum_params = bk_params->tlwe_params;
    const LweParams *extract_params = &accum_params->extracted_lweparams;
    const int32_t N = accum_params->N;
    const int32_t _2N = 2 * N;
    // Test polynomial 

    // Accumulator
    TLweSample *acc = new_TLweSample(accum_params);
  
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
    tfhe_blindRotate_FFT(acc, bk, bara, n, bk_params);
    //Many LUT Extraction
    for (size_t i = 0; i < 1<<thet; i++)
    {
      tLweExtractLweSampleIndex(&result[i],acc,i,extract_params,accum_params);
    }
    

    delete_TLweSample(acc);
}


/*
  result is lweSample array
*/
void tfhe_functional_bootstrap_woKS_FFT_BML(LweSample *result,
                                    const LweBootstrappingKeyFFT *bk,
                                    TLweSample * testvect,
                                    const LweSample *x, const int torus_base,
                                    int thet){

  const TGswParams *bk_params = bk->bk_params;
  const TLweParams *accum_params = bk->accum_params;
  const LweParams *in_params = bk->in_out_params;
  const int32_t N = accum_params->N;
  const int32_t Nx2 = 2 * N;
  const int32_t n = in_params->n;
  int32_t Nx2_2thet = Nx2/(1<<thet);

  // precision offset
  Torus32 prec_offset = modSwitchToTorus32(N/(2*torus_base), Nx2);
  

  int32_t *bara = new int32_t[N];

  //redefined_modulus switching
  int32_t barb = modSwitchFromTorus32(x->b + prec_offset, Nx2_2thet);

  barb = barb * (1<<thet);
 

  for (size_t i = 0; i < n; i++)
  {
    bara[i] = modSwitchFromTorus32(x->a[i], Nx2_2thet)* (1<<thet);

  }
  
  tfhe_functional_blindRotateAnd_Many_Extract_FFT(result, testvect, bk->bkFFT, barb, bara, n, bk_params, thet);

  delete[] bara;
}





TLweKeySwitchKey_BML* new_TLweKeySwitchKey_BML(int32_t N, int32_t k, int32_t t, int32_t basebit, const TLweParams * out_params,const int torus_base){
  TLweKeySwitchKey_BML * tlweKS = (TLweKeySwitchKey_BML *) malloc(sizeof(TLweKeySwitchKey_BML));
  tlweKS->N = N;
  tlweKS->k = k;
  tlweKS->t = t;
  tlweKS->basebit = basebit;
  tlweKS->out_params = out_params;
  int32_t base = 1<<basebit;

  tlweKS->ks0_raw = new_TLweSample_array((base - 1)*N*t, out_params);
  tlweKS->ks1_raw = new TLweSample*[N*t];
  tlweKS->ks = new TLweSample**[N];


  for (int32_t i = 0; i < N*t ; i++)
  {
    tlweKS->ks1_raw[i] = tlweKS->ks0_raw + (base - 1)*i;
  }
  for (int32_t p = 0; p < N; p++)
  {
    tlweKS->ks[p] = tlweKS->ks1_raw + t*p;
  }

  return tlweKS;
}

TLweKeySwitchKey* new_TLweKeySwitchKey_BML_no_rotated(int32_t N, int32_t k, int32_t t, int32_t basebit, const TLweParams * out_params,const int torus_base, const int thet){
  TLweKeySwitchKey * tlweKS = (TLweKeySwitchKey *) malloc(sizeof(TLweKeySwitchKey));
  tlweKS->N = N;
  tlweKS->k = k;
  tlweKS->t = t;
  tlweKS->basebit = basebit;
  tlweKS->out_params = out_params;
  int32_t base = 1<<basebit;

  tlweKS->ks0_raw = new_TLweSample_array((base - 1)*N*t*torus_base*(1<<thet), out_params);
  tlweKS->ks1_raw = new TLweSample*[N*t*torus_base*(1<<thet)];
  tlweKS->ks2_raw = new TLweSample**[N*torus_base*(1<<thet)];
  tlweKS->ks = new TLweSample***[N];


  for (int32_t i = 0; i < N*t*torus_base*(1<<thet) ; i++)
  {
    tlweKS->ks1_raw[i] = tlweKS->ks0_raw + (base - 1)*i;
  }
  for (int32_t p = 0; p < N*torus_base*(1<<thet); ++p)
	  tlweKS->ks2_raw[p] = tlweKS->ks1_raw + t*p;
  for (int32_t p = 0; p < N; p++)
  {
    tlweKS->ks[p] = tlweKS->ks2_raw + torus_base*(1<<thet)*p;
  }

  return tlweKS;
}

void tLweCreateKeySwitchKey_BML(TLweKeySwitchKey_BML* result, const LweKey* in_key, const TLweKey* out_key, const int torus_base, const int thet){
  const int32_t n = result->N;
  const int32_t t = result->t;
  const int32_t basebit = result->basebit;
  const int32_t base = 1<<basebit;
  const double alpha = out_key->params->alpha_min;


  for (int32_t i = 0; i < n; i++){
    for (size_t j = 0; j < t; j++){
      for (int32_t h = 1; h < base; h++){
        Torus32 mess = (in_key->key[i]*h)*(1<<(32-(j+1)*basebit));
        tLweSymEncryptZero(&result->ks[i][j][h - 1], alpha, out_key);
        for (size_t q = 0; q < n/torus_base; q = q + (1<<thet)){
          result->ks[i][j][h - 1].b->coefsT[q] += mess;
        }
      }
    }
  }
}

void tLweCreateKeySwitchKey_BML_no_rotated(TLweKeySwitchKey* result, const LweKey* in_key, const TLweKey* out_key, const int torus_base, const int thet){
  const int32_t n = result->N;
  const int32_t t = result->t;
  const int32_t basebit = result->basebit;
  const int32_t base = 1<<basebit;
  const double alpha = out_key->params->alpha_min;


  for (int32_t i = 0; i < n; i++){
    for(size_t e = 0; e < torus_base; e++){
      for(size_t k = 0; k < (1<<thet); k++){
        for (size_t j = 0; j < t; j++){
          for (int32_t h = 1; h < base; h++){
            Torus32 mess = (in_key->key[i]*h)*(1<<(32-(j+1)*basebit));
            tLweSymEncryptZero(&result->ks[i][e*(1<<thet) + k][j][h - 1], alpha, out_key);
            for (size_t q = e*(n/torus_base) + k; q < (e + 1)*(n/torus_base); q = q + (1<<thet)){
              result->ks[i][e*(1<<thet) + k][j][h - 1].b->coefsT[q] += mess;
            }
          }
        }
      }
    }
  }
}

typedef struct TFheFunctionalBootstrappingSecretKeySet_BML{
  TFheGateBootstrappingSecretKeySet* tfhe_keys;
  const TFheGateBootstrappingCloudKeySet* cloud_key;
  LweKey *extracted_key;
  TLweKeySwitchKey_BML * tlweKS;
  int log_torus_base;
} TFheFunctionalBootstrappingSecretKeySet_BML;



TFheFunctionalBootstrappingSecretKeySet_BML* new_random_functional_bootstrapping_secret_keyset_BML(TFheGateBootstrappingParameterSet* params, const int log_torus_base, const int thet){
  TFheFunctionalBootstrappingSecretKeySet_BML * result = (TFheFunctionalBootstrappingSecretKeySet_BML *) malloc(sizeof(TFheFunctionalBootstrappingSecretKeySet_BML));
  uint32_t seed[8];
  generate_random_bytes(32, (uint8_t * ) seed);
  result->log_torus_base = log_torus_base;
  const int torus_base = 1 << log_torus_base;
  tfhe_random_generator_setSeed(seed, 8);
  result->tfhe_keys = new_random_gate_bootstrapping_secret_keyset(params);
  result->cloud_key = &result->tfhe_keys->cloud;
  const TLweKey *tlwe_key = &result->tfhe_keys->tgsw_key->tlwe_key;
  result->extracted_key = new_LweKey(&result->tfhe_keys->tgsw_key->tlwe_params->extracted_lweparams);
  tLweExtractKey(result->extracted_key, tlwe_key); 
  result->tlweKS = new_TLweKeySwitchKey_BML(result->tfhe_keys->cloud.bkFFT->accum_params->N, result->tfhe_keys->cloud.bkFFT->accum_params->k, params->ks_t, params->ks_basebit, result->tfhe_keys->cloud.bkFFT->bk_params->tlwe_params, torus_base);
  tLweCreateKeySwitchKey_BML(result->tlweKS, result->extracted_key, tlwe_key, torus_base, thet);
  return result;
  
}

TFheFunctionalBootstrappingSecretKeySet* new_random_functional_bootstrapping_secret_keyset_BML_no_rotated(TFheGateBootstrappingParameterSet* params, const int log_torus_base, const int thet){
  TFheFunctionalBootstrappingSecretKeySet * result = (TFheFunctionalBootstrappingSecretKeySet *) malloc(sizeof(TFheFunctionalBootstrappingSecretKeySet));
  uint32_t seed[8];
  generate_random_bytes(32, (uint8_t * ) seed);
  result->log_torus_base = log_torus_base;
  const int torus_base = 1 << log_torus_base;
  tfhe_random_generator_setSeed(seed, 8);
  result->tfhe_keys = new_random_gate_bootstrapping_secret_keyset(params);
  result->cloud_key = &result->tfhe_keys->cloud;
  const TLweKey *tlwe_key = &result->tfhe_keys->tgsw_key->tlwe_key;
  result->extracted_key = new_LweKey(&result->tfhe_keys->tgsw_key->tlwe_params->extracted_lweparams);
  tLweExtractKey(result->extracted_key, tlwe_key);  
  result->tlweKS = new_TLweKeySwitchKey_BML_no_rotated(result->tfhe_keys->cloud.bkFFT->accum_params->N, result->tfhe_keys->cloud.bkFFT->accum_params->k, params->ks_t, params->ks_basebit, result->tfhe_keys->cloud.bkFFT->bk_params->tlwe_params, torus_base, thet);
  tLweCreateKeySwitchKey_BML_no_rotated(result->tlweKS, result->extracted_key, tlwe_key, torus_base, thet);
  return result;
  
}

/*b = (b + a.x^pow) mod(x^N + 1)*/
EXPORT void combina_b(Torus32 * b, Torus32 * a, const int N, const int pow){


  for (size_t i = 0; i < pow; i++)
  {
    b[i] -= a[N - pow + i];
  }
  for (size_t i = pow; i < N; i++)
  {
    b[i] += a[i - pow];
  }
}






void tLweKeySwitchTranslate_fromArray_BML_FULL(TLweSample* result,
  TLweSample*** ks, const TLweParams* params, 
  const TorusPolynomial* ai,
  const int32_t n, const int32_t t, const int32_t basebit, const int torus_base, const int thet){
  
  const int32_t base=1<<basebit;
  const int32_t prec_offset=1<<(32-(1+basebit*t));
  const int32_t mask=base-1;
  const int k = params->k;
  

  TorusPolynomial* temp = new_TorusPolynomial_array((k + 1)*(1<<thet)*torus_base, n);
  for (size_t i = 0; i < (k + 1)*(1<<thet)*torus_base; i++){
    torusPolynomialClear(&temp[i]);
  }
  
  for (int32_t i = 0; i < n; i++){
    for (int e = 0; e < (1<<thet)*torus_base; e++){
      const uint32_t aibar_1 = ai[i].coefsT[e] + prec_offset;
      for(int32_t j=0;j<t;j++){
        const uint32_t aij_1=(aibar_1>>(32-(j+1)*basebit)) & mask;
        if(aij_1 != 0){
          for(int q=0; q<=k; q++){
            for(int p = 0; p < n;p++){
              temp[q + (k + 1)*(e)].coefsT[p] -= ks[i][j][aij_1 - 1].a[q].coefsT[p];
            }
          }
        }
      }
    }
  }
  
  for(int j = 0; j < torus_base; j++){
    for(int i = (1<<thet) - 1; i>0; i--){
      for(int q=0; q<=k; q++){
        combina_b(temp[q + (k + 1)*(j*(1<<thet) + i - 1)].coefsT, temp[q + (k + 1)*(j*(1<<thet) + i)].coefsT, n, 1);
      }
    }
  }

  for(int i = torus_base - 1; i>0; i--){
    for(int q=0; q<=k; q++){
      combina_b(temp[q + (k + 1)*(i-1)*(1<<thet)].coefsT, temp[q + (k + 1)*i*(1<<thet)].coefsT, n, n/torus_base);
    }
  }

  for(size_t q=0; q<=k; q++){
    torusPolynomialCopy(&(result->a[q]), &temp[q]);
  }

  
  delete_TorusPolynomial_array((k + 1)*(1<<thet)*torus_base, temp);
}
  


void public_functional_tlweKeySwitch_BML_FULLThetB(TLweSample* result, const TLweKeySwitchKey_BML* ks, LweSample* samples, const int torus_base, const int thet){
  const TLweParams* params=ks->out_params;
  const int32_t n=ks->N;
  const int32_t basebit=ks->basebit;
  const int32_t t=ks->t;
  
  TorusPolynomial * f_a = new_TorusPolynomial_array(n, (1<<thet)*torus_base);

  int a1 = n/(torus_base*(1 << thet));
  int a2 = n/torus_base;

  for(int i = 0; i < n; i++){
    for(int j = 0; j <= params->k; j++){
      result->a[j].coefsT[i] = 0;
    }
  }

  for(size_t i = 0; i < n; i++){
    for(size_t j = 0; j < torus_base; j++){
      if(thet == 0){
        f_a[i].coefsT[j] = samples[j].a[i];
        }
      else{
        for(size_t k = 0; k < (1<<thet); k++){
          f_a[i].coefsT[(1<<thet)*j + k] = samples[j + k*torus_base ].a[i];
        }
      }
      
    }
  } 
 
  tLweKeySwitchTranslate_fromArray_BML_FULL(result, ks->ks, params, f_a, n, t, basebit, torus_base, thet);
  //result->current_variance = samples->current_variance + params->extracted_lweparams.n *t*n*ks->ks0_raw->current_variance + params->extracted_lweparams.n/(12*(1<< ks->basebit)*(1<< ks->basebit));

  
  for (size_t i = 0; i < torus_base; i++){
    for (size_t j = 0; j < a1; j++){
      for (size_t q = 0; q < (1 << thet); q++){
        result->b->coefsT[i*a2 + j*(1<<thet) + q] += samples[i + q*torus_base].b;//排序b
      }
    }
  }

  delete_TorusPolynomial_array(n, f_a);
}


void tLweKeySwitchTranslate_fromArray_BML_FULL_no_rotated(TLweSample* result,
    TLweSample**** ks, const TLweParams* params, 
    const TorusPolynomial* ai,
    const int32_t n, const int32_t t, const int32_t basebit, const int torus_base, const int thet){
    
    const int32_t base=1<<basebit;
    const int32_t prec_offset=1<<(32-(1+basebit*t));
    const int32_t mask=base-1;
    const int k = params->k;
    
    for (int32_t i = 0; i < n; i++){
      for (int e = 0; e < (1<<thet)*torus_base; e++){
        const uint32_t aibar_1 = ai[i].coefsT[e] + prec_offset;
        for(int32_t j=0;j<t;j++){
          const uint32_t aij_1=(aibar_1>>(32-(j+1)*basebit)) & mask;
          if(aij_1 != 0){
            for(int q=0; q<=k; q++){
              for(int p = 0; p < n;p++){
                result->a[q].coefsT[p] -= ks[i][e][j][aij_1 - 1].a[q].coefsT[p];                
              }
            }
          }
        }
      }
    }
  }

void public_functional_tlweKeySwitch_BML_FULLThetB_no_rotated(TLweSample* result, const TLweKeySwitchKey* ks, LweSample* samples, const int torus_base, const int thet){
  const TLweParams* params=ks->out_params;
  const int32_t n=ks->N;
  const int32_t basebit=ks->basebit;
  const int32_t t=ks->t;
  
  TorusPolynomial * f_a = new_TorusPolynomial_array(n, (1<<thet)*torus_base);

  int a1 = n/(torus_base*(1 << thet));
  int a2 = n/torus_base;

  for(int i = 0; i < n; i++){
    for(int j = 0; j <= params->k; j++){
      result->a[j].coefsT[i] = 0;
    }
  }

  for(size_t i = 0; i < n; i++){
    for(size_t j = 0; j < torus_base; j++){
      if(thet == 0){
        f_a[i].coefsT[j] = samples[j].a[i];
        }
      else{
        for(size_t k = 0; k < (1<<thet); k++){
          f_a[i].coefsT[(1<<thet)*j + k] = samples[j + k*torus_base ].a[i];
        }
      }
      
    }
  }

  tLweKeySwitchTranslate_fromArray_BML_FULL_no_rotated(result, ks->ks, params, f_a, n, t, basebit, torus_base, thet);
  
  for (size_t i = 0; i < torus_base; i++){
    for (size_t j = 0; j < a1; j++){
      for (size_t q = 0; q < (1 << thet); q++){
        result->b->coefsT[i*a2 + j*(1<<thet) + q] += samples[i + q*torus_base].b;//排序b
      }
    }
  }

  delete_TorusPolynomial_array(n, f_a);
}


void public_functional_tlweKeySwitch_BML_divisible(TLweSample* result, const TLweKeySwitchKey_BML* ks, LweSample* samples, const int torus_base, const int thet, int number_of_samples){
  int pad_times = ((1<<thet)*torus_base)/number_of_samples;

  for(size_t i = 1; i <= pad_times; i++){
    for(size_t j = 0; j < number_of_samples; j++){
      lweCopy(samples + (i * number_of_samples + j),samples + j, &(ks->out_params->extracted_lweparams));
    }
  }
  public_functional_tlweKeySwitch_BML_FULLThetB(result,ks,samples,torus_base,thet);
} 

void public_functional_tlweKeySwitch_BML_padding_then_divis(TLweSample* result, const TLweKeySwitchKey_BML* ks, LweSample* samples, const int torus_base, const int thet, int number_of_samples){
  int pad_times = (1<<((int)(log2(number_of_samples/torus_base) + 1)))*torus_base - number_of_samples;
  for(size_t i = 0; i < pad_times; i++){
    lweClear(samples + number_of_samples + i, &(ks->out_params->extracted_lweparams));
  }
  if((number_of_samples + pad_times) == (1<<thet)*torus_base)
    public_functional_tlweKeySwitch_BML_FULLThetB(result,ks,samples,torus_base,thet);
  else public_functional_tlweKeySwitch_BML_divisible(result, ks, samples, torus_base, thet, number_of_samples + pad_times);
}

void public_functional_tlweKeySwitch_BML_divisible_no_rotated(TLweSample* result, const TLweKeySwitchKey* ks, LweSample* samples, const int torus_base, const int thet, int number_of_samples){
  int pad_times = ((1<<thet)*torus_base)/number_of_samples;
  for(size_t i = 1; i <= pad_times; i++){
    for(size_t j = 0; j < number_of_samples; j++){
      lweCopy(samples + (i * number_of_samples + j),samples + j, &(ks->out_params->extracted_lweparams));
    }
  }
  public_functional_tlweKeySwitch_BML_FULLThetB_no_rotated(result,ks,samples,torus_base,thet);
}

void public_functional_tlweKeySwitch_BML_padding_then_divis_no_rotated(TLweSample* result, const TLweKeySwitchKey* ks, LweSample* samples, const int torus_base, const int thet, int number_of_samples){
  int pad_times = (1<<((int)(log2(number_of_samples/torus_base) + 1)))*torus_base - number_of_samples;
  for(size_t i = 0; i < pad_times; i++){
    lweClear(samples + number_of_samples + i, &(ks->out_params->extracted_lweparams));
  }
  if((number_of_samples + pad_times) == (1<<thet)*torus_base)
    public_functional_tlweKeySwitch_BML_FULLThetB_no_rotated(result,ks,samples,torus_base,thet);
  else public_functional_tlweKeySwitch_BML_divisible_no_rotated(result, ks, samples, torus_base, thet, number_of_samples + pad_times);
}
