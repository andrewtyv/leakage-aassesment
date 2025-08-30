/*
 * SimpleSerial V2 Template C code
 * Can be freely used to implement ChipWhisperer target binaries.
 *
 * Date: 14th March 2021
 */

/*
 * When debugging locally compile using `gcc -o debug-app.exe main.c network.c debug-source.c -DDEBUGGING=1`
 */
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "network.h"
#include "network_config.h"

#include "hal/hal.h"
#include "hal/stm32f3/stm32f3_hal.h"

#define SS_VER SS_VER_2_1

#include "simpleserial/simpleserial.h"

#ifdef DEBUGGING  // If debugging import windows for QueryPerformanceCounter() to measure overhead time with high resolution
#include <windows.h>

/// A Debugging test handle 
uint8_t test_handle(uint8_t cmd, uint8_t scmd, uint8_t len, uint8_t *buf)
{
  int arr[7] = {0, 1, 2, 3, 4, 5, 6};
  unsigned int s1[7] = {1, 1, 1, 3, 2, 1, 4};
  unsigned int s2[7] = {1, 1, 1, 3, 3, 1, 2};
  unsigned int length = 0;
  unsigned int num_of_neurons = 7; 
  while (num_of_neurons > 0) {
    num_of_neurons >>= 1;
    length++;
  }

  // fisher_yates_masked(arr, 7, s1, s2, length);

  // printf("Shuffled array: ");
  // for (int i = 0; i < 7; i++) {
  //   printf("%d", arr[i]);
  // }
  for (int i=6; i > 1; i--){
    unsigned int j = modulo_masked(i, s1, s2, length);
    printf("%d", j);
  }
  printf("\n");
  return 0;
}
#endif 

/// This function will handle the 'p' command send from the capture board.
uint8_t handle(uint8_t cmd, uint8_t scmd, uint8_t len, uint8_t *buf)
{
  int num_layers = NET_NUM_LAYERS;
  int *num_neurons_arr = NET_NUM_NEURONS;
  
  // Initialize the network with the pre-defined structure
  network net = init_network(num_layers, num_neurons_arr, net_config_layer_weights);

  // USED WHEN USING FIXED RANDOMLY GENERATED INPUTS INSTEAD OF 0.5 or 0.0
  // float new_inputs[7] = {
  //   -0.039964473300472925,
  //   0.3367311187095563,
  //   1.034674935250851746,
  //   0.7525740869461268,
  //   -1.5751881331414292,
  //   -0.6688709437380642,
  //   -0.9791366791093283
  // };
  // for (int i=0; i < net.layers[0].num_neurons; i++){
  //   net.layers[0].neurons[i].a = new_inputs[i];
  // }

  //Change the input of the first neuron in the first layer to the provided number
  //convert to float from a 4-byte buffer
  float input_value;
  memcpy(&input_value, buf, sizeof(float));
  //net.layers[0].neurons[0].a = input_value;
  ///net.layers[0].neurons[1].a = input_value;
  int n0 = net.layers[0].num_neurons;
  for (int i = 0; i < n0; i++) {
      net.layers[0].neurons[i].a = (i == 1 ? input_value : 0.5f);
  }
  
  //uint8_t rnd = buf[0];
  //float in_val = (float)rnd;
  //for (int i = 0; i < net.layers[0].num_neurons; i++) {
  //      net.layers[0].neurons[i].a = (i == 0 ? in_val : 0.5f);
  //} 


  #ifdef DEBUGGING
  #if 1
  LARGE_INTEGER frequency, start, middle, end;
  double overhead_time, forward_pass_time, overall_time;
  // Get the frequency of the performance counter
  QueryPerformanceFrequency(&frequency);
  QueryPerformanceCounter(&start);
  #endif
  #endif
  //scmd = 0; 
  if (scmd) {
    for (int i = 1; i < net.num_layers; i++) {
     
     net = shuffle_mul_indices_masked(net, i);
     //net = shuffle_mul_indices_deranged(net, i);
     //net = shuffle_mul_indices(net, i);
    }
    //net = shuffle_mul_indices_deranged(net, 1);
    //net = shuffle_mul_indices_masked(net, 1);
    //net = shuffle_mul_indices(net, 1);
  }

  #ifdef DEBUGGING
  #if 1
  //print_network(net);
  QueryPerformanceCounter(&middle);
  #endif
  #endif

  // Start Measurement
  trigger_high(); 
  if (scmd){ 
    net = forward_shuffled(net);
  }
  else {
    net = forward(net);
  }
  // Stop Measurement
  trigger_low();


  #ifdef DEBUGGING
  #if 1
  QueryPerformanceCounter(&end);
  overhead_time = (double)(middle.QuadPart - start.QuadPart) / frequency.QuadPart;
  overall_time = (double)(end.QuadPart - start.QuadPart) / frequency.QuadPart;
  forward_pass_time = (double)(end.QuadPart - middle.QuadPart) / frequency.QuadPart;
  double percentage = overhead_time / overall_time * 100;
  print_network(net);
  printf("Overall Time: %.16f\nForward Pass Time: %.16f\nOverhead Time: %.16f\nOverhead/Total percentage: %.16f\%\n", overall_time, forward_pass_time, overhead_time , percentage);
  
  #endif
  #endif
  
  //free dynamically allocated memory
  free_network(&net);
  
  simpleserial_put('r', len, buf);

  return 0;
}

int main(void) {
  srand(time(NULL));
  //Initialize network weights
  init_weights();
  // Setup the specific chipset.
  platform_init();
  // Setup serial communication line.
  init_uart();
  // Setup measurement trigger.
  trigger_setup();

  simpleserial_init();

  // Insert your handlers here.
  simpleserial_addcmd('p', 1*sizeof(float), handle);

#ifdef DEBUGGING
  simpleserial_addcmd('t', 16, test_handle);
#endif
  // What for the capture board to send commands and handle them.
  while (1)
    simpleserial_get();
}
