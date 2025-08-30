#include <stdio.h>
#include <stdlib.h>
#include <time.h>

typedef struct neuron_struct {
    int num_weights;
    float *weights;
    float bias;
    float z;
    float a;

    int *mul_indices; // the indices that dictate the order of multiplications
} neuron;

typedef struct layer_struct {
    int num_neurons;
    neuron *neurons;
} layer;

typedef struct network_struct {
    int num_layers;
    layer *layers;
} network;

//Utility functions
void print_network(network net);
void free_network(network *net);

//Network control
neuron create_neuron(void* weights, int num_out_weights, int layer_idx, int neuron_idx); //legacy - neuron create_neuron(int num_out_weights);
layer create_layer(int num_neurons);
network create_network(int num_layers);
network init_network(int num_layers, int *num_neurons, void* weights);  //legacy - network construct_network(int num_outputs, int num_layers, int *num_neurons);

network shuffle_mul_indices(network net, int layer_idx);
network shuffle_mul_indices_masked(network net, int layer_idx);
network shuffle_mul_indices_deranged(network net, int layer_idx);

network forward(network net);   //legacy - void forward(network net);
network forward_shuffled(network net);

//Random Shuffling
void swap(int *a, int *b);
void fisher_yates(int arr[], int size);
void fisher_yates_masked(int arr[], unsigned int size, unsigned int s1[], unsigned int s2[], unsigned int length);
unsigned int modulo_masked(int i, unsigned int s1[], unsigned int s2[], unsigned int length);



// LEGACY FUNCTIONS
//void forward_shuffled(network net);
// No Overhead
//void forward_shuffled_NO(network net, int**** random_indices);
// No Overhead, Activations At the End
//void forward_shuffled_NO_AAE(network net, int**** random_indices);
// No Overhead, Activations At the End, Random Dummy Operations
//void forward_shuffled_NO_AAE_RDO(network net, int**** random_indices, int ***random_dummy_operations);
//int ****generate_random_indices(network net);
//int* get_random_indices(int size);







