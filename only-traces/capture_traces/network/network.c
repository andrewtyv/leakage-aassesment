#include "network.h"
#include <stdint.h>
#include <math.h>

unsigned int blakely(unsigned int a, unsigned int b, unsigned int n, unsigned int length) {
    unsigned int R = 0;
    for (int i = length - 1; i >= 0; i--) {
        unsigned int ai = (a >> i) & 1;
        R = 2 * R + ai * b;
        if (R >= n) {
            R = R - n;
        }
        if (R >= n) {
            R = R - n;
        }
    }
    return R;
}

unsigned int modulo_masked(int i, unsigned int s1[], unsigned int s2[], unsigned int length) {
    unsigned int r1 = rand();
    unsigned int r2 = rand();
    unsigned int tmp = (r1 * (s1[i] % (i + 1)) + r2 * (i + 1)) % (i + 1);
    // Compute (tmp * s2[i]) % (i + 1)
    unsigned int j = blakely(tmp, (s2[i] % (i + 1)), i + 1, length);
    return j;
}

void fisher_yates_masked(int arr[], unsigned int size, unsigned int s1[], unsigned int s2[], unsigned int length) {
    #ifdef DEBUGGING
    #if 0
    printf("arr: ");
    for (int i=0; i<size; i++)printf("%d", arr[i]);
    printf(" size: %d, s1: ", size);
    for (int i=0; i<size; i++)printf("%d", s1[i]);
    printf(" s2: ");
    for (int i=0; i<size; i++)printf("%d", s2[i]);
    printf(" length: %d Shuffled arr: ", length);
    #endif
    #endif
    for (int i = size - 1; i >= 2; i--) {
        unsigned int j = modulo_masked(i, s1, s2, length);
        unsigned int temp = arr[i];
        arr[i] = arr[j];
        arr[j] = temp;
    }
    unsigned int j = rand() % 2;
    unsigned int temp = arr[1];
    arr[1] = arr[j];
    arr[j] = temp;

    #ifdef DEBUGGING
    #if 0
    for (int i=0; i<size; i++)printf("%d", arr[i]);
    printf("\n");
    #endif
    #endif
}

void swap(int *a, int *b){
    int temp = *a;
    *a = *b;
    *b = temp;
}
/** 
* @brief
* Shuffles the array using the Fisher-Yates shuffle
*/
void fisher_yates(int arr[], int size){
    for (int i = size - 1; i > 0; i--) {
        int j = rand() % (i + 1);
        swap(&arr[i], &arr[j]);
    }
}

void fisher_yates_deranged(int arr[], int size){
    for (int i = size - 1; i > 0; i--) {
        int j = rand() % (i + 1);
        swap(&arr[i], &arr[j]);
    }

    if (size > 2){
        for (int i = 0; i < size; i++){
            if (arr[i] == i) {

                int swap_index = -1;
                do {
                    swap_index = rand() % size;
                } while (swap_index == i);

                swap(&arr[i], &arr[swap_index]);
            }
        }
    }
}



void free_network(network *net){
    // free the dynamically allocated fields inside the network struct
    for (int i=0; i < net->num_layers; i++){
        for(int j=0; j< net->layers[i].num_neurons; j++){
            if (net->layers[i].neurons[j].weights != NULL) free(net->layers[i].neurons[j].weights);
            if (net->layers[i].neurons[j].mul_indices != NULL) free(net->layers[i].neurons[j].mul_indices);
        }
        free(net->layers[i].neurons);
    }
    free(net->layers);
}

/*
* Prints out the values contained withing the network - number of layers, and then for each layer prints out each neuron
* for each neuron this prints out the a value, the z value, and the weight values (an array of size equal to the number of neurons in the previous layer - except for the first layer - the input layer ) 
*/
void print_network(network net){
    printf("\n");
    printf("-----------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
    printf("Network - num_layers = %d\n", net.num_layers);
    for (int i = 0; i < net.num_layers; i++){
        printf("Layer %d:\n", i);
        for (int j = 0; j < net.layers[i].num_neurons; j++){
            printf("\tNeuron %d | a=%f z=%f\t| ", j, net.layers[i].neurons[j].a,  net.layers[i].neurons[j].z );
            if (i >= 1){
                printf("Mul Indices: ");
                for (int k = 0; k < net.layers[i - 1].num_neurons; k++){
                    printf("%d", net.layers[i].neurons[j].mul_indices[k]);
                }
                printf("\tWeights: {");
                for (int k = 0; k < net.layers[i - 1].num_neurons; k++){
                    printf("%f", net.layers[i].neurons[j].weights[k]);
                    if (k < net.layers[i - 1].num_neurons - 1)
                    printf(", ");
                    else
                    printf("}");
                }
            }
            printf("\n");
        }
    }
    printf("-----------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
}


neuron create_neuron(void* weights, int num_in_weights, int layer_idx, int neuron_idx){
    neuron new_neuron;
    new_neuron.a = 0.5;
    new_neuron.z = 0.0;
    new_neuron.bias = 0.0;
    new_neuron.num_weights = num_in_weights;

    if (num_in_weights > 0) {
        new_neuron.weights = (float*) malloc(num_in_weights * sizeof(float));
        new_neuron.mul_indices = (int*) malloc(num_in_weights * sizeof(int));
    } else {
        new_neuron.weights = NULL;
        new_neuron.mul_indices = NULL;
    }
    if (weights != NULL && num_in_weights > 0){
        //TODO: Dont question it... it works.
        float (*layer_weights)[num_in_weights] = ((float (*)[num_in_weights])((float**)weights)[layer_idx]);
        for (int i=0; i<num_in_weights; i++){
            new_neuron.weights[i] = layer_weights[neuron_idx][i];
            new_neuron.mul_indices[i] = i;
        }
    }
    return new_neuron;
}

layer create_layer(int num_neurons){
    layer lay;
    lay.num_neurons = num_neurons;
    lay.neurons = (neuron*) malloc(num_neurons * sizeof(neuron));
    return lay;
}

network create_network(int num_layers){
    network net;
    net.num_layers = num_layers;
    net.layers = (layer*) malloc(num_layers * sizeof(layer));
    return net;
}

network init_network(int num_layers, int *num_neurons, void* weights) {
    network net = create_network(num_layers);
    int curr_layer_idx, curr_neuron_idx;
    for (curr_layer_idx = 0; curr_layer_idx < num_layers; curr_layer_idx++){
        net.layers[ curr_layer_idx ] = create_layer(num_neurons[ curr_layer_idx ]);
    }
    // create neurons for the first (input) layer - they dont have weights
    for (curr_neuron_idx = 0; curr_neuron_idx < net.layers[0].num_neurons; curr_neuron_idx++){
        net.layers[0].neurons[ curr_neuron_idx ] = create_neuron(NULL, 0, 0, curr_neuron_idx);
    }
    // For each following layer create neurons with number of weights eqaual to the number of neurons in the previous layer
    for (curr_layer_idx = 1; curr_layer_idx < num_layers; curr_layer_idx++){
        int prev_layer_idx = curr_layer_idx - 1;
        for (curr_neuron_idx = 0; curr_neuron_idx <net.layers[ curr_layer_idx ].num_neurons; curr_neuron_idx++){
            net.layers[ curr_layer_idx ].neurons[ curr_neuron_idx ] = create_neuron(weights, net.layers[ prev_layer_idx ].num_neurons, curr_layer_idx, curr_neuron_idx );
        }
    }
    return net;
}


network shuffle_mul_indices_masked(network net, int layer_idx) {
   unsigned int s1[] = {1,1,1,3,2,1,4}; //2,3,1,4,1,2,3,1};
    unsigned int  s2[] = {1,1,1,3,3,1,2}; //3,1,2,1,4,2,1,3};

    unsigned int length = 0;
    unsigned int num_of_neurons = net.layers[ layer_idx ].num_neurons; 
    while (num_of_neurons > 0) {
        num_of_neurons >>= 1;
        length++;
    }
    if (layer_idx > 0 && layer_idx < net.num_layers) {
        for (int i = 0; i < net.layers[ layer_idx ].num_neurons; i++){
            fisher_yates_masked(net.layers[ layer_idx ].neurons[ i ].mul_indices, net.layers[ layer_idx ].neurons[ i ].num_weights, s1, s2, length);
        }
    }
    return net;
}

network shuffle_mul_indices(network net, int layer_idx) {
    if (layer_idx > 0 && layer_idx < net.num_layers) {
        for (int i = 0; i < net.layers[ layer_idx ].num_neurons; i++){
            fisher_yates(net.layers[ layer_idx ].neurons[ i ].mul_indices, net.layers[ layer_idx ].neurons[ i ].num_weights);
        }
    }
    return net;
}

network shuffle_mul_indices_deranged(network net, int layer_idx) {
    if (layer_idx > 0 && layer_idx < net.num_layers) {
        for (int i = 0; i < net.layers[ layer_idx ].num_neurons; i++){
            fisher_yates_deranged(net.layers[ layer_idx ].neurons[ i ].mul_indices, net.layers[ layer_idx ].neurons[ i ].num_weights);
        }
    }
    return net;
}

network forward(network net){
    volatile int curr_layer_idx, curr_neuron_idx, prev_layer_neuron_idx;
    // for each layer
    for (curr_layer_idx=1; curr_layer_idx < net.num_layers; curr_layer_idx++){
        
        int prev_layer_idx = curr_layer_idx - 1;
        // for each neuron in this layer
        for (curr_neuron_idx=0; curr_neuron_idx < net.layers[ curr_layer_idx ].num_neurons; curr_neuron_idx++){   
            net.layers[ curr_layer_idx ].neurons[ curr_neuron_idx ].z = net.layers[ curr_layer_idx ].neurons[ curr_neuron_idx ].bias;

            // for all neurons on the previous layer
            for (prev_layer_neuron_idx = 0; prev_layer_neuron_idx <net.layers[ prev_layer_idx ].num_neurons; prev_layer_neuron_idx++){
                net.layers[ curr_layer_idx ].neurons[ curr_neuron_idx ].z =
                    net.layers[ curr_layer_idx ].neurons[ curr_neuron_idx ].z
                    +
                    (
                        (net.layers[ curr_layer_idx ].neurons[ curr_neuron_idx ].weights[ prev_layer_neuron_idx ])
                        *
                        (net.layers[ prev_layer_idx ].neurons[ prev_layer_neuron_idx ].a)
                    );
                // We are looking for THIS MULTIPLICATION
            }
            //get a values
            net.layers[curr_layer_idx].neurons[ curr_neuron_idx ].a = net.layers[curr_layer_idx].neurons[ curr_neuron_idx ].z;
            //apply relu
            if(curr_layer_idx < net.num_layers-1){
                if((net.layers[curr_layer_idx].neurons[ curr_neuron_idx ].z) < 0)
                {
                    net.layers[curr_layer_idx].neurons[ curr_neuron_idx ].a = 0;
                }
                else
                {
                    net.layers[curr_layer_idx].neurons[ curr_neuron_idx ].a = net.layers[curr_layer_idx].neurons[ curr_neuron_idx ].z;
                }
            }
            //apply sigmoid to the last layer
            else{
                net.layers[curr_layer_idx].neurons[ curr_neuron_idx ].a = 1/(1+exp(-net.layers[curr_layer_idx].neurons[ curr_neuron_idx ].z));
            }
        }
    }
    return net;
}

network forward_shuffled(network net) {
    volatile int curr_layer_idx, curr_neuron_idx, prev_layer_neuron_idx;
    // for each layer
    for (curr_layer_idx=1; curr_layer_idx < net.num_layers; curr_layer_idx++){
        
        int prev_layer_idx = curr_layer_idx - 1;
        // for each neuron in this layer
        for (curr_neuron_idx=0; curr_neuron_idx < net.layers[ curr_layer_idx ].num_neurons; curr_neuron_idx++){   
            net.layers[ curr_layer_idx ].neurons[ curr_neuron_idx ].z = net.layers[ curr_layer_idx ].neurons[ curr_neuron_idx ].bias;

            // for all neurons on the previous layer
            for (prev_layer_neuron_idx = 0; prev_layer_neuron_idx < net.layers[ prev_layer_idx ].num_neurons; prev_layer_neuron_idx++){
                int mul_index = net.layers[ curr_layer_idx ].neurons[ curr_neuron_idx ].mul_indices[ prev_layer_neuron_idx ]; // CHANGE from forward - added this line

                net.layers[ curr_layer_idx ].neurons[ curr_neuron_idx ].z =
                    net.layers[ curr_layer_idx ].neurons[ curr_neuron_idx ].z
                    +
                    (
                        (net.layers[ curr_layer_idx ].neurons[ curr_neuron_idx ].weights[ mul_index ]) // CHANGE from forward - .weights[ prev_layer_neuron_idx ] -> .weights[ mul_index ]
                        *
                        (net.layers[ prev_layer_idx ].neurons[ mul_index ].a) // CHANGE from forward - .neurons[ prev_layer_neuron_idx ].a -> .neurons[ mul_index ].a
                    );
                // We are looking for THIS MULTIPLICATION
            }
            //get a values
            net.layers[curr_layer_idx].neurons[ curr_neuron_idx ].a = net.layers[curr_layer_idx].neurons[ curr_neuron_idx ].z;
            //apply relu
            if(curr_layer_idx < net.num_layers - 1){
                if((net.layers[curr_layer_idx].neurons[ curr_neuron_idx ].z) < 0)
                {
                    net.layers[curr_layer_idx].neurons[ curr_neuron_idx ].a = 0;
                }
                else
                {
                    net.layers[curr_layer_idx].neurons[ curr_neuron_idx ].a = net.layers[curr_layer_idx].neurons[ curr_neuron_idx ].z;
                }
            }
            //apply sigmoid to the last layer
            else{
                //for (int i = 0; i < 15; i++) a = a * a;
                net.layers[curr_layer_idx].neurons[ curr_neuron_idx ].a = 1/(1+exp(-net.layers[curr_layer_idx].neurons[ curr_neuron_idx ].z));
            }
        }
    }
    return net;
}