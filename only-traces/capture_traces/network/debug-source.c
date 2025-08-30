#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include "main.h"



// Put your handler definitions here.
// uint8_t set_key(uint8_t cmd, uint8_t scmd, uint8_t len, uint8_t* buf);

uint8_t TRIGGER_SETUP = 0;
uint8_t UART_SETUP = 0;
uint8_t PLATFORM_SETUP = 0;
uint8_t SS_SETUP = 0;

void platform_init() {
    printf("Initiated platform!\n");
    PLATFORM_SETUP = 1;
}
void init_uart() {
    if (PLATFORM_SETUP == 0) {
        printf("Tried to setup UART without Platform setup");
        exit(-1);
    }
    printf("Initiated UArt!\n");
    UART_SETUP = 1;
}
void trigger_setup() {
    if (UART_SETUP == 0) {
        printf("Tried to setup Trigger without UART setup");
        exit(-1);
    }
    printf("Trigger Setup!\n");
    TRIGGER_SETUP = 1;
}
void trigger_high() {
    if (TRIGGER_SETUP == 0) {
        printf("Tried to set trigger to high without trigger setup");
        exit(-1);
    }
    printf("Trigger put at High!\n");
}
void trigger_low() {
    if (TRIGGER_SETUP == 0) {
        printf("Tried to set trigger to low without trigger setup");
        exit(-1);
    }
    printf("Trigger put at Low!\n");
}

void simpleserial_init() {
    if (TRIGGER_SETUP == 0) {
        printf("Tried to setup SimpleSerial without Trigger setup");
        exit(-1);
    }
    printf("Initiated SimpleSerial!\n");
    SS_SETUP = 1;
}
void simpleserial_put(uint8_t cmd, uint8_t len, uint8_t* buf) {
    printf("Put '%c' on SS with len %u", cmd, len);
    printf("Data: [\n");

    for (int i = 0; i < len; ++i) {
        printf("%x, ", buf[i]);
        if (i != 0 && i != (len - 1) && i % 3 == 0)
            printf("\n");
    }

    printf("\n]\n");
}
void simpleserial_addcmd(uint8_t cmd, uint8_t len, uint8_t (*f)(uint8_t cmd, uint8_t scmd, uint8_t len, uint8_t *buf)) {
    printf("Listening for '%c' on SS with preferred len %u\n", cmd, len);
}
void simpleserial_get(uint8_t cmd, uint8_t len, uint8_t (*f)(uint8_t cmd, uint8_t scmd, uint8_t len, uint8_t *buf)) {
    // Put your debug code here
    printf("Debugging started!\n");
    
    uint8_t cmd2 = 'p';
    uint8_t scmd2 = 0x00;
    uint8_t len2 = 4; 
    float value = 0.657f;
    uint8_t buffer[4];  // 4-byte buffer

    // Copy the float's bytes into the buffer using pointer casting
    uint8_t *float_as_bytes = (uint8_t*)&value;
    for (int i = 0; i < 4; i++) {
        buffer[i] = float_as_bytes[i];
    }

    //for (int i = 0; i < 25; i++)
    //test_handle(cmd2, scmd2, len2, buffer);
    handle(cmd2, scmd2, len2, buffer);

    printf("Debugging ended!\n");

    exit(0);
}