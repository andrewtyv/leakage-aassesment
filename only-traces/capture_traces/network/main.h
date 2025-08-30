#include <stdint.h>

uint8_t handle(uint8_t cmd, uint8_t scmd, uint8_t len, uint8_t *buf);
#ifdef DEBUGGING
uint8_t test_handle(uint8_t cmd, uint8_t scmd, uint8_t len, uint8_t *buf);
#endif