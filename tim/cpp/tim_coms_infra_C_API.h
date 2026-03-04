#pragma once

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

int64_t tim_chksum_c(double* field, size_t field_size, double* mask_val);

#ifdef __cplusplus
}
#endif
