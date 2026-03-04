#pragma once

#include <cstdint>
#include <cstddef>

namespace TIM {
    int64_t checksum(double* field, size_t field_size, double* mask_val);
}
