#include "BatchSet.h"

__m128i M128iWrapper::masks[129];

template class BuiltinUIntWrapper<uint8_t>;
template class BuiltinUIntWrapper<uint16_t>;
template class BuiltinUIntWrapper<uint32_t>;
template class BuiltinUIntWrapper<uint64_t>;
template <typename UInt>
UInt BuiltinUIntWrapper<UInt>::masks[sizeof(UInt) * 8 + 1];
