#include "mechanism.H"
const int rmap[NUM_REACTIONS] = {
  11,  24,  85,  98,  141, 189, 195, 223, 226, 235, 199, 6,   7,   8,   9,
  10,  18,  20,  21,  30,  31,  52,  73,  128, 129, 155, 183, 207, 232, 236,
  0,   1,   2,   3,   4,   5,   12,  13,  14,  15,  16,  17,  19,  22,  23,
  25,  26,  27,  28,  29,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,
  42,  43,  44,  45,  46,  47,  48,  49,  50,  51,  53,  54,  55,  56,  57,
  58,  59,  60,  61,  62,  63,  64,  65,  66,  67,  68,  69,  70,  71,  72,
  74,  75,  76,  77,  78,  79,  80,  81,  82,  83,  84,  86,  87,  88,  89,
  90,  91,  92,  93,  94,  95,  96,  97,  99,  100, 101, 102, 103, 104, 105,
  106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120,
  121, 122, 123, 124, 125, 126, 127, 130, 131, 132, 133, 134, 135, 136, 137,
  138, 139, 140, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153,
  154, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169,
  170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 184, 185,
  186, 187, 188, 190, 191, 192, 193, 194, 196, 197, 198, 200, 201, 202, 203,
  204, 205, 206, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219,
  220, 221, 222, 224, 225, 227, 228, 229, 230, 231, 233, 234, 237, 238, 239,
  240, 241, 242};

// Returns 0-based map of reaction order
void
GET_RMAP(int* _rmap)
{
  for (int j = 0; j < NUM_REACTIONS; ++j) {
    _rmap[j] = rmap[j];
  }
}

// Returns a count of gas species in a gas reaction, and their indices
// and stoichiometric coefficients. (Eq 50)
void
CKINU(const int i, int& nspec, int ki[], int nu[])
{
  const int ns[NUM_GAS_REACTIONS] = {
    4, 4, 4, 4, 3, 3, 2, 2, 3, 3, 3, 3, 4, 3, 4, 4, 4, 4, 2, 3, 3, 3, 3, 3, 2,
    4, 4, 4, 4, 4, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4,
    4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 3, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 2, 3, 3, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 3, 4,
    4, 4, 4, 4, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5,
    4, 4, 4, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 3, 3, 4, 4, 3,
    4, 4, 3, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 2, 4,
    4, 3, 3, 4, 4, 4, 5, 3, 3, 4, 3, 3, 4, 4, 4, 4, 4, 4};
  const int kiv[NUM_GAS_REACTIONS * 5] = {
    1,  8,  3,  4,  0,  2,  3,  1,  4,  0,  2,  3,  1,  4,  0,  2,  4,  1,  5,
    0,  4,  5,  3,  0,  0,  4,  5,  3,  0,  0,  2,  1,  0,  0,  0,  3,  8,  0,
    0,  0,  1,  3,  4,  0,  0,  5,  1,  4,  0,  0,  5,  1,  4,  0,  0,  1,  8,
    6,  0,  0,  1,  6,  2,  8,  0,  1,  6,  4,  0,  0,  1,  6,  5,  3,  0,  6,
    3,  8,  4,  0,  6,  4,  5,  8,  0,  6,  4,  5,  8,  0,  1,  2,  0,  0,  0,
    1,  8,  4,  0,  0,  1,  3,  4,  0,  0,  1,  4,  5,  0,  0,  6,  7,  8,  0,
    0,  6,  7,  8,  0,  0,  7,  4,  0,  0,  0,  1,  7,  5,  4,  0,  1,  7,  2,
    6,  0,  7,  3,  6,  4,  0,  7,  4,  5,  6,  0,  7,  4,  5,  6,  0,  0,  9,
    0,  0,  0,  9,  3,  18, 0,  0,  9,  20, 19, 3,  0,  9,  8,  18, 3,  0,  9,
    4,  1,  18, 0,  9,  18, 0,  3,  0,  1,  10, 2,  9,  0,  10, 3,  1,  18, 0,
    10, 3,  9,  4,  0,  10, 4,  1,  21, 0,  10, 4,  5,  9,  0,  10, 8,  21, 3,
    0,  10, 8,  18, 4,  0,  10, 2,  0,  0,  0,  10, 1,  0,  0,  0,  9,  10, 1,
    0,  0,  10, 18, 1,  19, 0,  10, 18, 0,  4,  0,  10, 20, 19, 4,  0,  10, 20,
    21, 18, 0,  25, 10, 11, 20, 0,  19, 10, 21, 0,  0,  11, 1,  10, 0,  0,  1,
    11, 2,  10, 0,  11, 3,  10, 4,  0,  11, 3,  1,  21, 0,  11, 3,  1,  21, 0,
    11, 4,  5,  10, 0,  6,  11, 22, 4,  0,  6,  11, 12, 8,  0,  11, 8,  22, 3,
    0,  11, 8,  21, 4,  0,  9,  11, 1,  0,  0,  11, 18, 13, 4,  0,  11, 18, 5,
    0,  0,  11, 18, 5,  0,  0,  11, 20, 22, 18, 0,  11, 20, 5,  19, 0,  21, 11,
    12, 18, 0,  25, 11, 12, 20, 0,  10, 9,  11, 0,  0,  11, 10, 12, 0,  0,  10,
    11, 9,  12, 0,  14, 1,  13, 0,  0,  1,  14, 2,  13, 0,  14, 3,  11, 18, 0,
    14, 3,  13, 4,  0,  14, 4,  5,  13, 0,  14, 18, 19, 11, 0,  14, 20, 25, 13,
    0,  14, 10, 11, 13, 0,  14, 11, 12, 13, 0,  10, 11, 1,  14, 0,  11, 2,  14,
    0,  0,  6,  14, 7,  13, 0,  11, 16, 0,  0,  0,  16, 2,  14, 0,  0,  16, 1,
    15, 0,  0,  1,  16, 2,  15, 0,  16, 4,  5,  15, 0,  16, 11, 15, 12, 0,  14,
    16, 15, 0,  0,  16, 10, 15, 11, 0,  16, 3,  5,  14, 0,  16, 3,  15, 4,  0,
    16, 18, 21, 15, 0,  16, 20, 25, 15, 0,  16, 20, 26, 15, 0,  15, 1,  14, 0,
    0,  1,  15, 2,  14, 0,  1,  15, 2,  17, 0,  1,  15, 10, 12, 0,  15, 3,  14,
    4,  0,  15, 3,  21, 11, 0,  15, 3,  1,  11, 18, 15, 4,  5,  14, 0,  15, 4,
    17, 5,  0,  15, 4,  21, 12, 0,  6,  15, 7,  14, 0,  6,  15, 16, 8,  0,  15,
    11, 14, 12, 0,  15, 11, 17, 12, 0,  15, 10, 14, 11, 0,  14, 15, 16, 13, 0,
    15, 0,  12, 0,  0,  14, 15, 13, 0,  0,  11, 2,  17, 0,  0,  17, 1,  13, 0,
    0,  1,  17, 2,  13, 0,  1,  17, 1,  14, 0,  17, 3,  13, 4,  0,  17, 3,  11,
    18, 0,  17, 4,  5,  13, 0,  17, 4,  1,  11, 18, 17, 6,  11, 18, 4,  17, 6,
    7,  13, 0,  17, 8,  11, 20, 0,  17, 11, 12, 13, 0,  22, 23, 0,  0,  0,  23,
    1,  21, 0,  0,  1,  23, 11, 4,  0,  1,  23, 2,  21, 0,  23, 3,  21, 4,  0,
    23, 3,  21, 4,  0,  23, 4,  5,  21, 0,  23, 6,  7,  21, 0,  23, 8,  21, 6,
    0,  23, 11, 15, 4,  0,  23, 11, 17, 5,  0,  23, 11, 21, 12, 0,  23, 20, 21,
    25, 0,  24, 11, 4,  0,  0,  1,  24, 2,  23, 0,  1,  24, 2,  22, 0,  24, 3,
    23, 4,  0,  24, 3,  22, 4,  0,  24, 4,  5,  23, 0,  24, 4,  22, 5,  0,  11,
    24, 23, 12, 0,  11, 24, 22, 12, 0,  10, 24, 23, 11, 0,  10, 24, 22, 11, 0,
    6,  24, 7,  23, 0,  6,  24, 22, 7,  0,  21, 23, 24, 18, 0,  12, 1,  11, 0,
    0,  1,  12, 2,  11, 0,  12, 4,  5,  11, 0,  6,  12, 7,  11, 0,  12, 3,  11,
    4,  0,  11, 12, 2,  15, 0,  13, 1,  0,  0,  0,  1,  13, 2,  0,  0,  13, 4,
    5,  0,  0,  13, 18, 21, 0,  0,  10, 13, 0,  11, 0,  11, 13, 0,  12, 0,  13,
    3,  1,  19, 0,  13, 3,  10, 18, 0,  13, 3,  0,  4,  0,  13, 8,  6,  0,  0,
    13, 0,  14, 0,  0,  21, 1,  18, 0,  0,  21, 3,  18, 4,  0,  1,  21, 2,  18,
    0,  21, 4,  5,  18, 0,  21, 8,  6,  18, 0,  21, 9,  10, 18, 0,  21, 9,  1,
    19, 0,  21, 10, 11, 18, 0,  21, 18, 19, 4,  0,  21, 20, 25, 18, 0,  21, 5,
    19, 0,  0,  27, 1,  18, 0,  0,  1,  27, 1,  21, 0,  1,  27, 10, 4,  0,  27,
    3,  18, 4,  0,  27, 4,  1,  25, 0,  27, 8,  20, 4,  0,  18, 4,  25, 0,  0,
    6,  18, 20, 4,  0,  2,  20, 1,  25, 0,  2,  20, 1,  26, 0,  1,  20, 18, 4,
    0,  20, 3,  18, 8,  0,  18, 3,  20, 0,  0,  20, 18, 8,  0,  0,  6,  20, 25,
    8,  0,  6,  20, 26, 8,  0,  19, 0,  3,  0,  0,  1,  19, 0,  4,  0,  1,  19,
    0,  4,  0,  19, 3,  18, 0,  0,  19, 3,  0,  8,  0,  19, 4,  6,  0,  0,  19,
    18, 0,  20, 0,  9,  19, 0,  18, 0,  22, 1,  21, 0,  0,  1,  22, 2,  21, 0,
    1,  22, 11, 4,  0,  22, 3,  21, 4,  0,  22, 4,  5,  21, 0,  22, 11, 21, 12,
    0,  22, 6,  7,  21, 0,  22, 8,  21, 6,  0,  22, 18, 21, 0,  0,  22, 20, 21,
    25, 0,  25, 3,  20, 4,  0,  1,  25, 21, 4,  0,  1,  25, 5,  18, 0,  25, 4,
    5,  20, 0,  25, 20, 29, 18, 0,  25, 5,  18, 20, 0,  26, 25, 0,  0,  0,  26,
    3,  20, 4,  0,  26, 4,  5,  20, 0,  20, 3,  28, 0,  0,  20, 18, 28, 0,  0,
    1,  28, 20, 4,  0,  28, 3,  20, 8,  0,  28, 4,  6,  20, 0,  6,  28, 20, 8,
    4,  28, 18, 8,  0,  0,  28, 20, 8,  0,  0,  6,  28, 29, 8,  0,  20, 4,  29,
    0,  0,  6,  18, 29, 0,  0,  1,  29, 2,  28, 0,  1,  29, 5,  20, 0,  1,  29,
    25, 4,  0,  29, 4,  5,  28, 0,  29, 11, 12, 28, 0,  29, 11, 12, 28, 0};
  const int nuv[NUM_GAS_REACTIONS * 5] = {
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -2, 1,  1, 0, 0, -2, 1,  1, 0, 0, -1, 2,  0, 0, 0, -2, 1,  0, 0, 0,
    -1, -1, 1, 0, 0, -1, 1,  1, 0, 0, -1, 1,  1, 0, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 2, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -2, 1,  0, 0, 0, -2, -1, 2, 0, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -2, 1,  1, 0, 0, -2, 1,  1, 0, 0,
    -1, 2,  0, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, 2,  0, 0, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -2, 1,  1, 0, 0,
    -2, 2,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 2, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -2, 1,  1, 0, 0, -2, 1,  1, 0, 0,
    -1, -1, 1, 1, 0, -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -2, 1,  1, 0, 0,
    -1, -1, 1, 1, 0, -2, 1,  0, 0, 0, -1, 1,  1, 0, 0, -1, 1,  1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 2, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, 1,  1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 1, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -2, 1,  2, 0, 0, -2, 1,  1, 0, 0,
    -2, 1,  1, 0, 0, -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 1,
    -1, -1, 1, 1, 1, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, 1,  0, 0, 0, -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, 1,  1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -2, 1,  1, 0, 0,
    -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -2, 1,  1, 0, 0, -1, 1,  1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0,
    -2, 2,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, 1,  1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 2, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, 1,  1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 2, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -2, 1,  1, 1, 0, -1, 1,  0, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -2, 1,  1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 1,
    -1, 1,  1, 0, 0, -2, 2,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0};
  if (i < 1) {
    // Return max num species per reaction
    nspec = 5;
  } else {
    if (i > NUM_GAS_REACTIONS) {
      nspec = -1;
    } else {
      nspec = ns[i - 1];
      for (int j = 0; j < nspec; ++j) {
        ki[j] = kiv[(i - 1) * 5 + j] + 1;
        nu[j] = nuv[(i - 1) * 5 + j];
      }
    }
  }
}

// Returns the progress rates of each reactions
// Given P, T, and mole fractions
void
CKKFKR(
  const amrex::Real P,
  const amrex::Real T,
  const amrex::Real x[],
  amrex::Real q_f[],
  amrex::Real q_r[])
{
  amrex::Real c[30]; // temporary storage
  amrex::Real PORT =
    1e6 * P / (8.31446261815324e+07 * T); // 1e6 * P/RT so c goes to SI units

  // Compute conversion, see Eq 10
  for (int id = 0; id < 30; ++id) {
    c[id] = x[id] * PORT;
  }

  // convert to chemkin units
  progressRateFR(q_f, q_r, c, T);

  // convert to chemkin units
  for (int id = 0; id < 243; ++id) {
    q_f[id] *= 1.0e-6;
    q_r[id] *= 1.0e-6;
  }
}

// compute the progress rate for each reaction
// USES progressRate : todo switch to GPU
void
progressRateFR(
  amrex::Real* q_f, amrex::Real* q_r, amrex::Real* sc, amrex::Real T)
{
  const amrex::Real invT = 1.0 / T;
  const amrex::Real logT = log(T);
  // compute the Gibbs free energy
  amrex::Real g_RT[30];
  gibbs(g_RT, T);

  amrex::Real sc_qss[1];
  comp_qfqr(q_f, q_r, sc, sc_qss, T, invT, logT);
}

// save atomic weights into array
void
atomicWeight(amrex::Real* awt)
{
  awt[0] = 15.999000; // O
  awt[1] = 1.008000;  // H
  awt[2] = 14.007000; // N
}

// get atomic weight for all elements
void
CKAWT(amrex::Real* awt)
{
  atomicWeight(awt);
}

// Returns the elemental composition
// of the speciesi (mdim is num of elements)
void
CKNCF(int* ncf)
{
  int kd = 3;
  // Zero ncf
  for (int id = 0; id < kd * 30; ++id) {
    ncf[id] = 0;
  }

  // N2
  ncf[0 * kd + 2] = 2; // N

  // H
  ncf[1 * kd + 1] = 1; // H

  // H2
  ncf[2 * kd + 1] = 2; // H

  // O
  ncf[3 * kd + 0] = 1; // O

  // OH
  ncf[4 * kd + 1] = 1; // H
  ncf[4 * kd + 0] = 1; // O

  // H2O
  ncf[5 * kd + 1] = 2; // H
  ncf[5 * kd + 0] = 1; // O

  // HO2
  ncf[6 * kd + 1] = 1; // H
  ncf[6 * kd + 0] = 2; // O

  // H2O2
  ncf[7 * kd + 1] = 2; // H
  ncf[7 * kd + 0] = 2; // O

  // O2
  ncf[8 * kd + 0] = 2; // O

  // N
  ncf[9 * kd + 2] = 1; // N

  // NH
  ncf[10 * kd + 1] = 1; // H
  ncf[10 * kd + 2] = 1; // N

  // NH2
  ncf[11 * kd + 1] = 2; // H
  ncf[11 * kd + 2] = 1; // N

  // NH3
  ncf[12 * kd + 1] = 3; // H
  ncf[12 * kd + 2] = 1; // N

  // NNH
  ncf[13 * kd + 1] = 1; // H
  ncf[13 * kd + 2] = 2; // N

  // N2H2
  ncf[14 * kd + 1] = 2; // H
  ncf[14 * kd + 2] = 2; // N

  // N2H3
  ncf[15 * kd + 1] = 3; // H
  ncf[15 * kd + 2] = 2; // N

  // N2H4
  ncf[16 * kd + 1] = 4; // H
  ncf[16 * kd + 2] = 2; // N

  // H2NN
  ncf[17 * kd + 1] = 2; // H
  ncf[17 * kd + 2] = 2; // N

  // NO
  ncf[18 * kd + 2] = 1; // N
  ncf[18 * kd + 0] = 1; // O

  // N2O
  ncf[19 * kd + 2] = 2; // N
  ncf[19 * kd + 0] = 1; // O

  // NO2
  ncf[20 * kd + 2] = 1; // N
  ncf[20 * kd + 0] = 2; // O

  // HNO
  ncf[21 * kd + 1] = 1; // H
  ncf[21 * kd + 2] = 1; // N
  ncf[21 * kd + 0] = 1; // O

  // H2NO
  ncf[22 * kd + 1] = 2; // H
  ncf[22 * kd + 2] = 1; // N
  ncf[22 * kd + 0] = 1; // O

  // HNOH
  ncf[23 * kd + 1] = 2; // H
  ncf[23 * kd + 2] = 1; // N
  ncf[23 * kd + 0] = 1; // O

  // NH2OH
  ncf[24 * kd + 1] = 3; // H
  ncf[24 * kd + 2] = 1; // N
  ncf[24 * kd + 0] = 1; // O

  // HONO
  ncf[25 * kd + 1] = 1; // H
  ncf[25 * kd + 2] = 1; // N
  ncf[25 * kd + 0] = 2; // O

  // HNO2
  ncf[26 * kd + 1] = 1; // H
  ncf[26 * kd + 2] = 1; // N
  ncf[26 * kd + 0] = 2; // O

  // HON
  ncf[27 * kd + 1] = 1; // H
  ncf[27 * kd + 2] = 1; // N
  ncf[27 * kd + 0] = 1; // O

  // NO3
  ncf[28 * kd + 2] = 1; // N
  ncf[28 * kd + 0] = 3; // O

  // HNO3
  ncf[29 * kd + 1] = 1; // H
  ncf[29 * kd + 2] = 1; // N
  ncf[29 * kd + 0] = 3; // O
}

// Returns the vector of strings of element names
void
CKSYME_STR(amrex::Vector<std::string>& ename)
{
  ename.resize(3);
  ename[0] = "O";
  ename[1] = "H";
  ename[2] = "N";
}

// Returns the vector of strings of species names
void
CKSYMS_STR(amrex::Vector<std::string>& kname)
{
  kname.resize(30);
  kname[0] = "N2";
  kname[1] = "H";
  kname[2] = "H2";
  kname[3] = "O";
  kname[4] = "OH";
  kname[5] = "H2O";
  kname[6] = "HO2";
  kname[7] = "H2O2";
  kname[8] = "O2";
  kname[9] = "N";
  kname[10] = "NH";
  kname[11] = "NH2";
  kname[12] = "NH3";
  kname[13] = "NNH";
  kname[14] = "N2H2";
  kname[15] = "N2H3";
  kname[16] = "N2H4";
  kname[17] = "H2NN";
  kname[18] = "NO";
  kname[19] = "N2O";
  kname[20] = "NO2";
  kname[21] = "HNO";
  kname[22] = "H2NO";
  kname[23] = "HNOH";
  kname[24] = "NH2OH";
  kname[25] = "HONO";
  kname[26] = "HNO2";
  kname[27] = "HON";
  kname[28] = "NO3";
  kname[29] = "HNO3";
}

// compute the sparsity pattern of the chemistry Jacobian
void
SPARSITY_INFO(int* nJdata, const int* consP, int NCELLS)
{
  amrex::GpuArray<amrex::Real, 961> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 30> conc = {0.0};
  for (int n = 0; n < 30; n++) {
    conc[n] = 1.0 / 30.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 31; k++) {
    for (int l = 0; l < 31; l++) {
      if (Jac[31 * k + l] != 0.0) {
        nJdata_tmp = nJdata_tmp + 1;
      }
    }
  }

  *nJdata = NCELLS * nJdata_tmp;
}

// compute the sparsity pattern of the system Jacobian
void
SPARSITY_INFO_SYST(int* nJdata, const int* consP, int NCELLS)
{
  amrex::GpuArray<amrex::Real, 961> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 30> conc = {0.0};
  for (int n = 0; n < 30; n++) {
    conc[n] = 1.0 / 30.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 31; k++) {
    for (int l = 0; l < 31; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[31 * k + l] != 0.0) {
          nJdata_tmp = nJdata_tmp + 1;
        }
      }
    }
  }

  *nJdata = NCELLS * nJdata_tmp;
}

// compute the sparsity pattern of the simplified (for preconditioning) system
// Jacobian
void
SPARSITY_INFO_SYST_SIMPLIFIED(int* nJdata, const int* consP)
{
  amrex::GpuArray<amrex::Real, 961> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 30> conc = {0.0};
  for (int n = 0; n < 30; n++) {
    conc[n] = 1.0 / 30.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 31; k++) {
    for (int l = 0; l < 31; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[31 * k + l] != 0.0) {
          nJdata_tmp = nJdata_tmp + 1;
        }
      }
    }
  }

  nJdata[0] = nJdata_tmp;
}

// compute the sparsity pattern of the chemistry Jacobian in CSC format -- base
// 0
void
SPARSITY_PREPROC_CSC(int* rowVals, int* colPtrs, const int* consP, int NCELLS)
{
  amrex::GpuArray<amrex::Real, 961> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 30> conc = {0.0};
  for (int n = 0; n < 30; n++) {
    conc[n] = 1.0 / 30.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int nc = 0; nc < NCELLS; nc++) {
    int offset_row = nc * 31;
    int offset_col = nc * 31;
    for (int k = 0; k < 31; k++) {
      for (int l = 0; l < 31; l++) {
        if (Jac[31 * k + l] != 0.0) {
          rowVals[nJdata_tmp] = l + offset_row;
          nJdata_tmp = nJdata_tmp + 1;
        }
      }
      colPtrs[offset_col + (k + 1)] = nJdata_tmp;
    }
  }
}

// compute the sparsity pattern of the chemistry Jacobian in CSR format -- base
// 0
void
SPARSITY_PREPROC_CSR(
  int* colVals, int* rowPtrs, const int* consP, int NCELLS, int base)
{
  amrex::GpuArray<amrex::Real, 961> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 30> conc = {0.0};
  for (int n = 0; n < 30; n++) {
    conc[n] = 1.0 / 30.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtrs[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 31;
      for (int l = 0; l < 31; l++) {
        for (int k = 0; k < 31; k++) {
          if (Jac[31 * k + l] != 0.0) {
            colVals[nJdata_tmp - 1] = k + 1 + offset;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
        rowPtrs[offset + (l + 1)] = nJdata_tmp;
      }
    }
  } else {
    rowPtrs[0] = 0;
    int nJdata_tmp = 0;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 31;
      for (int l = 0; l < 31; l++) {
        for (int k = 0; k < 31; k++) {
          if (Jac[31 * k + l] != 0.0) {
            colVals[nJdata_tmp] = k + offset;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
        rowPtrs[offset + (l + 1)] = nJdata_tmp;
      }
    }
  }
}

// compute the sparsity pattern of the system Jacobian
// CSR format BASE is user choice
void
SPARSITY_PREPROC_SYST_CSR(
  int* colVals, int* rowPtr, const int* consP, int NCELLS, int base)
{
  amrex::GpuArray<amrex::Real, 961> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 30> conc = {0.0};
  for (int n = 0; n < 30; n++) {
    conc[n] = 1.0 / 30.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 31;
      for (int l = 0; l < 31; l++) {
        for (int k = 0; k < 31; k++) {
          if (k == l) {
            colVals[nJdata_tmp - 1] = l + 1 + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[31 * k + l] != 0.0) {
              colVals[nJdata_tmp - 1] = k + 1 + offset;
              nJdata_tmp = nJdata_tmp + 1;
            }
          }
        }
        rowPtr[offset + (l + 1)] = nJdata_tmp;
      }
    }
  } else {
    rowPtr[0] = 0;
    int nJdata_tmp = 0;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 31;
      for (int l = 0; l < 31; l++) {
        for (int k = 0; k < 31; k++) {
          if (k == l) {
            colVals[nJdata_tmp] = l + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[31 * k + l] != 0.0) {
              colVals[nJdata_tmp] = k + offset;
              nJdata_tmp = nJdata_tmp + 1;
            }
          }
        }
        rowPtr[offset + (l + 1)] = nJdata_tmp;
      }
    }
  }
}

// compute the sparsity pattern of the simplified (for precond) system Jacobian
// on CPU BASE 0
void
SPARSITY_PREPROC_SYST_SIMPLIFIED_CSC(
  int* rowVals, int* colPtrs, int* indx, const int* consP)
{
  amrex::GpuArray<amrex::Real, 961> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 30> conc = {0.0};
  for (int n = 0; n < 30; n++) {
    conc[n] = 1.0 / 30.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int k = 0; k < 31; k++) {
    for (int l = 0; l < 31; l++) {
      if (k == l) {
        rowVals[nJdata_tmp] = l;
        indx[nJdata_tmp] = 31 * k + l;
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[31 * k + l] != 0.0) {
          rowVals[nJdata_tmp] = l;
          indx[nJdata_tmp] = 31 * k + l;
          nJdata_tmp = nJdata_tmp + 1;
        }
      }
    }
    colPtrs[k + 1] = nJdata_tmp;
  }
}

// compute the sparsity pattern of the simplified (for precond) system Jacobian
// CSR format BASE is under choice
void
SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(
  int* colVals, int* rowPtr, const int* consP, int base)
{
  amrex::GpuArray<amrex::Real, 961> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 30> conc = {0.0};
  for (int n = 0; n < 30; n++) {
    conc[n] = 1.0 / 30.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int l = 0; l < 31; l++) {
      for (int k = 0; k < 31; k++) {
        if (k == l) {
          colVals[nJdata_tmp - 1] = l + 1;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[31 * k + l] != 0.0) {
            colVals[nJdata_tmp - 1] = k + 1;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
      }
      rowPtr[l + 1] = nJdata_tmp;
    }
  } else {
    rowPtr[0] = 0;
    int nJdata_tmp = 0;
    for (int l = 0; l < 31; l++) {
      for (int k = 0; k < 31; k++) {
        if (k == l) {
          colVals[nJdata_tmp] = l;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[31 * k + l] != 0.0) {
            colVals[nJdata_tmp] = k;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
      }
      rowPtr[l + 1] = nJdata_tmp;
    }
  }
}
