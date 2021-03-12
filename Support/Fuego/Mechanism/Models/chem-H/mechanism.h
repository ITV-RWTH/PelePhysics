#ifndef MECHANISM_h
#define MECHANISM_h
#include <AMReX_Gpu.H>
#include <AMReX_REAL.H>

#if 0
/* Elements
0  O
1  H
2  N
*/
#endif

/* Species */
#define H2_ID 0
#define H_ID 1
#define O_ID 2
#define O2_ID 3
#define OH_ID 4
#define H2O_ID 5
#define HO2_ID 6
#define H2O2_ID 7
#define N2_ID 8

#define NUM_ELEMENTS 3
#define NUM_SPECIES 9
#define NUM_REACTIONS 27

#define NUM_FIT 4
/* Transport function declarations  */


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetLENIMC(int* LENIMC ) {
    *LENIMC = 38;};


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetLENRMC(int* LENRMC ) {
    *LENRMC = 1854;};


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetNO(int* NO ) {
    *NO = 4;};


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetKK(int* KK ) {
    *KK = 9;};


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetNLITE(int* NLITE ) {
    *NLITE = 2;};


/*Patm in ergs/cm3 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetPATM(amrex::Real* PATM) {
    *PATM =   0.1013250000000000E+07;};


/*the molecular weights in g/mol */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetWT(amrex::Real* WT ) {
    WT[0] = 2.01594000E+00;
    WT[1] = 1.00797000E+00;
    WT[2] = 1.59994000E+01;
    WT[3] = 3.19988000E+01;
    WT[4] = 1.70073700E+01;
    WT[5] = 1.80153400E+01;
    WT[6] = 3.30067700E+01;
    WT[7] = 3.40147400E+01;
    WT[8] = 2.80134000E+01;
};


/*the lennard-jones potential well depth eps/kb in K */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetEPS(amrex::Real* EPS ) {
    EPS[0] = 3.80000000E+01;
    EPS[1] = 1.45000000E+02;
    EPS[2] = 8.00000000E+01;
    EPS[3] = 1.07400000E+02;
    EPS[4] = 8.00000000E+01;
    EPS[5] = 5.72400000E+02;
    EPS[6] = 1.07400000E+02;
    EPS[7] = 1.07400000E+02;
    EPS[8] = 9.75300000E+01;
};


/*the lennard-jones collision diameter in Angstroms */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetSIG(amrex::Real* SIG ) {
    SIG[0] = 2.92000000E+00;
    SIG[1] = 2.05000000E+00;
    SIG[2] = 2.75000000E+00;
    SIG[3] = 3.45800000E+00;
    SIG[4] = 2.75000000E+00;
    SIG[5] = 2.60500000E+00;
    SIG[6] = 3.45800000E+00;
    SIG[7] = 3.45800000E+00;
    SIG[8] = 3.62100000E+00;
};


/*the dipole moment in Debye */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetDIP(amrex::Real* DIP ) {
    DIP[0] = 0.00000000E+00;
    DIP[1] = 0.00000000E+00;
    DIP[2] = 0.00000000E+00;
    DIP[3] = 0.00000000E+00;
    DIP[4] = 0.00000000E+00;
    DIP[5] = 1.84400000E+00;
    DIP[6] = 0.00000000E+00;
    DIP[7] = 0.00000000E+00;
    DIP[8] = 0.00000000E+00;
};


/*the polarizability in cubic Angstroms */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetPOL(amrex::Real* POL ) {
    POL[0] = 7.90000000E-01;
    POL[1] = 0.00000000E+00;
    POL[2] = 0.00000000E+00;
    POL[3] = 1.60000000E+00;
    POL[4] = 0.00000000E+00;
    POL[5] = 0.00000000E+00;
    POL[6] = 0.00000000E+00;
    POL[7] = 0.00000000E+00;
    POL[8] = 1.76000000E+00;
};


/*the rotational relaxation collision number at 298 K */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetZROT(amrex::Real* ZROT ) {
    ZROT[0] = 2.80000000E+02;
    ZROT[1] = 0.00000000E+00;
    ZROT[2] = 0.00000000E+00;
    ZROT[3] = 3.80000000E+00;
    ZROT[4] = 0.00000000E+00;
    ZROT[5] = 4.00000000E+00;
    ZROT[6] = 1.00000000E+00;
    ZROT[7] = 3.80000000E+00;
    ZROT[8] = 4.00000000E+00;
};


/*0: monoatomic, 1: linear, 2: nonlinear */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetNLIN(int* NLIN) {
    NLIN[0] = 1;
    NLIN[1] = 0;
    NLIN[2] = 0;
    NLIN[3] = 1;
    NLIN[4] = 1;
    NLIN[5] = 2;
    NLIN[6] = 2;
    NLIN[7] = 2;
    NLIN[8] = 1;
};


/*Poly fits for the viscosities, dim NO*KK */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetCOFETA(amrex::Real* COFETA) {
    COFETA[0] = -1.37549435E+01;
    COFETA[1] = 9.65530587E-01;
    COFETA[2] = -4.45720114E-02;
    COFETA[3] = 2.05871810E-03;
    COFETA[4] = -1.98744496E+01;
    COFETA[5] = 3.41660514E+00;
    COFETA[6] = -3.63206306E-01;
    COFETA[7] = 1.58671021E-02;
    COFETA[8] = -1.48001581E+01;
    COFETA[9] = 1.79491990E+00;
    COFETA[10] = -1.54008440E-01;
    COFETA[11] = 6.86719439E-03;
    COFETA[12] = -1.68118868E+01;
    COFETA[13] = 2.52362554E+00;
    COFETA[14] = -2.49309128E-01;
    COFETA[15] = 1.10211025E-02;
    COFETA[16] = -1.47696103E+01;
    COFETA[17] = 1.79491990E+00;
    COFETA[18] = -1.54008440E-01;
    COFETA[19] = 6.86719439E-03;
    COFETA[20] = -1.17770937E+01;
    COFETA[21] = -8.26742721E-01;
    COFETA[22] = 3.39009079E-01;
    COFETA[23] = -2.00674327E-02;
    COFETA[24] = -1.67963797E+01;
    COFETA[25] = 2.52362554E+00;
    COFETA[26] = -2.49309128E-01;
    COFETA[27] = 1.10211025E-02;
    COFETA[28] = -1.67813391E+01;
    COFETA[29] = 2.52362554E+00;
    COFETA[30] = -2.49309128E-01;
    COFETA[31] = 1.10211025E-02;
    COFETA[32] = -1.62526779E+01;
    COFETA[33] = 2.24839597E+00;
    COFETA[34] = -2.13428438E-01;
    COFETA[35] = 9.46192413E-03;
};


/*Poly fits for the conductivities, dim NO*KK */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetCOFLAM(amrex::Real* COFLAM) {
    COFLAM[0] = 1.15899058E+01;
    COFLAM[1] = -1.52427727E+00;
    COFLAM[2] = 2.72840752E-01;
    COFLAM[3] = -1.03392618E-02;
    COFLAM[4] = -3.24539191E-01;
    COFLAM[5] = 3.41660514E+00;
    COFLAM[6] = -3.63206306E-01;
    COFLAM[7] = 1.58671021E-02;
    COFLAM[8] = 1.98513952E+00;
    COFLAM[9] = 1.79491990E+00;
    COFLAM[10] = -1.54008440E-01;
    COFLAM[11] = 6.86719439E-03;
    COFLAM[12] = -3.01284291E+00;
    COFLAM[13] = 3.37554994E+00;
    COFLAM[14] = -3.43353119E-01;
    COFLAM[15] = 1.51043444E-02;
    COFLAM[16] = 1.53490799E+01;
    COFLAM[17] = -3.77958145E+00;
    COFLAM[18] = 6.13516524E-01;
    COFLAM[19] = -2.72295753E-02;
    COFLAM[20] = 2.28195645E+01;
    COFLAM[21] = -8.72278946E+00;
    COFLAM[22] = 1.49300487E+00;
    COFLAM[23] = -7.41524047E-02;
    COFLAM[24] = 5.56023763E-01;
    COFLAM[25] = 1.59073590E+00;
    COFLAM[26] = -5.28053839E-02;
    COFLAM[27] = 4.07601571E-04;
    COFLAM[28] = 6.27051982E-01;
    COFLAM[29] = 1.43139617E+00;
    COFLAM[30] = 1.80509282E-03;
    COFLAM[31] = -3.55624900E-03;
    COFLAM[32] = 1.15507063E+01;
    COFLAM[33] = -2.91452379E+00;
    COFLAM[34] = 5.55043580E-01;
    COFLAM[35] = -2.75172461E-02;
};


/*Poly fits for the diffusion coefficients, dim NO*KK*KK */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetCOFD(amrex::Real* COFD) {
    COFD[0] = -1.02395222E+01;
    COFD[1] = 2.15403244E+00;
    COFD[2] = -6.97480266E-02;
    COFD[3] = 3.23666871E-03;
    COFD[4] = -1.11808682E+01;
    COFD[5] = 2.66936727E+00;
    COFD[6] = -1.34411514E-01;
    COFD[7] = 5.92957488E-03;
    COFD[8] = -1.06250182E+01;
    COFD[9] = 2.15849701E+00;
    COFD[10] = -6.53886401E-02;
    COFD[11] = 2.81453370E-03;
    COFD[12] = -1.15797750E+01;
    COFD[13] = 2.43235504E+00;
    COFD[14] = -1.02890179E-01;
    COFD[15] = 4.52903603E-03;
    COFD[16] = -1.06283453E+01;
    COFD[17] = 2.15849701E+00;
    COFD[18] = -6.53886401E-02;
    COFD[19] = 2.81453370E-03;
    COFD[20] = -1.68758926E+01;
    COFD[21] = 4.49460303E+00;
    COFD[22] = -3.64766132E-01;
    COFD[23] = 1.56457153E-02;
    COFD[24] = -1.15806808E+01;
    COFD[25] = 2.43235504E+00;
    COFD[26] = -1.02890179E-01;
    COFD[27] = 4.52903603E-03;
    COFD[28] = -1.15815344E+01;
    COFD[29] = 2.43235504E+00;
    COFD[30] = -1.02890179E-01;
    COFD[31] = 4.52903603E-03;
    COFD[32] = -1.13253458E+01;
    COFD[33] = 2.31195095E+00;
    COFD[34] = -8.63988037E-02;
    COFD[35] = 3.77573452E-03;
    COFD[36] = -1.11808682E+01;
    COFD[37] = 2.66936727E+00;
    COFD[38] = -1.34411514E-01;
    COFD[39] = 5.92957488E-03;
    COFD[40] = -1.43693056E+01;
    COFD[41] = 4.03992999E+00;
    COFD[42] = -3.08044800E-01;
    COFD[43] = 1.32757775E-02;
    COFD[44] = -1.31860117E+01;
    COFD[45] = 3.38003453E+00;
    COFD[46] = -2.25783856E-01;
    COFD[47] = 9.85028660E-03;
    COFD[48] = -1.43712864E+01;
    COFD[49] = 3.70920439E+00;
    COFD[50] = -2.67274113E-01;
    COFD[51] = 1.15967481E-02;
    COFD[52] = -1.31877711E+01;
    COFD[53] = 3.38003453E+00;
    COFD[54] = -2.25783856E-01;
    COFD[55] = 9.85028660E-03;
    COFD[56] = -1.93611051E+01;
    COFD[57] = 5.51579726E+00;
    COFD[58] = -4.76061961E-01;
    COFD[59] = 1.96329391E-02;
    COFD[60] = -1.43717529E+01;
    COFD[61] = 3.70920439E+00;
    COFD[62] = -2.67274113E-01;
    COFD[63] = 1.15967481E-02;
    COFD[64] = -1.43721922E+01;
    COFD[65] = 3.70920439E+00;
    COFD[66] = -2.67274113E-01;
    COFD[67] = 1.15967481E-02;
    COFD[68] = -1.40298830E+01;
    COFD[69] = 3.55837688E+00;
    COFD[70] = -2.47785790E-01;
    COFD[71] = 1.07555332E-02;
    COFD[72] = -1.06250182E+01;
    COFD[73] = 2.15849701E+00;
    COFD[74] = -6.53886401E-02;
    COFD[75] = 2.81453370E-03;
    COFD[76] = -1.31860117E+01;
    COFD[77] = 3.38003453E+00;
    COFD[78] = -2.25783856E-01;
    COFD[79] = 9.85028660E-03;
    COFD[80] = -1.29877365E+01;
    COFD[81] = 2.80841511E+00;
    COFD[82] = -1.52629888E-01;
    COFD[83] = 6.72604927E-03;
    COFD[84] = -1.40864894E+01;
    COFD[85] = 3.07458927E+00;
    COFD[86] = -1.86899591E-01;
    COFD[87] = 8.19829781E-03;
    COFD[88] = -1.30027772E+01;
    COFD[89] = 2.80841511E+00;
    COFD[90] = -1.52629888E-01;
    COFD[91] = 6.72604927E-03;
    COFD[92] = -1.91096797E+01;
    COFD[93] = 5.02608697E+00;
    COFD[94] = -4.26959993E-01;
    COFD[95] = 1.80709910E-02;
    COFD[96] = -1.40916052E+01;
    COFD[97] = 3.07458927E+00;
    COFD[98] = -1.86899591E-01;
    COFD[99] = 8.19829781E-03;
    COFD[100] = -1.40964661E+01;
    COFD[101] = 3.07458927E+00;
    COFD[102] = -1.86899591E-01;
    COFD[103] = 8.19829781E-03;
    COFD[104] = -1.38756407E+01;
    COFD[105] = 2.98558426E+00;
    COFD[106] = -1.75507216E-01;
    COFD[107] = 7.71173691E-03;
    COFD[108] = -1.15797750E+01;
    COFD[109] = 2.43235504E+00;
    COFD[110] = -1.02890179E-01;
    COFD[111] = 4.52903603E-03;
    COFD[112] = -1.43712864E+01;
    COFD[113] = 3.70920439E+00;
    COFD[114] = -2.67274113E-01;
    COFD[115] = 1.15967481E-02;
    COFD[116] = -1.40864894E+01;
    COFD[117] = 3.07458927E+00;
    COFD[118] = -1.86899591E-01;
    COFD[119] = 8.19829781E-03;
    COFD[120] = -1.53110708E+01;
    COFD[121] = 3.37317428E+00;
    COFD[122] = -2.24900439E-01;
    COFD[123] = 9.81228151E-03;
    COFD[124] = -1.41066459E+01;
    COFD[125] = 3.07458927E+00;
    COFD[126] = -1.86899591E-01;
    COFD[127] = 8.19829781E-03;
    COFD[128] = -2.10640014E+01;
    COFD[129] = 5.50980695E+00;
    COFD[130] = -4.78335488E-01;
    COFD[131] = 1.98515434E-02;
    COFD[132] = -1.53187643E+01;
    COFD[133] = 3.37317428E+00;
    COFD[134] = -2.24900439E-01;
    COFD[135] = 9.81228151E-03;
    COFD[136] = -1.53261114E+01;
    COFD[137] = 3.37317428E+00;
    COFD[138] = -2.24900439E-01;
    COFD[139] = 9.81228151E-03;
    COFD[140] = -1.50096240E+01;
    COFD[141] = 3.25515933E+00;
    COFD[142] = -2.09710110E-01;
    COFD[143] = 9.15941830E-03;
    COFD[144] = -1.06283453E+01;
    COFD[145] = 2.15849701E+00;
    COFD[146] = -6.53886401E-02;
    COFD[147] = 2.81453370E-03;
    COFD[148] = -1.31877711E+01;
    COFD[149] = 3.38003453E+00;
    COFD[150] = -2.25783856E-01;
    COFD[151] = 9.85028660E-03;
    COFD[152] = -1.30027772E+01;
    COFD[153] = 2.80841511E+00;
    COFD[154] = -1.52629888E-01;
    COFD[155] = 6.72604927E-03;
    COFD[156] = -1.41066459E+01;
    COFD[157] = 3.07458927E+00;
    COFD[158] = -1.86899591E-01;
    COFD[159] = 8.19829781E-03;
    COFD[160] = -1.30182843E+01;
    COFD[161] = 2.80841511E+00;
    COFD[162] = -1.52629888E-01;
    COFD[163] = 6.72604927E-03;
    COFD[164] = -1.91256261E+01;
    COFD[165] = 5.02608697E+00;
    COFD[166] = -4.26959993E-01;
    COFD[167] = 1.80709910E-02;
    COFD[168] = -1.41119732E+01;
    COFD[169] = 3.07458927E+00;
    COFD[170] = -1.86899591E-01;
    COFD[171] = 8.19829781E-03;
    COFD[172] = -1.41170372E+01;
    COFD[173] = 3.07458927E+00;
    COFD[174] = -1.86899591E-01;
    COFD[175] = 8.19829781E-03;
    COFD[176] = -1.38948667E+01;
    COFD[177] = 2.98558426E+00;
    COFD[178] = -1.75507216E-01;
    COFD[179] = 7.71173691E-03;
    COFD[180] = -1.68758926E+01;
    COFD[181] = 4.49460303E+00;
    COFD[182] = -3.64766132E-01;
    COFD[183] = 1.56457153E-02;
    COFD[184] = -1.93611051E+01;
    COFD[185] = 5.51579726E+00;
    COFD[186] = -4.76061961E-01;
    COFD[187] = 1.96329391E-02;
    COFD[188] = -1.91096797E+01;
    COFD[189] = 5.02608697E+00;
    COFD[190] = -4.26959993E-01;
    COFD[191] = 1.80709910E-02;
    COFD[192] = -2.10640014E+01;
    COFD[193] = 5.50980695E+00;
    COFD[194] = -4.78335488E-01;
    COFD[195] = 1.98515434E-02;
    COFD[196] = -1.91256261E+01;
    COFD[197] = 5.02608697E+00;
    COFD[198] = -4.26959993E-01;
    COFD[199] = 1.80709910E-02;
    COFD[200] = -1.31492641E+01;
    COFD[201] = 1.48004311E+00;
    COFD[202] = 1.60499553E-01;
    COFD[203] = -1.19765679E-02;
    COFD[204] = -2.04177482E+01;
    COFD[205] = 5.31457079E+00;
    COFD[206] = -4.58216496E-01;
    COFD[207] = 1.91825910E-02;
    COFD[208] = -2.04230073E+01;
    COFD[209] = 5.31457079E+00;
    COFD[210] = -4.58216496E-01;
    COFD[211] = 1.91825910E-02;
    COFD[212] = -2.08123325E+01;
    COFD[213] = 5.42470154E+00;
    COFD[214] = -4.69700416E-01;
    COFD[215] = 1.95706904E-02;
    COFD[216] = -1.15806808E+01;
    COFD[217] = 2.43235504E+00;
    COFD[218] = -1.02890179E-01;
    COFD[219] = 4.52903603E-03;
    COFD[220] = -1.43717529E+01;
    COFD[221] = 3.70920439E+00;
    COFD[222] = -2.67274113E-01;
    COFD[223] = 1.15967481E-02;
    COFD[224] = -1.40916052E+01;
    COFD[225] = 3.07458927E+00;
    COFD[226] = -1.86899591E-01;
    COFD[227] = 8.19829781E-03;
    COFD[228] = -1.53187643E+01;
    COFD[229] = 3.37317428E+00;
    COFD[230] = -2.24900439E-01;
    COFD[231] = 9.81228151E-03;
    COFD[232] = -1.41119732E+01;
    COFD[233] = 3.07458927E+00;
    COFD[234] = -1.86899591E-01;
    COFD[235] = 8.19829781E-03;
    COFD[236] = -2.04177482E+01;
    COFD[237] = 5.31457079E+00;
    COFD[238] = -4.58216496E-01;
    COFD[239] = 1.91825910E-02;
    COFD[240] = -1.53265780E+01;
    COFD[241] = 3.37317428E+00;
    COFD[242] = -2.24900439E-01;
    COFD[243] = 9.81228151E-03;
    COFD[244] = -1.53340417E+01;
    COFD[245] = 3.37317428E+00;
    COFD[246] = -2.24900439E-01;
    COFD[247] = 9.81228151E-03;
    COFD[248] = -1.50168028E+01;
    COFD[249] = 3.25515933E+00;
    COFD[250] = -2.09710110E-01;
    COFD[251] = 9.15941830E-03;
    COFD[252] = -1.15815344E+01;
    COFD[253] = 2.43235504E+00;
    COFD[254] = -1.02890179E-01;
    COFD[255] = 4.52903603E-03;
    COFD[256] = -1.43721922E+01;
    COFD[257] = 3.70920439E+00;
    COFD[258] = -2.67274113E-01;
    COFD[259] = 1.15967481E-02;
    COFD[260] = -1.40964661E+01;
    COFD[261] = 3.07458927E+00;
    COFD[262] = -1.86899591E-01;
    COFD[263] = 8.19829781E-03;
    COFD[264] = -1.53261114E+01;
    COFD[265] = 3.37317428E+00;
    COFD[266] = -2.24900439E-01;
    COFD[267] = 9.81228151E-03;
    COFD[268] = -1.41170372E+01;
    COFD[269] = 3.07458927E+00;
    COFD[270] = -1.86899591E-01;
    COFD[271] = 8.19829781E-03;
    COFD[272] = -2.04230073E+01;
    COFD[273] = 5.31457079E+00;
    COFD[274] = -4.58216496E-01;
    COFD[275] = 1.91825910E-02;
    COFD[276] = -1.53340417E+01;
    COFD[277] = 3.37317428E+00;
    COFD[278] = -2.24900439E-01;
    COFD[279] = 9.81228151E-03;
    COFD[280] = -1.53416186E+01;
    COFD[281] = 3.37317428E+00;
    COFD[282] = -2.24900439E-01;
    COFD[283] = 9.81228151E-03;
    COFD[284] = -1.50236516E+01;
    COFD[285] = 3.25515933E+00;
    COFD[286] = -2.09710110E-01;
    COFD[287] = 9.15941830E-03;
    COFD[288] = -1.13253458E+01;
    COFD[289] = 2.31195095E+00;
    COFD[290] = -8.63988037E-02;
    COFD[291] = 3.77573452E-03;
    COFD[292] = -1.40298830E+01;
    COFD[293] = 3.55837688E+00;
    COFD[294] = -2.47785790E-01;
    COFD[295] = 1.07555332E-02;
    COFD[296] = -1.38756407E+01;
    COFD[297] = 2.98558426E+00;
    COFD[298] = -1.75507216E-01;
    COFD[299] = 7.71173691E-03;
    COFD[300] = -1.50096240E+01;
    COFD[301] = 3.25515933E+00;
    COFD[302] = -2.09710110E-01;
    COFD[303] = 9.15941830E-03;
    COFD[304] = -1.38948667E+01;
    COFD[305] = 2.98558426E+00;
    COFD[306] = -1.75507216E-01;
    COFD[307] = 7.71173691E-03;
    COFD[308] = -2.08123325E+01;
    COFD[309] = 5.42470154E+00;
    COFD[310] = -4.69700416E-01;
    COFD[311] = 1.95706904E-02;
    COFD[312] = -1.50168028E+01;
    COFD[313] = 3.25515933E+00;
    COFD[314] = -2.09710110E-01;
    COFD[315] = 9.15941830E-03;
    COFD[316] = -1.50236516E+01;
    COFD[317] = 3.25515933E+00;
    COFD[318] = -2.09710110E-01;
    COFD[319] = 9.15941830E-03;
    COFD[320] = -1.47639290E+01;
    COFD[321] = 3.15955654E+00;
    COFD[322] = -1.97590757E-01;
    COFD[323] = 8.64692156E-03;
};


/*List of specs with small weight, dim NLITE */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetKTDIF(int* KTDIF) {
    KTDIF[0] = 0;
    KTDIF[1] = 1;
};


/*Poly fits for thermal diff ratios, dim NO*NLITE*KK */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetCOFTD(amrex::Real* COFTD) {
    COFTD[0] = 0.00000000E+00;
    COFTD[1] = 0.00000000E+00;
    COFTD[2] = 0.00000000E+00;
    COFTD[3] = 0.00000000E+00;
    COFTD[4] = -1.52534742E-01;
    COFTD[5] = -5.46404022E-05;
    COFTD[6] = 2.93412470E-08;
    COFTD[7] = -4.87091914E-12;
    COFTD[8] = 4.15583337E-01;
    COFTD[9] = 1.09738399E-05;
    COFTD[10] = -3.96021963E-09;
    COFTD[11] = 1.14414443E-12;
    COFTD[12] = 4.42739084E-01;
    COFTD[13] = 7.11770818E-05;
    COFTD[14] = -3.84768062E-08;
    COFTD[15] = 6.86323437E-12;
    COFTD[16] = 4.21932443E-01;
    COFTD[17] = 1.11414935E-05;
    COFTD[18] = -4.02072219E-09;
    COFTD[19] = 1.16162418E-12;
    COFTD[20] = 6.02028221E-02;
    COFTD[21] = 5.61561867E-04;
    COFTD[22] = -2.55372862E-07;
    COFTD[23] = 3.63389913E-11;
    COFTD[24] = 4.44452569E-01;
    COFTD[25] = 7.14525507E-05;
    COFTD[26] = -3.86257187E-08;
    COFTD[27] = 6.88979640E-12;
    COFTD[28] = 4.46070183E-01;
    COFTD[29] = 7.17126069E-05;
    COFTD[30] = -3.87662996E-08;
    COFTD[31] = 6.91487226E-12;
    COFTD[32] = 4.45261966E-01;
    COFTD[33] = 4.94697174E-05;
    COFTD[34] = -2.63023442E-08;
    COFTD[35] = 4.90306217E-12;
    COFTD[36] = 1.52534742E-01;
    COFTD[37] = 5.46404022E-05;
    COFTD[38] = -2.93412470E-08;
    COFTD[39] = 4.87091914E-12;
    COFTD[40] = 0.00000000E+00;
    COFTD[41] = 0.00000000E+00;
    COFTD[42] = 0.00000000E+00;
    COFTD[43] = 0.00000000E+00;
    COFTD[44] = 2.70010150E-01;
    COFTD[45] = 3.61555093E-04;
    COFTD[46] = -1.80744752E-07;
    COFTD[47] = 2.75321248E-11;
    COFTD[48] = 2.20482843E-01;
    COFTD[49] = 4.80164288E-04;
    COFTD[50] = -2.32927944E-07;
    COFTD[51] = 3.46470436E-11;
    COFTD[52] = 2.72041664E-01;
    COFTD[53] = 3.64275376E-04;
    COFTD[54] = -1.82104647E-07;
    COFTD[55] = 2.77392722E-11;
    COFTD[56] = -1.41883744E-01;
    COFTD[57] = 7.66558810E-04;
    COFTD[58] = -3.06550003E-07;
    COFTD[59] = 4.02959502E-11;
    COFTD[60] = 2.20907853E-01;
    COFTD[61] = 4.81089870E-04;
    COFTD[62] = -2.33376944E-07;
    COFTD[63] = 3.47138305E-11;
    COFTD[64] = 2.21308399E-01;
    COFTD[65] = 4.81962174E-04;
    COFTD[66] = -2.33800100E-07;
    COFTD[67] = 3.47767730E-11;
    COFTD[68] = 2.40744421E-01;
    COFTD[69] = 4.45343451E-04;
    COFTD[70] = -2.18173874E-07;
    COFTD[71] = 3.26958506E-11;
};
#endif

/* End of file  */
