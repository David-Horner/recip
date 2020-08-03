#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <inttypes.h>
#include <stdio.h>
#include <string.h>
#include <fenv.h>
#include <float.h>
#include <algorithm>
#include <vector>
#include <iostream>

#if !defined(__STDC_IEC_559__) && !defined(__APPLE__)
# error Need IEEE 754 FP
#endif

using namespace std;

union f32_union {
  float f;
  uint32_t i;

  f32_union(uint32_t i) : i(i) {}
  f32_union(unsigned long i) : i(uint32_t(i)) {}
  f32_union(float f) : f(f) {}
};


using std::string;
using std::vector;

#define S 23               // significand bits in binary32
#define E 8                // exponent bits in binary32
#define B ((1UL<<(E-1))-1) // binary32 exponent bias

int ip = 0, op= 0;  // ip is input index width, op is output estimate width
int ipN;            // number of elements in LUT based on ip

#define ipMax   10
#define N (1UL<<ipMax)        // max number of LUT entries

uint32_t rsqrt_lut[N];
uint32_t recip_lut[N];

#define RIGHT 1
#define LEFT  2
#define OTHER 3
#define RERR  4
#define LERR  5
#define OERR  6
#define IDXN  RIGHT|LEFT|OTHER|RERR|LERR|OERR


double rsqrt_lut_errf[N][IDXN];   // error val per LUT entry
double recip_lut_errf[N][IDXN];

// input arguments settings
  int verilogn = 0;  // print verilog LUT definition
  int testn    = 0;  // 1 is short test, 2 is long test
  int detailn  = 0;  // prints LUT entry values/errors
  int summary  = 1;  // default 1 line summary of inputs and informatives


#define rsqrt_lut_idx(sig, exp) (((sig) >> (S-ip+1)) | ((~exp & 1UL) << (ip-1)))
#define recip_lut_idx(sig) ((sig) >> (S-ip))

uint32_t estimate_rsqrt_sig(uint32_t idx)
{

  // P-bit index corresponds to {exp[0], sig[S-1:S-(ip-1)]}
  uint32_t exp = (idx >> (ip-1)) ? B-1 : B-2; // 1 bit from exp -> [0.25, 1.0)
  uint32_t sig = (idx & ((1UL<<(ip-1))-1)) << (S-(ip-1)); // ip-1 bits from sig

  // sqrt(leftmost point on interval)
  // (If P is increased substantially, need to increase precision beyond double.)
  f32_union in = (exp << S) | sig;
  double left = sqrt(double(in.f));

  // sqrt(rightmost point on interval)
  f32_union in1 = in.i + (1UL<<(S-ip+1));
  double right = sqrt(         (double(in1.f)     ));

  // Naively search the space of 2^P output values for the one that minimizes
  // the maximum error on the interval.  Since the function is monotonic,
  // evaluating the error on the extremes of the interval suffices.
  // (This could obviously be done more efficiently, but 2^P is small.)
  double best_error = INFINITY;
  f32_union best = 0.0f;
  f32_union base = B << S; // [1.0, 2.0)
  for (f32_union cand = base; cand.i < base.i + (1UL<<S); cand.i += 1UL<<(S-op)) {
    double lerr = fabs(1.0 - double(cand.f) * left);
    double rerr = fabs(1.0 - double(cand.f) * right);
    double error = max(lerr,rerr);
    if (error < best_error) {
      best_error = error;
      best = cand;
  //printf(" lerr %f in.f %f rerr %f in1.f %f\n", lerr, in.f, rerr, in1.f);
      rsqrt_lut_errf[idx][LERR]  = lerr;  // error val per LUT entry
      rsqrt_lut_errf[idx][LEFT]  = in.f;
      rsqrt_lut_errf[idx][RERR]  = rerr;
      rsqrt_lut_errf[idx][RIGHT] = in1.f;
    }
  }

  // Return P MSBs of mantissa
  return (best.i >> (S-op)) & ((1UL<<op)-1);
}

uint32_t estimate_recip_sig(uint32_t idx)
{

  // P-bit index corresponds to sig[S-1:S-ip]
  uint32_t sig = idx << (S-ip);
  uint32_t exp = B-1; // [0.5, 1.0)

  // Leftmost point on interval
  f32_union in = (exp << S) | sig;
  double left = in.f;

  // Rightmost point on interval
  f32_union in1 = in.i + (1UL<<(S-ip));
  double right =          (double(in1.f)     );
if (right != double(in1.f)) printf(" right != in1.f %A %A\n", right, double(in1.f));

  // Naively search the space of 2^P output values for the one that minimizes
  // the maximum error on the interval.  Since the function is monotonic,
  // evaluating the error on the extremes of the interval suffices.
  // (This could obviously be done more efficiently, but 2^P is small.)
  double best_error = INFINITY;
  f32_union best = 0.0f;
  f32_union base = B << S; // [1.0, 2.0)
  for (f32_union cand = base; cand.i < base.i + (1UL<<S); cand.i += 1UL<<(S-op)) {
    double lerr = fabs(1.0 - double(cand.f) * left);
    double rerr = fabs(1.0 - double(cand.f) * right);
    double error = max(lerr,rerr);
    if (error < best_error) {
      best_error = error;
      best = cand;
      recip_lut_errf[idx][LERR]  = lerr;  // error val per LUT entry
      recip_lut_errf[idx][LEFT]  = in.f;
      recip_lut_errf[idx][RERR]  = rerr;
      recip_lut_errf[idx][RIGHT] = in1.f;
    }
  }

  // Return P MSBs of mantissa
  return (best.i >> (S-op)) & ((1UL<<op)-1);
}

float rsqrt(float a)
{
  f32_union in = a;

  bool sign = in.i >> (S+E);
  uint32_t exp = (in.i >> S) & ((1UL<<E)-1);
  uint32_t sig = in.i & ((1UL<<S)-1);

  if (exp == ((1UL<<E)-1) && sig == 0) {
    // inf => zero of same sign
    return copysignf(0, a);
  } else if (exp == ((1UL<<E)-1)) {
    // NaN => canonical NaN
    if (!(sig >> (S-1))) // raise invalid on sNaN
      feraiseexcept(FE_INVALID);
    return NAN;
  } else if (exp == 0 && sig == 0) {
    // zero => inf of same sign; raise divide-by-zero
    feraiseexcept(FE_DIVBYZERO);
    return copysignf(INFINITY, a);
  } else if (sign) {
    // nonzero negative => NaN; raise invalid
    feraiseexcept(FE_INVALID);
    return NAN;
  } else if (exp == 0) {
    // normalize the subnormal
    while ((sig & (1UL<<(S-1))) == 0)
      exp--, sig <<= 1;
    sig = (sig << 1) & ((1UL<<S)-1);
  }

  uint32_t out_sig = rsqrt_lut[rsqrt_lut_idx(sig, exp)] << (S-op);
  uint32_t out_exp = (3 * B + ~exp) / 2;
  f32_union res = (out_exp << S) | out_sig;
  return res.f;
}

float recip(float a)
{
  f32_union in = a;

  bool sign = in.i >> (S+E);
  uint32_t exp = (in.i >> S) & ((1UL<<E)-1);
  uint32_t sig = in.i & ((1UL<<S)-1);

  if (exp == ((1UL<<E)-1) && sig == 0) {
    // inf => zero of same sign
    return copysignf(0, a);
  } else if (exp == ((1UL<<E)-1)) {
    // NaN => canonical NaN
    if (!(sig >> (S-1))) // raise invalid on sNaN
      feraiseexcept(FE_INVALID);
    return NAN;
  } else if (exp == 0 && sig == 0) {
    // zero => inf of same sign; raise divide-by-zero
    feraiseexcept(FE_DIVBYZERO);
    return copysignf(INFINITY, a);
  } else if (exp == 0) {
    // normalize the subnormal
    while ((sig & (1UL<<(S-1))) == 0)
      exp--, sig <<= 1;
    sig = (sig << 1) & ((1UL<<S)-1);

    if (exp != 0 && exp != (uint32_t)-1) {
      // overflow to inf or max value of same sign, depending on sign and
      // rounding mode
      feraiseexcept(FE_INEXACT | FE_OVERFLOW);

      if (fegetround() == FE_TOWARDZERO ||
          (fegetround() == FE_DOWNWARD && !sign) ||
          (fegetround() == FE_UPWARD && sign))
        return copysignf(FLT_MAX, a);
      else
        return copysignf(INFINITY, a);
    }
  }

  uint32_t out_exp = 2 * B + ~exp;
  uint32_t out_sig = recip_lut[recip_lut_idx(sig)] << (S-op);

  if (out_exp == 0 || out_exp == (uint32_t)-1) {
    // the result is subnormal, but don't raise the underflow exception,
    // because there's no additional loss of precision.
    out_sig = (out_sig >> 1) | (1UL << (S-1));
    if (out_exp == (uint32_t)-1) {
      out_sig >>= 1;
      out_exp = 0;
    }
  }

  f32_union res = (out_exp << S) | out_sig;
  return res.f;
}

void populate_luts()
{
  for (size_t i = 0; i < ipN; i++) {
    rsqrt_lut[i] = estimate_rsqrt_sig(i);
    recip_lut[i] = estimate_recip_sig(i);
  }
}

void verilog()
{
  printf("module RSqrt%dx%dLUT (input [%d:0] in, output reg [%d:0] out);\n", ip, op, ip-1, op-1);
  printf("  // in[%d] corresponds to exp[0]\n", ip-1);
  printf("  // in[%d:0] corresponds to sig[S-1:S-%d]\n", ip-2, ip-2);
  printf("  // out[%d:0] corresponds to sig[S-1:S-%d]\n", op-1, op-1);
  printf("  always @(*)\n");
  printf("    case (in)\n");
  for (size_t i = 0; i < ipN; i++)
  printf("      %zu: out = %" PRIu32 ",;\n", i, rsqrt_lut[i]);
  printf("    endcase\n");
  printf("endmodule\n");

  printf("module Recip%dx%dLUT (input [%d:0] in, output reg [%d:0] out);\n", ip, op, ip-1, op-1);
  printf("  // in[%d:0]  corresponds to sig[S-1:S-%d]\n", ip-1, ip-1);
  printf("  // out[%d:0] corresponds to sig[S-1:S-%d]\n", op-1, op-1);
  printf("  always @(*)\n");
  printf("    case (in)\n");
  for (size_t i = 0; i < ipN; i++)
    printf("      %zu: out = %" PRIu32 ";\n", i, recip_lut[i]);
  printf("    endcase\n");
  printf("endmodule\n");
}

void detail(int detailn)
{

int lshift =  max(0,int(op - ip));
int rshift =  max(0,int(ip - op));

  printf("Recip%dx%dLUT (input [%d:0] in, output reg [%d:0] out);\n", ip, op, ip-1, op-1);
  printf(" in[%d:0]  corresponds to sig[S-1:S-%d]\n", ip-1, ip-1);
  printf(" out[%d:0] corresponds to sig[S-1:S-%d]\n", op-1, op-1);
  printf(" biased : ((ipN-1) - in) << (op - ip) // or >> if neg\n");
  printf(" base bias %d  left-shift %d right-shift %d\n",
           int(((ipN - 1) << lshift)>>rshift), lshift, rshift);

  for (size_t i = 0; i < ipN; i++) {
    int bias = ((((ipN - 1) - i) << lshift)>>rshift);

    printf(" %zu: out = %" PRIu32 " biased %d; lerr %G rerr %G larg %G rarg %G\n", i,
               recip_lut[i],
               int(bias - recip_lut[i]),
               recip_lut_errf[i][LERR],
               recip_lut_errf[i][RERR],
               recip_lut_errf[i][LEFT],
               recip_lut_errf[i][RIGHT]
               );
  }

  printf("\n");

  for (size_t i = 0; i < ipN; i++)

    printf(" %zu: out = %X lerr %A rerr %A larg %A rarg %A\n", i,
               recip_lut[i],
               recip_lut_errf[i][LERR],
               recip_lut_errf[i][RERR],
               recip_lut_errf[i][LEFT],
               recip_lut_errf[i][RIGHT]
               );
  printf("\n");

  printf("RSqrt%dx%dLUT (input [%d:0] in, output reg [%d:0] out);\n", ip, op, ip-1, op-1);
  printf("  // in[%d] corresponds to exp[0]\n", ip-1);
  printf("  // in[%d:0] corresponds to sig[S-1:S-%d]\n", ip-2, ip-2);
  printf("  // out[%d:0] corresponds to sig[S-1:S-%d]\n", op-1, op-1);
  printf("  // biased : ((ipN-1) - in) << (op - ip)\n");
  for (size_t i = 0; i < ipN; i++) {
    int bias = ((((ipN - 1) - i) << lshift)>>rshift);

    printf(" %zu: out %" PRIu32 " biased %d; lerr %G rerr %G larg %G rarg %G\n", i,
               rsqrt_lut[i],
               int(bias - rsqrt_lut[i]),
               rsqrt_lut_errf[i][LERR],
               rsqrt_lut_errf[i][RERR],
               rsqrt_lut_errf[i][LEFT],
               rsqrt_lut_errf[i][RIGHT]
               );
  }
  printf("\n");
  for (size_t i = 0; i < ipN; i++)
        printf(" %zu: out = %X lerr %A rerr %A larg %A rarg %A\n", i,
               rsqrt_lut[i],
               rsqrt_lut_errf[i][LERR],
               rsqrt_lut_errf[i][RERR],
               rsqrt_lut_errf[i][LEFT],
               rsqrt_lut_errf[i][RIGHT]
               );
  printf("\n");

}

void test(int testn)
{

  uint32_t first   = 0;
  uint32_t firstsq = 0;
  uint32_t last    = 0x7f7fffff;

  if (testn == 1) {          // fast test that covers LUT
       first   = 0x3F000000;
       firstsq = 0x3E800000;
       last    = 0x3F800000;
  }

  if (detailn == 1)
    for (uint32_t i = 0; i < ipN; i++) {
      recip_lut_errf[i][OTHER] = 0.0;
      recip_lut_errf[i][OERR]  = 0.0;
      rsqrt_lut_errf[i][OTHER] = 0.0;
      rsqrt_lut_errf[i][OERR]  = 0.0;
    }

  double max_error = -1.0;
  double max_errf  = -1.0;

  int idx =-1;

  for (uint32_t i = first; i <= last; i++) {
    f32_union r = i;
    double rcp = recip(r.f);
    double error = 1.0 - rcp * double(r.f);

    if (!isfinite(rcp)) {
      assert(!isfinite(1.0f / r.f));
    } else {
      if (max_error < fabs(error)) {
         max_error = fabs(error);
         max_errf = r.f;
      }
    }
  }
  printf("max recip %dx%d error at %G: %G or 2^%g\n", ip, op, max_errf, max_error, log2(max_error));


  max_error = -1.0;
  max_errf  = -1.0;

  for (uint32_t i = firstsq; i <= last; i++) {
    f32_union r = i;
    double error = 1.0 - double(rsqrt(r.f)) * sqrt(double(r.f));
    if (max_error < fabs(error)) {
       max_error = fabs(error);
       max_errf = r.f;
    }
  }
  printf("max rsqrt %dx%d error at %G: %G or  2^%g\n", ip, op, max_errf, max_error, log2(max_error));

}

int main(int argc, char** argv)
{

  // handle arguements


  if (argc > 1) {
    vector<string> args(argv + 0, argv + (argc));
    for(vector<string>::iterator it = args.begin(); it != args.end(); ++it) {
      if(*it == "--verilog")
        verilogn = 1;
      if(*it == "--test")
        testn = 1;
      if(*it == "--test-long")
        testn = 2;
      if(*it == "--nosum")
        summary = 0;
      if(*it == "--detail")
        detailn = 1;
   }
  } else {
    printf(" No arguments were provided\n");
    printf(" parameters seperated by spaces are\n");
    printf(" optional index-size(default 7 in bits), estimate size (both in bits )\n");
    printf(" these can be followed by\n --verilog to print out LUT table definition \n");
    printf(" --test to test minimal range that covers the LUT table (default) -or- \n");
    printf(" --test-long which tests all single float non-NaN values\n");
    printf(" --detail prints LUT entry values/errors\n");
    return (-1);
  }

  if (argc>1) ip= strtoul(argv[1], NULL, 10);
  if (ip>ipMax || ip <= 0) {
    if (summary > 0) printf("index of estimates, ip=%d, out of range reset to default\n",ip);
    ip = 0;
  }

  if (argc>2) op= strtoul(argv[2], NULL, 10);
  if (op>20 || op <= 0) {
    if (summary > 0) printf("estimate width, op=%d, out of range reset to default\n",op);
    op = 0;
  }

  if (ip==0) ip = 7; // default to 7x7
  if (op==0) op = ip;

  ipN = (1UL<<ip);

if (verilogn == 0 && testn==0) testn = 1;

  if (summary > 0)
    printf ("ip %d op %d LUT #bits %d verilog %d  test/test-long %d \n",
                ip, op, int((1UL <<ip) * op), verilogn, testn);

  populate_luts();

  if (detailn == 1)  detail(detailn);

  if (verilogn == 1) verilog();
  else if (testn==0) testn = 1;

  if (testn != 0) test(testn);
}
