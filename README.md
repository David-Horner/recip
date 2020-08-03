# recip
Reciprical float and rsqrt estimation docs and code / discussion

recip.cc - Andrew's code for generating and testing hard coded n-by-n Look Up Table (LUT)approximations.
 
revised version accepts arguments for LUT size :

    parameters seperated by spaces are
    optional index-size(default 7 in bits), estimate size (both in bits 
     these can be followed by\n --verilog to print out LUT table definition 
     --test to test minimal range that covers the LUT table (default) -or-
     --test-long which tests all single float non-NaN values.
     
Also added is a --nosum directive that suppresses print of informitive and summary parameter info.

and --detail prints LUT entries with details on float values and errors per entry
