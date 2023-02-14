## A family of programs that compute eigenvalues by restriction to a modular curve.

Six programs. Hecke Restriction programs (BPHR, GQHR, BPTp2HR, GQTp2HR, TwHR, GritMonHR) using mod arithmetic, and helpers MakeWHJF, FindBPEcombo. 

* Implements speed up theorems.  
* T2 versions of GQ and BP, but not Tw. 
* BPHR and BPTp2HR: do Hecke restriction of Borcherds products.  
* GQHR and GQTp2HR: do Hecke restriction of Gritsenko Quotients.  
* TwHR: do Hecke restriction of a quotient where the numerator can include T2(Grit*Grit) and (Tweak*Tweak). 
* GritMonHR: do Hecke restriction of a quotient where the numerator and denominator are monomials, potentially of high degree.
* GQHRbadprime: do Hecke restriction of a Gritsenko Quotient at a bad prime (only T(p)).
* Also helper program MakeWHJF that makes and saves a weakly holomorphic Jacobi form that can be read by BPHR and BPTp2HR.
* Helper FindBPEcombo will try to express a weight 0 form as a combo of BPE weight 0 forms.  Use Mathematica program to make input of the format that MakeWHJF outputs. This program is used to find GritMon representations.

## Requirements

* flint, and implicitly anything that flint requires.
