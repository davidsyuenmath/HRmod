## A family of programs that compute eigenvalues by restriction to a modular curve.

Seven programs. Hecke Restriction programs (BPHR, GQHR, BPTp2HR, GQTp2HR, TwHR, GritMonHR, GQHRbadprime) using mod arithmetic, and helpers MakeWHJF, FindBPEcombo. 

* Implements speed up theorems.  
* Tp2 versions of GQ and BP, but not Tw. 
* BPHR and BPTp2HR: do Hecke restriction of Borcherds products.  
* GQHR and GQTp2HR: do Hecke restriction of Gritsenko Quotients.  
* TwHR: do Hecke restriction of a quotient where the numerator can include T2(Grit*Grit) and (Tweak*Tweak). 
* GritMonHR: do Hecke restriction of a quotient where the numerator and denominator are monomials, potentially of high degree.
* GQHRbadprime: do Hecke restriction of a Gritsenko Quotient at a bad prime (only T(p)) that simply divides the level.
* Also helper program MakeWHJF that makes and saves a weakly holomorphic Jacobi form that can be read by BPHR and BPTp2HR.
* Helper FindBPEcombo will try to express a weight 0 form as a combo of BPE weight 0 forms.  Use Mathematica program to make input of the format that MakeWHJF outputs. This program is used to find GritMon representations.

## Requirements

* flint, and implicitly anything that flint requires.
* The versions that worked with flint2.5 are in the branch flint2.5-version.
* The current (2023-06-01) version of flint is flint3.0 and is **not backwards compatible** with the version that was compiled with flint2.5. 
We are currently working to rewrite all the programs to be compatible with flint3.0.  Currently only GQHR and GQTp2HR have been updated to be compatible with flint3.0

## Notes

* These programs were used in papers of Cris Poor, Jerry Shurman, David S. Yuen
  - *On the paramodularity of typical abelian surfaces, Journal of Algebra and Number Theory*, 13 (2019) 1145-1195. (Armand Brumer, Ariel Pacetti, Cris Poor, Gonzalo Tornaria, John Voight, David S. Yuen)
  - *Eigenvalues and congruences for the weight 3 paramodular nonlifts of levels 61, 73, and 79*, preprint (Cris Poor, Jerry Shurman, David S. Yuen)
  - Future papers
