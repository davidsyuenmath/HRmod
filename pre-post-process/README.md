## GQHR
### Format of input files
* The GritQuo input file describes the Gritsenko quotient.
* The HR input file describes the restriction modular curve and the numbers needed to calculation each eigenvalue.
* To be written

### How to post-process output file
* To be written

## GQHRTp2
### Format of GritQuo input file
* The GritQuo input file is the same as that used for GQHR.
* The HR input file is slightly different than that used for GQHR.
* By theorem (see paper "On the paramodularity of typical abelian surfaces"), we have a bound on the T_1(p^2)-eigenvalue in weight k:

    |a_{1,p^2}| <= p^{2k-6}(1+p)(1+p^2)p.
    
* For simplicity, denote the righthand side as boundTp2 = p^{2k-6}(1+p)(1+p^2)p.
* Find a prime _modn_ and a number _root_ such that:

    modn > 2 * boundTp2
    
    modn > 2 * p^2 * boundTp2 _if weight k=2_
    
    root^(p^2) = 1 mod modn
    
    root^p != 1 mod modn
    
* __IMPORTANT__: Note the increase in size of boundTp2 when the weight is 2.
* The format of the "HR" is as follows

>weight
>level
>a b c
>d
>numHeckeToBeDone
>p modn root root^p
>p modn root root^p
>...

* Example:
3
61
122 11 1
2
3
2 61 11 60
3 271 258 242
5 1601 1419 442
-----(the following is included at the end of the file)-------
weight
level
Restriction matrix
Desired up to dot product with restriction matrix
numHecke
list of  heckeL  modn  L2root Lroot
where L2root^(heckeL * heckeL) modn = 1
and L2root^heckeL = Lroot.
This redundancy is so that you can be sure that Lroot != 1.


### How to post-process output file
* Same as post-processing output file of GQHR _except when weight is 2_
* When weight is 2, you must do the following. To be written.
