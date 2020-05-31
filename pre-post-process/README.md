## GQHR
### Format of input files
* The GritQuo input file describes the Gritsenko quotient.
```
weight

level

numJFs
JFlist (entrynum  type  weight  level  power of eta  length  numHeckeOps HeckeOps d1 ... dlength)

number of numerator terms
numertermlist (coefficient  num JFs  JF indices)

number of denominator terms
denomtermlist (coefficient  num JFs  JF indices)
```

* The HR input file describes the restriction modular curve and the numbers needed to calculation each eigenvalue.
* By theorem (see paper "On the paramodularity of typical abelian surfaces"), we have a bound on the T_1(p^2)-eigenvalue in weight k:

    |a_p| <= p^{k-3}(1+p)(1+p^2).
    
* For simplicity, denote the righthand side as boundTp = p^{k-3}(1+p)(1+p^2).
* For each p where you wish to compute T(p), find a prime _modn_ and a number _root_ such that:

    modn > 2 * boundTp
    root^p = 1 mod modn (and root != 1)
    
```
weight
level

a b c (restriction matrix)
uptoDot

number of Hecke operators
p  modn  root
p  modn  root
...
```
* Remark: A requirement is that L^2 | (modn-1).  A useful Mathematica command is PrimitiveRoot.

### How to post-process output file
* Set the Mathematica function:
```
pList = {};Clear[myInfo,q];
domeSomething := (
   myInfo[L] = {modn, origqs, totqs};
   AppendTo[pList, L];
   );
```
* Then read in the output file.
* Then compute the eigenvalues.
```
resList = {}; Clear[TpvalHash];
For[i = 1, i <= Length[pList], i++,
 pp = pList[[i]];
 {nnn, origf, Tpf} = myInfo[pp];
 Tpval = Mod[Coefficient[Tpf, q^(targetDot*pp)]*
    PowerMod[Coefficient[origf, q^targetDot], -1, nnn],
   nnn, -nnn/2
   ];
 TpvalHash[pp] = Tpval;
 AppendTo[resList, {pp, Tpval, nnn}];
 ]
```
* The eigenvalues are now in the list resList and also in the hash TpvalHash

## GQHRTp2
### Format of GritQuo input file
* The GritQuo input file is the same as that used for GQHR.
* The HR input file is slightly different than that used for GQHR.
* By theorem (see paper "On the paramodularity of typical abelian surfaces"), we have a bound on the T_1(p^2)-eigenvalue in weight k:

    |a_{1,p^2}| <= p^{2k-6}(1+p)(1+p^2)p.
    
* For simplicity, denote the righthand side as boundTp2 = p^{2k-6}(1+p)(1+p^2)p.
* For each p where you wish to compute T_1(p^2), find a prime _modn_ and a number _root_ such that:

    modn > 2 * boundTp2
    
    modn > 2 * p^2 * boundTp2 _if weight k=2_
    
    root^(p^2) = 1 mod modn
    
    root^p != 1 mod modn
    
* __IMPORTANT__: Note the increase in size of boundTp2 when the weight is 2.
* The format of the "HR" is as follows

```
weight
level
a b c
d
numHeckeToBeDone
p modn root root^p
p modn root root^p
...
```

* Example:
```
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
```

* Remark: A requirement is that L^2 | (modn-1).  A useful Mathematica command is PrimitiveRoot.

### How to post-process output file
* Same as post-processing output file of GQHR _except when weight is 2_.  So do the following:
```
resListTp2 = {}; Clear[Tp2valHash];
For[i = 1, i <= Length[pList2], i++,
 pp = pList2[[i]];
 {nnn, origf, Tpf} = myInfoTp2[pp];
 If[weight == 2,
  Tpval = Mod[pp^2*Coefficient[Tpf, q^(targetDot*pp*pp)]*
      PowerMod[Coefficient[origf, q^targetDot], -1, nnn],
     nnn, -nnn/2
     ]/pp^2,
  Tpval = Mod[Coefficient[Tpf, q^(targetDot*pp*pp)]*
     PowerMod[Coefficient[origf, q^targetDot], -1, nnn],
    nnn, -nnn/2
    ]
  ];
 Tp2valHash[pp] = Tpval;
 AppendTo[resListTp2, {pp, Tpval, nnn}]
 ];
 ```

## BPHR
* To be written

## BPHRTp2
* To be written

## TwHR
* To be written
