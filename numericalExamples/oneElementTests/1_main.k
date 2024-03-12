$# LS-DYNA Keyword file created by LS-PrePost(R) V4.7.20 -24Aug2020
$# Created on May-6-2022 (08:53:11)
*KEYWORD
*TITLE
$#                                                                         title
LS-DYNA keyword deck by LS-PrePost
*CONTROL_IMPLICIT_GENERAL
$#  imflag       dt0    imform      nsbs       igs     cnstn      form    zero_v
         1       1.0         2         1         1         0         0         0
*CONTROL_IMPLICIT_SOLUTION
$#  nsolvr    ilimit    maxref     dctol     ectol     rctol     lstol    abstol
        12         1        10      0.011.0000E-201.0000E-20       0.9-1.0000E-4
$#   dnorm    diverg     istif   nlprint    nlnorm   d3itctl     cpchk     
         1         2         1         3         2         0         0
$#  arcctl    arcdir    arclen    arcmth    arcdmp    arcpsi    arcalf    arctim
         0         0       0.0         1         2         0         0         0
$#   lsmtd     lsdir      irad      srad      awgt      sred    
         2         2       0.0       0.0       0.0       0.0
*CONTROL_TERMINATION
$#  endtim    endcyc     dtmin    endeng    endmas     nosol     
      10.0         0       1.0       0.01.000000E8         0
*DATABASE_ELOUT
$#      dt    binary      lcur     ioopt   option1   option2   option3   option4
       0.0         0         0         1         0         0         0         0
*DATABASE_GLSTAT
$#      dt    binary      lcur     ioopt     
       0.0         0         0         1
*DATABASE_MATSUM
$#      dt    binary      lcur     ioopt     
       0.0         0         0         1
*DATABASE_NCFORC
$#      dt    binary      lcur     ioopt     
       0.0         0         0         1
*DATABASE_NODFOR
$#      dt    binary      lcur     ioopt     
       1.0         0         0         1
*DATABASE_NODOUT
$#      dt    binary      lcur     ioopt   option1   option2       
       1.0         0         0         1       0.0         0
*DATABASE_RCFORC
$#      dt    binary      lcur     ioopt     
       0.0         0         0         1
*DATABASE_SECFORC
$#      dt    binary      lcur     ioopt     
       0.0         0         0         1
*DATABASE_BINARY_D3PLOT
$#      dt      lcdt      beam     npltc    psetid      
       1.0         0         0         0         0
$#   ioopt      rate    cutoff    window      type      pset    
         0       0.0       0.0       0.0         0         0
*DATABASE_EXTENT_BINARY
$#   neiph     neips    maxint    strflg    sigflg    epsflg    rltflg    engflg
         6         0         3         0         1         1         1         1
$#  cmpflg    ieverp    beamip     dcomp      shge     stssz    n3thdt   ialemat
         0         0         0         1         1         1         2         1
$# nintsld   pkp_sen      sclp     hydro     msscl     therm    intout    nodout
         0         0       1.0         0         0         0           
$#    dtdt    resplt     neipb     quadr     cubic     
         0         0         0         0         0
*BOUNDARY_PRESCRIBED_MOTION_SET
$#    nsid       dof       vad      lcid        sf       vid     death     birth
         1         2         2         1       1.0         01.00000E28       0.0
*DEFINE_CURVE
$#    lcid      sidr       sfa       sfo      offa      offo    dattyp     lcint
         1         0       1.0       1.0       0.0       0.0         0         0
$#                a1                  o1  
                 0.0                 0.0
                10.0                0.01
*PART
$#                                                                         title
boxsolid
$#     pid     secid       mid     eosid      hgid      grav    adpopt      tmid
         1         1         1         0         0         0         0         0
$
*INCLUDE
2_OET-mesh.inc
$3_BC-tension.inc
3_BC-tensionShear.inc
$
*INCLUDE
3_control-UEL.inc
$3_control-ELFORM2_MAT24.inc
*END
