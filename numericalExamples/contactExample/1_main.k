$# LS-DYNA Keyword file created by LS-PrePost(R) V4.7.20 -24Aug2020
$# Created on Apr-29-2022 (12:04:40)
*KEYWORD
$
$ Based on [https://www.dynaexamples.com/implicit/basic-examples/contact-ii]
$
$ Units: m, s, kg, N
$
*INCLUDE
3_control.inc
$
*DATABASE_NODFOR
$#      dt    binary      lcur     ioopt     
&dt                3         0         1
*DATABASE_NODAL_FORCE_GROUP
$#    nsid       cid   
         1         0
*DATABASE_EXTENT_BINARY
         6         0         3         0         1         1         1         1
         0         0         0         1         1         1         2         1
         0         0       1.0         0         0         0                    
         0         0         0         0         0
*PART
$#                                                                         title
elastic ball
$#     pid     secid       mid     eosid      hgid      grav    adpopt      tmid
         1         1         1         0         0         0         0         0
*PART
$#                                                                         title
elastic plate
$#     pid     secid       mid     eosid      hgid      grav    adpopt      tmid
         2         1         1         0         0         0         0         0
*INCLUDE
3_control-UEL.inc
$3_control-ELFORM2_MAT24.inc
*END
