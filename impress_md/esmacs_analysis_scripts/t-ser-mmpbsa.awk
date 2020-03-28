#------------------------------------------------------------------------------
# time series of energies: 
# bonded
# vdw (VDWAALS + 1-4 VDW)
# elec (EEL + 1-4 EEL)
# epb
# enp from esurf (not ECAVITY + EDISPER) 
BEGIN {
} 
{  
  p1[NR]=$1; p2[NR]=$2; p3[NR]=$3; p4[NR]=$4; p5[NR]=$5; n=NR;
} 

END {
    for(i=1;i<=n;i++) {
       enp=p5[i]/72*54.2+0.92;
       etot=p1[i]+p2[i]+p3[i]+p4[i]+enp;
       printf ("%12.3f %12.3f %12.3f %12.3f %12.3f %12.3f\n",p1[i],p2[i],p3[i],p4[i],enp,etot);
    }

}
