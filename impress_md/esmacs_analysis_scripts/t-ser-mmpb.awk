#------------------------------------------------------------------------------
# time series of energies: 
# bonded
# vdw (VDWAALS + 1-4 VDW)
# elec (EEL + 1-4 EEL)
# epb
# enp (ECAVITY + EDISPER) 
BEGIN { 
} 
{  
  p1[NR]=$1; p2[NR]=$2; p3[NR]=$3; n=NR;
} 

END {
    nstep=n/3;
    for(i=1;i<=nstep;i++) {
       nn=3*(i-1);

       ebond=p1[nn+1]+p2[nn+1]+p3[nn+1];
       evdw1=p1[nn+2];
       eeel1=p2[nn+2];
       epb=p3[nn+2];
       evdw2=p1[nn+3];
       eeel2=p2[nn+3];

       evdw=evdw1+evdw2;
       eeel=eeel1+eeel2;
       printf ("%12.3f %12.3f %12.3f %12.3f\n",ebond,evdw,eeel,epb);
    }

}
