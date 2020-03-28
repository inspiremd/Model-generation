#------------------------------------------------------------------------------
# Average and rms of deviation
BEGIN {
}
{
   v[NR]=$1; n=NR;
}
END {
   sum=0.0;
   sum2=0.0;
   flu=0.0;
   for(i=1;i<=n;i++) {
      sum+=v[i];
      sum2+=v[i]*v[i];
   }
   sum=sum/n;
   sum2=sum2/n;
   flu=sum2-sum*sum;
   if (flu < 0) flu=0;
   flu=sqrt(flu);
   # print out
#   printf("%s%12.5f  %12.5f\n","average and fluctuation:  ",sum,flu);
#   print TITLE
   printf("%12.5f  %12.5f\n",sum,flu);
}

