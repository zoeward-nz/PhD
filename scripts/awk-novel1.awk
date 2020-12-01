# only consider tx that are of gffcompare classes i,x,s,u,y
# http://ccb.jhu.edu/software/stringtie/gffcompare.shtml#transfrag-class-codes
{if ($3=="i" || $3=="x" || $3=="s" || $3=="u" || $3=="y") print "\""$1"\"";}
