/* define SAS library */
libname ANOVA 'C:\Users\dchoi2\Box\EDPS_971\Simulation studies\ANOVA';

/* Importing data from a .csv file */
/*
proc import datafile = 'C:\Users\dchoi2\Box\EDPS_971\Simulation studies\ANOVA\final_R.csv'
 dbms = csv
 out = a
 replace;
 getnames = yes;
 run;
*/

proc import datafile = 'C:\Users\dchoi2\Box\EDPS_971\Simulation studies\ANOVA\final_SAS.csv'
 dbms = csv
 out = a
 replace;
 getnames = yes;
 run;


/* Conduct ANOVA and compute measures of effect size */ 
PROC GLM DATA=a; 
CLASS mtype msize mlevel nsize; 
MODEL chisq = mtype|msize|mlevel|nsize /effectsize; 
run;

PROC GLM DATA=a; 
CLASS mtype msize mlevel nsize; 
MODEL rmr = mtype|msize|mlevel|nsize /effectsize; 
run;

PROC GLM DATA=a; 
CLASS mtype msize mlevel nsize; 
MODEL srmr = mtype|msize|mlevel|nsize /effectsize; 
run;

PROC GLM DATA=a; 
CLASS mtype msize mlevel nsize; 
MODEL gfi = mtype|msize|mlevel|nsize /effectsize; 
run;

PROC GLM DATA=a; 
CLASS mtype msize mlevel nsize; 
MODEL agfi = mtype|msize|mlevel|nsize /effectsize; 
run;

PROC GLM DATA=a; 
CLASS mtype msize mlevel nsize; 
MODEL rmsea = mtype|msize|mlevel|nsize /effectsize; 
run;

PROC GLM DATA=a; 
CLASS mtype msize mlevel nsize; 
MODEL cfi = mtype|msize|mlevel|nsize /effectsize; 
run;

PROC GLM DATA=a; 
CLASS mtype msize mlevel nsize; 
MODEL nfi = mtype|msize|mlevel|nsize /effectsize; 
run;

PROC GLM DATA=a; 
CLASS mtype msize mlevel nsize; 
MODEL nnfi = mtype|msize|mlevel|nsize /effectsize; 
run;

PROC GLM DATA=a; 
CLASS mtype msize mlevel nsize; 
MODEL gamma = mtype|msize|mlevel|nsize /effectsize; 
run;
quit;

*chisq	rmr	srmr	gfi	agfi	rmsea	cfi	nnfi	nfi	gamma


/* descriptive analysis for the 9 descriptive model fit indices */

DATA D2; SET SEM.SEM_FITS; 
PROC SORT data=d2; BY MODEL MTHD N; 
PROC MEANS N MEAN STD MAX MIN; 
	BY MODEL MTHD N;
	VAR FIT DF NVAR CHISQUAR RMR SRMR GFI AGFI RMSEAEST COMPFITI BB_NONOR BB_NORMD BOL_RHO1
		BOL_DEL2 CENTRALI NPARM; 
RUN; 





proc glm outstat=summary;
   class Rep Current Time Number;
   model MuscleWeight = Rep Current|Time|Number;
   contrast 'Time in Current 3'
      Time 1 0 0 -1 Current*Time 0 0 0 0 0 0 0 0 1 0 0 -1,
      Time 0 1 0 -1 Current*Time 0 0 0 0 0 0 0 0 0 1 0 -1,
      Time 0 0 1 -1 Current*Time 0 0 0 0 0 0 0 0 0 0 1 -1;
   contrast 'Current 1 versus 2' Current 1 -1;
   lsmeans Current*Time / slice=Current;
run;
