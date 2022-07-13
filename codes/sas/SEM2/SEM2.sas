/*******************************************************************/ 
/* this program conducts Monte Carlo simulation of SEM, and outputs*/ 
/* parameter estimates and fit indices to a SAS system file named */ 
/* ’SEM_FITS’. 2 models are simulated in this program: true, */
/* & misspecified models. */
/*******************************************************************/

*OPTIONS MLOGIC MLOGICNEST MPRINT MPRINTNEST SYMBOLGEN MERROR SERROR;
*OPTIONS NONOTES;

%LET MC= 500; * # of Monte Carlo replications for each cell condition; 

libname SEM 'C:\Users\dchoi2\Box\EDPS_971\Simulation studies\montecarlo\SEM2';

* -- to direct the SAS log to a disk file; 
PROC PRINTTO LOG='C:\Users\dchoi2\Box\EDPS_971\Simulation studies\montecarlo\SEM2\LOGFILE.TMP';

DATA A (TYPE=CORR); 
_TYPE_='CORR';
INPUT _NAME_$ X1 X2 Y1 Y2 Y3 Y4; 
CARDS;
X1  1.0000000  0.5368755 -0.3791639 -0.3839751 -0.4586165 -0.4715445
X2  0.5368755  1.0000000 -0.2908055 -0.2944954 -0.3517428 -0.3616581
Y1 -0.3791639 -0.2908055  1.0000000  0.6206536  0.5579119  0.4849042
Y2 -0.3839751 -0.2944954  0.6206536  1.0000000  0.4775941  0.6172980
Y3 -0.4586165 -0.3517428  0.5579119  0.4775941  1.0000000  0.7139541
Y4 -0.4715445 -0.3616581  0.4849042  0.6172980  0.7139541  1.0000000
;

/* 
*SAS code;
X1 1.000000 . . . . . 
X2 .536875 1.000000 . . . . 
Y1 -.379164 -.290805 1.000000 . . . 
Y2 -.414038 -.317552 .669246 1.000000 . . 
Y3 -.375884 -.288290 .569164 .482667 1.000000 .
Y4 -.388047 -.297619 .456315 .529719 .672364 1.000000

*SEM2;
X1  1.0000000  0.5368755 -0.3791639 -0.3839751 -0.4586165 -0.4715445
X2  0.5368755  1.0000000 -0.2908055 -0.2944954 -0.3517428 -0.3616581
Y1 -0.3791639 -0.2908055  1.0000000  0.6206536  0.5579119  0.4849042
Y2 -0.3839751 -0.2944954  0.6206536  1.0000000  0.4775941  0.6172980
Y3 -0.4586165 -0.3517428  0.5579119  0.4775941  1.0000000  0.7139541
Y4 -0.4715445 -0.3616581  0.4849042  0.6172980  0.7139541  1.0000000
*/

* obtain factor pattern matrix for later data generation;
PROC FACTOR N=6 OUTSTAT=FACOUT NOPRINT; RUN; 
DATA PATTERN; SET FACOUT; 
	IF _TYPE_='PATTERN'; 
	DROP _TYPE_ _NAME_;
RUN;

%MACRO SEM_MC; * start of monte carlo simulation macro ’SEM_MC’; 
* do-loop for 3 conditions for true and misspecified models;
%DO MOD = 1 %TO 3;
	%IF &MOD=1 %THEN %DO; %LET MODEL=TRU; %END;
	%IF &MOD=2 %THEN %DO; %LET MODEL=LV1; %END; 
	%IF &MOD=3 %THEN %DO; %LET MODEL=LV2; %END;

* do-loop for 2 estimation procedures; 
%DO A = 1 %TO 1;
	%IF &A=1 %THEN %DO; %LET METHOD=MAX; %END; 
	%IF &A=2 %THEN %DO; %LET METHOD=GLS; %END;


* do-loop for 4 sample size conditions; 
%DO B = 1 %TO 4;
	%IF &B=1 %THEN %DO; %LET SMPLN=100; %END; 
	%IF &B=2 %THEN %DO; %LET SMPLN=200; %END; 
	%IF &B=3 %THEN %DO; %LET SMPLN=500; %END; 
	%IF &B=4 %THEN %DO; %LET SMPLN=1000; %END;

* do-loop for the number of replications in each cell;
%DO C=1 %TO &MC;


* use SAS PROC IML for data generation;
* use the factor pattern matrix; 
* diagonal matrix containing variances for 6 variables;

PROC IML;
	USE PATTERN;
	READ ALL VAR _NUM_ INTO F; 
	F=F`;
									
*SAS CODE
VAR={10	0	0	0	0	0,
	0	4.250	0	0	0	0,
	0	0	12.27	0	0	0,
	0	0	0	9.2868	0	0,
	0	0	0	0	12.9047	0,
	0	0	0	0	0	9.807807};

*SEM2;
VAR={10	0	0	0	0	0,
	0	4.25000	0	0	0	0,
	0	0	12.2700	0	0	0,
	0	0	0	6.730000	0	0,
	0	0	0	0	14.724175	0,
	0	0	0	0	0	11.281582};


STD=SQRT(VAR);				* matrix containing stds for the 6 variables;
X=RANNOR(J(&SMPLN,6,&C));	* generate 6 random normal variables; 
XT=X`;						* transpose the data matrix for multiplication;
XTCORR=F*XT;
							* transform uncorrelated variables to correlated ones;
							* transform the scale of the variables (from std=1 to std=specified above);
XTSTD=STD*XTCORR;							* transpose the data matrix back;
XY=XTSTD`; 

* create SAS data set ’DAT’;
CREATE DAT FROM XY[COLNAME={X1 X2 Y1 Y2 Y3 Y4}]; 
APPEND FROM XY;
QUIT;

* implement the true model; 
* output the model fitting results to data set ’SEMOUT’;

%IF &MOD=1 %THEN %DO;

PROC CALIS DATA=DAT METHOD=&METHOD COV OUTFIT=SEMOUT NOPRINT; *NOPRINT;
	LINEQS
		X1 = 		   FK 				+	EX1,
		X2 = LX2(.50)  FK				+	EX2,
		Y1 =		   FE1 				+	EY1,
		Y2 = LY2(.75)  FE1 				+	EY2, 
		Y3 =		 				FE2 +	EY3,
		Y4 = 				LY4(.90)FE2 +	EY4, 
		FE1= GA1(-.6)  FK				+	DE1,
		FE2= GA2(-.435)FK+ 	BE1(.60)FE1 +	DE2; 

STD
	FK EX1 EX2 EY1 EY2 EY3 EY4 DE1 DE2 = 
	VFK(7) VEX1(5) VEX2(4) VEY1(4.75) VEY2(2.5) VEY3(4.5) VEY4(3) VDE1(3) VDE2(2.5);
COV 
	EY1 EY3 = C_EY13(1.16), EY2 EY4 = C_EY24(1.1); 
RUN;

DATA SEMOUT; SET SEMOUT; 			* keep only the relevant results; 
KEEP FitIndex FitValue;
RUN; 

PROC TRANSPOSE DATA=SEMOUT OUT=SEMOUT; ID FitIndex ; RUN;			* transpose the SAS data set ’SEMOUT’ to ’NEWFITS’;


	DATA NEWFITS; SET SEMOUT;
	MTHD		= "&METHOD";
	MODEL		= "&MODEL";
	N			= Number_of_Observations;
	NVAR		= Number_of_Variables;
	NPARM		= Number_of_Parameters;
	FIT			= Fit_Function;
	DF		 	= Chi_Square_DF;
	CHISQUAR	= Chi_Square;
	P_CHISQ		= Pr___Chi_Square;
	CHISQNUL	= Baseline_Model_Chi_Square;
	RMR			= Root_Mean_Square_Residual__RMR_;
	SRMR		= Standardized_RMR__SRMR_;
	GFI			= Goodness_of_Fit_Index__GFI_;
	AGFI		= Adjusted_GFI__AGFI_;
	RMSEAEST	= RMSEA_Estimate;
	CENTRALI	= McDonald_Centrality;
	COMPFITI	= Bentler_Comparative_Fit_Index;
	BB_NONOR	= Bentler_Bonett_Non_normed_Index;
	BB_NORMD	= Bentler_Bonett_NFI;
	BOL_RHO1	= Bollen_Normed_Index_Rho1;
	BOL_DEL2 	= Bollen_Non_normed_Index_Delta2;

	KEEP MTHD MODEL N NVAR NPARM FIT DF CHISQUAR P_CHISQ CHISQNUL RMR SRMR GFI AGFI 
		RMSEAEST COMPFITI BB_NONOR BB_NORMD BOL_RHO1 BOL_DEL2 CENTRALI;		* keep only the relevant results; 

PROC APPEND BASE=SEM.SEM_FITS FORCE;
RUN;

%END;							
* end of implementing the true model; 
* implement misspecified model 1;
			
%IF &MOD=2 %THEN %DO;

PROC CALIS DATA=DAT METHOD=&METHOD COV OUTFIT=SEMOUT NOPRINT; *NOPRINT;
	LINEQS
		X1 = 		   FK 				+	EX1,
		X2 = LX2(.50)  FK				+	EX2,
		Y1 =		   FE1 				+	EY1,
		Y2 = LY2(.75)  FE1 				+	EY2, 
		Y3 =		 				FE2 +	EY3,
		Y4 = 				LY4(.90)FE2 +	EY4, 
		FE1= GA1(-.6)  FK				+	DE1,
		FE2= 0		   FK+ 	BE1(.60)FE1 +	DE2; 	* misspecification: GA2 fixed;

STD
	FK EX1 EX2 EY1 EY2 EY3 EY4 DE1 DE2 = 
	VFK(7) VEX1(5) VEX2(4) VEY1(4.75) VEY2(2.5) VEY3(4.5) VEY4(3) VDE1(3) VDE2(2.5);
COV 
	EY1 EY3 = C_EY13(1.16), EY2 EY4 = C_EY24(1.1); 
RUN;

DATA SEMOUT; SET SEMOUT; 			* keep only the relevant results; 
KEEP FitIndex FitValue;
RUN; 

PROC TRANSPOSE DATA=SEMOUT OUT=SEMOUT; ID FitIndex ; run;

	DATA NEWFITS; SET SEMOUT;
	MTHD		= "&METHOD";
	MODEL		= "&MODEL";
	N			= Number_of_Observations;
	NVAR		= Number_of_Variables;
	NPARM		= Number_of_Parameters;
	FIT			= Fit_Function;
	DF		 	= Chi_Square_DF;
	CHISQUAR	= Chi_Square;
	P_CHISQ		= Pr___Chi_Square;
	CHISQNUL	= Baseline_Model_Chi_Square;
	RMR			= Root_Mean_Square_Residual__RMR_;
	SRMR		= Standardized_RMR__SRMR_;
	GFI			= Goodness_of_Fit_Index__GFI_;
	AGFI		= Adjusted_GFI__AGFI_;
	RMSEAEST	= RMSEA_Estimate;
	CENTRALI	= McDonald_Centrality;
	COMPFITI	= Bentler_Comparative_Fit_Index;
	BB_NONOR	= Bentler_Bonett_Non_normed_Index;
	BB_NORMD	= Bentler_Bonett_NFI;
	BOL_RHO1	= Bollen_Normed_Index_Rho1;
	BOL_DEL2 	= Bollen_Non_normed_Index_Delta2;

	KEEP MTHD MODEL N NVAR NPARM FIT DF CHISQUAR P_CHISQ CHISQNUL RMR SRMR GFI AGFI 
		RMSEAEST COMPFITI BB_NONOR BB_NORMD BOL_RHO1 BOL_DEL2 CENTRALI;		* keep only the relevant results; 

PROC APPEND BASE=SEM.SEM_FITS FORCE;
RUN;

%END;
* end of implementing the misspecified model 1; 
* implement misspecified model 2;
	
%IF &MOD=3 %THEN %DO;

PROC CALIS DATA=DAT METHOD=&METHOD COV OUTFIT=SEMOUT NOPRINT; *NOPRINT;
	LINEQS
		X1 = 		   FK 				+	EX1,
		X2 = LX2	   FK				+	EX2,
		Y1 =		   FE1 				+	EY1,
		Y2 = LX2	   FE1 				+	EY2, 
		Y3 =		 				FE2 +	EY3,
		Y4 = 				LY4(.90)FE2 +	EY4, 
		FE1= GA1(-.6)  FK				+	DE1,
		FE2= 0		   FK+ 	BE1(.60)FE1 +	DE2; 	* misspecification: GA2 fixed;
													* constrained: LY2 = LX2;
STD
	FK EX1 EX2 EY1 EY2 EY3 EY4 DE1 DE2 = 
	VFK(7) VEX1(5) VEX2(4) VEY1(4.75) VEY2(2.5) VEY3(4.5) VEY4(3) VDE1(3) VDE2(2.5);
*COV 
	EY1 EY3 = C_EY13(1.16), EY2 EY4 = C_EY24(1.1); * misspecification: error covariances fixed to zeros;
RUN;

DATA SEMOUT; SET SEMOUT; 			* keep only the relevant results; 
KEEP FitIndex FitValue;
RUN; 

PROC TRANSPOSE DATA=SEMOUT OUT=SEMOUT; ID FitIndex ; run;

	DATA NEWFITS; SET SEMOUT;
	MTHD		= "&METHOD";
	MODEL		= "&MODEL";
	N			= Number_of_Observations;
	NVAR		= Number_of_Variables;
	NPARM		= Number_of_Parameters;
	FIT			= Fit_Function;
	DF		 	= Chi_Square_DF;
	CHISQUAR	= Chi_Square;
	P_CHISQ		= Pr___Chi_Square;
	CHISQNUL	= Baseline_Model_Chi_Square;
	RMR			= Root_Mean_Square_Residual__RMR_;
	SRMR		= Standardized_RMR__SRMR_;
	GFI			= Goodness_of_Fit_Index__GFI_;
	AGFI		= Adjusted_GFI__AGFI_;
	RMSEAEST	= RMSEA_Estimate;
	CENTRALI	= McDonald_Centrality;
	COMPFITI	= Bentler_Comparative_Fit_Index;
	BB_NONOR	= Bentler_Bonett_Non_normed_Index;
	BB_NORMD	= Bentler_Bonett_NFI;
	BOL_RHO1	= Bollen_Normed_Index_Rho1;
	BOL_DEL2 	= Bollen_Non_normed_Index_Delta2;

	KEEP MTHD MODEL N NVAR NPARM FIT DF CHISQUAR P_CHISQ CHISQNUL RMR SRMR GFI AGFI 
		RMSEAEST COMPFITI BB_NONOR BB_NORMD BOL_RHO1 BOL_DEL2 CENTRALI;		* keep only the relevant results; 

PROC APPEND BASE=SEM.SEM_FITS FORCE;
RUN;

	
%END;		
* end of implementing the misspecified model 2; 
	
								* close the do-loop for replications in each cell;
%END;
								* close the do-loop for sample size conditions;
%END;
							* close the do-loop for estimation procedure;
%END;
							* close the do-loop for model specification conditions;
%END;

%MEND SEM_MC; 					* end of simulation macro ’SEM_MC’; 
%SEM_MC;
								* running the macro ’SEM_MC’;
RUN;

PROC PRINTTO 
PRINT=PRINT; * direct output to SAS Output window; 

RUN;


* descriptive analysis for the 9 descriptive model fit indices;
*libname SEM 'C:\Users\dchoi2\Box\EDPS_971\Simulation studies\montecarlo\SEM2';

DATA D2; SET SEM.SEM_FITS; 
PROC SORT data=d2; BY MODEL MTHD N; 
PROC MEANS N MEAN STD MAX MIN; 
	BY MODEL MTHD N;
	VAR FIT DF NVAR CHISQUAR RMR SRMR GFI AGFI RMSEAEST COMPFITI BB_NONOR BB_NORMD BOL_RHO1
		BOL_DEL2 CENTRALI NPARM; 
RUN; 

/*
PROC IML;
* define the 8 population matrices;
LX = {1.00, 0.50}; 
LY = {1.00 0.00, 0.95 0.00,0.00 1.00, 0.00 0.90}; 
GA = {-0.60, -0.25}; 
PH = {7.00};
PS = {5.00 0.00, 0.00 4.00}; 
TD = {3.00 0.00, 0.00 2.50};
TE = {4.75 0.00 1.60 0.00, 0.00 2.50 0.00 0.30, 1.60 0.00 4.50 0.00, 0.00 0.30 0.00 3.00}; 
B = {0.00 0.00, 0.60 0.00}; 
I = {1 0, 0 1};

* obtain the four quadrants in Equation 2; 
COVY = LY*(INV(I-B))*(GA*PH*GA`+PS)*(INV(I-B`))*LY`+TE; 
COVX = LX*PH*LX` + TD;
COVYX = LY*(INV(I-B))*GA*PH*LX`; 
COVXY = LX*PH*GA`*(INV(I-B`))*LY`; 

* put 4 quadrants together as the population covariance matrix; 
COV =  (COVYX || COVY)// (COVX || COVXY); 
PROC PRINT COV; 
RUN;
QUIT;
*/

