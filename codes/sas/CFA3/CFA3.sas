/*******************************************************************/ 
/* this program conducts Monte Carlo simulation of SEM, and outputs*/ 
/* parameter estimates and fit indices to a SAS system file named */ 
/* ’SEM_FITS’. 2 models are simulated in this program: true, */
/* & misspecified models. */
/*******************************************************************/

*OPTIONS MLOGIC MLOGICNEST MPRINT MPRINTNEST SYMBOLGEN MERROR SERROR;
*OPTIONS NONOTES;

%LET MC= 500; * # of Monte Carlo replications for each cell condition; 

libname SEM 'C:\Users\dchoi2\Box\EDPS_971\Simulation studies\montecarlo\CFA3';

* -- to direct the SAS log to a disk file; 
PROC PRINTTO LOG='C:\Users\dchoi2\Box\EDPS_971\Simulation studies\montecarlo\CFA3\LOGFILE.TMP';

DATA A (TYPE=CORR); 
_TYPE_='CORR';
INPUT _NAME_$ X1 X2 X3 X4;
CARDS;

X1 1 		 0.6918166 0.53568   0.53568
X2 0.6918166 1 		   0.6432823 0.6432823
X3 0.53568 	 0.6432823 1 		 0.64
X4 0.53568 	 0.6432823 0.64 	 1

;

* obtain factor pattern matrix for later data generation;
PROC FACTOR N=4 OUTSTAT=FACOUT NOPRINT; RUN; 
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

VAR={1 0 0 0,
	0 1.560857 0 0,
	0 0 1 0,
	0 0 0 1};

STD=SQRT(VAR);				* matrix containing stds for the 6 variables;
X=RANNOR(J(&SMPLN,4,&C));		* generate 6 random normal variables; *rep # was put to be used as seed;
XT=X`;						* transpose the data matrix for multiplication;
XTCORR=F*XT;
							* transform uncorrelated variables to correlated ones;
							* transform the scale of the variables (from std=1 to std=specified above);
XTSTD=STD*XTCORR;							* transpose the data matrix back;
XY=XTSTD`; 

* create SAS data set ’DAT’;
CREATE DAT FROM XY[COLNAME={X1 X2 X3 X4}]; 
APPEND FROM XY;
QUIT;

* implement the true model; 
* output the model fitting results to data set ’SEMOUT’;

%IF &MOD=1 %THEN %DO;

PROC CALIS DATA=DAT METHOD=&METHOD COV OUTFIT=SEMOUT NOPRINT; *NOPRINT;
	LINEQS
		X1 = LX1(.8)FK1 				+	EX1,
		X2 = LX1(.8)FK1	+	LX5(.335)FK2+	EX2,
		X3 = LX3(.8)FK2					+	EX3,
		X4 = LX4(.8)FK2					+	EX4;

STD
	FK1 FK2 EX1 EX2 EX3 EX4  = 
	1 1 VEX1(.36) VEX2(.36) VEX3(.36) VEX4(.36);
COV 
	FK1 FK2 = C_FK12(.837); 
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
		X1 = LX1(.8)FK1 				+	EX1,
		X2 = LX1(.8)FK1	+	LX5(.335)FK2+	EX2,
		X3 = LX3(.8)FK2					+	EX3,
		X4 = LX4(.8)FK2					+	EX4;

STD
	FK1 FK2 EX1 EX2 EX3 EX4  = 
	1 1 VEX1(.36) VEX2(.36) VEX3(.36) VEX4(.36);
COV 
	FK1 FK2 = 1; 							*covariance is constrained to be 1;
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
		X1 = LX1(.8)FK1 				+	EX1,
		X2 = LX1(.8)FK1	+	0		FK2	+	EX2,
		X3 = LX3(.8)FK2					+	EX3,
		X4 = LX4(.8)FK2					+	EX4;	* misspecification: LX5 fixed;

STD
	FK1 FK2 EX1 EX2 EX3 EX4  = 
	1 1 VEX1(.36) VEX2(.36) VEX3(.36) VEX4(.36);
COV 
	FK1 FK2 = 1; 							*covariance is constrained to be 1;
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
*libname SEM 'C:\Users\dchoi2\Box\EDPS_971\Simulation studies\montecarlo\CFA3';

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

