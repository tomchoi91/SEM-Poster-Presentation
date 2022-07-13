/*******************************************************************/ 
/* this program conducts Monte Carlo simulation of SEM, and outputs*/ 
/* parameter estimates and fit indices to a SAS system file named */ 
/* ’SEM_FITS’. 2 models are simulated in this program: true, */
/* & misspecified models. */
/*******************************************************************/

*OPTIONS MLOGIC MLOGICNEST MPRINT MPRINTNEST SYMBOLGEN MERROR SERROR;
*OPTIONS NONOTES;

%LET MC= 500; * # of Monte Carlo replications for each cell condition; 

libname SEM 'C:\Users\dchoi2\Box\EDPS_971\Simulation studies\montecarlo\CFA1';

* -- to direct the SAS log to a disk file; 
PROC PRINTTO LOG='C:\Users\dchoi2\Box\EDPS_971\Simulation studies\montecarlo\CFA1\LOGFILE.TMP';

DATA A (TYPE=CORR); 
_TYPE_='CORR';
INPUT _NAME_$ X1 X2 X3 X4 X5 X6 X7 X8 X9 X10 X11 X12 X13 X14 X15; 
CARDS;

X1	1 		  0.5093248 0.5457052 0.6104977 0.5820855 0.2829582 0.2829582 0.3031695 0.5659721 0.3233808 0.4414148 0.4414148 0.4729445 0.5044741 0.5044741
X2	0.5093248 1 		0.525	  0.5721967 0.56	  0.245 	0.245	  0.2625	0.3522545 0.28	    0.196	  0.196	    0.21	  0.224	    0.224
X3	0.5457052 0.525	    1		  0.6130679 0.6		  0.2625	0.2625	  0.28125	0.3774155 0.3		0.21	  0.21	    0.225	  0.24	    0.24
X4	0.6104977 0.5721967 0.6130679 1 		0.6539391 0.4904543 0.4904543 0.5254868 0.6043425 0.5605192 0.2561262 0.2561262 0.2744209 0.2927156 0.2927156
X5	0.5820855 0.56 	    0.6		  0.6539391 1		  0.28	    0.28	  0.3		0.4025766 0.32	    0.224	  0.224	    0.24	  0.256	    0.256
X6	0.2829582 0.245	    0.2625	  0.4904543 0.28	  1		    0.49	  0.525	    0.5232016 0.56	    0.147	  0.147	    0.1575	  0.168	    0.168
X7	0.2829582 0.245	    0.2625	  0.4904543 0.28	  0.49	    1		  0.525	    0.5232016 0.56	    0.147	  0.147	    0.1575	  0.168	    0.168
X8	0.3031695 0.2625	0.28125   0.5254868 0.3		  0.525	    0.525	  1		    0.5605731 0.6		0.1575	  0.1575	0.16875	  0.18	    0.18
X9	0.5659721 0.3522545 0.3774155 0.6043425 0.4025766 0.5232016 0.5232016 0.5605731 1		  0.5979446 0.4869401 0.4869401 0.5217215 0.5565029 0.5565029
X10 0.3233808 0.28	    0.3	      0.5605192 0.32	  0.56	    0.56	  0.6		0.5979446 1		    0.168	  0.168	    0.18	  0.192	    0.192
X11 0.4414148 0.196	    0.21	  0.2561262 0.224	  0.147	    0.147	  0.1575	0.4869401 0.168	    1		  0.49	    0.525	  0.56	    0.56
X12 0.4414148 0.196	    0.21	  0.2561262 0.224	  0.147	    0.147	  0.1575	0.4869401 0.168	    0.49	  1		    0.525	  0.56	    0.56
X13 0.4729445 0.21	    0.225	  0.2744209 0.24	  0.1575	0.1575	  0.16875	0.5217215 0.18	    0.525	  0.525	    1		  0.6		0.6
X14 0.5044741 0.224	    0.24	  0.2927156 0.256	  0.168	    0.168	  0.18	    0.5565029 0.192	    0.56	  0.56	    0.6		  1		    0.64
X15 0.5044741 0.224	    0.24	  0.2927156 0.256	  0.168	    0.168	  0.18	    0.5565029 0.192	    0.56	  0.56	    0.6		  0.64	    1

;

* obtain factor pattern matrix for later data generation;
PROC FACTOR N=15 OUTSTAT=FACOUT NOPRINT; RUN; 
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

VAR={1.53 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ,
	0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 ,
	0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 ,
	0 0 0 1.65 0 0 0 0 0 0 0 0 0 0 0 ,
	0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 ,
	0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 ,
	0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 ,
	0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 ,
	0 0 0 0 0 0 0 0 1.826 0 0 0 0 0 0 ,
	0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 ,
	0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 ,
	0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 ,
	0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 ,
	0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 ,
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 1};


STD=SQRT(VAR);				* matrix containing stds for the 6 variables;
X=RANNOR(J(&SMPLN,15,&C));	* generate 6 random normal variables; *rep # was used as seed;
XT=X`;						* transpose the data matrix for multiplication;
XTCORR=F*XT;
							* transform uncorrelated variables to correlated ones;
							* transform the scale of the variables (from std=1 to std=specified above);
XTSTD=STD*XTCORR;							* transpose the data matrix back;
XY=XTSTD`; 

* create SAS data set ’DAT’;
CREATE DAT FROM XY[COLNAME={X1 X2 X3 X4 X5 X6 X7 X8 X9 X10 X11 X12 X13 X14 X15}]; 
APPEND FROM XY;
QUIT;

* implement the true model; 
* output the model fitting results to data set ’SEMOUT’;

%IF &MOD=1 %THEN %DO;

PROC CALIS DATA=DAT METHOD=&METHOD COV OUTFIT=SEMOUT NOPRINT; *NOPRINT;
	LINEQS
		X1 = LX1(.70)	FK1 +	LX16(.50)FK3+	EX1,
		X2 = LX2(.70)	FK1					+	EX2,
		X3 = LX3(.75)	FK1					+	EX3,
		X4 = LX4(.80)	FK1	+	LX17(.50)FK2+	EX4,
		X5 = LX5(.80)	FK1			 		+	EX5,
		X6 = LX6(.70)	FK2 				+	EX6,
		X7 = LX7(.70)	FK2 				+	EX7,
		X8 = LX8(.75) 	FK2 				+	EX8,
		X9 = LX9(.80)	FK2 +	LX18(.70)FK3+	EX9,
		X10 = LX10(.80)	FK2					+	EX10,
		X11 = LX11(.70)	FK3					+	EX11,
		X12 = LX12(.70)	FK3					+	EX12,
		X13 = LX13(.75)	FK3					+	EX13,
		X14 = LX14(.80)	FK3					+	EX14,
		X15 = LX15(.80)	FK3					+	EX15;

STD
	FK1 FK2 FK3 EX1 EX2 EX3 EX4 EX5 EX6 EX7 EX8 EX9 EX10 EX11 EX12 EX13 EX14 EX15  = 
	1 1 1 VEX1(.51) VEX2(.51) VEX3(.4375) VEX4(.36) VEX5(.36) VEX6(.51) VEX7(.51) VEX8(.4375) VEX9(.36) VEX10(.36) VEX11(.51) VEX12(.51) VEX13(.4375) VEX14(.36) VEX15(.36);
COV 
	FK1 FK2 = C_FK12(.50), FK1 FK3 = C_FK13(.40), FK2 FK3 = C_FK23(.30); 
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
		X1 = LX1(.70)	FK1 +	0		FK3	+	EX1,
		X2 = LX2(.70)	FK1					+	EX2,
		X3 = LX3(.75)	FK1					+	EX3,
		X4 = LX4(.80)	FK1	+	LX17(.50)FK2+	EX4,
		X5 = LX5(.80)	FK1			 		+	EX5,
		X6 = LX6(.70)	FK2 				+	EX6,
		X7 = LX7(.70)	FK2 				+	EX7,
		X8 = LX8(.75) 	FK2 				+	EX8,
		X9 = LX9(.80)	FK2 +	LX18(.70)FK3+	EX9,
		X10 = LX10(.80)	FK2					+	EX10,
		X11 = LX11(.70)	FK3					+	EX11,
		X12 = LX12(.70)	FK3					+	EX12,
		X13 = LX13(.75)	FK3					+	EX13,
		X14 = LX14(.80)	FK3					+	EX14,
		X15 = LX15(.80)	FK3					+	EX15;	* misspecification: LX16 fixed;

STD
	FK1 FK2 FK3 EX1 EX2 EX3 EX4 EX5 EX6 EX7 EX8 EX9 EX10 EX11 EX12 EX13 EX14 EX15  = 
	1 1 1 VEX1(.51) VEX2(.51) VEX3(.4375) VEX4(.36) VEX5(.36) VEX6(.51) VEX7(.51) VEX8(.4375) VEX9(.36) VEX10(.36) VEX11(.51) VEX12(.51) VEX13(.4375) VEX14(.36) VEX15(.36);
COV 
	FK1 FK2 = C_FK12(.50), FK1 FK3 = C_FK13(.40), FK2 FK3 = C_FK23(.30); 
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
		X1 = LX1(.70)	FK1 +	0		FK3	+	EX1,
		X2 = LX2(.70)	FK1					+	EX2,
		X3 = LX3(.75)	FK1					+	EX3,
		X4 = LX4(.80)	FK1	+	0		FK2 +	EX4,
		X5 = LX5(.80)	FK1			 		+	EX5,
		X6 = LX6(.70)	FK2 				+	EX6,
		X7 = LX7(.70)	FK2 				+	EX7,
		X8 = LX8(.75) 	FK2 				+	EX8,
		X9 = LX9(.80)	FK2 +	LX18(.70)FK3+	EX9,
		X10 = LX10(.80)	FK2					+	EX10,
		X11 = LX11(.70)	FK3					+	EX11,
		X12 = LX12(.70)	FK3					+	EX12,
		X13 = LX13(.75)	FK3					+	EX13,
		X14 = LX14(.80)	FK3					+	EX14,
		X15 = LX15(.80)	FK3					+	EX15;	* misspecification: LX16 LX17 fixed;

STD
	FK1 FK2 FK3 EX1 EX2 EX3 EX4 EX5 EX6 EX7 EX8 EX9 EX10 EX11 EX12 EX13 EX14 EX15  = 
	1 1 1 VEX1(.51) VEX2(.51) VEX3(.4375) VEX4(.36) VEX5(.36) VEX6(.51) VEX7(.51) VEX8(.4375) VEX9(.36) VEX10(.36) VEX11(.51) VEX12(.51) VEX13(.4375) VEX14(.36) VEX15(.36);
COV 
	FK1 FK2 = C_FK12(.50), FK1 FK3 = C_FK13(.40), FK2 FK3 = C_FK23(.30); 
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
*libname SEM 'C:\Users\dchoi2\Box\EDPS_971\Simulation studies\montecarlo\CFA1';

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

