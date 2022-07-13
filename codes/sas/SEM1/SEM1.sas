/*******************************************************************/ 
/* this program conducts Monte Carlo simulation of SEM, and outputs*/ 
/* parameter estimates and fit indices to a SAS system file named */ 
/* ’SEM_FITS’. 2 models are simulated in this program: true, */
/* & misspecified models. */
/*******************************************************************/

*OPTIONS MLOGIC MLOGICNEST MPRINT MPRINTNEST SYMBOLGEN MERROR SERROR;
*OPTIONS NONOTES;

%LET MC= 500; * # of Monte Carlo replications for each cell condition; 

libname SEM 'C:\Users\dchoi2\Box\EDPS_971\Simulation studies\montecarlo\SEM1';

* -- to direct the SAS log to a disk file; 
PROC PRINTTO LOG='C:\Users\dchoi2\Box\EDPS_971\Simulation studies\montecarlo\SEM1\LOGFILE.TMP';

DATA A (TYPE=CORR); 
_TYPE_='CORR';
INPUT _NAME_$ X1 X2 X3 X4 X5 X6 Y1 Y2 Y3 Y4 Y5 Y6; 
CARDS;
X1	1			0.7674419	0.760813	0.6244444	0.6499423	0.7076567	0.7165681	0.6937219	0.6899983	0.5487462	0.5970102	0.650363
X2	0.7674419	1			0.8480283	0.7517513	0.7824476	0.7346499	0.7439012	0.7585539	0.766399	0.6606204	0.718724	0.6751707
X3	0.760813	0.8480283	1			0.7882247	0.8204103	0.7283043	0.7374757	0.7679327	0.7805722	0.6926723	0.753595	0.6693389
X4	0.6244444	0.7517513	0.7882247	1			0.8184328	0.5977626	0.6052901	0.6819668	0.7081134	0.6910027	0.7517785	0.5493661
X5	0.6499423	0.7824476	0.8204103	0.8184328	1			0.622171	0.6300059	0.7098136	0.7370278	0.7192184	0.7824759	0.5717984
X6	0.7076567	0.7346499	0.7283043	0.5977626	0.622171	1			0.6859499	0.6640799	0.6605154	0.5252989	0.5715006	0.6225737
Y1	0.7165681	0.7439012	0.7374757	0.6052901	0.6300059	0.6859499	1			0.7973243	0.790735	0.6189543	0.6733933	0.7544619
Y2	0.6937219	0.7585539	0.7679327	0.6819668	0.7098136	0.6640799	0.7973243	1			0.8231293	0.7111343	0.7736808	0.723658
Y3	0.6899983	0.766399	0.7805722	0.7081134	0.7370278	0.6605154	0.790735	0.8231293	1			0.7420744	0.8073422	0.7176774
Y4	0.5487462	0.6606204	0.6926723	0.6910027	0.7192184	0.5252989	0.6189543	0.7111343	0.7420744	1			0.8042626	0.5617679
Y5	0.5970102	0.718724	0.753595	0.7517785	0.7824759	0.5715006	0.6733933	0.7736808	0.8073422	0.8042626	1			0.6111771
Y6	0.650363	0.6751707	0.6693389	0.5493661	0.5717984	0.6225737	0.7544619	0.723658	0.7176774	0.5617679	0.6111771	1
;


* obtain factor pattern matrix for later data generation;
PROC FACTOR N=12 OUTSTAT=FACOUT NOPRINT; RUN; 
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

VAR={115.05 0 0 0 0 0 0 0 0 0 0 0 ,
	0 185.2875 0 0 0 0 0 0 0 0 0 0 ,
	0 0 296 0 0 0 0 0 0 0 0 0 ,
	0 0 0 210.6 0 0 0 0 0 0 0 0 ,
	0 0 0 0 135 0 0 0 0 0 0 0 ,
	0 0 0 0 0 155 0 0 0 0 0 0 ,
	0 0 0 0 0 0 148.1605 0 0 0 0 0 ,
	0 0 0 0 0 0 0 215.4305 0 0 0 0 ,
	0 0 0 0 0 0 0 0 306.0005 0 0 0 ,
	0 0 0 0 0 0 0 0 0 153.4 0 0 ,
	0 0 0 0 0 0 0 0 0 0 160 0 ,
	0 0 0 0 0 0 0 0 0 0 0 222.05};


STD=SQRT(VAR);				* matrix containing stds for the 6 variables;
X=RANNOR(J(&SMPLN,12,&C));	* generate 6 random normal variables; *rep # was used as seed;
XT=X`;						* transpose the data matrix for multiplication;
XTCORR=F*XT;
							* transform uncorrelated variables to correlated ones;
							* transform the scale of the variables (from std=1 to std=specified above);
XTSTD=STD*XTCORR;							* transpose the data matrix back;
XY=XTSTD`; 

* create SAS data set ’DAT’;
CREATE DAT FROM XY[COLNAME={X1 X2 X3 X4 X5 X6 Y1 Y2 Y3 Y4 Y5 Y6}]; 
APPEND FROM XY;
QUIT;

* implement the true model; 
* output the model fitting results to data set ’SEMOUT’;

%IF &MOD=1 %THEN %DO;

PROC CALIS DATA=DAT METHOD=&METHOD COV OUTFIT=SEMOUT NOPRINT;  *NOPRINT;
	LINEQS
		X1 = LX1(.90)FK1 				+	EX1,
		X2 = LX2(.80)FK1+	LX4(.45)FK2 +	EX2,
		X3 = LX3(.80)FK1+	LX5(.80)FK2 +	EX3,
		X4 = 				LX6(1.2)FK2	+	EX4,
		X5 =				 		FK2	+	EX5,
		X6 = 		 FK1				+	EX6,
		Y1 = LY1(.90)FE1 				+	EY1,
		Y2 = LY2(.70)FE1 +	LY4(.45)FE2	+	EY2, 
		Y3 = LY3(.70)FE1 +	LY5(.70)FE2	+	EY3,
		Y4 = 				LY6(.90)FE2	+	EY4, 
		Y5 = 						FE2	+	EY5,
		Y6 = 		 FE1				+	EY6,
		FE1= GA1(1.1)FK1				+	DE1,
		FE2= GA2(1.0)FK2				+   DE2; 

STD
	FK1 FK2 EX1 EX2 EX3 EX4 EX5 EX6 EY1 EY2 EY3 EY4 EY5 EY6 DE1 DE2 = 
	VFK1(105) VFK2(115) VEX1(30) VEX2(30) VEX3(40) VEX4(45) VEX5(20) VEX6(50) VEY1(25) VEY2(40) VEY3(50) VEY4(40) VEY5(20) VEY6(70) VDE1(25) VDE2(25);
COV 
	FK1 FK2 = C_FK12(90), DE1 DE2 = C_DE12(16.20); 
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
		X1 = LX1	FK1 				+	EX1,
		X2 = LX21	FK1	+	0		FK2	+	EX2,
		X3 = LX31	FK1	+	LX32	FK2	+	EX3,
		X4 = 				LX4		FK2	+	EX4,
		X5 =				 		FK2	+	EX5,
		X6 = 		FK1					+	EX6,
		Y1 = LY1	FE1 				+	EY1,
		Y2 = LY21 	FE1 +	0		FE2	+	EY2, 
		Y3 = LY31	FE1 +	LY32	FE2	+	EY3,
		Y4 = 				LY4 	FE2	+	EY4, 
		Y5 = 						FE2	+	EY5,
		Y6 = 		FE1					+	EY6,
		FE1= GA1 	FK1					+	DE1,
		FE2= GA2 	FK2					+   DE2; 
													* misspecification: LX22 LY22 fixed;

STD
	FK1 FK2 EX1 EX2 EX3 EX4 EX5 EX6 EY1 EY2 EY3 EY4 EY5 EY6 DE1 DE2 = 
	VFK1 VFK2 VEX1 VEX2 VEX3 VEX4 VEX5 VEX6 VEY1 VEY2 VEY3 VEY4 VEY5 VEY6 VDE1 VDE2;
COV 
	FK1 FK2 = C_FK12, DE1 DE2 = C_DE12; 
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
		X1 = LX1	FK1 				+	EX1,
		X2 = LX21	FK1	+	0		FK2	+	EX2,
		X3 = LX31	FK1	+	LX32	FK2	+	EX3,
		X4 = 				LX4		FK2	+	EX4,
		X5 =				 		FK2	+	EX5,
		X6 = 		FK1					+	EX6,
		Y1 = LY1	FE1 				+	EY1,
		Y2 = LY21 	FE1 +	0		FE2	+	EY2, 
		Y3 = LY31	FE1 +	LY32	FE2	+	EY3,
		Y4 = 				LY4 	FE2	+	EY4, 
		Y5 = 						FE2	+	EY5,
		Y6 = 		FE1					+	EY6,
		FE1= GA1 	FK1					+	DE1,
		FE2= GA2 	FK2					+   DE2; 
													* misspecification: LX22 LY22  fixed;

STD
	FK1 FK2 EX1 EX2 EX3 EX4 EX5 EX6 EY1 EY2 EY3 EY4 EY5 EY6 DE1 DE2 = 
	VFK1 VFK2 VEX1 VEX2 VEX3 VEX4 VEX5 VEX6 VEY1 VEY2 VEY3 VEY4 VEY5 VEY6 VDE1 VDE2;
COV 
	FK1 FK2 = C_FK12; 
	*DE1 DE2 = C_DE12; 								* misspecification: error covariances fixed to zeros;
	
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
*libname SEM 'C:\Users\dchoi2\Box\EDPS_971\Simulation studies\montecarlo\SEM1';

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

