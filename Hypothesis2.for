      MODULE GlobalVars
	IMPLICIT NONE
	SAVE

      INTEGER MinPYr,MaxPYr,MaxAge,MaxSpawnArea,MaxBox,MinSim,MaxSim
      INTEGER MaxFpar
      PARAMETER (MinPYr=0,MaxPYr=700,MaxAge = 50,MaxSpawnArea = 10,
     +   MaxBox = 2,MinSim = 0,MaxSim = 2000,MaxFpar = 50)

     	INTEGER, PARAMETER :: MaxClosedAreas= 1

!   MinPYr          - First possible year in the model
!   MaxPYr          - Max possible year in the model
!   MaxAge          - Max age in the model
!   MaxSpawnArea    - Maximum number of spawning times/spaces
!   MaxBox          - Maximum number of spatial boxes (not spawning times/spaces)     

!     Major parameters
      INTEGER :: Nsim,Ntime,Nspawnarea,Nage,Nbox
      !INTEGER :: MPASize
      INTEGER :: EnvPreDD,GlobalDD,FindVariance
      INTEGER :: BurnInTime
      REAL(8) :: Fpar
	REAL(8) :: rho,Mpar,hpar,r0
      REAL(8) :: Linf,kpar,t0par
      REAL(8) :: alpha,beta
      REAL(8) :: a50mat,a95mat

      REAL(8),DIMENSION(MaxAge) :: BetaMeans
      REAL(8),DIMENSION(Maxage) :: BetaVar
      REAL(8) :: YoungBetaVar,OldBetavar
      REAL(8) :: YoungBetaMean,OldBetaMean
      !INTEGER :: ConstantBeta !this is no longer needed YoungBetaVar and OldBetaVar take care of this. 
      INTEGER :: OntoEggs
      INTEGER :: NOntoGroup
      INTEGER,DIMENSION(MaxAge) :: OntoGroupsVec

      INTEGER :: DoubleNormal
      REAL(8) :: a50sel,a95sel,Sfull,leftvar,rightvar

! Major Parameter definitions
!   Fpar        - Fishing mortality (Instantaneous)
!   Nsim        - Number of simulations
!   Ntime       - Number of years in the projection
!   Nspawnarea       - Number of spawning times/spaces in the model (boxes are unrelated to spawning time/area)
!   Nbox        - Number of boxes (spatial areas regardless of spawning pattern)
!   MPAsize     - Proportion of total area within an MPA in a 2 box model (which is a special case of this model)
!   EnvPreDD    - 1 if environmental variation happens before density dependent mortality; 0 otherwise
!   BurnInTime  - Number of years for determining effective SSB0 and unfished numbers at age when EnvPreDD = 1
!   GlobalDD    - 1 if density dependent mortality is Global; 0 if DD mortality is local
!   rho         - autocorrelation in environment between spawning times or spaces (areas); 0 is uncorrelated (full portfolio effect), 1 is no portfolio effect
!   Mpar        - Natural mortality
!   hpar        - Steepness
!   r0          - R0, unfished number of recruits
!   Linf        - Asymptotic length in von bert curve
!   kpar        - parameter of von bert curve
!   t0par       - parameter of von bert curve
!   alpha       - parameter of allometric weight-length relationship
!   beta        - parameter of allometric weight-length relationship
!   a50mat      - age at which 50% of animals are mature
!   a95mat      - age at which 95% of animals are mature
!   YoungBetaVar - variance of the beta function distributing the age 0s over spawning time/space
!   OldBetaVar   - variance of the beta function distribution the maxage inidividuals over spawning time/space
!   YoungBetaMean   - Mean of the beta function distribution the youngest individuals over spawning time/space
!   OldBetaMean   - Mean of the beta function distribution for the oldest inidividuals over spawning time/space
!   BetaVar      - vector of Variances of the beta functions (over age) used to allocate numbers at age to spawning times/spaces
!   ConstantBeta     -Is the variance of the beta functions constant or does it increase with age? (1 if constant; 0 otherwise)
!   OntoEggs    - 1 if adults move ontogenetically to spawning time/space and deposit eggs/larvae; 0 if only age-related, but not in big ontogenetic groups 
!   DoubleNormal    - 1 if double normal selectivity is used instead of logistic (asymptotic) selectivity
!   a50sel          - Age at which 50% of individuals are selected to fishery under logistic selectivity
!   a95sel          - Age at which 95% of individuals are selected to fishery under logistic selectivity
!   Sfull           - Age at full selectivity under double normal selectivity
!   leftvar         - variance of left side of double normal curve
!   rightvar        - variance of right side of double normal curve

      REAL(8),DIMENSION(MaxAge) :: Selectivity
      REAL(8),DIMENSION(MaxAge) :: LengthAtAge
      REAL(8),DIMENSION(MaxAge) :: WeightAtAge
      REAL(8),DIMENSION(MaxAge) :: MaturityAtAge
      REAL(8),DIMENSION(MaxAge) :: FecundityAtAge

!  Selectivity      - Selectivity at age (either double normal if DoubleNormal = 1 or logistic if DoubleNormal = 0)
!  LengthAtAge      - Length at age (Von Bertalanffy growth curve)
!  WeightAtAge      - WeightAtAge (Allometric)
!  MaturityAtAge    - Maturity at age (Logistic ogive)

      !Parameters for making random deviate files
	INTEGER, PARAMETER :: MaxSigmaR= 20
	INTEGER, PARAMETER :: Maxrho= 5
      REAL(8) :: MySigmaR 
      INTEGER :: NsigmaR
!      INTEGER :: Nrho
!	REAL(8) :: Crho
	REAL(8),DIMENSION(MinPYr:MaxPYr,MaxSpawnArea,MinSim:MaxSim)::MyError
	REAL(8),DIMENSION(MinPYr:MaxPYr,MaxSpawnArea,MinSim:MaxSim)::CorrError
	REAL(8),DIMENSION(MaxSigmaR)::SigmaRvec
!	REAL(8),DIMENSION(Maxrho)::rhovec
	INTEGER :: ISEEDO
	INTEGER:: SEED1
	INTEGER:: SEED2

C	NsigmaR			- Number of sigmas in SigmaRvec
C	Nrho			- Number of rho values in rhovec
C	MyError			- Non-correlated random errors
C	CorrError		- Correlated random errors
C	SigmaRvec		- A vector of SigmaRs for which to create correlated random devs
C	rhovec			- A vector of rhos specifying the amount of correlation in random devs
C	ISEED0, etc.	- Seeds (from file) for generating random numbers
! For unfished conditions:
      REAL(8),DIMENSION(MaxAge,MaxSpawnArea,MaxBox) :: N0
      REAL(8),DIMENSION(Maxage,MaxBox) :: NOneArea0
      REAL(8),DIMENSION(MaxAge,MaxSpawnArea,MaxBox) :: Ninf
      REAL(8),DIMENSION(Maxage,MaxBox) :: NOneAreaEqm
      REAL(8),DIMENSION(MaxBox) :: SSB0vec
      REAL(8),DIMENSION(MaxSpawnArea,MaxBox) :: SSBbyArea0
      REAL(8),DIMENSION(MaxBox) :: B0vec
      REAL(8),DIMENSION(MaxBox) :: R0vec
      REAL(8),DIMENSION(MaxSpawnArea,MaxBox) :: Qvec0
      REAL(8),DIMENSION(MaxSpawnArea,MaxBox) :: R0byArea
      REAL(8),DIMENSION(MaxBox) :: Rinf,Einf,Binf


!  N0                - Unfished Numbers at age, specific to areas and boxes
!  NOneArea0         - Unfished Numbers at age summed over spawning areas and boxes
!  SSB0vec           - Unfished Spawning stock biomass 
!  SSBbyArea0        - Unfished spawning stock biomass by spawning area (equal to EggsbyArea if EnvPreDD = 0)
!  B0vec             - Unfished total biomass
!  R0vec             - Unfished recruitment
!  Rinf,Einf,Binf    - Recruitment, eggs (SSB), and biomass at fished equilibrium

!     Beta function parameters:
      REAL(8),DIMENSION(Maxage,MaxSpawnArea) :: P
      REAL(8),DIMENSION(Maxage) :: CheckP

!     Population dynamics over time
      REAL(8),DIMENSION(MaxAge,MaxSpawnArea,MaxBox,MinPYr:MaxPYr) :: N
      REAL(8),DIMENSION(MaxAge,MaxBox,MinPYr:MaxPYr) :: NCompact
      REAL(8),DIMENSION(MaxBox,MinPYr:MaxPYr) :: SSB
      REAL(8),DIMENSION(MaxBox,MinPYr:MaxPYr) :: Eggs
      REAL(8),DIMENSION(MaxBox,MinPYr:MaxPYr) :: B
      REAL(8),DIMENSION(MaxSpawnArea,MaxBox,MinPYr:MaxPYr) :: SSBbyArea
      REAL(8),DIMENSION(MaxSpawnArea,MaxBox,MinPYr:MaxPYr) :: EggsbyArea
	REAL(8),DIMENSION(MaxSpawnArea,MaxBox,MinPYr:MaxPYr) :: BbyArea
      REAL(8),DIMENSION(MaxSpawnArea,MaxBox,MinPYr:MaxPYr) :: Recruits
      REAL(8),DIMENSION(MaxBox) :: PreRecruits
      REAL(8),DIMENSION(MaxSpawnArea,MaxBox,MinPYr:MaxPYr) :: Qvec
      REAL(8),DIMENSION(MaxBox,MinPYr:MaxPYr,MinSim:MaxSim)::TotRecruits
      REAL(8),DIMENSION(Maxage,MaxBox,MinPYr:MaxPYr) :: CatchAtAge
      REAL(8),DIMENSION(MinPYr:MaxPYr,MinSim:MaxSim):: TotalCatchBiomass
      REAL(8),DIMENSION(MinPYr:MaxPYr,MinSim:MaxSim)::NonSpatialRecruits
      REAL(8),DIMENSION(MaxSigmaR) :: MeanRecruitCV
      
!     SigmaR-related variables
      REAL(8) :: sigmaR
      REAL(8) :: TargetCV
      REAL(8) :: LowTolerance,HighTolerance
      INTEGER :: MaxIterations
      REAL(8) :: MinDiff
      REAL(8),DIMENSION(MaxSigmaR) :: SigmaRVecNew
      REAL(8) :: BestSigmaR

!  sigmaR       - This is the sigmaR value that actually gets used in a given run
!  TargetCV     - The recruitment CV that is representative of nature; when FindVariance = 1 trying to find the sigmaR that produces the TargetCV
!  LowTolerance    - The amount by which it is acceptable precision for an objective function (LowTolerance used for fmultiply and FMSY objective functions)
!  HighTolerance    - The amount by which it is acceptable precision for an objective function (HighTolerance used for hpar adn MeanDays objective functions)


!     Fishing mortality-related stuff and MPA stuff
      INTEGER :: NFpar
      REAL(8),DIMENSION(MaxFPar) :: Fvec
      REAL(8),DIMENSION(MaxFPar) ::  AvgCumCatch
      INTEGER :: NClosedAreas
      INTEGER :: MPAYear
      INTEGER,DIMENSION(MaxClosedAreas) :: MinClosed,MaxClosed
      INTEGER,DIMENSION(MaxBox,MinPYr:MaxPYr) :: Closed

!     Date and time output variables
      CHARACTER,DIMENSION(8) :: date
      CHARACTER,DIMENSION(10) :: time
      CHARACTER,DIMENSION(5) :: zone
      INTEGER,DIMENSION(8) :: values

!  Permanent movement-related parameters
      INTEGER :: TwoBoxLarvalMovement
      REAL(8) :: MPASize
      REAL(8) :: LMoveIn, LMoveOut
      INTEGER :: PostDDLmove
      REAL(8),DIMENSION(MaxBox) :: WasMoved
      INTEGER,DIMENSION(MinPYr:MaxPYr) :: ClosedArea,OpenArea
      INTEGER :: MoveBurnInTime
      REAL(8) :: SaveSigmaR
      REAL(8) :: SaveMoveBITime
      INTEGER :: DetSim

!  TwoBoxLarvalMovement  - Yes if doing a 2 box model with an exchange rate of larvae between boxes; LMoveIn is the exchange rate.
!                        - in the two box larval movement model the MPA can be a disproportionate size.
!  MoveIn                - The rate of movement between a reserve and open area in a 2 box model
!  MoveOut               - MoveIn, but adjusted for size of reserve, so there is really only 1 rate of movement
!  MoveBurnInTime        - Number of years for burning in deterministic larval movement dynamics (runs for 1 simulation)
!  DetSim      - !simulation number to assign to deterministic movement burn in simulations
      
      !Variables added when adding MAE to the program
      INTEGER :: PEorMAE
      REAL(8),DIMENSION(MinPYr:MaxPYr,MinSim:MaxSim) :: DaysError
      REAL(8) :: SigmaD
      REAL(8) :: MeanDays
      REAL(8), DIMENSION(1:Maxage) :: d50
      INTEGER :: Berkeley
      REAL(8),DIMENSION(1:MaxAge) :: SbyAge
      REAL(8) :: DiffPar
      REAL(8),DIMENSION(1:MaxBox,MinPYr:MaxPYr) :: Sstarve
      REAL(8),DIMENSION(1:MaxBox,MinPYr:MaxPYr) ::UnfishedSstarve
      INTEGER :: MAEctl
      REAL(8),DIMENSION(1:MaxBox) :: NobodySstarves
      REAL(8) :: BestDay
      INTEGER :: NDays
      REAL(8),DIMENSION(MaxSigmaR) :: DaysVec
      REAL(8),DIMENSION(MaxSigmaR) :: DaysVecNew
      INTEGER :: FindFmultiply
      INTEGER :: Nfmultiply
      REAL(8),DIMENSION(MaxSigmaR) :: FmultiplyVec
      REAL(8) :: fmultiply
      REAL(8),DIMENSION(MaxSigmaR) :: AvgR0
      REAL(8),DIMENSION(MaxSigmaR) :: FmultVecNew
      REAL(8),DIMENSION(MaxSigmaR) :: SaveFmultiplyVec
      REAL(8) :: BestFmult
      REAL(8) :: BestR0
      REAL(8) :: MinDiffFmult
      !REAL(8) :: DetFmultMain



!   PEorMAE               - 1 = do portfolio effects model; 2=do maternal age effects model
!   DaysError             - Array with random deviates chosen from the distribution of number of days of starvation ~N(log(MeanDays),SigmaD)
!   SigmaD                - standard deviation of distribution of number of days of starvation
!   MeanDays              - mean number of days of starvation (log(MeanDays) is mean of distribution of # days of starvation, which is then exp(days-sigmaD^2/2)
!   d50                   - days of starvation at which 50% larval mortality occurs for a particular maternal age
!   Berkeley              - 1= use the relationship between maternal age and d50 defined in Berkeley (2005) (no other option currently exists)
!   SbyAge                - larval survival by age, this is folded into one survival proportion (Sstarve or Unfished Sstarve) for eggs integrated over age.
!   DiffPar               - Controls the steepness of the logistic curves relating #days starvation to probability of survival at each maternal age (one curve for each maternal age)
!   Sstarve               - The one survival proportion integrated over ages that is applied to eggs if EnvPreDD = 1 and to recruits if EnvPreDD = 0
!   UnfishedSstarve       - Same as Sstarve, but for MAEctl = 1 (uses the unfished age structure to calculate the value)
!   MAEctl                - Do MAE control test: use the unfished age structure, regardless of the true age structure.
!   NobodySstarves        - Dummy argument for the recruitment function when no MAE (control or not) is used.
!   BestDay               - The number of days of starvation that produces a RecruitCV of TargetCV at unfished biomass
!   NDays                 - Length of the vector of trial #days of starvations that are read in from hypothesis2.dat
!   DaysVec               - Starter several trials for # of days of starvation that result in a RecruitCV at unfished equal to the TargetCV
!   FindFmultiply         - Find a multiplier of fecundity that makes it so that effective R0 is the real R0.
!   FmultiplyVec          - vector of potential fmultiply values that should span the reasonable possibilities for an appropriate fmultiply
!   fmultiply             - The fecundity multiplier that makes effective B0 = B0 (applied to fecundity for envPreDD = 1 adn to recruits otherwise)
!   AvgR0                 - The average R0 over years and sims (post burn-in period) used to find the correct fmultiply.
!   FmultVecNew           - New vector of fmultiply values to try out in the next round
!   BestFmult             - This is the fmultiply value that produced an R0 closest to the true r0

! Params related to stock recruit curve generation for MAE runs
      INTEGER :: SRsims
      REAL(8),DIMENSION(MaxAge,MaxBox) :: FishedNperRecruit
      REAL(8),DIMENSION(MaxAge,MaxBox) :: UnfishedNperRecruit  
      REAL(8),DIMENSION(MaxFpar) :: DRecSR,DSBPRsr,DSSBsr
      REAL(8) :: SBPR0
      REAL(8) :: Dfmultiply
      INTEGER :: SRsimsNFpar
      REAL(8),DIMENSION(MaxFpar) :: SRFvec
      REAL(8) :: MaxSSBsrsims
      REAL(8),DIMENSION(1:100) :: SSBbinVec
      INTEGER :: Nbins
      REAL(8)::SSB0nospace
      REAL(8) :: SaveFmultiply
      ! REAL(8) :: DetermInDays

!   SRsims                - Number of simulations for generating the stock-recruit curve for mae runs (if 0 don't generate stock-recruit curve info)
!   FishedNperRecruit     - what it looks like!
!   UnfishedNperRecruit   - what it looks like!

!    Dfmultiply           - Deterministic fmultiply values
!    MaxSSBsrsims         - Maximum SSB value for SR sims section; necessary info for binning the SSB data
!    DetermInDays         - Deterministic transformation of MeanDays - a trick so that it ends up being the right thing in the end

      REAL(8) :: TotMeanRecruits,TotMeanCVRecruits
      REAL(8),DIMENSION(MinPYr:MaxPyr,0:MaxSim) :: TotSSB
      INTEGER :: ISEED4,ISEED5
      INTEGER :: randominc
      REAL(8),DIMENSION(4,1000) :: FmultDaysCombo
      INTEGER :: itrue
      INTEGER :: Findh
      REAL(8),DIMENSION(500) :: HparVec
      INTEGER :: Nhpar
      REAL(8),DIMENSION(30) :: FindhFvec
      INTEGER :: NhFpar
      REAL(8),DIMENSION(30) :: AvgSSB
      REAL(8),DIMENSION(30) :: AvgDepl
      REAL(8) :: BestDepl,BesthFpar,MinDiffDepl
      REAL(8),DIMENSION(30) :: FindhFvecNew
      REAL(8),DIMENSION(500) :: AvgR
      REAL(8),DIMENSION(500) :: AvgSteepness
      REAL(8) :: hparOriginal
      REAL(8) :: MinDiffh,Besth
      REAL(8),DIMENSION(30) :: hparVecNew
      REAL(8),DIMENSION(30) :: FindhFvecOriginal
      REAL(8),DIMENSION(2,5000) :: HfmultCombo
      INTEGER :: IMeanDays
      REAL(8) :: BestSteep
      REAL(8),DIMENSION(4,1000) :: HFparCombo
      INTEGER :: iHFparCombo
      INTEGER :: md
      REAL(8),DIMENSION(4,1000) :: ThreeWayCombo
      REAL(8) :: hFpar
      REAL(8) :: ReportBestCV


!     TotMeanRecruits       - The average number of recruits at a particular F post-burn in periods that is written out to a data file
!     TotMeanCVRecruits     - The CV of Recruits (average interannual CV) at a particular F post burn-in-periods that is written out to a data file
!     FmultDaysCombo        - A matrix detailing the MeanDays and associated fmultiply for which AvgR0 = r0 AND CVRecruits = TargetCV (added because the program was using the last fmultiply instead of the one associated with the correct MeanDays)
!     Findh                 - Yes = find the h for which R/R0 = 0.6 when SSB/SSB0 = 0.2 (for MAE, h is not steepness)
!     HparVec               - A vector of potential values for h that makes effective steepness = 0.6
!     Nhpar                 - Length of the HparVec
!     FindhFvec             - A vector of Fs used to find the F at which SSB/SSB0 = 0.2
!     NhFpar                - The length of FindhFvec to be read in from Hypothesis2.dat
!     AvgSSB                - Average SSB for each of several Fs (in CalcAvgDepletion(); Used in DoFindH and also coudl be useful output when finding MSY)
!     AvgDepl               - Average Depletion for each of several Fs (in CalcAvgDepletion(); Used in DoFindH and also coudl be useful output when finding MSY)
!     BestDepl              - Depletion level closest to 0.2
!     BesthFpar             - Fpar associated with depletion level closest to 0.2
!     MinDiffDepl           - Minimum difference between BestDepl and 0.2
!     AvgR                  - Average R for a given F
!     AvgSteepness          - Average EFFECTIVE steepness for a given R
!     hparOriginal          - hpar all over again - save it because h changes.
!     MinDiffh              - Minimum diff between Besth and hparOriginal
!     Besth                 - the h that makes steepness closest to hparOriginal
!     HfmultCombo           - A matrix to align the best fmultiply value associated with each h (this is for each MeanDays - the info for the BEST h and its associated fmultiply is passed along to FmultDaysCombo()
!     IMeanDays             - Counts the number of MeanDays that are tried out.
!     BestSteep             - Steepness closest to R/R0 = 0.6 when Depl = 0.2
!     HFparCombo            - a matrix that records BestHFpar for a particular fmultiply, hpar, and MeanDays combo and reports all of these together.
!     md                    - Counts the number of mean days iterations that are done and uses to assign ThreeWayCombo, which records corresponding hpar, fmultiply,and hFpar for each MeanDays value.
!     hFpar                 - This is the fpar associated with the best h that makes depletion = 0.2
!     ReportBestCV          - Same as BestCV, but global - don't want to mess up BestCV itself.

!Variables for finding MSY
       REAL(8) :: FindMSY
       REAL(8),DIMENSION(MaxFpar) :: FindMSYVec,FindMSYvecNew
       REAL(8),DIMENSION(MaxFpar) :: AvgNegCumCatch
       REAL(8),DIMENSION(MaxFpar) :: AvgCatch
       REAL(8) :: MinCatches
       REAL(8) :: BestFMSY
       REAL(8) :: BestSSBMSY
       REAL(8) :: BestDeplMSY
       REAL(8) :: BestMSY

       REAL(8) :: FMSY,MSY,SSBMSY,DeplMSY
       REAL(8),DIMENSION(MinPYr:MaxPYr,MinSim:MaxSim) ::  TotBiomass


!     FindMSY               - Switch for DoMAEInnards so that EQMPOP is run again (otherwise ends on EQMPOP for a high F from the main run)
!     FindMSYVec            - Vector of trial F values that could be Fmsy
!     AvgNegCumCatch        - Average negative cumulative catch (cumulative catch over years, averaged over sims: the objfun for finding MSY
!     AvgCatch              - catch averaged over sims and time (post burn in periods)
!     MinCatches            - Best objective function for FMSY (minimum avgerage neg. cum catches)
!     BestFMSY,BestSSBMSY,BestDeplMSY,BestMSY - Values closest to MSY for that round of trial Fs in the end these end up being the MSY values.
!     FMSY,MSY,SSBMSY,DeplMSY  - these are the final MSY-related values

      END MODULE GlobalVars
!======================================================================
! Input files
! =============
! HYPOTHESIS2.DAT (Unit 13) Main input file

! Output files
! =============
!! Outparams.OUT  (Unit 10) Parameter values at the end of a run (for checking to see if they still match the input file?)
! Recruits.OUT    (Unit 11) Sim, OutRecruits,InRecruits,Time,F,TotalRecruits 
! Eggs.OUT        (Unit 12) Sim,OutEggs(Year-1),InEggs(Year-1),Year-1,F,TotalEggs(Year-1) !Eggs operate a year behind
! EggsByArea.OUT  (Unit 14; 13 is input file) Sim,Year-1,F,TotEggsByArea(1:Nspawnarea)
! SSB.OUT         (UNIT 15) Sim,OutSSB(Year),InSSB(Year),Year,TotalSSB(Year)
! SSBbyArea.OUT   (UNIT 16) Sim,Year,F,TotSSBByArea(1:Nspawnarea) (same as eggs in a given year if EnvPreDD = 0; diff is environ. var. if EnvPreDD = 1)
! TotalCatch.OUT  (UNIT 17) 
! NCompact.OUT    (UNIT 18) Sim,Year,F,NCompact(1:Nage,box = 1,Year) (made for debugging - only includes box = 1, not box = 2 if there's an MPA)
! FProfile.OUT    (UNIT 19) Fpar, cumulative catch from BurnInTime to the last time period averaged over simulations
! TimeDateStamp.OUT  (UNIT 20) Start time and date and end time and date (see how long the simulation took)
! FindVarianceInfo.OUT  (UNIT 21) MinDiff,BestCV, and BestSigmaR for each set of sigmaRvecs used to find a sigmaR producing RecruitCV = 0.5
! RecruitsByAreaAndBox.OUT (UNIT 22) Sim,Year,F,Box,SpawnArea,Recruits
! NdetailedFparX.XX.OUT (UNIT 50+1 to 22+NFpar) Sim,Year,F,Box,SpawnArea,Age,N
!=====================================================================
      PROGRAM Hypothesis2
      
      USE GlobalVars
      IMPLICIT NONE
      
!      Variables
       INTEGER :: IYear,Isim,IsigmaR,iterations,IFpar,Iiterbug
       INTEGER :: IYr,Ibox,IClosed,Idays
       CHARACTER*4 TheF
       CHARACTER*18 MyNname
       REAL(8) :: RoundedCombo1,RoundedMeanDays
!       REAL*8 XNORM
!       EXTERNAL XNORM

      NobodySstarves = 1.0
      DetSim = 0 !simulation number to assign to deterministic movement burn in simulations
      OPEN(UNIT=20,FILE='TimeDateStamp.OUT')
      CALL mydate('Start time and date      ')
!     Output files
      OPEN(UNIT=11,FILE="RECRUITS.OUT")
      OPEN(UNIT=12,FILE="EGGS.OUT")
      OPEN(UNIT=14,FILE="EGGSBYAREA.OUT")
      OPEN(UNIT=15,FILE="SSB.OUT")
      OPEN(UNIT=16,FILE="SSBbyArea.OUT")
      OPEN(UNIT=17,FILE="TotalCatch.OUT")
      OPEN(UNIT=18,FILE="NCompact.OUT")
      OPEN(UNIT=19,FILE="FProfile.OUT")
      OPEN(UNIT=21,FILE="FindVarianceInfo.OUT")
      OPEN(UNIT=22,FILE="RecruitsByAreaAndBox.OUT")
      OPEN(UNIT=23,FILE="SRsims.OUT")
      OPEN(UNIT=24,FILE="DeterministicSRsims.OUT")
      OPEN(UNIT=25,FILE="DaysErrorSR.OUT")
      OPEN(UNIT=26,FILE="MeansAndCVs.OUT")
      OPEN(UNIT=27,FILE="DaysErrorForDebug.OUT")
      OPEN(UNIT=28,FILE="FindMSYInfo.OUT")
      OPEN(UNIT=29,FILE="FinalMSYInfo.OUT")
      OPEN(UNIT=30,FILE="TotalBiomass.OUT")
      OPEN(UNIT=31,FILE="CumInfobySim.OUT")
      OPEN(UNIT=32,FILE="CumInfo.OUT")
      ! Ndetailed files have to be opened after readspec.


!------------------------------------------------------------------------------------------------------------
!      Deterministic things that only need to happen once, regardless of FindVariance, Burn in situation, vector of fishing mortalities
       CALL READSPEC()

       DO 667 IFpar = 1,NFpar
         WRITE(TheF,'(F4.2)') Fvec(IFpar) 
         MyNname = 'NdetailedF'//TheF//'.OUT'
         OPEN(UNIT=50+IFpar,FILE=MyNname)
667    CONTINUE

!Moved this to UnfishedPop
       !Distribute R0 proportionally across the boxes
!       DO 7659 Ibox = 1,Nbox
!         R0vec(Ibox) = r0*(1/NBox)
!7659   CONTINUE

!Assign closed areas here, in an obvious spot, so that it is clear how this happens
       Closed = 1 !Area is not closed
       DO 7660 IYr = MPAYear,MinPYr+MoveBurnInTime +BurnInTime+Ntime
         DO 7660 IClosed = 1,NClosedAreas
           DO 7660 Ibox = MinClosed(IClosed),MaxClosed(IClosed)
             Closed(Ibox,IYr) = 0 !Area is closed
7660   CONTINUE
       !Temporary:
    !   ClosedArea = 1
    !   OpenArea = 1
    !   DO 110 IYr = MinPYr,(MinPYr+BurnInTime+Ntime)
    !   IF (TwoBoxLarvalMovement.EQ.1.AND.IYr.GE.MPAYear) THEN
    !     DO 111 Ibox = 1,Nbox
    !       IF (Closed(Ibox,MinPYr).EQ.0) THEN
    !         ClosedArea(IYr) = Ibox
    !       ELSE
    !         OpenArea(IYr) = Ibox
    !       ENDIF
!111      CONTINUE
!       ENDIF
!110    CONTINUE
       
       CALL BASICS()
       
       IF (PEorMAE.EQ.2.AND.Nspawnarea.GT.1) THEN
         IF (OntoEggs.EQ.0) THEN
           CALL BETASTUFF()
         ELSE 
          CALL BetaOntogenetic()
        ENDIF
      ELSE
        P = 1.0
      ENDIF
       
       CALL UNFISHEDPOP()
       
       IF (PEorMAE.EQ.1) THEN
         CALL EQMPOP()
         !CALL DetFmultMain()
         CALL CalcsByAge(MeanDays)
         CALL DetFmultMainTry2()
         !Run MoveBurnIn period for F = 0 here
         fmultiply = DFmultiply !double check variable name
         !DetermInDays = log(MeanDays)+(SigmaD**2)/2 ! don't know if you need this or not
         CALL FirstTime(MinPYr)
         DO 4201 Isim = 1,Nsim  !you don't really need the Isim loop: this is deterministic
           DO 4201 IYear = MinPYr+1,MinPYr + MoveBurnInTime 
             !Maternal age stuff
             !CALL CalcsByAge(MeanDays) !Already ran this function - it doesn't depend on Age or anything.
             IF (MAEctl.EQ.1) THEN
              CALL UnfishedCalcSstarve(IYear)
              CALL RECRUITMENT(IYear,Isim,UnfishedSstarve(1:Nbox,IYear))
             ELSEIF (MAEctl.EQ.0) THEN
               CALL CalcSstarve(IYear)
               CALL RECRUITMENT(IYear,Isim,Sstarve(1:Nbox,IYear))
              ELSEIF (MAEctl.GT.1) THEN
               CALL ControlCalcSstarve(IYear,MAEctl)
               CALL RECRUITMENT(IYear,Isim,Sstarve(1:Nbox,IYear))
              ENDIF
               CALL POPDY(IYear,Isim)
               !Do you want to write anything out? do that here.   
4201    CONTINUE
       ENDIF

      IF (FindVariance.EQ.0.AND.PEorMAE.EQ.1) THEN
        AvgR0 = 0.0
        SaveMoveBITime = MoveBurnInTime
        CALL GenStarveDistnBI(MeanDays)
        CALL DoFindFmult()
        MoveBurnInTime = SaveMoveBITime
      ENDIF
!-------------------------------------------------------------------------------------------------------------

!Find Variance Chunk for MAE runs:
       !------------------------------------------------------------------------------------
       !Find Variance Chunk for MAE
       
       IF (FindVariance.EQ.1.AND.PEorMAE.EQ.1) THEN
        Closed = 1 !Area is not closed: all areas are open for optimizations!!!
       itrue = 0
       Fpar = 0.0
       iterations = 0  !md is defined globally
       IMeanDays = 0
       iHFparCombo = 0
       SaveMoveBITime = MoveBurnInTime
 !      MoveBurnInTime = 0
       SaveFmultiplyVec = FmultiplyVec
       PRINT *,'Just entering FindVariance' 
         DO
           iterations = iterations+1
           DO 4000 Idays = 1,NDays 
            md = md + 1 !Same as IMeanDays
            MeanDays = DaysVec(Idays)
            PRINT *,'MeanDays in FindVariance loop is now ',MeanDays
            IMeanDays = IMeanDays+1 !Same as md

          !----------------------------------------------------------------------------------
          !Right here technically need to get proper Dfmultiply and do the burn-in period with it.
          CALL EQMPOP()
         !CALL DetFmultMain()
         CALL CalcsByAge(MeanDays)
         CALL DetFmultMainTry2()
         !Run MoveBurnIn period for F = 0 here
         fmultiply = DFmultiply !double check variable name
         !DetermInDays = log(MeanDays)+(SigmaD**2)/2 ! don't know if you need this or not
         CALL FirstTime(MinPYr)
         DO 5301 Isim = 1,Nsim  !you don't really need the Isim loop: this is deterministic
           DO 5301 IYear = MinPYr+1,MinPYr + MoveBurnInTime 
             !Maternal age stuff
             !CALL CalcsByAge(MeanDays) !Already ran this function - it doesn't depend on Age or anything.
             IF (MAEctl.EQ.1) THEN
              CALL UnfishedCalcSstarve(IYear)
              CALL RECRUITMENT(IYear,Isim,UnfishedSstarve(1:Nbox,IYear))
             ELSEIF (MAEctl.EQ.0) THEN
               CALL CalcSstarve(IYear)
               CALL RECRUITMENT(IYear,Isim,Sstarve(1:Nbox,IYear))
              ELSEIF (MAEctl.GT.1) THEN
               CALL ControlCalcSstarve(IYear,MAEctl)
               CALL RECRUITMENT(IYear,Isim,Sstarve(1:Nbox,IYear))
              ENDIF
               CALL POPDY(IYear,Isim)
               !Do you want to write anything out? do that here.   
5301    CONTINUE
            !---------------------------------------------------------
            !End finding Dfmultiply
            !--------------------------------------------------------
             CALL GenStarveDistnBI(MeanDays)
             !Within each DaysVec option find the correct fmultiply, then apply that fmultiply when running sims to find RecruitCV below
             AvgR0 = 0.0
          !Change so that CALL DoFindH() and DoFindFmult() will be inside of DoFindH()
              IF (FindH.EQ.1) THEN
                CALL DoFindH() !coming into being
              ELSE
               CALL DoFindFmult() !old and now called within DoFindH
              ENDIF
  !This is already initialized by the MoveBurnIn: different from earlier (check this out for PEorMAE = 2)          
 !            CALL EQMPOP()          
 !            CALL FirstTime(MinPYr) 

          
!---------------------------------NEED THIS SECTION TO RUN THINGS WITH CORRECT FMULTIPLY-----------------------          
            !set F back to 0 at end of DoFindH or here?
            CALL DoMAEInnards()
!Seeing if DoMAEInnards() works instead of code below (down to the !----)
!             DO 4001 Isim = 1,Nsim
!             DO 4001 IYear = (MinPYr+MoveBurnInTime+1),(MinPYr+
!     +                       MoveBurnInTime+BurnInTime+Ntime)
!               !Maternal age stuff
!               CALL CalcsByAge(exp(DaysError(IYear,Isim)-(SigmaD**2/2)))
!
!               IF (MAEctl.EQ.1) THEN
!                CALL UnfishedCalcSstarve(IYear)
!              CALL RECRUITMENT(IYear,Isim,UnfishedSstarve(1:Nbox,IYear))
!               ELSEIF (MAEctl.EQ.0) THEN
!                 CALL CalcSstarve(IYear)
!                 CALL RECRUITMENT(IYear,Isim,Sstarve(1:Nbox,IYear))
!              ELSEIF (MAEctl.GT.1) THEN
!                 CALL ControlCalcSstarve(IYear,MAEctl)
!                 CALL RECRUITMENT(IYear,Isim,Sstarve(1:Nbox,IYear))
!              ENDIF
!               
!               CALL POPDY(IYear,Isim)
!               !Do you want to write anything out? do that here.
!4001         CONTINUE
!---------------------------------NEED ABOVE TO RUN THINGS WITH CORRECT FMULTIPLY, but also need this after get correct MeanDays-----------------------

          !Calculate the CV of recruits based on NonSpatialRecruits(Year,Sim) - make a subroutine for calculating this. Use R.
          !Record for each MeanDays
           

           CALL CalcRecruitCV(Idays)
           
!DEBUGGING LINES------------------------------------------------------------------
     
!DEBUGGING LINES------------------------------------------------------------------


4000       CONTINUE
        
         !Find value closest to 0.5
         !Make a new sigmaRvec and do it all over again if needed.
           CALL DiffsAndNewDaysVec
        IF (iterations.GE.MaxIterations.OR.MinDiff.LE.LowTolerance) THEN
               MeanDays = BestDay
               fmultiply = -1.0
               !Go find the fmultiply value for this MeanDays
               
               !works when no FindH
               IF (FindH.EQ.0) THEN
               RoundedMeanDays = 0.0
               RoundedCombo1 = 0.0
               DO 5547 Iiterbug = 1,(MaxIterations*NDays)
                 !  RoundedMeanDays = nint(MeanDays*1000.0)*0.
                    RoundedMeanDays = nint(MeanDays*1000000.0)*1E-06
                   RoundedCombo1 = nint(FmultDaysCombo(1,Iiterbug)*
     +                                  1000000.0)*1E-06
                 IF (RoundedCombo1.EQ.RoundedMeanDays) THEN
                   fmultiply = FmultDaysCombo(2,Iiterbug)
                   IF (FindH.EQ.1) THEN
                     hpar = FmultDaysCombo(3,Iiterbug)
                   ENDIF
            !     ELSE 
            !       fmultiply = -1.0
                 ENDIF
5547           CONTINUE
          
               ELSE !when FIndH = 1
               DO 5549 Iiterbug = 1,(MaxIterations*NDays)
                 IF (MeanDays.EQ.ThreeWayCombo(4,Iiterbug)) THEN
                   hFpar = ThreeWayCombo(1,Iiterbug)
                   fmultiply = ThreeWayCombo(2,Iiterbug)
                   hpar = ThreeWayCombo(3,Iiterbug)
                 ENDIF !meanDays equals ThreeWayCombo...
5549           CONTINUE

               ENDIF !if FindH.EQ.0 else FindH = 1
         
              IF (FindH.EQ.1) THEN      
                WRITE(21,980) MinDiff,ReportBestCV,BestDay,fmultiply,
     +hpar,hFpar
              ELSE !If FIndH equals 0
               WRITE(21,981) MinDiff,ReportBestCV,BestDay,fmultiply,hpar
             ENDIF !IF FindH.eQ.1


               PRINT *,'Final MeanDays ', MeanDays ,'and MinDiff is',
     +           MinDiff
              PRINT *,'Final Fmultiply ', fmultiply 
              PRINT *,'Final hpar ', hpar
              IF (FindH.EQ.1) THEN
                PRINT *,'Final hFpar ',hFpar
              ENDIF !IF FindH.EQ.1
              IF (MinDiff.GT.LowTolerance) THEN
                PRINT *, 'LowTolerance for MeanDays not reached; 
     + tolerance was ',LowTolerance
              ENDIF !If MinDiff.GT.LowTolerance

              WRITE(21,699) MeanDays,fmultiply,hpar
        ELSE  !not IF (iterations.GE.MaxIterations.OR.MinDiff.LE.LowTolerance) 
          DaysVec = DaysVecNew
        ENDIF  !IF (iterations.GE.MaxIterations.OR.MinDiff.LE.LowTolerance) 
        IF (iterations.GE.MaxIterations.OR.MinDiff.LE.LowTolerance) THEN
             MoveBurnInTime =  SaveMoveBITime
            exit
        ENDIF
       !Before the end of this chunk, re-assign sigmaR for future chunks
       END DO  !End iterating over FmultiplyVecs
      
      !Reset Closed for F simulations - everything needs to be open during optimizations.
      !Assign closed areas here, in an obvious spot, so that it is clear how this happens
       Closed = 1 !Area is not closed
       DO 7661 IYr = MPAYear,MinPYr+MoveBurnInTime +BurnInTime+Ntime
         DO 7661 IClosed = 1,NClosedAreas
           DO 7661 Ibox = MinClosed(IClosed),MaxClosed(IClosed)
             Closed(Ibox,IYr) = 0 !Area is closed
7661   CONTINUE

      ENDIF !FindVariance = 1
       !End FindVariance for MAE Runs
       !-----------------------------------------------------------------------------------



       !------------------------------------------------------------------------------------
       !Find Variance Chunk
       IF (FindVariance.EQ.1.AND.PEorMAE.EQ.2) THEN
       Fpar = 0
       iterations = 0
       SaveMoveBITime = MoveBurnInTime
       MoveBurnInTime = 0

         DO
           iterations = iterations+1
           DO 4010 IsigmaR = 1,NsigmaR
             sigmaR = SigmaRVec(IsigmaR)

             CALL EQMPOP()
             CALL FirstTime(MinPYr)
             CALL CORRMODEL2(sigmaR,rho) !rho will be set in Hypothesis2.dat every time now that we have the seed thing worked out and don't need to compare to the global case; just gotta get TargetCV = 0.5

             DO 4011 Isim = 1,Nsim
             DO 4011 IYear = MinPYr +1,MinPYr+BurnInTime+Ntime
               CALL RECRUITMENT(IYear,Isim,NobodySstarves(1:Nbox))
               CALL POPDY(IYear,Isim)
               !Do you want to write anything out? do that here.
4011         CONTINUE
          !Calculate the CV of recruits based on NonSpatialRecruits(Year,Sim) - make a subroutine for calculating this. Use R.
          !Record for each sigmaR
           CALL CalcRecruitCV(IsigmaR)
           !CALL DiffsToleranceNewSigmaRvec

4010       CONTINUE
        
         !Find value closest to 0.5
         !Make a new sigmaRvec and do it all over again if needed.
           CALL DiffsAndNewSigmaRvec
        IF (iterations.GE.MaxIterations.OR.MinDiff.LE.LowTolerance) THEN
               sigmaR = BestSigmaR
               PRINT *,'Final sigmaR ', sigmaR ,'and MinDiff is',
     +           MinDiff
                 IF (MinDiff.GT.LowTolerance) THEN
                   PRINT *, 'LowTolerance for sigmaR not reached; 
     +               tolerance was ',LowTolerance
                 ENDIF
             ELSE
               SigmaRVec = SigmaRVecNew
           ENDIF
        IF (iterations.GE.MaxIterations.OR.MinDiff.LE.LowTolerance) exit
       !Before the end of this chunk, re-assign sigmaR for future chunks
       END DO
       MoveBurnInTime =  SaveMoveBITime
         ENDIF
       !End FindVariance
       !-----------------------------------------------------------------------------------

       IF (FindVariance.EQ.0.AND.PEorMAE.EQ.2) THEN
         sigmaR = SigmaRvec(1)
       ENDIF
       CALL CORRMODEL2(sigmaR,rho)

       IF (PEorMAE.EQ.1) THEN
          CALL GenStarveDistnBI(MeanDays) !I think this could be outside of Fpar.
 !         CALL GenStarveDistnBI(MeanDays) !Just for debugging purposes, otherwise get rid of this line!
       ENDIF
       
       !-----------------------------------------------------------------------------------
       !Run model over a range of fishing mortality levels
       DO 5400 IFpar = 1,NFpar
         Fpar = Fvec(IFpar)
         
       SaveSigmaR = sigmaR  !sigmaR gets set to 0.000001 for the MoveBurnIn Period and then back to what it was.
       CALL EQMPOP()
       CALL FirstTime(MinPYr)

       !-----------------------------------------------------------------------------
       !Deterministic Burn-In Period
       IF (MoveBurnInTime.GT.0) THEN
         sigmaR = 0.000000001
       ENDIF

       IF (MoveBurnInTime.GT.0.AND.PEorMAE.EQ.1) THEN
         SaveFmultiply = fmultiply
         !fmultiply = Dfmultiply !nope - this is not necessarily for the correct meandays; have to run the calculation over again!
         CALL CalcsbyAge(MeanDays)
         CALL DetFmultMainTry2()
       ENDIF

       Isim = DetSim
       !for the burn-in period use MeanDays for SbyAge and Sstarve is based on eqm population numbers
       IF (PEorMAE.EQ.1) THEN
         CALL CalcsbyAge(MeanDays)
 !        CALL DetFmultMainTry2() !Need this to be outside of Fpar loop
          IF (MAEctl.EQ.1) THEN
              CALL UnfishedCalcSstarve(MinPYr)
          ELSEIF (MAEctl.GT.1) THEN
              CALL ControlCalcSstarve(MinPYr+1,MAEctl) !1 year gets subtracted here but not in UnfishedCalcSstarve
          ENDIF
         ENDIF

       IF (PEorMAE.EQ.1) THEN
         fmultiply = Dfmultiply
       ENDIF
       DO 2001 IYear = MinPYr+1,MinPYr+MoveBurnInTime
        IF (PEorMAE.EQ.1) THEN
         !calculate Sstarve here, so you are taking into account the burn-in stuff.
         IF (MAEctl.EQ.1) THEN
          DO 4333 Ibox = 1,Nbox
           UnfishedSstarve(Ibox,IYear)= UnfishedSstarve(Ibox,MinPYr)
 4333     CONTINUE
         ELSEIF (MAEctl.EQ.0) THEN
          CALL CalcSstarve(IYear) !still use CalcsbyAge(MeanDays) from above because still deterministic
         ELSEIF (MAEctl.GT.1) THEN
          CALL ControlCalcSstarve(IYear,MAEctl)
         ENDIF !If MAEctl.EQ.1
        ENDIF !If PEorMAE.EQ.1

         !this does not create stochsticity because InDays is always MeanDays.
         IF (PEorMAE.EQ.2) THEN
           CALL RECRUITMENT(IYear,Isim,NobodySstarves(1:Nbox))
         ELSE
          IF (MAEctl.EQ.1) THEN
           CALL RECRUITMENT(IYear,Isim,UnfishedSstarve(1:Nbox,IYear))
          ELSE 
           CALL RECRUITMENT(IYear,Isim,Sstarve(1:Nbox,IYear))
          ENDIF
         ENDIF
         CALL POPDY(IYear,Isim)
         CALL SUMMARY2(IYear,Isim)
2001  CONTINUE
      !Deterministic Burn-In is over; return sigmaR to correct value
      sigmaR = SaveSigmaR
      IF (MoveBurnInTime.GT.0.AND.PEorMAE.EQ.1) THEN
        fmultiply = SaveFmultiply
      ENDIF
      !End Deterministic burn in for larval movement dynamics
      !----------------------------------------------------------------------------------------

      !----------------------------------------------------------------------------------------
      !Do stochastic runs
      !----------------------------------------------------------------------------------------
        DO 2000 Isim = 1,Nsim
          DO 2000 IYear=MinPYr+MoveBurnInTime+1,MinPYr+MoveBurnInTime+
     +BurnInTime+Ntime !Change years to start after larval burn in
      !IF PEorMAE.EQ.1 THEN 
            !- choose the days deviate (already exists)
            !- run SbyAgeCalcs(this year's deviate)
            IF (PEorMAE.EQ.1) THEN
             CALL CalcSbyAge(exp(DaysError(IYear,Isim)-(SigmaD**2)/2))
             !- run UnfishedCalcSstarve and CalcSstarve
             IF (MAEctl.EQ.1) THEN
               CALL UnfishedCalcSstarve(IYear)
             ELSEIF (MAEctl.EQ.0) THEN
              CALL CalcSstarve(IYear)
             ELSEIF (MAEctl.GT.1) THEN
              CALL ControlCalcSstarve(IYear,MAEctl)
             ENDIF
            !- go through same if statement sequence as above to call the recruitment function
            ENDIF
          IF (PEorMAE.EQ.2) THEN
           CALL RECRUITMENT(IYear,Isim,NobodySstarves(1:Nbox))
         ELSE
          IF (MAEctl.EQ.1) THEN
           CALL RECRUITMENT(IYear,Isim,UnfishedSstarve(1:Nbox,IYear))
          ELSE 
           CALL RECRUITMENT(IYear,Isim,Sstarve(1:Nbox,IYear))
          ENDIF
         ENDIF


           CALL POPDY(IYear,Isim)
           CALL SUMMARY2(IYear,Isim) !to read data out to files (could do this outside the year loop)
2000    CONTINUE
        
        !Output average cumulative catch
        CALL CumulativeCatch(IFpar) !Change so that calculations start after larval burn in and other burn in.

        !Add something here to process data to make things quicker/easier in R?
        !-----------------------------------------------------------------------------------
        !end stochastic runs
        !-----------------------------------------------------------------------------------
          CALL SummaryStats(IFpar) !Write out average and CVs
          CALL MoreSummaryStats(IFpar)
5400    CONTINUE !End loop over IFpar
!--------------------------------------------------------------------------------
!End profiling over vector of Fs
!--------------------------------------------------------------------------------

 !Find MSY
      CALL DoFindMSY()

 !Andre's SBPR re-parameterization exercise: might not use this in the end.
       CALL DoSRsims()
       !add here DoDeterministicSRsims() and put inside of that DetFmult (change name from DetFmult to FindDfmultiply); run this outside of an IFpar loop.
!!Do a loop over SBPRsims and F
!!Write things out as you go.
!      IF (SRsims.GT.0) THEN
!        DO 5412 IFpar = 1,NFpar
!         CALL EQMPOP()
!          DO 5411 Isim = 1,SRsims
!            !Calculate gamma for each year and sim
!      
!5411    CONTINUE
!5412   CONTINUE
!      ENDIF 


        
      
!------------------------------------------------------------
!       Close the output files
        CLOSE(11)
        CLOSE(12)
        !13 is the input file: Hypothesis2.dat
        CLOSE(14)
        CLOSE(15)
        CLOSE(16)
        CLOSE(17)
        CLOSE(18)
        CLOSE(19)
        CLOSE(21)
        CLOSE(22)
        CLOSE(23)
        CLOSE(24)
        CLOSE(25)
        CLOSE(26)
        CLOSE(27)
        CLOSE(28)
        CLOSE(29)
        CLOSE(30)
        CLOSE(31)
        CLOSE(32)
!-------------------------------------------------------------


!       End the program
        CALL mydate('End time and date        ')
        CLOSE(20)

!MeanDays,fmultiply,IYear,ISim,DaysError(IYear,ISim)        
877     FORMAT(F20.10,1x,F20.15,1x,I4,1x,I5,1x,F20.10)
699     FORMAT('Final MeanDays is ',F15.7,' and corresponding fmultiply
     + is ',F15.7,' and Besth is ',F15.12)
980   FORMAT('MinDiff MeanDays= ',F10.7,', Best CV= ',F10.7,
     + ', Best MeanDays= ',F10.7,' using fmultiply ',F10.7,' and hpar '
     +,F10.7,' and hfpar ',F10.7)    
981   FORMAT('MinDiff MeanDays= ',F10.7,', Best CV= ',F10.7,
     + ', Best MeanDays= ',F10.7,' using fmultiply ',F10.7,' and hpar '
     +,F10.7)  
        STOP
        END PROGRAM
       INCLUDE 'COMMON.FOR'




!================================================================
      SUBROUTINE ReadSpec()
!
! This subroutine reads in the specifications for the run
!      
	  USE GlobalVars
      IMPLICIT NONE
!
!     Local variables (Example if needed)
!      INTEGER Ispec,Isex,Itime,Iage,IblkNS,IblkEW,JblkNS,JblkEW
!      INTEGER Yr,Iyr,Iarea,II
      INTEGER :: IsigmaR,Irho,IFpar,Ibox,Idays,Ifmult,Ihpar
      REAL(8) :: hdist
!
!     Read in the specifications
      OPEN(UNIT=13,FILE='HYPOTHESIS2.DAT')
!     Basic Model info
      READ(13,'(40X,I2)') PEorMAE 
      READ(13,'(40X,I4)') Nsim 
      READ(13,'(40X,I4)') Ntime
      READ(13,'(40X,I4)') Nspawnarea
      READ(13,'(40X,I4)') Nage
      READ(13,'(40X,I4)') Nbox

!     Hypothesis 2 variables
      READ(13,*)
      READ(13,'(40X,I4)') EnvPreDD
      READ(13,'(40X,I4)') BurnInTime
      READ(13,'(40X,I4)') MoveBurnInTime
      READ(13,'(40X,I4)') GlobalDD
      READ(13,'(40X,I4)') FindVariance
	READ(13,'(40X,F5.3)') rho

!     Age structured model variables
      READ(13,*)
	READ(13,'(40X,F5.3)') Mpar
	READ(13,'(40X,F5.3)') hpar
	READ(13,'(40X,F5.3)') r0
	READ(13,'(40X,F5.3)') Linf
	READ(13,'(40X,F5.3)') kpar
	READ(13,'(40X,F5.3)') t0par
	READ(13,'(40X,F10.7)') alpha
	READ(13,'(40X,F5.3)') beta
	READ(13,'(40X,F5.3)') a50mat
	READ(13,'(40X,F5.3)') a95mat
 


 !    Characteristics of spawning time/space
      READ(13,*)
 	READ(13,'(40X,F7.5)') YoungBetaVar
 	READ(13,'(40X,F7.5)') OldBetaVar
    !  READ(13,'(40X,I4)') ConstantBeta
 	READ(13,'(40X,F7.5)') YoungBetaMean
      READ(13,'(40X,F7.5)') OldBetaMean
      READ(13,'(40X,I4)') OntoEggs
      READ(13,'(40X,I4)') NOntoGroup
      READ(13,'(40X,50(I2,1x))') OntoGroupsVec 

!     Fishing: Selectivity parameters
      READ(13,*)
      READ(13,'(40X,I4)') DoubleNormal
      READ(13,'(40X,F5.3)') a50sel
 	READ(13,'(40X,F5.3)') a95sel
 	READ(13,'(40X,F5.3)') Sfull
 	READ(13,'(40X,F8.3)') leftvar
 	READ(13,'(40X,F8.3)') rightvar

!     Specifications for 2-box larval movement parameters
      READ(13,*)
      READ(13,*)
      READ(13,'(40X,I4)') TwoBoxLarvalMovement
      READ(13,'(40X,F5.3)') MPAsize
 	READ(13,'(40X,F8.3)') LMoveOut
      READ(13,'(40X,I4)') PostDDLmove

!     Specifications for closed areas
      READ(13,*)
      READ(13,*)
      READ(13,'(40X,I2)') NClosedAreas
      READ(13,'(40X,I4)') MPAYear
      READ(13,'(40X,50(I2,1x))') (MinClosed(Ibox),Ibox=1,NClosedAreas)
      READ(13,'(40X,50(I2,1x))') (MaxClosed(Ibox),Ibox=1,NClosedAreas)

!     Means of spawning beta functions
!      READ(13,*)
!      READ(13,*)
!      READ(13,'(5X,7(F8.7,1x))') BetaMeans(1:7)
!      READ(13,'(5X,7(F8.7,1x))') BetaMeans(8:14)
!      READ(13,'(5X,7(F8.7,1x))') BetaMeans(15:21)
!      READ(13,'(5X,7(F8.7,1x))') BetaMeans(22:28)
!      READ(13,'(5X,7(F8.7,1x))') BetaMeans(29:35)

!      READ(13,'(5X,6(F10.9,1x))') BetaMeans(0:5)
!      READ(13,'(5X,6(F10.9,1x))') BetaMeans(6:11)
!      READ(13,'(5X,6(F10.9,1x))') BetaMeans(12:17)
!      READ(13,'(5X,6(F10.9,1x))') BetaMeans(18:23)
!      READ(13,'(5X,6(F10.9,1x))') BetaMeans(24:29)
!      READ(13,'(5X,6(F10.9,1x))') BetaMeans(30:35)

!      READ(13,'(5X,5(F10.9,1x))') BetaMeans(30:34)
!      READ(13,'(5X,5(F10.9,1x))') BetaMeans(35)

      READ(13,*)
      READ(13,*)
      READ(13,'(40X,I4)') NsigmaR
    !  READ(13,'(40X,I4)') Nrho
    !  READ(13,'(40X,F5.3)') Crho
      READ(13,'(40X,F5.3)') TargetCV
      READ(13,'(40X,F10.8)') LowTolerance
      READ(13,'(40X,F10.8)') HighTolerance 
      READ(13,'(40X,I4)') MaxIterations

      READ(13,*)
      READ(13,*)
      DO 71000 IsigmaR= 1,NsigmaR
        READ(13,*) SigmaRvec(IsigmaR)
71000 CONTINUE 
!
!      READ(13,*)
!      READ(13,*)
!      DO 72000 Irho= 1,Nrho
!       READ(13,*) rhovec(Irho)
!72000 CONTINUE

      READ(13,*)
      READ(13,'(40X,I4)') NFpar
      READ(13,*)
      DO 71011 IFpar= 1,NFpar
        READ(13,*) Fvec(IFpar)
71011 CONTINUE 
     
      READ(13,*)
      READ(13,*)

      READ(13,'(40X,I4)') MAEctl
      READ(13,'(40X,I4)') Findh
      READ(13,'(40X,F10.7)') hdist
      READ(13,'(40X,F20.15)') MeanDays
      READ(13,'(40X,F5.3)') SigmaD
      READ(13,'(40X,I2)') Berkeley
      READ(13,'(40X,F5.3)') DiffPar
      READ(13,'(40X,I10)') SRsims
      READ(13,'(40X,I10)') Nbins

      READ(13,*)

      READ(13,'(40X,I2)') NDays
      READ(13,*)
      DO 71012 Idays= 1,NDays
        READ(13,*) DaysVec(Idays)
71012 CONTINUE     
      
      READ(13,*)
      READ(13,'(40X,I2)') FindFmultiply
      READ(13,'(40X,I2)') NFmultiply
      DO 71013 Ifmult = 1,NFmultiply
        READ(13,*) FmultiplyVec(Ifmult)
71013 CONTINUE
      
      IF (FindFmultiply.EQ.0) THEN
       fmultiply = FmultiplyVec(1)
      ENDIF

      READ(13,*)
      READ(13,'(40X,I4)') SRsimsNFpar
      DO 71014 IFpar = 1,SRsimsNFpar
        READ(13,*) SRFvec(IFpar)
71014 CONTINUE

      READ(13,*)
      READ(13,'(40X,I4)') NhFpar
      READ(13,*)
      DO 71015 IFpar = 1,NhFpar
        READ(13,*) FindhFvec(IFpar)
71015 CONTINUE

      FindhFvecOriginal = FindhFvec  !Save this info for later

      IF (PEorMAE.EQ.2) THEN
      fmultiply = 1.0
      ENDIF

      !Make hparVec based on fact that h producing effective steepness will be larger than hpar.
      hparVec = 0.0
      hparVec(1) = 0.2000
      !hdist = 0.01
      Nhpar = nint((1.0-0.2)/hdist + 1) !needs to be calculated using math in the future.
      hparOriginal = hpar
      DO 5554 Ihpar = 2,Nhpar
         hparVec(Ihpar) = hparVec(Ihpar-1) + hdist
5554  CONTINUE

      CLOSE(13)
      RETURN
      END SUBROUTINE ReadSpec
!=================================================================



!==================================================================
! Using COMMON.FOR so need to include a FUNK function
      REAL*8 FUNCTION FUNK(X)
!
! This function is the objective function that projects forward to see how close depletion is
! to what is specified as initial depletion or final depletion
	USE GlobalVars
      IMPLICIT NONE
      
!     Global variables
      REAL*8 X(100)
!
!     Local variables
!      REAL*8 Kvals(MaxSpec),CatPred(MaxSpec),Tmp,TheEffort(10),Penalty1

!
!     Reset FUNK
      FUNK = 1.0E20 
      FUNK = 2*X(1)

      RETURN
      END FUNCTION FUNK
!
! ===========================================================================
! 

! Age structured functions
!==========================================================
      SUBROUTINE Basics
    
      USE GlobalVars
      IMPLICIT NONE
!   Local variables
      INTEGER :: Iage

       DO 100 Iage = 1,Nage
C     Find Selectivity
        IF (DoubleNormal.EQ.0) THEN
          Selectivity(Iage) = 1/(1+EXP(-(Log(19.0)*(Iage-a50sel))/
     +     (a95sel-a50sel)))
        ELSE
          IF (Iage.LE.Sfull) THEN
            Selectivity(Iage) = EXP(-((Iage-Sfull)**2)/leftvar)
          ELSE
            Selectivity(Iage) = EXP(-((Iage-Sfull)**2)/rightvar)
          ENDIF
        ENDIF

C     Find Length at age (von bert)
        LengthAtAge(Iage) = Linf*(1-EXP(-kpar*(Iage-t0par)))
        LengthAtAge(Iage) = max(0.0,LengthAtAge(Iage))

C     Weight at age (allometric)
       WeightAtAge(Iage) = max(0.0,alpha*(LengthAtAge(Iage)**beta))

C     Maturity at age (Logistic)
        MaturityAtAge(Iage) = 1/(1+exp(-(log(19.0))*(Iage-a50mat)/
     +   (a95mat-a50mat)))

C     Fecundity at age
        FecundityAtAge(Iage) = MaturityAtAge(Iage)*WeightAtAge(Iage)
100    CONTINUE

	RETURN
	END SUBROUTINE Basics
!==========================================================

!==========================================================
      SUBROUTINE EqmPop
    
      USE GlobalVars
      IMPLICIT NONE
!   Local variables
      INTEGER :: Iage,Ibox

C     Local variables
!      REAL(8),DIMENSION(Nbox) :: Rinf  !Moved to a global definition
      REAL(8),DIMENSION(Nbox) :: Q1,Q1Unfished
      REAL(8),DIMENSION(Nbox) :: Q2,Q2Unfished
      REAL(8),DIMENSION(Nbox) :: Q3
      REAL(8),DIMENSION(Nbox) :: Q4
!      REAL(8),DIMENSION(Nage,Nbox) :: FishedNperRecruit
!      REAL(8),DIMENSION(Nage,Nbox) :: UnfishedNperRecruit      
      REAL(8) :: Wrpar,Fecrpar

      Wrpar = WeightAtAge(1)
      Fecrpar = FecundityAtAge(1)
      Q2 = 0
      Q2Unfished = 0
      Q1 = 0
      Q1Unfished = 0
      FishedNperRecruit = 0
      UnfishedNperRecruit = 0
!     This helps with the nasty part of Q2
      DO 203 Ibox = 1,Nbox
        FishedNperRecruit(1,Ibox) = 1.0
        UnfishedNperRecruit(1,Ibox) = 1.0
        DO 204 Iage = 2,Nage-1
          FishedNperRecruit(Iage,Ibox) = FishedNperRecruit(Iage-1,Ibox)*
     +     EXP(-(Closed(Ibox,MinPYr)*Fpar*Selectivity(Iage-1) + Mpar))
          UnfishedNperRecruit(Iage,Ibox) = 
     +      UnfishedNperRecruit(Iage-1,Ibox)*EXP(-Mpar)
204     CONTINUE
       FishedNperRecruit(Nage,Ibox) = FishedNperRecruit(Nage-1,Ibox)*
     +   EXP(-(Closed(Ibox,MinPYr)*Fpar*Selectivity(Nage-1)+Mpar))/
     +   (1-EXP(-(Closed(Ibox,MinPYr)*Fpar*Selectivity(Nage)+Mpar)))
       UnfishedNperRecruit(Nage,Ibox) = UnfishedNperRecruit(Nage-1,Ibox)
     +   *EXP(-Mpar)/(1-EXP(-Mpar))
       
        DO 205 Iage = 1,Nage  !I think maybe this should be Nage-1: come back to this!
          Q2(Ibox) = Q2(Ibox) + 0.5*(FecundityAtAge(Iage)*
     +     FishedNperRecruit(Iage,Ibox))
          Q2Unfished(Ibox) = Q2Unfished(Ibox) + 
     +     0.5*(FecundityAtAge(iage)*
     +      UnfishedNperRecruit(Iage,Ibox))
          Q1(Ibox) = Q1(Ibox) + WeightAtAge(Iage)*
     +      FishedNperRecruit(Iage,Ibox)
! Q1Unfished is not used for anything.
 !         Q1Unfished(Ibox) = Q1Unfished(Ibox) + WeightAtAge(Iage)*
 !    +      UnfishedNperRecruit(Iage,Ibox)

205     CONTINUE
! Calculate Rinf, Binf, Einf, etc. directly without Q1,Q3,etc.
      Rinf(Ibox) = (((4*hpar*Q2(Ibox))/Q2Unfished(Ibox))-(1-hpar))*
     +  (1/((5*hpar)-1))*((R0vec(Ibox)*Q2Unfished(Ibox))/Q2(Ibox))
      Rinf(Ibox) = max(0.0000001,Rinf(Ibox))
203    CONTINUE

!     Any Post-DD larval movement is incorporated here into Rinf(Ibox)
        IF (TwoBoxLarvalMovement.EQ.1.AND.PostDDLmove.EQ.1) THEN
          CALL TwoBoxMovement(Rinf,MinPYr)
          Rinf = WasMoved
        ENDIF

      DO 344 Ibox = 1,Nbox
!     Any pre-DD larval movement is incorporated into Q2, I think.
      Einf(Ibox) = Rinf(Ibox)*Q2(Ibox)
      Binf(Ibox) = Rinf(Ibox)*Q1(Ibox)
      NOneAreaEqm(1:Nage,Ibox) = Rinf(Ibox)*
     + FishedNperRecruit(1:Nage,Ibox)
!203    CONTINUE
344   CONTINUE

      	RETURN
	END SUBROUTINE EqmPop
!==========================================================

!==========================================================
      SUBROUTINE UnfishedPop
    
      USE GlobalVars
      IMPLICIT NONE
!   Local variables
      INTEGER :: Iage,Ibox,Ispawn,Ibin
      REAL(8) :: RealNbox
      !REAL(8),DIMENSION(Nbins) :: SSBbinVec
      !Distribute R0 proportionally across the boxes
      R0vec = 0

       IF (Nbox.EQ.2) THEN
       DO 8759 Ibox = 1,Nbox
        IF (Ibox.EQ.MinClosed(1)) THEN
        !R0vec(Ibox) = r0/RealNbox
          R0vec(Ibox) = r0*MPASize
        ELSE
          R0vec(Ibox) = r0*(1-MPASize)
        ENDIF
8759   CONTINUE
       ELSE
        R0vec = r0
       ENDIF
      !Assumes that R0vec has been defined
      DO 200 Ibox = 1,Nbox
        NOneArea0(1,Ibox) = R0vec(Ibox)
        DO 201 Iage = 2,Nage-1
          NOneArea0(Iage,Ibox) = NOneArea0(Iage-1,Ibox)*EXP(-Mpar)
201     CONTINUE
        NOneArea0(Nage,Ibox) = NOneArea0(Iage-1,Ibox)*EXP(-Mpar)/
     +    (1-exp(-Mpar))
        
200   CONTINUE

C     Calculate Spawning stock biomass and total biomass at unfished
      SSB0vec = 0
      B0vec = 0
      N0 = 0
      SSBbyArea0 = 0
      DO 202 Ibox = 1,Nbox
        DO 202 Iage = 1,Nage
        SSB0vec(Ibox) = SSB0vec(Ibox) + NOneArea0(Iage,Ibox)*
     +    FecundityAtAge(Iage)*0.5 !equal sex ratios

        B0vec(Ibox) = B0vec(Ibox) + NOneArea0(Iage,Ibox)*
     +   WeightAtAge(Iage)

          DO 202 Ispawn = 1,Nspawnarea
   !  #-----------------------------------------------------
               N0(Iage,Ispawn,Ibox) = NOneArea0(Iage,Ibox)*
     +           P(Iage,Ispawn)
              SSBbyArea0(Ispawn,Ibox) = SSBbyArea0(Ispawn,Ibox) + 
     +          0.5*N0(Iage,Ispawn,Ibox)*FecundityAtAge(Iage)
202   CONTINUE
    
      DO 203 Ibox = 1,Nbox
          DO 203 Ispawn = 1,Nspawnarea
            Qvec0(Ispawn,Ibox) = SSBbyArea0(Ispawn,Ibox)/SSB0Vec(Ibox)
            R0byArea(Ispawn,Ibox) = R0vec(Ibox)*Qvec0(Ispawn,Ibox) !For local density dependent mortality scenarios
203   CONTINUE


      SSB0nospace = 0.0
      !Bin SSB values:
      DO 8899 Ibox = 1,Nbox
        SSB0nospace = SSB0nospace + SSB0vec(Ibox)
8899  CONTINUE
      !Bin the SSB data:
      !Make bins:
      SSBbinVec(1) = 0.0
      DO 5413 Ibin = 2,Nbins
       SSBbinVec(Ibin) = SSBbinVec(Ibin-1) + 2*SSB0nospace/Nbins
5413  CONTINUE      
      
       
	RETURN
	END SUBROUTINE UnfishedPop
!==========================================================

!==========================================================
      SUBROUTINE BetaStuff
    
      USE GlobalVars
      IMPLICIT NONE

!   Local variables
      INTEGER :: Iage,Ispawn
      REAL(8),DIMENSION(1:Nage,Nspawnarea) :: Pprelim
      REAL(8),DIMENSION(1:Nage) :: alphaBfunc,betaBfunc,SumPrelim
      REAL(8),DIMENSION(1:Nspawnarea) :: SpawnVec

!    Figure out vector of beta variances
      BetaVar(1) = YoungBetaVar
      DO 7832 Iage = 2,Nage
        BetaVar(Iage) = BetaVar(Iage-1) + (OldBetaVar-YoungBetaVar)/
     +     real(Nage)
7832  CONTINUE

!    Figure out vector of beta means (do this in betastuff as well,but in terms of age)
!      BetaMeans(0) = YoungBetaMean
      BetaMeans(1) = YoungBetaMean
      DO 7838 Iage = 2,Nage
        BetaMeans(Iage) = BetaMeans(Iage-1) + 
     +   (OldBetaMean-YoungBetaMean)/real(Nage)
7838  CONTINUE


      SpawnVec(1) = 1.0
      DO 657 ISpawn = 2,Nspawnarea
        SpawnVec(Ispawn) = 1.0 + SpawnVec(ISpawn-1)
657   CONTINUE

 !     SpawnVec(1) = 0.999/REAL(Nspawnarea)
 !     DO 657 ISpawn = 2,Nspawnarea
 !       SpawnVec(Ispawn) = 1.0/REAL(Nspawnarea) + SpawnVec(ISpawn-1)
!657   CONTINUE
      
!      DO 658 ISpawn = 1,Nspawnarea
!        SpawnVec(Ispawn) = SpawnVec(ISpawn)/REAL(Nspawnarea)
!658   CONTINUE
      !Maybe we want something like this
      !SpawnVec = SpawnVec/REAL(Nspawnarea)

      !Calculate betameans without having to read it in all the time
      !read in FirstBetaMean and lastbetamean?


!     Calculate alpha and beta for the Beta functions from means and variance
      SumPrelim = 0
      Pprelim = 0
      P = 0
      CheckP = 0
      DO 300 Iage = 1,Nage
        alphaBfunc(Iage) = (BetaMeans(Iage)**2-BetaMeans(Iage)**3-
     +   BetaMeans(Iage)*BetaVar(Iage))
        alphaBfunc(Iage) = alphaBfunc(Iage)/BetaVar(Iage)

      betaBfunc(Iage) = (BetaMeans(Iage)-2*BetaMeans(Iage)**2+
     + BetaMeans(Iage)**3-BetaVar(Iage)+BetaMeans(Iage)*BetaVar(Iage))/
     + BetaVar(Iage)

        DO 301 Ispawn = 1,Nspawnarea
          Pprelim(Iage,Ispawn) = (SpawnVec(Ispawn)**(alphaBfunc(Iage)-1)
     +     *(Nspawnarea-SpawnVec(Ispawn)+1)**(betaBfunc(Iage)-1))
          SumPrelim(Iage) = SumPrelim(Iage) + Pprelim(Iage,Ispawn)
301        CONTINUE
          DO 302 Ispawn = 1,Nspawnarea
           P(Iage,Ispawn) = Pprelim(Iage,Ispawn)/SumPrelim(Iage)
           CheckP(Iage) = CheckP(Iage) + P(Iage,Ispawn)
302       CONTINUE
300   CONTINUE
      
      !Incorporate P into population dynamics in other functions
	RETURN
	END SUBROUTINE BetaStuff
!==========================================================

!For generating random errors for various sigmaRs and rhos:
C-------------------------------------------------------------------------------------------	
      SUBROUTINE CORRMODEL2(ThesigmaR,Therho)

      USE GlobalVars
C	USE mpamodelsubs
      IMPLICIT NONE

	INTEGER :: Icount
      REAL(8),INTENT(IN) :: ThesigmaR,Therho
	CHARACTER(len = 4):: rhochar
	CHARACTER(len = 4):: sigmaRchar
	CHARACTER(len = 29) :: BigChar
C	REAL*8 myrho,mysigmaR
      REAL(8) :: myrho
C	Open the random seed file
      OPEN(UNIT=2,FILE="RANDOM.NUM")
C	UNIT = 3 is INPUT.DAT, opened in READSPEC()
C      OPEN(UNIT=4,FILE="RANDOMDEVS.OUT")
      READ(2,*) SEED1,SEED2
	!ISEEDO = SEED1
C	ISEED0 = 2+2
C	Open the input files
C      CALL ReadSpecCORR()
C	CALL ReadSigmas()
C	CALL ReadRhos()

	Icount = 3
!	DO 12000 ThesigmaR = 1,NsigmaR
C	Find the random errors
      !Use the same initialization for each sigmaR by resetting ISEEDO to ISEED1 (ISEEDO gets reset to 1 when GenNormalError is run)
      ISEEDO = SEED1
	CALL GenNormalErrorBI(ThesigmaR)
!	  DO 12000 Therho = 1,Nrho
		Icount = Icount + 1
C	Find the correlated random errors
	CALL GenCorrError(Therho)

C	Write out results
C	Open file

C	Convert numbers to characters
	myrho = Therho
	mysigmaR = ThesigmaR

	Write(rhochar,'(F4.2)') myrho
	Write(sigmaRchar,'(F4.2)') mysigmaR
	BigChar = 'CorrDevsroe' // rhochar // 'SigmaR' //
     +	 sigmaRchar //'.OUT'
      OPEN(UNIT=Icount,FILE=BigChar)
	CALL SUMMARY2CORR(Icount)
	CLOSE(Icount)
12000 CONTINUE



C	Make the control comparison file (sigmaR = 0.6 and rho = 1)
      ISEEDO = SEED1
	!CALL GenNormalError(CsigmaR)
      !CALL GenCorrError(Crho)

	!Write(rhochar,'(F4.2)') Crho
	!Write(sigmaRchar,'(F4.2)') CsigmaR
	!BigChar = 'CorrDevsroe' // rhochar // 'SigmaR' //
      !+	 sigmaRchar //'.OUT'
      !OPEN(UNIT=10000,FILE=BigChar)
	!CALL SUMMARY2CORR(10000)
	!CLOSE(10000)

      RETURN
      END SUBROUTINE CORRMODEL2
C---------------------------------------------------------------------------------------------

CC ===========================================================================
CC
      SUBROUTINE GenNormalError(insigmaR)
C
C This subroutine generates the density ratio survey errors that will be needed for the
C current simulation

C	Call this for the years for which density ratio sampling happens and only
C	when SuperSim = 1
C
      USE GlobalVars
      IMPLICIT NONE
C
C     Local variables
      INTEGER Itime,Isim,Iarea,IsigmaR
	REAL(8), INTENT(IN) :: insigmaR
      REAL*8 XNORM
      EXTERNAL XNORM
C	SigmaR = 0.6
C	Nsim = 20
C	Ntime = 10
C	Nspawnarea = 15
C	
C	DO 20000 IsigmaR = 1,NsigmaR
C	SigmaR = SigmaRvec(IsigmaR)
	MyError = 0
C     Generate the recruitment error terms
      DO 10000 Iarea = 1,Nspawnarea
      ! DO 10000 Itime = 1,Ntime
          !could make this time loop more efficient if this process takes too long
       DO 10000 Itime = MinPYr,MaxPYr
	   DO 10000 Isim = 1,Nsim
           MyError(Itime,Iarea,Isim) = XNORM(2,0.0d0,insigmaR,ISEEDO)
10000 CONTINUE  
	     
C20000	CONTINUE     
C
C	Calculate correlated error terms
	


      RETURN
      END SUBROUTINE GenNormalError
C============================================================================
	SUBROUTINE GenCorrError(inrho)
C	Calculate correlated random errors
      USE GlobalVars
      IMPLICIT NONE

C	Local Variables
      INTEGER Itime,Isim,Iarea,IsigmaR
	REAL(8), INTENT(IN) :: inrho


	CorrError = 0
    !   DO 10000 Itime = 1,Ntime
    !could make this time loop more efficient if this process takes too long
       DO 10000 Itime = MinPYr,MaxPYr
	   DO 10000 Isim = 1,Nsim
		 CorrError(Itime,1,Isim) = MyError(Itime,1,Isim)
           DO 10000 Iarea = 2,Nspawnarea
           CorrError(Itime,Iarea,Isim) = 
     +    	 inrho*CorrError(Itime,Iarea-1,Isim)+ 
     +         SQRT(1-inrho*inrho)*MyError(Itime,Iarea,Isim)
10000 CONTINUE  

	END SUBROUTINE GenCorrError
C ===========================================================================

C============================================================================
	SUBROUTINE SUMMARY2CORR(Icount)
      USE GlobalVars
	IMPLICIT NONE

	INTEGER, INTENT(IN):: Icount
	INTEGER Iarea,Isim,Itime


      
	DO 14000 Isim = 1,Nsim
	  DO 14000 Itime = 1,Ntime
	    DO 14000 Iarea = 1,Nspawnarea
C	WRITE(4,'(F20.17,3(1x,I4))') CorrError(Itime,Iarea,Isim),
C     +	Iarea,Isim,Itime
	WRITE(Icount,600) CorrError(Itime,Iarea,Isim),Iarea,Isim,
     +	Itime
C	WRITE(*,*) 2
14000 CONTINUE
C601   FORMAT(1x,I4,10(3(1x,F8.2)))
600   FORMAT(1x,F30.17,',',I4,',',I4,',',I4)
	RETURN
	END SUBROUTINE SUMMARY2CORR

C ============================================================================

!==========================================================
      SUBROUTINE FirstTime(MinYr)
    
      USE GlobalVars
      IMPLICIT NONE

!   Local variables
      INTEGER :: Iage,Ispawn,Ibox
      INTEGER,INTENT(IN) :: MinYr
      REAL(8),DIMENSION(Nage,Nbox) :: Ncheck

    ! Setup initial numbers at age and spawning stock biomass
      SSB = 0
      SSBbyArea = 0
      Eggs = 0
      EggsbyArea = 0
      B = 0
      BbyArea = 0
!       IF (EnvPreDD.EQ.0.OR.BurnInTime.GT.0) THEN
        DO 400 Iage = 1,Nage
          DO 400 Ibox = 1,Nbox
            NCompact(Iage,Ibox,MinYr) = NOneAreaEqm(Iage,Ibox)
            !Spawning stock biomass, biomass, eggs, what else? Step through Rcode
            SSB(Ibox,MinYr) = SSB(Ibox,MinYr) + 
     +         0.5*NCompact(Iage,Ibox,MinYr)*FecundityAtAge(Iage)
            B(Ibox,MinYr) = B(Ibox,MinYr) + NCompact(Iage,Ibox,MinYr)*
     +        WeightAtAge(Iage)
            !Add in First time recruits here - is equal to Rinf.
            TotRecruits(Ibox,MinYr,:) = Rinf(Ibox)
            DO 400 Ispawn = 1,Nspawnarea
          N(Iage,Ispawn,Ibox,MinYr) = NOneAreaEqm(Iage,Ibox)*
     +      P(Iage,Ispawn)
         SSBbyArea(Ispawn,Ibox,MinYr) = SSBbyArea(Ispawn,Ibox,MinYr) + 
     +     0.5*N(Iage,Ispawn,Ibox,MinYr)*FecundityAtAge(Iage)
400     CONTINUE
!      ENDIF
      NonSpatialRecruits(MinYr,:) = 0
      DO 444 Ibox = 1,Nbox
        NonSpatialRecruits(MinYr,:) = NonSpatialRecruits(MinYr,:) + 
     +    TotRecruits(Ibox,MinYr,:)
444   CONTINUE

      !Qvec is the proportion of eggs in each area prior to Density Dependent mortality
 !     DO 411 Ispawn = 1,NspawnArea
 !       DO 411 Ibox = 1,Nbox
 !     Qvec(Ispawn,Ibox,MinYr) = SSBbyArea(Ispawn,Ibox,MinYr)/
 !    + SSB(Ibox,MinYr)
!411   CONTINUE

      RETURN
      END SUBROUTINE FirstTime
!==========================================================

!==========================================================
      SUBROUTINE Recruitment(Year,Sim,TheSstarve)
!Stuff that goes into the recruitment function in various ways
  !      Eggs = SSB !For this case.
  !      EggsbyArea = SSBbyArea

  !      !Qvec is the proportion of eggs in each area prior to Density Dependent mortality
  !    DO 411 Ispawn = 1,NspawnArea
  !      DO 411 Ibox = 1,Nbox
  !    Qvec(Ispawn,Ibox,MinYr) = EggsbyArea(Ispawn,Ibox,MinYr)/
  !   + Eggs(Ibox,MinYr)
!411   CONTINUE
      USE GlobalVars
      IMPLICIT NONE

!   Local variables
      INTEGER,INTENT(IN) :: Year,Sim
      INTEGER :: Ibox,Ispawn
      REAL(8),DIMENSION(1:Nbox) :: TheSstarve
      
      TotRecruits(:,Year,Sim) = 0
      IF (EnvPreDD.EQ.0.AND.GLOBALDD.EQ.0) THEN
        !Eggs and SSB are equivalent in this scenario (where environmental var. applies after DD)
        DO 111 Ibox = 1,Nbox
         Eggs(Ibox,Year-1) = SSB(Ibox,Year-1) 
         DO 111 Ispawn = 1,Nspawnarea
           EggsbyArea(Ispawn,Ibox,Year-1)=SSBbyArea(Ispawn,Ibox,Year-1)
          !Apply local density dependent mortality, then environmental variation to each spawning time/space
          Recruits(Ispawn,Ibox,Year) = R0byArea(Ispawn,Ibox)*
     +     (4*hpar*(EggsbyArea(Ispawn,Ibox,Year-1)/
     +     SSBbyArea0(Ispawn,Ibox)))/((1-hpar)+(5*hpar-1)*
     +     (EggsbyArea(Ispawn,Ibox,Year-1)/
     +     SSBbyArea0(Ispawn,Ibox)))
         Recruits(Ispawn,Ibox,Year) = Recruits(Ispawn,Ibox,Year)*
     +    exp(CorrError(Year,Ispawn,Sim)-(sigmaR**2)/2)*TheSstarve(Ibox)
     +*fmultiply
         TotRecruits(Ibox,Year,Sim) = TotRecruits(Ibox,Year,Sim) + 
     +     Recruits(Ispawn,Ibox,Year)
 111   CONTINUE
       ELSEIF (EnvPreDD.EQ.0.AND.GLOBALDD.EQ.1) THEN
          PreRecruits = 0
          TotRecruits(:,Year,Sim) = 0
          DO 112 Ibox = 1,Nbox
         Eggs(Ibox,Year-1) = SSB(Ibox,Year-1) 
          !Apply global density dependent mortality (env. var. is later)
          PreRecruits(Ibox) = R0vec(Ibox)*
     +     (4*hpar*(Eggs(Ibox,Year-1)/
     +     SSB0vec(Ibox)))/((1-hpar)+(5*hpar-1)*
     +     (Eggs(Ibox,Year-1)/
     +     SSB0vec(Ibox)))
           DO 112 Ispawn = 1,Nspawnarea
            !Split out pre-recruits into spawning times/spaces
            Qvec(Ispawn,Ibox,Year-1) = SSBbyArea(Ispawn,Ibox,Year-1)/
     +        SSB(Ibox,Year-1)
            Recruits(Ispawn,Ibox,Year) = PreRecruits(Ibox)*
     +        Qvec(Ispawn,Ibox,Year-1)
            Recruits(Ispawn,Ibox,Year) = Recruits(Ispawn,Ibox,Year)*
     +    exp(CorrError(Year,Ispawn,Sim)-(sigmaR**2)/2)*TheSstarve(Ibox)
     +*fmultiply
            TotRecruits(Ibox,Year,Sim) = TotRecruits(Ibox,Year,Sim) + 
     +        Recruits(Ispawn,Ibox,Year)
 112   CONTINUE       
      ELSEIF (EnvPreDD.EQ.1.AND.GLOBALDD.EQ.1) THEN
      !Define eggs locally by applying CorrError to SSBbyArea
      Eggs(1:Nbox,Year-1) = 0
      DO 212 Ibox = 1,Nbox
        DO 213 Ispawn = 1,Nspawnarea
      EggsbyArea(Ispawn,Ibox,Year-1) = SSBbyArea(Ispawn,Ibox,Year-1)*
     + exp(CorrError(Year,Ispawn,Sim)-(sigmaR**2)/2)*TheSstarve(Ibox)
     +*fmultiply
      Eggs(Ibox,Year-1) = Eggs(Ibox,Year-1)+
     + EggsbyArea(Ispawn,Ibox,Year-1)
213     CONTINUE
       !This will change either here or in initialization stages
       !SSB0vec will be EffAvgSSB0vec (right??)
       !Note that Recruits in each spawning area are not defined here.
       TotRecruits(Ibox,Year,Sim) = R0vec(Ibox)*
     +     (4*hpar*(Eggs(Ibox,Year-1)/
     +     SSB0vec(Ibox)))/((1-hpar)+(5*hpar-1)*
     +     (Eggs(Ibox,Year-1)/
     +     SSB0vec(Ibox)))
212   CONTINUE
      ELSEIF (EnvPreDD.EQ.1.AND.GLOBALDD.EQ.0) THEN
      TotRecruits(:,Year,Sim) = 0
      DO 214 Ibox = 1,Nbox
        DO 214 Ispawn = 1,Nspawnarea
      EggsbyArea(Ispawn,Ibox,Year-1) = SSBbyArea(Ispawn,Ibox,Year-1)*
     + exp(CorrError(Year,Ispawn,Sim)-(sigmaR**2)/2)*TheSstarve(Ibox)
     +*fmultiply
      Eggs(Ibox,Year-1) = Eggs(Ibox,Year-1)+
     + EggsbyArea(Ispawn,Ibox,Year-1)

        Recruits(Ispawn,Ibox,Year) = R0byArea(Ispawn,Ibox)*
     +     (4*hpar*(EggsbyArea(Ispawn,Ibox,Year-1)/
     +     SSBbyArea0(Ispawn,Ibox)))/((1-hpar)+(5*hpar-1)*
     +     (EggsbyArea(Ispawn,Ibox,Year-1)/
     +     SSBbyArea0(Ispawn,Ibox)))
        TotRecruits(Ibox,Year,Sim) = TotRecruits(Ibox,Year,Sim) + 
     +        Recruits(Ispawn,Ibox,Year)
214     CONTINUE

      !Still need to add in EnvPreDD = 1 here; local DD ( think I did that right above.
      ENDIF !this will be an elseif eventually, but start with one case

      NonSpatialRecruits(Year,Sim) = 0
      DO 100 Ibox = 1,Nbox
        NonSpatialRecruits(Year,Sim) = NonSpatialRecruits(Year,Sim) + 
     +   TotRecruits(Ibox,Year,Sim)
100   CONTINUE

      !Add an option to move larvae in the two-box model
      IF (TwoBoxLarvalMovement.EQ.1) THEN
        CALL TwoBoxMovement(TotRecruits(:,Year,Sim),Year)
        TotRecruits(:,Year,Sim) = WasMoved        
      ENDIF

      RETURN
      END SUBROUTINE Recruitment
!==========================================================
!==============================================================
      SUBROUTINE PopDy(Year,Sim)
      USE GlobalVars
      IMPLICIT NONE

!   Local variables
      INTEGER, INTENT(IN) :: Year,Sim
      INTEGER :: Ibox, Iage, Ispawn
      
      SSBbyArea(:,:,Year) = 0.0
      SSB(:,Year) = 0.0
      B(:,Year) = 0.0
      TotalCatchBiomass(Year,Sim) = 0.0
      TotBiomass(Year,Sim) = 0.0
      DO 500 Ibox = 1,Nbox
        NCompact(1,Ibox,Year) = TotRecruits(Ibox,Year,Sim)
        DO 501 Iage = 2,Nage-1
          NCompact(Iage,Ibox,Year) = NCompact(Iage-1,Ibox,Year-1)*
     +      exp(-(Mpar+Closed(Ibox,Year-1)*Fpar*Selectivity(Iage-1)))
501     CONTINUE
        NCompact(Nage,Ibox,Year) = NCompact(Nage-1,Ibox,Year-1)*
     +    exp(-(Mpar+Closed(Ibox,Year-1)*Fpar*Selectivity(Nage-1))) + 
     +    NCompact(Nage,Ibox,Year-1)*exp(-(Mpar+Closed(Ibox,Year-1)*
     +       Fpar*Selectivity(Nage)))

          !Now define SSB and B
          DO 500 Iage = 1,Nage
            !Calculate spawning stock biomass and total biomass
            SSB(Ibox,Year) = SSB(Ibox,Year) + 
     +        0.5*NCompact(Iage,Ibox,Year)*FecundityAtAge(Iage)
            B(Ibox,Year) = B(Ibox,Year) + NCompact(Iage,Ibox,Year)*
     +         WeightAtAge(Iage)
            !Calculate catches
            CatchAtAge(Iage,Ibox,Year) = (Closed(Ibox,Year)*Fpar*
     +        Selectivity(Iage)/
     +        (Closed(Ibox,Year)*Fpar*Selectivity(Iage)+Mpar))*
     +        NCompact(Iage,Ibox,Year)*
     +        (1-exp(-(Mpar+Closed(Ibox,Year)*Fpar*Selectivity(Iage))))
           
           TotalCatchBiomass(Year,Sim) = TotalCatchBiomass(Year,Sim) + 
     +       CatchAtAge(Iage,Ibox,Year)*WeightAtAge(Iage)
           
            ! Now define spatial dynamics
            DO 500 Ispawn = 1,Nspawnarea
              N(Iage,Ispawn,Ibox,Year) = NCompact(Iage,Ibox,Year)*
     +          P(Iage,Ispawn)
          SSBbyArea(Ispawn,Ibox,Year) = SSBbyArea(Ispawn,Ibox,Year) + 
     +     0.5*N(Iage,Ispawn,Ibox,Year)*FecundityAtAge(Iage)
400     CONTINUE
500   CONTINUE
            TotSSB(Year,Sim) = 0.0
            DO 511 Ibox = 1,Nbox
              TotSSB(Year,Sim) = TotSSB(Year,Sim) + SSB(Ibox,Year)
              TotBiomass(Year,Sim) = TotBiomass(Year,Sim) + B(Ibox,Year)
511         CONTINUE
      RETURN
      END SUBROUTINE PopDy
!=================================================================

!=================================================================
      SUBROUTINE Summary2(Year,Sim)
      USE GlobalVars
      IMPLICIT NONE

!   Local variables
      INTEGER, INTENT(IN) :: Year,Sim
      INTEGER :: Ibox,Ispawn,Iage,i,Ibin
      REAL(8) :: NonSpatialEggs,NonSpatialSSB
      REAL(8) :: NonSpatialB
      REAL(8),DIMENSION(Nspawnarea) :: TotEggsbyArea,TotSSBbyArea
      REAL(8) :: SSB0forBins,SSBbinned
 !     CHARACTER*4 NsimChar
 !     CHARACTER*4 NtimeChar
 !     CHARACTER*1 NboxChar
 !     CHARACTER*2 NspawnareaChar
 !     CHARACTER*3 BurnInTimeChar
 !     CHARACTER*2 NFparChar
 !     CHARACTER*3 MoveBurnInChar
 !     CHARACTER*2 NageChar

      !Preliminary format stuff:
            !------------------------------------------------
      !Make some characters so that you can write out header info to output files
  !    Write(NsimChar,'(I4)') Nsim
  !    Write(NtimeChar,'(I4)') Ntime
  !    Write(NboxChar,'(I1)') Nbox
  !    Write(NspawnareaChar,'(I2)') Nspawnarea
  !    Write(BurnInTimeChar,'(I3)') BurnInTime
  !    Write(NFparChar,'(I2)') NFpar
  !    Write(MoveBurnInChar,'(I3)') MoveBurnInTime
  !    Write(NageChar,'(I2)') Nage

      !-----------------------------------------------
 

      !Preliminary calcs for writing things out
      !-----------------------------------------------
  
      NonSpatialEggs = 0
      NonSpatialSSB = 0
      NonSpatialB = 0
      DO 100 Ibox = 1,Nbox
!        Switched to defining NonSpatialRecruits in Recruitment function for use in FindVariance chunk.
!        NonSpatialRecruits = NonSpatialRecruits + 
!     +   TotRecruits(Ibox,Year,Sim)
        NonSpatialEggs = NonSpatialEggs + Eggs(Ibox,Year-1)
        NonSpatialSSB = NonSpatialSSB + SSB(Ibox,Year)
        NonSpatialB = NonSpatialB + B(Ibox,Year)
100   CONTINUE

      TotEggsbyArea = 0
      TotSSBbyArea = 0
      DO 101 Ispawn = 1,Nspawnarea
        DO 101 Ibox = 1,Nbox
      TotEggsbyArea(Ispawn) = TotEggsbyArea(Ispawn)+
     + EggsbyArea(Ispawn,Ibox,Year-1)

      TotSSBbyArea(Ispawn) = TotSSBbyArea(Ispawn) +
     +  SSBbyArea(Ispawn,Ibox,Year)
101   CONTINUE
      

      !NonSpatialSSB
      DO 5415 Ibin = 1,(Nbins-1)
       IF (NonSpatialSSB.GE.SSBbinVec(Ibin).AND.
     +NonSpatialSSB.LT.SSBbinVec(Ibin+1)) THEN
            SSBBinned = SSBbinVec(Ibin)
       ENDIF
5415  CONTINUE
      IF (NonSpatialSSB.GE.SSBbinVec(Nbins)) THEN
        SSBBinned = SSBbinVec(Nbins)
      ENDIF
      !---------------------------------------------------------
      !End Preliminary Calcs
      
      ! Write out recruits
      IF (Sim.EQ.DetSim.AND.Year.EQ.1.AND.Fpar.EQ.0.0) THEN
        WRITE(11,707) Nsim,BurnInTime,Ntime,Nbox,NSpawnArea,NFpar,
     +MoveBurnInTime,Nage
        WRITE(11,*) 'Sim,Year,Fpar,TotalRecruits,RecruitsByBox'
      ENDIF
      WRITE(11,701) Sim,Year,Fpar,NonSpatialRecruits(Year,Sim),
     + (TotRecruits(Ibox,Year,Sim),Ibox = 1,Nbox)

      ! Write out total eggs and EggsbyArea
      IF (Sim.EQ.DetSim.AND.Year.EQ.1.AND.Fpar.EQ.0.0) THEN
        WRITE(12,707) Nsim,BurnInTime,Ntime,Nbox,NSpawnArea,NFpar,
     +MoveBurnInTime,Nage
        WRITE(12,*) 'Sim,Year,Fpar,NonSpatialEggs,EggsByBox'
      ENDIF
      WRITE(12,701) Sim,Year-1,Fpar,NonSpatialEggs,
     + (Eggs(Ibox,Year-1),Ibox = 1,Nbox)

      IF (Sim.EQ.DetSim.AND.Year.EQ.1.AND.Fpar.EQ.0.0) THEN
        WRITE(14,707) Nsim,BurnInTime,Ntime,Nbox,NSpawnArea,NFpar,
     +MoveBurnInTime,Nage
        WRITE(14,*) 'Sim,Year,Fpar,Ispawn,TotEggsbyarea'
      ENDIF
      DO 5548 Ispawn = 1,Nspawnarea
        WRITE(14,702) Sim,Year-1,Fpar,Ispawn,
     +   TotEggsbyArea(Ispawn)
5548  CONTINUE
      !Write out SSB information (SSB and SSBbyArea)
      IF (Sim.EQ.DetSim.AND.Year.EQ.1.AND.Fpar.EQ.0.0) THEN
        WRITE(15,707) Nsim,BurnInTime,Ntime,Nbox,NSpawnArea,NFpar,
     +MoveBurnInTime,Nage
        WRITE(15,*) 'Sim,Year,Fpar,NonSpatialSSB,SSBbinned,SSBbyBox'
      ENDIF
      WRITE(15,701) Sim,Year,Fpar,NonSpatialSSB,SSBbinned,
     +(SSB(Ibox,Year),Ibox=1,Nbox)

      IF (Sim.EQ.DetSim.AND.Year.EQ.1.AND.Fpar.EQ.0.0) THEN
        WRITE(16,707) Nsim,BurnInTime,Ntime,Nbox,NSpawnArea,NFpar,
     +MoveBurnInTime,Nage
        WRITE(16,*) 'Sim,Year,Fpar,Ispawn,TotalSSBbySpawnArea'
      ENDIF
      DO 5549 Ispawn = 1,Nspawnarea
        WRITE(16,702) Sim,Year,Fpar,Ispawn,TotSSBbyArea(Ispawn)
5549  CONTINUE
      !Write out total catches (right now not writing out CatchAtAge, but could if needed)
      IF (Sim.EQ.DetSim.AND.Year.EQ.1.AND.Fpar.EQ.0.0) THEN
        WRITE(17,707) Nsim,BurnInTime,Ntime,Nbox,NSpawnArea,NFpar,
     +MoveBurnInTime,Nage
        WRITE(17,*) 'Sim,Year,Fpar,TotalCatchBiomass'
      ENDIF
      WRITE(17,703) Sim,Year,Fpar,TotalCatchBiomass(Year,Sim)

      
      !Write out total biomass
      IF (Sim.EQ.DetSim.AND.Year.EQ.1.AND.Fpar.EQ.0.0) THEN
        WRITE(17,707) Nsim,BurnInTime,Ntime,Nbox,NSpawnArea,NFpar,
     +MoveBurnInTime,Nage
        WRITE(17,*) 'Sim,Year,Fpar,TotBiomass,BiomassbyBox'
      ENDIF
      WRITE(30,708) Sim,Year,Fpar,TotBiomass(Year,Sim),
     +(B(Ibox,Year),Ibox=1,Nbox)



      !Write out Numbers at age (not by space) for each F for box 1 (used for debugging - change for actual model output)
      IF (Sim.EQ.DetSim.AND.Year.EQ.1.AND.Fpar.EQ.0.0) THEN
        WRITE(18,707) Nsim,BurnInTime,Ntime,Nbox,NSpawnArea,NFpar,
     +MoveBurnInTime,Nage
        WRITE(18,*) 'Sim,Year,Fpar,Iage,NCompactByAge'
      ENDIF
      DO 5550 Iage = 1,Nage
       WRITE(18,704) Sim,Year,Fpar,Iage,NCompact(Iage,1,Year)
5550  CONTINUE
      ! RecruitsByAreaAndBox.OUT (UNIT 22) Sim,Year,F,Box,SpawnArea,Recruits
      IF (Sim.EQ.DetSim.AND.Year.EQ.1.AND.Fpar.EQ.0.0) THEN
        WRITE(22,707) Nsim,BurnInTime,Ntime,Nbox,NSpawnArea,NFpar,
     +MoveBurnInTime,Nage
        WRITE(22,*)'Sim,Year,Fpar,Ibox,Ispawnarea,Recruits'
      ENDIF
      DO 5544 Ibox = 1,Nbox
        DO 5544 Ispawn = 1,NSpawnArea
      WRITE(22,705) Sim,Year,Fpar,Ibox,Ispawn,Recruits(Ispawn,Ibox,Year)
5544  CONTINUE

      !Write out detailed N info.
      IF (NFpar.EQ.7) THEN
      DO 5545 i = 1,NFpar
        IF (Sim.EQ.DetSim.AND.Year.EQ.1.AND.Fpar.EQ.Fvec(i)) THEN
      WRITE(50+i,707) Nsim,BurnInTime,Ntime,Nbox,NSpawnArea,NFpar,
     +MoveBurnInTime,Nage
      WRITE(50+i,*) 'Sim,Year,Fpar,Ibox,Ispawnarea,Iage,N'
        ENDIF
        !Write out NDetailed for just last year of deterministic burn-in
      IF (Sim.EQ.DetSim.AND.Year.EQ.MoveBurnInTime) THEN
       DO 5546 Ibox = 1,Nbox
       DO 5546 Ispawn = 1,NSpawnArea
         DO 5546 Iage = 1,Nage
          IF (Fpar.EQ.Fvec(i)) THEN
         WRITE(50+i,706) Sim,Year,Fpar,Ibox,Ispawn,Iage,
     +   N(Iage,Ispawn,Ibox,Year)
         ENDIF
5546   CONTINUE
      ENDIF
5545  CONTINUE
      ENDIF


      RETURN
700   FORMAT(I4,2(1x,F18.16),1x,I4,1x,F5.3,1x,F18.16)
701   FORMAT(2(I4,','),F5.3,51(',',F18.6))
702   FORMAT(2(I4,','),F5.3,',',I2,',',10(F18.6,','))
703   FORMAT(2(I4,','),F5.3,',',F18.6)
704   FORMAT(2(I4,','),F5.3,',',I2,100(',',F18.6))
705   FORMAT(2(I4,','),F5.3,',',2(I2,','),F18.6)
706   FORMAT(2(I4,','),F5.3,',',3(I2,','),F18.6)
707   FORMAT("Nsim",1x,I4,1x,"BurnInTime",1x,I4,1x,"Ntime",1x,I4,1x,
     +"Nbox",1x,I2,1x,"NSpawnArea",1x,I2,1x,"NFpar",1x,I4,1x,"MoveBI",
     +1x,I4,1x,"Nage ",1x,I2) 

!      WRITE(30,703) Sim,Year,Fpar,TotBiomass(Year,Sim),
!     +(B(Ibox,Year),Ibox=1,Nbox)
708   FORMAT(2(I4,','),F5.3,',',3(F18.6,','))
600   FORMAT(1x,F30.17,',',I4,',',I4,',',I4)
      END SUBROUTINE Summary2
!=================================================================
!=================================================================
      SUBROUTINE CalcRecruitCV(SigmaRindex)
      USE GlobalVars
      IMPLICIT NONE

      !Local variables
    !  INTEGER,INTENT(IN) :: Year,Sim
      INTEGER,INTENT(IN) :: SigmaRIndex
      INTEGER :: Isim,Iyear
      REAL(8),DIMENSION(Nsim) :: SumRecruits,MeanRecruits
      REAL(8),DIMENSION(Nsim) :: SumSDCalcs,SDRecruits,CVRecruits
      REAL(8),DIMENSION((MinPYr+MoveBurnInTime+BurnInTime+1):
     +(MinPYr+MoveBurninTime+BurnInTime+Ntime),
     +  Nsim) :: SDCalcs
      REAL(8) :: sumCVRecruits

      sumCVRecruits = 0
      SDRecruits = 0
      CVRecruits =0
      DO 322 Isim = 1,Nsim
      !Calc interannual CV for each sim then average the sims
      SumRecruits(Isim) = 0
      MeanRecruits(Isim) = 0
     
      
      DO 320 IYear = (MinPYr+MoveBurnInTime +BurnInTime+1),(MinPYr+
     +MoveBurnInTime +BurnInTime+Ntime)
       SumRecruits(ISim) = SumRecruits(ISim) + 
     +  NonSpatialRecruits(IYear,ISim)
320   CONTINUE
       MeanRecruits(ISim) = SumRecruits(ISim)/(Ntime)

       SumSDCalcs(ISim) = 0
       DO 321 IYear = (MinPYr+MoveBurnInTime +BurnInTime+1),(MinPYr+
     +MoveBurnInTime +BurnInTime+Ntime)
         SDCalcs(IYear,ISim) = (NonSpatialRecruits(IYear,ISim)-
     +     MeanRecruits(ISim))**2
         SumSDCalcs(ISim) = SumSDCalcs(ISim) + SDCalcs(IYear,ISim)
321    CONTINUE
         SDRecruits(ISim) = sqrt(SumSDCalcs(ISim)/(Ntime))
         CVRecruits(Isim) = SDRecruits(Isim)/MeanRecruits(Isim)
         sumCVRecruits = sumCVRecruits + CVRecruits(ISim)
322   CONTINUE

      MeanRecruitCV(SigmaRIndex) = sumCVRecruits/Nsim
      RETURN
      END SUBROUTINE CalcRecruitCV
!=================================================================

!=================================================================
      SUBROUTINE DiffsAndNewSigmaRvec
      USE GlobalVars
      IMPLICIT NONE

      !Local variables
      INTEGER :: IsigmaR,count,Iinc
      REAL(8),DIMENSION(NsigmaR) :: Diffs
      REAL(8) :: BestCV
      REAL(8) :: Increment,NewInc

      !Calculate differences between Recruit CVs and target recruit CV
      Diffs = 0
      DO 100 IsigmaR = 1,NsigmaR
        Diffs(IsigmaR) = abs(MeanRecruitCV(IsigmaR)-TargetCV)
100   CONTINUE

      !If none of the sigmaRs is not within the desired tolerance level, make a new SigmaRvec
      
      !MinDiff = min(Diffs)
      MinDiff = 1000
      BestCV = -1
      BestSigmaR = -1
      count = 0
      DO 101 IsigmaR =1,NsigmaR
        IF (Diffs(IsigmaR).EQ.MinDiff) THEN
          count = count+1
        ELSEIF (Diffs(IsigmaR).LT.MinDiff) THEN
          MinDiff = Diffs(IsigmaR)
          BestCV = MeanRecruitCV(IsigmaR)
          BestSigmaR = SigmaRVec(IsigmaR)
        ENDIF
101   CONTINUE

      IF (MinDiff.LE.LowTolerance) THEN
        sigmaR = BestSigmaR
        !FindVariance = 0
      ELSE
        !ReAssign SigmaRVec
        !Keep FindVariance = 1
        !Need a DO WHILE or DO UNTIL or something like that. Can FORTRAN do that?

        !Increment between SigmaR values
        Increment = abs(SigmaRVec(2) - SigmaRVec(1))
        NewInc = 2*Increment/(NSigmaR+2)
        SigmaRVecNew(1) = BestSigmaR-Increment+NewInc
        DO 200 Iinc = 2,NSigmaR
          SigmaRVecNew(Iinc) = SigmaRVecNew(Iinc-1) + NewInc
200     CONTINUE
      ENDIF

      WRITE(21,980) MinDiff,BestCV,BestSigmaR
980   FORMAT('Minimum Diff= ',F10.7,', Best CV= ',F10.7,
     +  ', Best SigmaR= ',F10.7)
      RETURN
      END SUBROUTINE DiffsAndNewSigmaRvec
!=================================================================

!=================================================================
      SUBROUTINE CumulativeCatch(IFvec)
      USE GlobalVars
      IMPLICIT NONE

      !Local variables
       INTEGER,INTENT(IN) :: IFvec
       INTEGER :: IYr,Isim
       REAL(8),DIMENSION(Nsim) :: CumCatch
       REAL(8) :: TotSum

       CumCatch = 0
       TotSum = 0
       DO 222 Isim = 1,Nsim
         DO 222 IYr = MinPYr+BurnInTime+1, MinPYr+BurnInTime+Ntime
           CumCatch(ISim) = CumCatch(ISim) + 
     +       TotalCatchBiomass(IYr,ISim)
           TotSum = TotSum + TotalCatchBiomass(IYr,Isim)
222      CONTINUE
         AvgCumCatch(IFvec) = TotSum/Nsim !Is this right?

      WRITE(19,778) Fvec(IFvec),AvgCumCatch(IFvec)

778   FORMAT(F5.3,',',F30.8)
      RETURN
      END SUBROUTINE CumulativeCatch
!=================================================================


!===============================================================
       SUBROUTINE mydate(mystring)
       IMPLICIT NONE

       ! this is based on a website example: http://download.oracle.com/docs/cd/E19205-01/819-5259/aetcf/index.html
       integer date_time(8)
       character*10 b(3)
       character*25 mystring

       call date_and_time(b(1), b(2), b(3), date_time)

      ! WRITE(20,*) 'End time and date'
        WRITE(20,*) mystring
    !   WRITE(20,300) date_time(5),date_time(6),date_time(7)
       WRITE(20,301) date_time(2),date_time(3),date_time(1),
     +  date_time(5),date_time(6),date_time(7)
       WRITE(*,301) date_time(2),date_time(3),date_time(1),
     +  date_time(5),date_time(6),date_time(7)

     
!       print *,'date_time    array values:'
!       print *,'year=',date_time(1)
!       print *,'month_of_year=',date_time(2)
!       print *,'day_of_month=',date_time(3)
!       print *,'time difference in minutes=',date_time(4)
!       print *,'hour of day=',date_time(5)
!       print *,'minutes of hour=',date_time(6)
!       print *,'seconds of minute=',date_time(7)
!       print *,'milliseconds of second=',date_time(8)
!       print *, 'DATE=',b(1)
!       print *, 'TIME=',b(2)
!       print *, 'ZONE=',b(3)

 300   FORMAT('Time: ',I2,':',I2,':',I2)
 301   FORMAT('Date:',I2,'/',I2,'/',I4,', Time: ',I2,':',I2,':',I2)

       END SUBROUTINE mydate
!=================================================================

!!=================================================================
!      SUBROUTINE WriteSigmaRStuff
!      USE GlobalVars
!      IMPLICIT NONE
!
!      WRITE(*,980) MinDiff,BestCV,BestSigmaR
!980   FORMAT('Minimum Diff= ',F10.7,', Best CV= ',F10.7,
!     +  ', Best SigmaR= ',F10.7)
!
!      RETURN
!      END SUBROUTINE WriteSigmaRStuff
!!================================================================

!==========================================================
      SUBROUTINE BetaOntogenetic
    
      USE GlobalVars
      IMPLICIT NONE

!   Local variables
      INTEGER :: Iage,Ispawn,Igroups
      REAL(8),DIMENSION(1:NOntoGroup,Nspawnarea) :: Pprelim
      REAL(8),DIMENSION(1:NOntoGroup) :: alphaBfunc,betaBfunc,SumPrelim
      REAL(8),DIMENSION(1:Nspawnarea) :: SpawnVec
      REAL(8),DIMENSION(1:Nage,2) :: AgeMatrix 
      INTEGER :: theage,thegroup
      REAL(8),DIMENSION(1:Nage) :: CheckAgain


      !OntoGroupsVec



      !try this out for arranging ontogenetic groups and matching with age:
!      AgeMatrix(1,1) = 1      
!      DO 3099 Iage = 1,Nage
!        AgeMatrix(Iage,1) = Iage
!3099  CONTINUE

! Figure out a matrix that maps ages to groups, with groups numbered 1 to NOntoGroup
      DO 3098 Igroups = 1,NOntoGroup-1
        DO 3098 Iage = OntoGroupsVec(Igroups),OntoGroupsVec(Igroups+1)
          AgeMatrix(Iage,1) = Iage
          AgeMatrix(Iage,2) = Igroups      
3098  CONTINUE

      DO 3099 Iage =  OntoGroupsVec(NOntoGroup),NAge
         AgeMatrix(Iage,1) = iage
         AgeMatrix(iage,2) = NOntoGroup
3099  CONTINUE
     

!!    Figure out vector of beta variances
!      BetaVar(0) = YoungBetaVar
!      DO 7832 Iage = 1,Nage
!        BetaVar(Iage) = BetaVar(Iage-1) + (OldBetaVar-YoungBetaVar)/
!     +     real(Nage)
!7832  CONTINUE

!    Figure out vector of beta variances
      BetaVar(1) = YoungBetaVar
      DO 7832 Igroups = 2,NOntoGroup
        BetaVar(Igroups) = BetaVar(Igroups-1)+(OldBetaVar-YoungBetaVar)/
     +     (real(NOntoGroup)-1.0)
7832  CONTINUE

!    Figure out vector of beta means (do this in betastuff as well,but in terms of age)
      BetaMeans(1) = YoungBetaMean
      DO 7838 Igroups = 2,NOntoGroup
        BetaMeans(Igroups) = BetaMeans(Igroups-1) + 
     +   (OldBetaMean-YoungBetaMean)/(real(NOntoGroup)-1.0)
7838  CONTINUE

      SpawnVec(1) = 1.0
      DO 657 ISpawn = 2,Nspawnarea
        SpawnVec(Ispawn) = 1.0 + SpawnVec(ISpawn-1)
657   CONTINUE

!     Calculate alpha and beta for the Beta functions from means and variance
      SumPrelim = 0
      Pprelim = 0
      P = 0
      CheckP = 0
      DO 300 Igroups = 1,NOntoGroup
        alphaBfunc(Igroups) = (BetaMeans(Igroups)**2-BetaMeans(Igroups)
     +   **3-BetaMeans(Igroups)*BetaVar(Igroups))
        alphaBfunc(Igroups) = alphaBfunc(Igroups)/BetaVar(Igroups)

      betaBfunc(Igroups) = (BetaMeans(Igroups)-2*BetaMeans(Igroups)**2+
     + BetaMeans(Igroups)**3-BetaVar(Igroups)+BetaMeans(Igroups)*
     + BetaVar(Igroups))/BetaVar(Igroups)

        DO 301 Ispawn = 1,Nspawnarea
          Pprelim(Igroups,Ispawn) = (SpawnVec(Ispawn)**
     +     (alphaBfunc(Igroups)-1)
     +     *(Nspawnarea-SpawnVec(Ispawn)+1)**(betaBfunc(Igroups)-1))
          SumPrelim(Igroups) = SumPrelim(Igroups) + Pprelim(Igroups,
     +      Ispawn)
301        CONTINUE
          DO 302 Ispawn = 1,Nspawnarea
           Pprelim(Igroups,Ispawn) = Pprelim(Igroups,Ispawn)/
     +       SumPrelim(Igroups)
           CheckP(Igroups) = CheckP(Igroups) + Pprelim(Igroups,Ispawn)
302       CONTINUE
300   CONTINUE
      
      !Incorporate P into population dynamics in other functions
      CheckAgain = 0
      DO 4500 Ispawn = 1,Nspawnarea
        DO 4500 Iage = 1,Nage
          theage = INT(AgeMatrix(Iage,1))
          thegroup = INT(AgeMatrix(Iage,2))
	    P(theage,Ispawn) = Pprelim(thegroup,Ispawn)
          CheckAgain(Iage) = CheckAgain(Iage) + P(theage,Ispawn)
4500  CONTINUE

      RETURN
	END SUBROUTINE BetaOntogenetic
!==========================================================

!====================================================================================
!Add a larval movement model by adding alternative recruitment and eqmpop subroutines
!Going for the two box model for now (???) How easy is the extension?
!Just a matter of different definitions of movement matrices??
!====================================================================================

!====================================================================================
      SUBROUTINE TwoBoxMovement(ToMove,IYr)
      !This function assumes that there are 2 boxes and that the MPA is in box 1 if there is an MPA
      USE GlobalVars
      IMPLICIT NONE

      INTEGER :: Ibox,IYr
      REAL(8),DIMENSION(Nbox),INTENT(IN) :: ToMove

      !Make constraint:
      !LMoveIn, LMoveOut: LMoveIn is defined
      WasMoved(1:2) = 0.0
      !LMoveOut = LMoveIn*(1-MPASize)/MPASize
      LMoveIn = LMoveOut*MPASize/(1-MPASize)
!      DO 100 Ibox = 1,Nbox
        !think about which box the MPA is! could be either in this setup!
!        IF (Ibox.EQ.OpenArea) THEN
!          !this is the open area
          WasMoved(2) = ToMove(2) - 
     +      ToMove(2)*LMoveIn + 
     +      ToMove(1)*LMoveOut
!        ELSEIF (Ibox.EQ.ClosedArea) THEN
          !this is the closed area
          WasMoved(1) = ToMove(1)-
     +     ToMove(1)*LMoveOut +
     +      ToMove(2)*LMoveIn
 !       ENDIF
!100   CONTINUE
      
      RETURN
	END SUBROUTINE TwoBoxMovement
!====================================================================================

CC ===========================================================================
CC
      SUBROUTINE GenNormalErrorBI(insigmaR)
C
      USE GlobalVars
      IMPLICIT NONE
C
C     Local variables
      INTEGER Itime,Isim,Iarea,IsigmaR
	REAL(8), INTENT(IN) :: insigmaR
      REAL*8 XNORM
      EXTERNAL XNORM

	MyError = 0
C     Generate the recruitment error terms
      DO 10000 Itime = MinPYr+MoveBurnInTime+1,MinPYr+MoveBurnInTime+
     +BurnInTime
       DO 10000 Iarea = 1,Nspawnarea
	  DO 10000 Isim = 1,Nsim
           MyError(Itime,Iarea,Isim) = XNORM(2,0.0d0,insigmaR,ISEEDO)
10000 CONTINUE  
      
      ISEEDO = SEED1
      DO 10001 Itime = MinPYr+MoveBurnInTime+BurnInTime+1,MinPYr+
     +MoveBurnInTime+BurnInTime+Ntime
       DO 10001 Iarea = 1,Nspawnarea
	  DO 10001 Isim = 1,Nsim
           MyError(Itime,Iarea,Isim) = XNORM(2,0.0d0,insigmaR,ISEEDO)
10001 CONTINUE 	     

      RETURN
      END SUBROUTINE GenNormalErrorBI
C============================================================================

C=====================================================================
#Generate errors for distribution of number of days of starvation for MAE

CC ===========================================================================
CC
      SUBROUTINE GenStarveDistnBI(InMeanDays)
C
      USE GlobalVars
      IMPLICIT NONE
C
C     Local variables
      INTEGER Itime,Isim,Iarea,Idays
	REAL(8), INTENT(IN) :: InMeanDays
      REAL*8 XNORM
      EXTERNAL XNORM

	DaysError = 0.0
    !  INTEGER :: ISEED4,ISEED5
      randominc = randominc + 1
      
      IF (randominc.EQ.1) THEN
        OPEN(UNIT=2,FILE="RANDOM.NUM")
        READ(2,*) ISEED4,ISEED5
        CLOSE(2)
      ENDIF

      ISEEDO = ISEED4
C     Generate the recruitment error terms
      DO 999 Itime = MinPYr,MinPYr+MoveBurnInTime
        DO 999 Isim = 1,Nsim
        DaysError(Itime,Isim) = log(InMeanDays)+(SigmaD**2)/2
999   CONTINUE
      DO 10000 Itime = MinPYr+MoveBurnInTime+1,MinPYr+MoveBurnInTime+
     +BurnInTime
	  DO 10000 Isim = 1,Nsim
           DaysError(Itime,Isim)=XNORM(2,log(InMeanDays),SigmaD,ISEEDO)
10000 CONTINUE  
      
      ISEEDO = ISEED4
      DO 10001 Itime = MinPYr+MoveBurnInTime+BurnInTime+1,MinPYr+
     +MoveBurnInTime+BurnInTime+Ntime
	  DO 10001 Isim = 1,Nsim
           DaysError(Itime,Isim)=XNORM(2,log(InMeanDays),SigmaD,ISEEDO)
10001 CONTINUE 	     

      RETURN
      END SUBROUTINE GenStarveDistnBI
C============================================================================

C=========================================================================
      SUBROUTINE CalcSbyAge(InDays)
C
      USE GlobalVars
      IMPLICIT NONE

      !Local Variables
      INTEGER :: Iage
      REAL(8),INTENT(IN) :: InDays

      IF (Berkeley.EQ.1) THEN
       DO 5000 Iage = 1,Nage
        d50(Iage) = max(0.001,-15.23 + 28.79*(1-exp(-0.23*Iage))) !Berkeley's original
5000   CONTINUE
      ELSEIF (Berkeley.EQ.2) THEN
        DO 5002 Iage = 1,Nage
        d50(Iage) = max(0.001,0.6534*Iage + 3.1395) !Linear
5002   CONTINUE
      ELSEIF (Berkeley.EQ.3) THEN
        DO 5003 Iage = 1,Nage
        d50(Iage) = max(0.001,4.5573*exp(0.0717*Iage)) !Exponential
5003   CONTINUE
      ELSEIF (Berkeley.EQ.4) THEN
        DO 5004 Iage = 1,Nage
        d50(Iage) = max(0.001,6.838*log(real(Iage))-5.5584) !Logarithmic
5004   CONTINUE
      ENDIF

      DO 5001 Iage = 1,Nage
       IF (MAEctl.EQ.0) THEN
        SbyAge(Iage) = (1+exp(-log(19.0)*(Indays-d50(Iage))/
     +   (-DiffPar)))**(-1)
       ELSEIF (MAEctl.GT.1) THEN
        SbyAge(Iage) = (1+exp(-log(19.0)*(Indays-d50(MAEctl))/
     +   (-DiffPar)))**(-1)
       ENDIF
5001  CONTINUE

      RETURN
      END SUBROUTINE CalcSbyAge
C==========================================================================

C==========================================================================
      SUBROUTINE CalcSstarve(Year)
C
      USE GlobalVars
      IMPLICIT NONE

      !Local Variables
      INTEGER :: Iage, Year, Ibox
      REAL(8),DIMENSION(1:Nage,1:Nbox) :: NumeratorBits
      REAL(8),DIMENSION(1:Nbox) :: Numerator,Denominator
      Numerator = 0.0
      Denominator = 0.0

      DO 6001 Ibox = 1,Nbox
       DO 6000 Iage = 1,Nage
      NumeratorBits(Iage,Ibox) = SbyAge(Iage)*NCompact(Iage,Ibox,Year-1)
     +*FecundityAtAge(Iage)*0.5
        Numerator(Ibox) = Numerator(Ibox) + NumeratorBits(Iage,Ibox)

       Denominator(Ibox)=Denominator(Ibox) + NCompact(Iage,Ibox,Year-1)*
     +FecundityAtAge(Iage)*0.5
6000   CONTINUE

      Sstarve(Ibox,Year) = Numerator(Ibox)/Denominator(Ibox) 
6001  CONTINUE
      RETURN
      END SUBROUTINE CalcSstarve
C==========================================================================


C==========================================================================
      SUBROUTINE UnfishedCalcSstarve(Year)
C
      USE GlobalVars
      IMPLICIT NONE

      !Local Variables
      INTEGER :: Iage, Year, Ibox
      REAL(8),DIMENSION(1:Nage,1:Nbox) :: NumeratorBits
      REAL(8),DIMENSION(1:Nbox) :: Numerator,Denominator
      Numerator = 0.0
      Denominator = 0.0

      DO 6001 Ibox = 1,Nbox
       DO 6000 Iage = 1,Nage
        NumeratorBits(Iage,Ibox) = SbyAge(Iage)*NOneArea0(Iage,Ibox)
     +*FecundityAtAge(Iage)*0.5
        Numerator(Ibox) = Numerator(Ibox) + NumeratorBits(Iage,Ibox)

        Denominator(Ibox)=Denominator(Ibox) + NOneArea0(Iage,Ibox)*
     +FecundityAtAge(Iage)*0.5
6000  CONTINUE

      UnfishedSstarve(Ibox,Year) = Numerator(Ibox)/Denominator(Ibox) 
6001  CONTINUE

      RETURN
      END SUBROUTINE UnfishedCalcSstarve
C==========================================================================

!=================================================================
      SUBROUTINE DiffsAndNewDaysvec
      USE GlobalVars
      IMPLICIT NONE

      !Local variables
      INTEGER :: IDays,count,Iinc
      REAL(8),DIMENSION(NDays) :: Diffs
      REAL(8) :: BestCV
      REAL(8) :: Increment,NewInc

      !Calculate differences between Recruit CVs and target recruit CV
      Diffs = 0
      DO 100 IDays = 1,NDays
        Diffs(IDays) = abs(MeanRecruitCV(IDays)-TargetCV)
100   CONTINUE

      !If none of the Dayss is not within the desired tolerance level, make a new Daysvec
   !   PRINT *,'FmultiplyVec is ',FmultiplyVec
   !   PRINT *,'Diffs is ',Diffs

      !MinDiff = min(Diffs)
      MinDiff = 1000
      BestCV = -1
      BestDay = -1
      count = 0
      DO 101 IDays =1,NDays
        IF (Diffs(IDays).EQ.MinDiff) THEN
          count = count+1
        ELSEIF (Diffs(IDays).LT.MinDiff) THEN
          MinDiff = Diffs(IDays)
          BestCV = MeanRecruitCV(IDays)
          BestDay = DaysVec(IDays)
        ENDIF
101   CONTINUE

      IF (MinDiff.LE.LowTolerance) THEN
        MeanDays = BestDay
        !FindVariance = 0
      ELSE
        !ReAssign DaysVec
        !Keep FindVariance = 1
        !Need a DO WHILE or DO UNTIL or something like that. Can FORTRAN do that?

        !Increment between Days values
        Increment = abs(DaysVec(2) - DaysVec(1))
        NewInc = 2*Increment/(NDays+2)
        DaysVecNew(1) = BestDay-Increment+NewInc
        DO 200 Iinc = 2,NDays
          DaysVecNew(Iinc) = DaysVecNew(Iinc-1) + NewInc
200     CONTINUE
      ENDIF

      ReportBestCV = BestCV
      RETURN
      END SUBROUTINE DiffsAndNewDaysVec
!=================================================================



 !=================================================================
      SUBROUTINE DiffsAndNewFmultVec
      USE GlobalVars
      IMPLICIT NONE

      !Local variables
      INTEGER :: Ifmult,count,Iinc
      REAL(8),DIMENSION(NFmultiply) :: Diffs
      REAL(8) :: Increment,NewInc

      !Calculate differences between Recruit CVs and target recruit CV
      Diffs = 0.0
      DO 100 Ifmult = 1,NFmultiply
        Diffs(Ifmult) = abs(AvgR0(Ifmult)-r0)
100   CONTINUE

      !If none of the Dayss is not within the desired tolerance level, make a new FmultiplyVec
      
      !MinDiff = min(Diffs)
      MinDiffFmult = 1000
      BestR0 = -1
      BestFmult = -1
      count = 0
      DO 101 Ifmult =1,NFmultiply
        IF (Diffs(Ifmult).EQ.MinDiffFmult) THEN
          count = count+1
        ELSEIF (Diffs(Ifmult).LT.MinDiffFmult) THEN
          MinDiffFmult = Diffs(Ifmult)
          BestR0 = AvgR0(Ifmult)
          BestFmult = FmultiplyVec(Ifmult)
        ENDIF
101   CONTINUE
          
          PRINT *, 'The numbers below are the best for the above 
     +FmultiplyVec, but not necessarily adopted'
          PRINT *,'MinDiffFmult is ',MinDiffFmult
          PRINT *,'BestR0 is ',BestR0
          PRINT *,'BestFmult is ',BestFmult
        !  PAUSE

      IF (MinDiffFmult.LE.HighTolerance) THEN
        fmultiply = BestFmult
        !FindVariance = 0
      ELSE
        !ReAssign FmultiplyVec
        !Keep FindVariance = 1
        !Need a DO WHILE or DO UNTIL or something like that. Can FORTRAN do that?

        !Increment between Days values
        Increment = abs(FmultiplyVec(2) - FmultiplyVec(1))
        IF (BestFmult.NE.FmultiplyVec(1).AND.BestFmult.NE.FmultiplyVec
     +(NFmultiply)) THEN
         PRINT *,'BestFmult is ',BestFmult,' and is NOT first or last in
     + FmultiplyVec'
    !     NewInc = 2*Increment/(NFmultiply)
        NewInc = 2*Increment/(NFmultiply-1) ! This covers all possibilities, but is repetitive - I don't think it's needed, especially not with bookends.
         FmultVecNew(1) = BestFmult-Increment
         FmultVecNew(1) = max(1.0,FmultVecNew(1))
!         DO 200 Iinc = 2,NFmultiply
!           FmultVecNew(Iinc) = FmultVecNew(Iinc-1) + NewInc
!200      CONTINUE
        ELSEIF (BestFmult.EQ.1.0) THEN
    !      NewInc = 2*Increment/(NFmultiply)
           NewInc = 2*Increment/(NFmultiply-1) ! This covers all possibilities, but is repetitive - I don't think it's needed, especially not with bookends.
           FmultVecNew(1) = 1.0
        ELSEIF (BestFmult.EQ.FmultiplyVec(1)) THEN
         PRINT *,'BestFmult is ',BestFmult,' and IS FIRST in
     + FmultiplyVec'
          !Bookend Stuff
          NewInc = Increment
          FmultVecNew(1) = max(1.0,BestFmult-(NewInc*Nfmultiply)+
     +2*NewInc)
!          DO 290 Iinc = 2,NFmultiply
!           FmultVecNew(Iinc) = FmultVecNew(Iinc-1) + NewInc
!290      CONTINUE
        ELSEIF (BestFmult.EQ.FmultiplyVec(Nfmultiply)) THEN
         PRINT *,'BestFmult is ',BestFmult,' and IS LAST in
     + FmultiplyVec'
         !Bookend stuf
         NewInc = Increment
         FmultVecNew(1) = BestFmult-NewInc
        ELSE
          PRINT *,'Could not assign FmultVecNew - time to debug'
          PAUSE
        ENDIF
        DO 200 Iinc = 2,NFmultiply
          FmultVecNew(Iinc) = FmultVecNew(Iinc-1) + NewInc
200     CONTINUE
      ENDIF

!      Do 201 Iinc = 1,NFmultiply
!        FmultvecNew(Iinc) = max(0.001,FmultvecNew(Iinc))
!201   CONTINUE

      WRITE(21,980) MinDiffFmult,BestR0,BestFmult,MeanDays,hpar
980   FORMAT('MinDiffFmult= ',F15.10,', Best R0= ',F10.5,
     +  ', Best fmultiply = ',F10.7,' This is for MeanDays= ',F10.7,
     +' and hpar is ',F10.7)
      RETURN
      END SUBROUTINE DiffsAndNewFmultVec
!=================================================================

!=================================================================
      SUBROUTINE CalcAvgR0(FmultIndex)  !Calculates AvgR (and happens to calc AvgR0 when F = 0
      USE GlobalVars
      IMPLICIT NONE

      !Local variables
      INTEGER,INTENT(IN) :: FmultIndex
      INTEGER :: Isim,IYear,Ibox,i
      REAL(8) :: SumR0

      i=0
      !AvgR0 = 0
      SumR0 = 0.0
      DO 322 Isim = 1,Nsim
       DO 322 IYear = (MinPYr+MoveBurnInTime+BurnInTime+1),(MinPYr+
     +MoveBurnInTime+BurnInTime+Ntime)
         i= i+1
         DO 322 Ibox = 1,Nbox
             SumR0 = SumR0 + TotRecruits(Ibox,IYear,ISim)
322   CONTINUE
      
      IF (Fpar.EQ.0.0) THEN
        AvgR0(FmultIndex) = SumR0/i
      ELSE
        AvgR(FmultIndex) = SumR0/i
      ENDIF
        AvgSteepness(FmultIndex) = AvgR(FmultIndex)/r0
      RETURN
      END SUBROUTINE CalcAvgR0
!=================================================================

!Find Fmultiply Chunk for MAE runs:
       !------------------------------------------------------------------------------------

       SUBROUTINE DoFindFmult
       USE GlobalVars
       IMPLICIT NONE

       !Local Variables
       INTEGER :: Isim,IYear,Ifmult,iterbug


       IF (FindFmultiply.NE.1.AND.PEorMAE.EQ.1) THEN
       fmultiply = FmultiplyVec(1)
       ENDIF
      
       IF (FindFmultiply.EQ.1.AND.PEorMAE.EQ.1) THEN
        FmultiplyVec = SaveFmultiplyVec
       Fpar = 0.0
       iterbug = 0
      
       !SaveMoveBITime = MoveBurnInTime
       !MoveBurnInTime = 0

       !RIGHT HERE ADD A LOOP FOR MOVEBURNINTIME THAT's DONE WITH DetFmultiply - only needs to happen once. OR this could just happen once at the beginning of the program before FindVariance OR FindFmultiply,
       !then no need to come back to MoveBurnInTime time period at all.  Put this just after CALL UNFISHEDPOP and CALL EQMPOP.

         DO
           iterbug = iterbug+1
           !debugging code
!           PRINT *,'FmultiplyVec is'
!           DO 5456 Ifmult = 1,Nfmultiply
!             PRINT *,FmultiplyVec(Ifmult)
!5456       CONTINUE
           !end debugging code
           DO 4000 Ifmult = 1,Nfmultiply
             fmultiply = FmultiplyVec(Ifmult)

       !      CALL GenStarveDistnBI(MeanDays)
    !These are only relevant if MoveBurnInTime = 0
              IF (MoveBurnInTime.EQ.0) THEN
               CALL EQMPOP()
               CALL FirstTime(MinPYr)
              ENDIF

          
             DO 4001 Isim = 1,Nsim
             DO 4001 IYear = MinPYr+MoveBurnInTime+1,MinPYr+
     +                 MoveBurnInTime+BurnInTime+Ntime
               !Maternal age stuff
               CALL CalcsByAge(exp(DaysError(IYear,Isim)-(SigmaD**2/2)))

               IF (MAEctl.EQ.1) THEN
                CALL UnfishedCalcSstarve(IYear)
              CALL RECRUITMENT(IYear,Isim,UnfishedSstarve(1:Nbox,IYear))
               ELSEIF (MAEctl.EQ.0) THEN
                 CALL CalcSstarve(IYear)
                 CALL RECRUITMENT(IYear,Isim,Sstarve(1:Nbox,IYear))
               ELSEIF (MAEctl.GT.1) THEN
                 CALL ControlCalcSstarve(IYear,MAEctl)
                 CALL RECRUITMENT(IYear,Isim,Sstarve(1:Nbox,IYear))
               ENDIF
               
               CALL POPDY(IYear,Isim)
               !Do you want to write anything out? do that here.
4001         CONTINUE
          !Calculate the CV of recruits based on NonSpatialRecruits(Year,Sim) - make a subroutine for calculating this. Use R.
          !Record for each MeanDays
           !CALL CalcRecruitCV(Idays)
           CALL CalcAvgR0(Ifmult)
           
           !DEBUGGING CODE--------------------------------------------------------------
           !Right here write out DaysError info so that you can compare to what it is for Fpar loop.
           !For debugging purposes:
           WRITE(27,*) 'MeanDays,Fmultiply,Year,Sim,DaysError'
           DO 4040 ISim = 2,3
            DO 4040 IYear=(MinPYr+MoveBurnInTime+BurnInTime+1),(MinPYr+
     +MoveBurnInTime+BurnInTime+Ntime)
            WRITE(27,877) MeanDays,fmultiply,IYear,ISim,
     +DaysError(IYear,ISim)
4040       CONTINUE  
         !DEBUGGING CODE--------------------------------------------------------------

4000       CONTINUE
        

    !       PRINT *,'iterbug is ',iterbug,' just before 
    ! +DiffsAndNewFmultVec'
         !Find value closest to 0.5
         !Make a new sigmaRvec and do it all over again if needed.
           CALL DiffsAndNewFmultVec
      IF (iterbug.GE.MaxIterations.OR.MinDiffFmult.LE.HighTolerance) 
     !THEN
               itrue = itrue + 1
               fmultiply = BestFmult
               IF (MinDiffFmult.GT.5.0) THEN
                 PRINT *,'About to pause in DoFindFmult'
                 PRINT *,'hpar is ',hpar,' MeanDays is ',MeanDays
                 PRINT *, 'MinDiffFmult is ',MinDiffFmult,' 
     +and MaxIterations has been reached; BestFmult is ',fmultiply
   !              PAUSE  !This pause you used a lot when debugging
               ENDIF
               IF (FindH.EQ.0) THEN
                FmultDaysCombo(1,itrue) = MeanDays
                FmultDaysCombo(2,itrue) = BestFmult
               ELSE
              ! FmultDaysCombo(3,itrue) = hpar
                HfmultCombo(1,itrue) = hpar
                HfmultCombo(2,itrue) = BestFmult
               ENDIF
               ! HfmultCombo(3,itrue) = MeanDays !this will start over for each MeanDays anyway.
         !  PRINT *,'Right now MeanDays is', MeanDays
         !  PRINT *,' and hpar is ',hpar
         !  PRINT *,'Final Fmultiply for this MeanDays is ', fmultiply ,
   !  +'and MinDiffFmult is',
   !  +           MinDiffFmult, 'and AvgR0 is ',BestR0
                 IF (MinDiffFmult.GT.HighTolerance) THEN
                   PRINT *, 'HighTolerance for Fmultiply not reached; 
     +               tolerance was ',HighTolerance
                 ENDIF
             ELSE
               FmultiplyVec = FmultVecNew
           ENDIF
      IF (iterbug.GE.MaxIterations.OR.MinDiffFmult.LE.HighTolerance) 
     !exit
       !Before the end of this chunk, re-assign sigmaR for future chunks
       END DO
    !   MoveBurnInTime =  SaveMoveBITime
         ENDIF
       
       FmultiplyVec = SaveFmultiplyVec

!MeanDays,fmultiply,IYear,ISim,DaysError(IYear,ISim)        
877     FORMAT(F20.10,1x,F20.15,1x,I4,1x,I5,1x,F20.10)

       RETURN
       END SUBROUTINE DoFindFmult
 !-----------------------------------------------------------------------------------

C==========================================================================
      SUBROUTINE ControlCalcSstarve(Year,CtlMatAge)
C
      USE GlobalVars
      IMPLICIT NONE

      !Local Variables
      INTEGER,INTENT(IN) :: CtlMatAge
      INTEGER :: Iage, Year, Ibox
      REAL(8),DIMENSION(1:Nage,1:Nbox) :: NumeratorBits
      REAL(8),DIMENSION(1:Nbox) :: Numerator,Denominator
      Numerator = 0.0
      Denominator = 0.0

      DO 6001 Ibox = 1,Nbox
       DO 6000 Iage = 1,Nage
      NumeratorBits(Iage,Ibox) = SbyAge(CtlMatAge)*
     +NCompact(Iage,Ibox,Year-1)
     +*FecundityAtAge(Iage)*0.5
        Numerator(Ibox) = Numerator(Ibox) + NumeratorBits(Iage,Ibox)

       Denominator(Ibox)=Denominator(Ibox) + NCompact(Iage,Ibox,Year-1)*
     +FecundityAtAge(Iage)*0.5
6000   CONTINUE

      Sstarve(Ibox,Year) = Numerator(Ibox)/Denominator(Ibox) 
6001  CONTINUE
      RETURN
      END SUBROUTINE ControlCalcSstarve
C==========================================================================


!CC ===========================================================================
!CC
!      SUBROUTINE GenStarveDistnSims(InMeanDays)
!C
!      USE GlobalVars
!      IMPLICIT NONE
!C
!C     Local variables
!      INTEGER Itime,Isim,Iarea,Idays
!	REAL(8), INTENT(IN) :: InMeanDays
!      REAL*8 XNORM
!      EXTERNAL XNORM
!
!	DaysErrorSims = 0
!	DO 10000 Isim = 1,SRsims
!        DaysErrorSims(Isim)=XNORM(2,log(InMeanDays),SigmaD,ISEEDO)
!10000 CONTINUE  
!   	     
!
!      RETURN
!      END SUBROUTINE GenStarveDistnSims
!C============================================================================



!-------------------------------------------------------------
      SUBROUTINE DetFmult()

      USE GlobalVars
      IMPLICIT NONE

      !Local variables
      INTEGER :: Iage,Ibox,IFpar
      REAL(8) :: SBPR0fmult,Part1,Part2,Part3


      !---------------------------------------------------------
      !Do calculations for deterministic MeanDays SR curve
      CALL CalcSbyAge(MeanDays)
      !-------------------------------------------------------------------------
      !When Fpar = 0.0 find SBPRFmult(0) which incorporates SbyAge, but not Dfmultiply
      CALL UnfishedCalcSstarve(1)
        Fpar = 0.0
          DO 8007 Iage = 1,Nage
           DO 8007 Ibox = 1,Nbox
          SBPR0fmult=SBPR0fmult+0.5*SbyAge(Iage)*
     +UnfishedNperRecruit(Iage,Ibox)*
     +FecundityAtAge(Iage)
8007    CONTINUE

       Part1 = r0*(5*hpar-1)*SBPR0fmult
       Part2 = Part1**2+16*r0*r0*hpar*SBPR0fmult*UnfishedSstarve(1,1)*
     +(1-hpar)*SBPR0
       Part3 = 8*r0*hpar*SBPR0fmult*UnfishedSstarve(1,1)

       Dfmultiply = (Part1+sqrt(Part2))/Part3 
   !-----------------------------------------------------------
      WRITE(24,624) SRsimsNFpar,Dfmultiply,MeanDays
      WRITE(24,623)

      DSBPRsr = 0.0
      DO 5433 IFpar = 1,SRsimsNFpar
       Fpar = SRFvec(IFpar)

      CALL EQMPOP()
      CALL FirstTime(MinPYr)

       IF (MAEctl.EQ.1) THEN
         CALL UnfishedCalcSstarve(1)
       ELSEIF (MAEctl.EQ.0) THEN
         CALL CalcSstarve(1)
       ELSEIF (MAEctl.GT.1) THEN
         CALL ControlCalcSstarve(1,MAEctl)
       ENDIF

      !Find SBPR
          DO 6007 Iage = 1,Nage
           DO 6007 Ibox = 1,Nbox
          DSBPRsr(IFpar)=DSBPRsr(IFpar)+0.5*Dfmultiply*SbyAge(Iage)*
     +FishedNperRecruit(Iage,Ibox)*
     +FecundityAtAge(Iage)
6007    CONTINUE
      
       DRecSR(IFpar) = MAX(r0*((4*hpar*DSBPRsr(IFpar)*
     +Sstarve(1,1)*Dfmultiply)-(1-hpar)*SBPR0)/((5*hpar-1)*
     +DSBPRsr(IFpar)),0.0)

      DSSBsr(IFpar) = DRecSR(IFpar)*DSBPRsr(IFpar)

      !CALL DetFmult(IFpar)

       WRITE(24,621) SRFvec(IFpar),DRecSR(IFpar),DSSBsr(IFpar),
     +DSBPRsr(IFPar)

5433  CONTINUE

621   FORMAT(F7.5,1x,2(F15.4,1x),F15.8)
623   FORMAT('F, Recruits, SSB, SBPR')
624   FORMAT('SRsimsNFpar',1x,I3,1x,
     +'Dfmultiply',1x,F10.7,1x,'MeanDays',1x,F10.4)

      RETURN
      END SUBROUTINE DetFmult
      
!------------------------------------------------------------

!--------------------------------------------------------------------------------
      SUBROUTINE DoSRsims()

      USE GlobalVars
      IMPLICIT NONE

      !Local variables
      INTEGER :: Isim,IFpar,Iage,Ibox,Ibin
      REAL(8) :: DaysErrorSR
      REAL(8) :: SBPR
      REAL(8) :: RecSRsims,SSBsrsims
      REAL(8) :: SSBsrSimsBinned
      REAL(8) :: SSB0
      REAL(8),DIMENSION(0:Nage) :: TheErrorSR
      INTEGER :: TheSim,TheAge
      REAL(8) :: XNORM
      EXTERNAL XNORM

      !Do a loop over SBPRsims and F
      !Write things out as you go.
      	DaysErrorSR = 0.0
      MaxSSBsrsims = 0.0
      !Calculate SBPR0

      SBPR0 = 0.0
      DO 6003 Iage= 1,Nage
        DO 6003 Ibox = 1,Nbox
          SBPR0 = SBPR0+0.5*UnfishedNperRecruit(Iage,Ibox)*
     +    FecundityAtAge(Iage)
6003  CONTINUE

      !Find SSB0
      SSB0 = 0.0
      DO 6033 Ibox = 1,Nbox
        SSB0 = SSB0 + SSB0vec(Ibox)
6033  CONTINUE

      !Bin the SSB data:
      !Make bins:
      SSBbinVec(1) = 0.0
      DO 5413 Ibin = 2,Nbins
       SSBbinVec(Ibin) = SSBbinVec(Ibin-1) + 2*SSB0/Nbins
5413  CONTINUE

      !Make and write out DaysErrorSR for the random deviates (in logspace)
      WRITE(23,622) SRsimsNFpar,SRsims,Nbins,fmultiply,MeanDays
      WRITE(23,625)

      IF (SRsims.GT.0) THEN
      ISEEDO = SEED1
      WRITE(25,628)
      DO 6001 Isim = 1,SRsims
       DO 6001 Iage = 0,Nage
         DaysErrorSR=XNORM(2,log(MeanDays),SigmaD,ISEEDO)
         WRITE(25,627) Isim,Iage,DaysErrorSR
6001  CONTINUE
  !    CLOSE(25)
  !    OPEN(25)  
        SBPR = 0
        DO 5412 IFpar = 1,SRsimsNFpar
          Fpar = SRFvec(IFpar)
          !rewind file 25 all the way
          REWIND(25)
          READ(25,*)
      !    READ(25,*)
         CALL EQMPOP()
         CALL FirstTime(MinPYr)
          DO 5411 Isim = 1,SRsims
         !Find SBPR
         SBPR = 0.0
         RecSRsims = 0.0
         SSBsrsims = 0.0
          DO 6022 Iage = 0,Nage
            READ(25,*) TheSim,TheAge,TheErrorSR(Iage)
6022      CONTINUE
          DO 6002 Iage = 1,Nage
             !READ(25,*) TheSim,TheAge,TheErrorSR(Iage)
             !Calculate gamma (Sstarve) for each sim
             CALL CalcSbyAge(exp(TheErrorSR(Iage)-(SigmaD**2)/2))
           DO 6002 Ibox = 1,Nbox
          SBPR=SBPR+0.5*fmultiply*SbyAge(Iage)*
     +FishedNperRecruit(Iage,Ibox)*
     +FecundityAtAge(Iage)
6002    CONTINUE
     
           !Calculate gamma (Sstarve) for each sim
           CALL CalcSbyAge(exp(TheErrorSR(0)-(SigmaD**2)/2))

             IF (MAEctl.EQ.1) THEN
               CALL UnfishedCalcSstarve(1)
             ELSEIF (MAEctl.EQ.0) THEN
              CALL CalcSstarve(1)
             ELSEIF (MAEctl.GT.1) THEN
              CALL ControlCalcSstarve(1,MAEctl)
             ENDIF

       RecSRsims = MAX(r0*((4*hpar*SBPR*
     +Sstarve(1,1)*fmultiply)-(1-hpar)*SBPR0)/((5*hpar-1)*
     +SBPR),0.0)

       SSBsrsims = RecSRsims*SBPR

       !Figure out the maximum value for SSBsrsims
       !Instead of doing this define upper bin to be 2*B0 or something.
!       IF (SSBsrsims.GT.MaxSSBsrsims) THEN
!         MaxSSBsrsims = SSBsrsims
!       ENDIF


         DO 5415 Ibin = 1,(Nbins-1)
      IF (SSBsrsims.GE.SSBbinVec(Ibin).AND.
     +SSBsrsims.LT.SSBbinVec(Ibin+1)) THEN
            SSBsrSimsBinned = SSBbinVec(Ibin)
          ENDIF
5415     CONTINUE
         IF (SSBsrsims.GE.SSBbinVec(Nbins)) THEN
           SSBsrSimsBinned = SSBbinVec(Nbins)
         ENDIF

      !Write out F, RecSRsims, SBPR, and SSBsrsims
      !loop over F and sims to do this

        WRITE(23,626) SRFvec(IFpar),Isim,RecSRsims,
     +SSBsrsims,SBPR,SSBsrSimsBinned
      !-----------------------------------------------------------

5411    CONTINUE
5412   CONTINUE
      ENDIF 

      !---------------------------------------------------------
      !Do calculations for deterministic MeanDays SR curve
       CALL DetFmult()

621   FORMAT(F7.5,1x,2(F15.4,1x),F15.8)
622   FORMAT('SRsimsNFpar',1x,I3,1x,'SRsims',1x,I4,1x,'Nbins',1x,I6,1x,
     +'fmultiply',1x,F10.7,1x,'MeanDays',1x,F10.4)
623   FORMAT('F, Recruits, SSB, SBPR')
624   FORMAT('SRsimsNFpar',1x,I3)
625   FORMAT('F, Sim, Recruits, SSB, SBPR, SSBbinned')
626   FORMAT(F7.5,1x,I4,1x,2(F15.4,1x),2(F15.8,1x))
627   FORMAT(I15,1x,I4,1x,F20.15)
628   FORMAT('Isim,Iage,DaysErrorSR')
      RETURN
      END SUBROUTINE DoSRsims

!----------------------------------------------------------------------------------

!-------------------------------------------------------------
      SUBROUTINE DetFmultMain()

      USE GlobalVars
      IMPLICIT NONE

      !Local variables
      INTEGER :: Iage,Ibox
      REAL(8) :: SBPR0fmult,Part1,Part2,Part3


      !---------------------------------------------------------
      !Do calculations for deterministic MeanDays SR curve
      CALL CalcSbyAge(MeanDays)
      !-------------------------------------------------------------------------
      !When Fpar = 0.0 find SBPRFmult(0) which incorporates SbyAge, but not Dfmultiply
      CALL UnfishedCalcSstarve(1)
        Fpar = 0.0
          DO 8007 Iage = 1,Nage
           DO 8007 Ibox = 1,Nbox
          SBPR0fmult=SBPR0fmult+0.5*SbyAge(Iage)*
     +UnfishedNperRecruit(Iage,Ibox)*
     +FecundityAtAge(Iage)
8007    CONTINUE

       Part1 = r0*(5*hpar-1)*SBPR0fmult
       Part2 = Part1**2+16*r0*r0*hpar*SBPR0fmult*UnfishedSstarve(1,1)*
     +(1-hpar)*SBPR0
       Part3 = 8*r0*hpar*SBPR0fmult*UnfishedSstarve(1,1)

       Dfmultiply = (Part1+sqrt(Part2))/Part3 
   !-----------------------------------------------------------
      RETURN
      END SUBROUTINE DetFmultMain
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!-------------------------------------------------------------
      SUBROUTINE DetFmultMainTry2()

      USE GlobalVars
      IMPLICIT NONE

      !Local variables
      INTEGER :: Iage,Ibox
      REAL(8) :: apar,bpar,cpar
      REAL(8) :: SSB0pieceMAE, SSB0pieceNoMAE
      REAL(8),DIMENSION(Nage) :: N0stuff
      !Don't know what the difference is between DetFmultMain and this, but can explore that later - let's just see if this works or if there's a mistake.
      
      N0stuff = 0.0
      DO 5554 Iage = 1,Nage
       DO 5554 Ibox = 1,Nbox
        N0stuff(Iage) = N0stuff(Iage) + NOneArea0(Iage,Ibox)
5554  CONTINUE

      !sum(N(a,t)*0.5*fec(a)*gamma(a)
      SSB0pieceMAE = 0.0
      SSB0pieceNoMAE = 0.0
      DO 5555 Iage = 1,Nage
        SSB0pieceMAE = SSB0pieceMAE+N0stuff(Iage)*FecundityAtAge(Iage)*
     +   SbyAge(Iage)
        SSB0pieceNoMAE = SSB0pieceNoMAE+N0stuff(Iage)*
     +FecundityAtAge(Iage)
5555  CONTINUE

      apar = 4*hpar*SSB0pieceMAE/SSB0pieceNoMAE
      bpar = -(5*hpar-1)
      cpar = -((1-hpar)*SSB0pieceNoMAE/SSB0pieceMAE)

      Dfmultiply = (-bpar + sqrt(bpar**2-4*apar*cpar))/(2*apar)
      IF (Dfmultiply.LT.0) THEN
        Dfmultiply = (-bpar - sqrt((bpar**2)-(4*apar*cpar)))/(2*apar)
      ENDIF


      RETURN
      END SUBROUTINE DetFmultMainTry2

!-------------------------------------------------------------------------
!-------------------------------------------------------------
      SUBROUTINE SummaryStats(IFpar)

      USE GlobalVars
      IMPLICIT NONE

      !Local variables
      INTEGER,INTENT(IN) :: IFpar
      INTEGER :: Isim,Iyear,Ibox
      REAL(8),DIMENSION(Nsim) :: SumRecruits,MeanRecruits
      REAL(8),DIMENSION(Nsim) :: SumSSB,MeanSSB,SDSSB,CVSSB
      REAL(8),DIMENSION(Nsim) :: SumCatch,MeanCatch,SDCatch,CVCatch
      REAL(8),DIMENSION(Nsim) :: SumSDCalcs,SDRecruits,CVRecruits
      REAL(8),DIMENSION(Nsim) :: SumSDstuffSSB,SumSDstuffCatch
      REAL(8),DIMENSION((MinPYr+MoveBurnInTime+BurnInTime+1):
     +(MinPYr+MoveBurninTime+BurnInTime+Ntime),
     +  Nsim) :: SDCalcs
      REAL(8) :: sumCVRecruits
      REAL(8) :: TotMeanCVSSB
      REAL(8) :: TotMeanSSB,TotMeanCatch,TotMeanCVCatch
      REAL(8),DIMENSION((MinPYr+MoveBurnInTime +BurnInTime+1):(MinPYr+
     +MoveBurnInTime +BurnInTime+Ntime),Nsim) :: MyTotRecs
  !Variables associated with debugging stuff
      REAL(8) :: SumR0,AvgR0thing
      INTEGER :: i

      

      sumCVRecruits = 0.0
      SDRecruits = 0.0
      CVRecruits =0.0

       TotMeanRecruits = 0.0
       TotMeanCVRecruits = 0.0
       TotMeanSSB = 0.0
       TotMeanCVSSB = 0.0
       TotMeanCatch = 0.0
       TotMeanCVCatch = 0.0

      !Repeat CalcCVRecruits:
      DO 322 Isim = 1,Nsim
      !Calc interannual CV for each sim then average the sims
      SumRecruits(Isim) = 0
      MeanRecruits(Isim) = 0
      SumSSB(Isim) = 0
      MeanSSB(Isim) = 0
      SumCatch(Isim) = 0
      MeanCatch(Isim) = 0
      
      DO 320 IYear = (MinPYr+MoveBurnInTime +BurnInTime+1),(MinPYr+
     +MoveBurnInTime +BurnInTime+Ntime)
 !      SumRecruits(ISim) = SumRecruits(ISim) + 
 !    +  NonSpatialRecruits(IYear,ISim)
 !      SumSSB(Isim) = SumSSB(Isim) + TotSSB(IYear,Isim)
 !      SumCatch = SumCatch + TotalCatchBiomass(IYear,Isim)

        !Test to see if NonSpatialRecruits ends up being the same as TotRecruits:
 !       MyTotRecs(IYear,Isim) = 0.0
 !       DO 677 Ibox = 1,Nbox
 !         MyTotRecs(IYear,Isim) = MyTotRecs(IYear,ISim) +
 !    +TotRecruits(Ibox,IYear,Isim)
!677     CONTINUE

       MeanRecruits(Isim) = MeanRecruits(Isim)+NonSpatialRecruits
     +(IYear,Isim)/Ntime
       MeanSSB(Isim) = MeanSSB(Isim) + TotSSB(IYear,Isim)/Ntime
       MeanCatch(Isim)=MeanCatch(Isim)+TotalCatchBiomass(Iyear,Isim)
     +/Ntime
320   CONTINUE

       SumSDCalcs(ISim) = 0
       SumSDstuffSSB(ISim) = 0
       SumSDstuffCatch(Isim) = 0
       DO 321 IYear = (MinPYr+MoveBurnInTime +BurnInTime+1),(MinPYr+
     +MoveBurnInTime +BurnInTime+Ntime)
         SDCalcs(IYear,ISim) = (NonSpatialRecruits(IYear,ISim)-
     +     MeanRecruits(ISim))**2
         SumSDCalcs(ISim) = SumSDCalcs(ISim) + SDCalcs(IYear,ISim)

         SumSDstuffSSB(Isim) = SumSDstuffSSB(Isim) + 
     +(TotSSB(IYear,Isim)-MeanSSB(Isim))**2 

         SumSDstuffCatch(Isim) = SumSDStuffCatch(Isim) +
     +(TotalCatchBiomass(Iyear,Isim)-MeanCatch(Isim))**2
321    CONTINUE

         SDRecruits(ISim) = sqrt(SumSDCalcs(ISim)/(Ntime))
         CVRecruits(Isim) = SDRecruits(Isim)/MeanRecruits(Isim)
         sumCVRecruits = sumCVRecruits + CVRecruits(ISim)
         
         SDSSB(Isim) = sqrt(SumSDstuffSSB(Isim)/Ntime)
         CVSSB(Isim) = SDSSB(Isim)/MeanSSB(Isim)

         IF (MeanCatch(Isim).GT.0) THEN
          SDCatch(Isim) = sqrt(SumSDstuffCatch(Isim)/Ntime)
          CVCatch(Isim) = SDCatch(Isim)/MeanCatch(Isim)
         ELSE
           SDCatch(Isim) = 0.0
           CVCatch(Isim) = 0.0
         ENDIF

         TotMeanRecruits = TotMeanRecruits + MeanRecruits(Isim)/Nsim
         TotMeanCVRecruits = TotMeanCVRecruits + CVRecruits(Isim)/Nsim
         TotMeanSSB = TotMeanSSB + MeanSSB(Isim)/Nsim
         TotMeanCVSSB = TotMeanCVSSB + CVSSB(Isim)/Nsim
         TotMeanCatch = TotMeanCatch + MeanCatch(Isim)/Nsim
         TotMeanCVCatch = TotMeanCVCatch + CVCatch(Isim)/Nsim

322   CONTINUE
      
      !Write to a file
      !Header of file has info on what years are included in the average, NFpar, true SSB0 (?), true R0 (?)
      IF (IFpar.EQ.1) THEN
        WRITE(26,781) NFpar
        WRITE(26,*) 'F, CVRecruits,MeanRecruits,CVSSB,MeanSSB,
     +CVCatch,MeanCatch'
      ENDIF
      WRITE(26,782) Fvec(IFpar),TotMeanCVRecruits,TotMeanRecruits,
     +TotMeanCVSSB,TotMeanSSB,TotMeanCVCatch,TotMeanCatch

      !--------------------------------------------------
      !FOR DEBUGGING, OW COMMENT OUT!
      i=0
      !AvgR0 = 0
      SumR0 = 0.0
      DO 3229 Isim = 1,Nsim
       DO 3229 IYear = (MinPYr+MoveBurnInTime+BurnInTime+1),(MinPYr+
     +MoveBurnInTime+BurnInTime+Ntime)
         i= i+1
         DO 3229 Ibox = 1,Nbox
             SumR0 = SumR0 + TotRecruits(Ibox,IYear,ISim)
3229   CONTINUE

      AvgR0thing = SumR0/i

      !---------------------------------------------------

781   FORMAT('NFpar',I4)
782   FORMAT(F10.7,1x,3(F20.10,1x,F20.6,1x))

      RETURN
      END SUBROUTINE SummaryStats

!----------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
      SUBROUTINE DoMAEInnards()
        USE GlobalVars
       IMPLICIT NONE

       !Local variables
       INTEGER :: Isim,IYear

!      IF (MoveBurnInTime.EQ.0.OR.Fpar.NE.0.0) THEN  !Why did you say Fpar.NE.0? what was that for? was it just efficiency?
       IF (MoveBurnInTime.EQ.0.OR.FindMSY.EQ.1) THEN
          CALL EQMPOP()
         CALL FirstTime(MinPYr)
      ENDIF

  !       IF (Fpar.NE.0.0) THEN
         SaveFmultiply = fmultiply
         IF (PEorMAE.EQ.1) THEN
           fmultiply = DFmultiply !Apply deterministic fmultiply
         ENDIF
         DO 4201 Isim = 1,Nsim  !you don't really need the Isim loop: this is deterministic
           DO 4201 IYear = MinPYr+1,MinPYr + MoveBurnInTime 
             !Maternal age stuff
             !CALL CalcsByAge(MeanDays) !Already ran this function - it doesn't depend on Age or anything.
             IF (PEorMAE.EQ.2) THEN
               Sstarve(1:Nbox,IYear) = 1.0
               fmultiply = 1.0
               CALL RECRUITMENT(IYear,Isim,Sstarve(1:Nbox,IYear))
             ELSEIF (MAEctl.EQ.1) THEN
              CALL UnfishedCalcSstarve(IYear)
              CALL RECRUITMENT(IYear,Isim,UnfishedSstarve(1:Nbox,IYear))
             ELSEIF (MAEctl.EQ.0) THEN
               CALL CalcSstarve(IYear)
               CALL RECRUITMENT(IYear,Isim,Sstarve(1:Nbox,IYear))
              ELSEIF (MAEctl.GT.1) THEN
               CALL ControlCalcSstarve(IYear,MAEctl)
               CALL RECRUITMENT(IYear,Isim,Sstarve(1:Nbox,IYear))
              ENDIF
               CALL POPDY(IYear,Isim)
               !Do you want to write anything out? do that here.   
4201    CONTINUE
        fmultiply = SaveFmultiply
  !      ENDIF

          
             DO 4001 Isim = 1,Nsim
             DO 4001 IYear = MinPYr+MoveBurnInTime+1,MinPYr+
     +                 MoveBurnInTime+BurnInTime+Ntime
               !Maternal age stuff
               CALL CalcsByAge(exp(DaysError(IYear,Isim)-(SigmaD**2/2)))
               IF (PEorMAE.EQ.2) THEN
                 Sstarve(1:Nbox,IYear) = 1.0
                 fmultiply = 1.0
                 CALL RECRUITMENT(IYear,Isim,Sstarve(1:Nbox,IYear))
               ELSEIF (MAEctl.EQ.1) THEN
                CALL UnfishedCalcSstarve(IYear)
              CALL RECRUITMENT(IYear,Isim,UnfishedSstarve(1:Nbox,IYear))
               ELSEIF (MAEctl.EQ.0) THEN
                 CALL CalcSstarve(IYear)
                 CALL RECRUITMENT(IYear,Isim,Sstarve(1:Nbox,IYear))
               ELSEIF (MAEctl.GT.1) THEN
                 CALL ControlCalcSstarve(IYear,MAEctl)
                 CALL RECRUITMENT(IYear,Isim,Sstarve(1:Nbox,IYear))
               ENDIF
               
               CALL POPDY(IYear,Isim)
               !Do you want to write anything out? do that here.
4001         CONTINUE


      RETURN
      END SUBROUTINE DoMAEInnards
!--------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------
      SUBROUTINE DoFindH
       USE GlobalVars
       IMPLICIT NONE

       !Local Variables
       INTEGER :: Isim,IYear,Ifmult,Ihpar,IFpar,buggysplat
       INTEGER :: ibuggy


       Fpar = 0.0
       AvgR = 0.0
       AvgR0 = 0.0
       buggysplat = 0
       itrue = 0  !This gets reset to 0 for each MeanDays, but accumulates across hpar values; also gets set to 0 at top of FindVariance


           DO 4000 Ihpar = 1,Nhpar
             hpar = hparVec(Ihpar)
             PRINT *,'Inside DoFindH, trying out hpar = ',hpar
       !      CALL GenStarveDistnBI(MeanDays)
    !These are only relevant if MoveBurnInTime = 0

    !Call fmultiply using this h and whatever MeanDays has been specified:
              PRINT *,'Now finding correct fmultiply with F = 0.0'
              CALL DoFindFmult() !Sets F to 0 inside
              !Find the Fpar that produces Depletion = 0.2 to work with steepness definition
              
              PRINT *,'Now finding the F that makes Depletion = 0.2,'
              PRINT *,'Current fmultiply is ',fmultiply
              CALL DoFindhFpar()

              !Do MAEInnards with correct Fpar, h, and fmultiply combo
              CALL DoMAEInnards()
              !Calculate R/R0 in vector form for each h value (which comes with an Fpar and fmultiply and MeanDays and should be recorded)
              CALL CalcAvgR0(Ihpar) !This calculates AvgR, not just AvgR0

4000       CONTINUE !end FindH loop
              
              !Find diffs between 0.6 (TargetSteepness) and AvgR/R0
              CALL DiffsAndNewHparVec()
              
              !Calculate a new hvec or whatever that is called.

              hpar = Besth
           PRINT *,'Conclusion of DoFindH:'
           PRINT *,'Final hpar is ', hpar,'for MeanDays is ',MeanDays
           PRINT *,'Average steepness is ',BestSteep
           PRINT *,'And corresponding Fpar is ',ThreeWayCombo(1,md)
           PRINT *,'And corresponding fmultiply is ',ThreeWayCombo(2,md)

                 !ID the BEST h and its associated fmultiply; can make FmultDaysCombo right here: MeanDays, BestH, associated bestFmultiply
                 fmultiply = -1.0
                 DO 5557 Ibuggy = 1,(MaxIterations*Nhpar)
                   IF (HfmultCombo(1,Ibuggy).EQ.Besth) THEN
                     fmultiply = HfmultCombo(2,Ibuggy)
                   ENDIF
5557             CONTINUE
                 !Here fill in FmultDaysCombo with correct BestH and associated fmultiply
                 FmultDaysCombo(1,IMeanDays) = MeanDays
                 FmultDaysCombo(2,IMeanDays) = fmultiply
                 FmultDaysCombo(3,IMeanDays) = Besth
                 FmultDaysCombo(4,IMeanDays) = ThreeWayCombo(1,md)

                 PRINT *,'FmultDaysCombo MeanDays is ',
     +FmultDaysCombo(1,IMeanDays)
                 PRINT *,'FmultDaysCombo fmultiply is ',
     +FmultDaysCombo(2,IMeanDays)
                 PRINT *,'FmultDaysCombo hpar is ',
     +FmultDaysCombo(3,IMeanDays)
                 PRINT *,'FmultDaysCombo Fpar is ',
     +FmultDaysCombo(4,IMeanDays)

                 IF (MinDiffh.GT.LowTolerance) THEN
                   PRINT *, 'LowTolerance for Besth was not reached; 
     +tolerance was ', LowTolerance
                 ENDIF

       !Reset Fpar to 0.0 so that you can continue figuring out what MeanDays should be so that CVRecruits = 0.5 under unfished conditions!
       Fpar = 0.0
       RETURN
       END SUBROUTINE DoFindH
!--------------------------------------------------------------------------------------

!=================================================================
      SUBROUTINE CalcAvgDepletion(FparIndex)
      USE GlobalVars
      IMPLICIT NONE

      !Local variables
      INTEGER,INTENT(IN) :: FparIndex
      INTEGER :: Isim,IYear,i
      REAL(8) :: SumSSB

      i=0
      !AvgR0 = 0
      SumSSB = 0.0
      AvgDepl(FparIndex) = 0.0
      DO 322 Isim = 1,Nsim
       DO 322 IYear = (MinPYr+MoveBurnInTime+BurnInTime+1),(MinPYr+
     +MoveBurnInTime+BurnInTime+Ntime)
         i= i+1
          SumSSB = SumSSB + TotSSB(IYear,ISim)
322   CONTINUE

      AvgSSB(FparIndex) = SumSSB/i
      AvgDepl(FparIndex) = AvgSSB(FparIndex)/SSB0nospace

      RETURN
      END SUBROUTINE CalcAvgDepletion
!=================================================================

!=================================================================
      SUBROUTINE DiffsAndNewFparVec()
      USE GlobalVars
      IMPLICIT NONE

      !Local variables
      !INTEGER,INTENT(IN) :: BMSYorH
      INTEGER :: IFpar,count,Iinc
      REAL(8),DIMENSION(NhFpar) :: Diffs
      REAL(8) :: Increment,NewInc

      !Calculate differences between 0.2 and AvgDepl(IFpar)
      Diffs = 0.0
      DO 100 IFpar = 1,NhFpar
        Diffs(IFpar) = abs(AvgDepl(IFpar)-0.2)
100   CONTINUE

      !If none of the Diffs is not within the desired tolerance level, make a new FindhFvec
      
      !MinDiff = min(Diffs)
      MinDiffDepl = 1000
      BestDepl = -1
      BesthFpar = -1
      count = 0
      DO 101 IFpar =1,NhFpar
        IF (Diffs(IFpar).EQ.MinDiffDepl) THEN
          count = count+1
        ELSEIF (Diffs(IFpar).LT.MinDiffDepl) THEN
          MinDiffDepl = Diffs(IFpar)
          BestDepl = AvgDepl(IFpar)
          BesthFpar = FindhFvec(IFpar)
        ENDIF
101   CONTINUE

      IF (MinDiffDepl.LE.LowTolerance) THEN
        Fpar = BesthFpar
        !FindVariance = 0
      ELSE
        !ReAssign FmultiplyVec
        !Keep FindVariance = 1
        !Need a DO WHILE or DO UNTIL or something like that. Can FORTRAN do that?

        !Increment between Fpars
        Increment = abs(FindhFVec(2) - FindhFvec(1))
        NewInc = 2*Increment/(NhFpar+2)
        FindhFvecNew(1) = BesthFpar-Increment+NewInc
        FindhFvecNew(1) = max(0.0,FindhFvecNew(1))
        DO 200 Iinc = 2,NhFpar
          FindhFvecNew(Iinc) = FindhFvecNew(Iinc-1) + NewInc
200     CONTINUE
      ENDIF

!      Do 201 Iinc = 1,NFmultiply
!        FmultvecNew(Iinc) = max(0.001,FmultvecNew(Iinc))
!201   CONTINUE

      WRITE(21,980) MinDiffDepl,BestDepl,BesthFpar,BestFmult,hpar
     +,MeanDays
980   FORMAT('Minimum Diff Depl= ',F20.7,', BestDepl= ',F10.7,'Best F '
     +,F10.7,
     +  ', using fmultiply = ',F10.7,'This is for hpar ',F10.7,
     +' and MeanDays= ',F10.7)
      RETURN
      END SUBROUTINE DiffsAndNewFparVec
!=================================================================


!-----------------------------------------------------------------------
      SUBROUTINE DoFindhFpar()
      USE GlobalVars

      !Local variables
      INTEGER :: IFpar,buggybuggy

              !Cycle through F values here
              !Need a DO loop that exits when MaxIterations or Tolerance are reached.
             buggybuggy = 0
             FindhFvec = FindhFvecOriginal
             DO
              buggybuggy = buggybuggy+1
              DO 5666 IFpar = 1,NhFpar  !Come back and change this to a second F vector.
                Fpar = FindhFvec(IFpar)
             !   CALL EQMPOP()
             !   CALL FirstTime(MinPYr)
                CALL DoMAEInnards()
                
                !Calculate AvgSSB/SSB0
                CALL CalcAvgDepletion(IFpar)



5666          CONTINUE !end FindhFpar loop

              !Calculate an ObjFun: (AvgSSB/SSB0 - 0.2)**2
              !Calculate min objfun and make a new vector around that

              CALL DiffsAndNewFparVec()
              !Go until reach tolerance or maxiterations OR use FIT function (will be pretty specific) with original vec. just between 0.06 and 0.20.
              !Make sure that Fpar = bestFpar at the end of this.
              
              IF (buggybuggy.GE.MaxIterations.OR.MinDiffDepl.LE.
     +LowTolerance) THEN
                 Fpar = BesthFpar
                 PRINT *,'Right now Fpar is ', Fpar
                 PRINT *,'And AvgDepl is ',BestDepl,' And hpar is ',hpar
                 IF (MinDiffDepl.GT.LowTolerance) THEN
                   PRINT *, 'LowTolerance for BestDepl was not reached; 
     +tolerance was ', LowTolerance
                 ENDIF
              ELSE
                FindhFvec = FindhFvecNew
              ENDIF
              IF (buggybuggy.GE.MaxIterations.OR.MinDiffDepl.LE.
     +LowTolerance) exit
            END DO

            iHFparCombo = iHFparCombo + 1 !set to 0 at beg. of FindVariance
            HFparCombo(1,iHFparCombo) = BesthFpar
            HFparCombo(2,iHFparCombo) = fmultiply
            HFparCombo(3,iHFparCombo) = hpar
            HFparcombo(4,iHFparCombo) = MeanDays

        RETURN
        END SUBROUTINE DoFindhFpar
!-----------------------------------------------------------------------


!=================================================================
      SUBROUTINE DiffsAndNewHparVec()
      USE GlobalVars
      IMPLICIT NONE

      !Local variables
      !INTEGER,INTENT(IN) :: BMSYorH
      INTEGER :: IHpar,count,Iinc,ii,jj
      REAL(8),DIMENSION(Nhpar) :: Diffs
      REAL(8) :: Increment,NewInc

      !Calculate differences between 0.2 and AvgDepl(IFpar)
      Diffs = 0.0
      DO 100 IHpar = 1,Nhpar
        Diffs(Ihpar) = abs(AvgSteepness(Ihpar)-hparOriginal)
100   CONTINUE

      !If none of the Diffs is not within the desired tolerance level, make a new FindhFvec
      
      !MinDiff = min(Diffs)
      MinDiffh = 1000
      Besth = -1
      !BesthFpar = -1
      count = 0
      DO 101 Ihpar =1,Nhpar
        IF (Diffs(Ihpar).EQ.MinDiffh) THEN
          count = count+1
        ELSEIF (Diffs(Ihpar).LT.MinDiffh) THEN
          MinDiffh = Diffs(Ihpar)
          BestSteep = AvgSteepness(Ihpar)
          Besth = hparVec(Ihpar)
        ENDIF
101   CONTINUE

 !     IF (MinDiffh.LE.LowTolerance) THEN
        hpar = Besth
        !If hpar = 0.4, then what was the corresponding fmultiply, and hFpar?
        !Pass it up to ThreeWayCombo: MeanDays, hpar,hFpar,fmultiply.
        DO 4434 ii = 1,Nhpar
          jj = ii+Nhpar*(md-1)
          IF (HFparCombo(3,jj).EQ.Besth) THEN
           ThreeWayCombo(1,md) = HFparCombo(1,jj)  !hFpar
           ThreeWayCombo(2,md) = HFparCombo(2,jj)  !fmultiply
           ThreeWayCombo(3,md) = HFparCombo(3,jj)  !hpar
           ThreeWaycombo(4,md) = HFparCombo(4,jj)  !MeanDays
          ENDIF
4434    CONTINUE


        !FindVariance = 0
 !     ELSE
        !ReAssign FmultiplyVec
        !Keep FindVariance = 1
        !Need a DO WHILE or DO UNTIL or something like that. Can FORTRAN do that?

        !Increment between Fpars
!        Increment = abs(hparVec(2) - hparVec(1))
!       ! NewInc = 2*Increment/(Nhpar+2)
!       NewInc = 2*Increment/(Nhpar-1)
!       ! hparVecNew(1) = Besth-Increment+NewInc
!        hparVecNew(1) = Besth-Increment
!        hparVecNew(1) = min(1.0,max(0.2,hparVecNew(1)))
!        DO 200 Iinc = 2,Nhpar
!          hparVecNew(Iinc) = min(1.0,hparVecNew(Iinc-1) + NewInc)
!         ! PRINT *,'hparVecNew(Iinc) is ',hparVecNew(Iinc)
!200     CONTINUE
!      ENDIF

!      Do 201 Iinc = 1,NFmultiply
!        FmultvecNew(Iinc) = max(0.001,FmultvecNew(Iinc))
!201   CONTINUE

      !Need to go get the Fpar and fmultiply that goes along with this h and meandays
      !From hfmultcombo before the stuff below is written out

      WRITE(21,980) MinDiffh,BestSteep,Besth,ThreeWayCombo(1,md),
     +ThreeWayCombo(2,md),MeanDays
980   FORMAT('MinimumDiffhpar= ',F10.7,' ,Best R/R0',F10.7,
     +', Best h= ',F10.7,' using Fpar ',F10.7,
     +  ', and fmultiply = ',F10.7,' and MeanDays= ',F10.7)
      RETURN
      END SUBROUTINE DiffsAndNewHparVec
!=================================================================

!=================================================================
      SUBROUTINE CalcMaxCatches(FparIndex)
      USE GlobalVars
      IMPLICIT NONE

      !Local variables
      INTEGER,INTENT(IN) :: FparIndex
      INTEGER :: Isim,IYear,i
      REAL(8):: SumCumCatches

      i=0
      SumCumCatches = 0.0
      !SumCatches over years within each simulation
      DO 322 Isim = 1,Nsim
       DO 322 IYear = (MinPYr+MoveBurnInTime+BurnInTime+1),(MinPYr+
     +MoveBurnInTime+BurnInTime+Ntime)
         i = i+1
         SumCumCatches=SumCumCatches+TotalCatchBiomass(IYear,Isim)
322   CONTINUE
      
      !Average the SumCatches over the simulations
      AvgNegCumCatch(FparIndex) = -SumCumCatches/Nsim
      AvgCatch(FparIndex) = SumCumCatches/i

      RETURN
      END SUBROUTINE CalcMaxCatches
!=================================================================

!-----------------------------------------------------------------------
      SUBROUTINE DoFindMSY()
      USE GlobalVars

      !Local variables
      INTEGER :: IFpar,buggybuggy
             FindMSY = 1
              !Cycle through F values here
              !Need a DO loop that exits when MaxIterations or Tolerance are reached.
             buggybuggy = 0
             !FindMSYvec = FindhFvecOriginal

             FindMSYvec = Fvec
             DO
              buggybuggy = buggybuggy+1
              DO 5666 IFpar = 1,NFpar  !Come back and change this to a second F vector.
                Fpar = FindMSYVec(IFpar)
                CALL DoMAEInnards()
                
                !Calculate AvgSSB/SSB0
                CALL CalcAvgDepletion(IFpar)
                CALL CalcMaxCatches(IFpar)


5666          CONTINUE !end FindMaxCatches loop

  
              !Calculate min objfun (AvgNegCumCatch) and make a new vector around that

              CALL DiffsAndNewFMSYVec()
              !Go until reach tolerance or maxiterations OR use FIT function (will be pretty specific) with original vec. just between 0.06 and 0.20.
              !Make sure that Fpar = bestFpar at the end of this.
              
              IF (buggybuggy.GE.MaxIterations.OR.Increment.LE.
     +LowTolerance) THEN
                 FMSY = BestFMSY
                 MSY = BestMSY
                 SSBMSY = BestSSBMSY
                 DeplMSY = BestDeplMSY

                 PRINT *,'FMSY,MSY,SSBMSY,DeplMSY is ', FMSY,MSY,SSBMSY,
     +DeplMSY
                 
                 IF (Increment.GT.LowTolerance) THEN
                   PRINT *, 'LowTolerance for FMSY was not reached; 
     +tolerance was ', LowTolerance
                 ENDIF
              ELSE
                FindMSYvec = FindMSYvecNew
              ENDIF
              IF (buggybuggy.GE.MaxIterations.OR.Increment.LE.
     +LowTolerance) exit
            END DO

        !Can return final MSY information to a file right here:
        WRITE(28,912) SSBMSY,FMSY,MSY,DeplMSY
        WRITE(29,913) SSBMSY,FMSY,MSY,DeplMSY

        FindMSY = 0
912     FORMAT('Final SSBMSY is ',F20.5,' Final Fmsy is ',F10.7,
     +' Final MSY is ',F20.5,' Final DeplMSY is ',F10.7)
      
913     FORMAT('SSBMSY ',F20.5,' FMSY ',F10.7,' MSY ',F20.7,' DeplMSY ',
     +F10.7)
        RETURN
        END SUBROUTINE DoFindMSY
!-----------------------------------------------------------------------


!=================================================================
      SUBROUTINE DiffsAndNewFMSYVec()
      USE GlobalVars
      IMPLICIT NONE

      !Local variables
      !INTEGER,INTENT(IN) :: BMSYorH
      INTEGER :: IFpar,count,Iinc
      REAL(8),DIMENSION(NFpar) :: Diffs
      REAL(8) :: Increment,NewInc

      !Write things out to check out whether this algorithm really works
      !601   FORMAT(' AvgNegCumCatch ',F30.10,' Fpar ',F20.17)
      DO 3445 IFpar = 1,NFpar
        WRITE(28,601) AvgNegCumCatch(IFpar),FindMSYvec(IFpar)
3445  CONTINUE
      
      !Calculate differences between 0.2 and AvgDepl(IFpar)
      Diffs = 0.0
      DO 100 IFpar = 1,NFpar
        Diffs(IFpar) = AvgNegCumCatch(IFpar)
100   CONTINUE

      !If none of the Diffs is not within the desired tolerance level, make a new FindhFvec
      
      !MinDiff = min(Diffs)
      MinCatches = 1.0
      BestMSY = -1.0
      BestFMSY = -1.0
      BestSSBMSY = -1.0
      BestDeplMSY = -1.0
      
      count = 0
      DO 101 IFpar =1,NFpar
        IF (Diffs(IFpar).EQ.MinCatches) THEN
          count = count+1
        ELSEIF (Diffs(IFpar).LT.MinCatches) THEN
          MinCatches = Diffs(IFpar)
          BestMSY = AvgCatch(IFpar) !double check on this
          BestFMSY = FindMSYvec(IFpar)
          BestSSBMSY = AvgSSB(IFpar)
          BestDeplMSY = AvgDepl(IFpar)
        ENDIF
101   CONTINUE

        !Increment between Fpars
        Increment = abs(FindMSYVec(2) - FindMSYvec(1))
      IF (Increment.LE.LowTolerance) THEN
        Fpar = BestFMSY
        !FindVariance = 0
      ELSE

        NewInc = 2*Increment/(NFpar-1)
        FindMSYvecNew(1) = BestFMSY-Increment+NewInc
        FindMSYvecNew(1) = max(0.0,FindMSYvecNew(1))
        DO 200 Iinc = 2,NFpar
          FindMSYvecNew(Iinc) = FindMSYvecNew(Iinc-1) + NewInc
200     CONTINUE
      ENDIF

! Write out intermediate finding MSY information to new file

      WRITE(28,600) MinCatches,BestFMSY,BestMSY,BestSSBMSY,BestDeplMSY

600   FORMAT(' MinCatches is ',F30.10,' BestFMSY is ',F20.17,' BestMSY 
     +is ',F20.5,' BestSSBMSY is ',F20.5,' BestDeplMSY is ',F10.7)

601   FORMAT(' AvgNegCumCatch ',F30.10,' Fpar ',F20.17)

      RETURN
      END SUBROUTINE DiffsAndNewFMSYVec
!=================================================================


!------------------------------------------------------------------
!-------------------------------------------------------------
      SUBROUTINE MoreSummaryStats(IFpar)

      USE GlobalVars
      IMPLICIT NONE

      !Local variables
      INTEGER,INTENT(IN) :: IFpar
      INTEGER :: Isim,Iyear,Ibox

      !Maybe local, maybe switch to global:
      REAL(8),DIMENSION(Nsim) :: CumCatch
      REAL(8),DIMENSION(Nsim) :: ExploitRate
      REAL(8) :: MeanCumCatch,AvgExploitRate
      INTEGER :: SumCount
      REAL(8) :: Thresh
      INTEGER,DIMENSION(Nsim) :: Times,MyCount
      REAL(8) :: ProbBelow

      !Calculate long term catches by simulation and F
      CumCatch = 0.0
      ExploitRate = 0.0
      AvgExploitRate = 0.0
      MeanCumcatch = 0.0
      DO 1002 Isim = 1,Nsim
        DO 1001 IYear = MinPYr+MoveBurnInTime+BurnInTime+1,MinPYr+
     +MoveBurnInTime+BurnInTime+Ntime
          CumCatch(Isim) = CumCatch(Isim)+TotalCatchBiomass(IYear,Isim)
          ExploitRate(Isim) = ExploitRate(Isim)+
     +                       (TotalCatchBiomass(IYear,Isim)/
     +                        TotBiomass(IYear,Isim))/real(Nsim)
1001  CONTINUE
          !average over simulations as well
          AvgExploitRate = AvgExploitRate + ExploitRate(Isim)/real(Nsim)
          MeanCumCatch = MeanCumCatch + CumCatch(Isim)/real(Nsim)
1002  CONTINUE

      
      !Probability of falling below a biomass threshold
      Times = 0
      Thresh = 0.2
      SumCount = 0
      MyCount = 0
      ProbBelow = 0.0
      DO 1003 Isim = 1,Nsim
        DO 1004 IYear = MinPYr+MoveBurnInTime+BurnInTime+1,MinPYr+
     +MoveBurnInTime+BurnInTime+Ntime
         IF (TotSSB(IYear,Isim).LT.(Thresh*SSB0nospace)) THEN        
           Times(Isim) = Times(Isim) + 1
         ENDIF
1004    CONTINUE
          IF (Times(Isim).GT.0) THEN
            MyCount(Isim) = 1
          ENDIF
            SumCount = SumCount + MyCount(Isim)
1003  CONTINUE

      ProbBelow = real(SumCount)/real(Nsim)


      !write out files
      IF (IFpar.EQ.1) THEN
       WRITE(31,801) NFpar,Nsim
       WRITE(31,*) 'Fpar, Sim, CumulativeCatch, AvgExploitationRate'
      ENDIF
      
      DO 4444 Isim = 1,Nsim
        WRITE(31,802) FVec(IFpar),Isim,CumCatch(Isim),
     +ExploitRate(Isim)
4444  CONTINUE

      IF (IFpar.EQ.1) THEN
        WRITE(32,800) NFpar
        WRITE(32,*) 'Fpar,AvgCumCatch,AvgExploitationRate,
     +ProbBelow20pct'
      ENDIF

      WRITE(32,803) FVec(IFpar),MeanCumCatch,AvgExploitRate,ProbBelow

800   FORMAT('NFpar ',I4)
801   FORMAT('NFpar ',I4,' Nsim ',I4)
802   FORMAT(F10.7,1x,I4,1x,F30.10,1x,F15.8)
803   FORMAT(F10.7,1x,F30.10,1x,F15.8,1x,F15.10)


      END SUBROUTINE MoreSummaryStats
!------------------------------------------------------------------




