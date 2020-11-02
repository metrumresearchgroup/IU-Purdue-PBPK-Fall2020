//Voriconazole PBPK model for a typical adult male ; modified for population simulation

[PARAM] 
//Tissue volumes (L); source: https://www.ncbi.nlm.nih.gov/pubmed/14506981
// Vad = 21.17 //adipose
// Vbo = 10.58 //bone
// Vbr = 1.5 //brain
// Vla_int = 0.47  //gut wall
// VguLumen = 0.35 //gut lumen
// Vhe = 0.33 //heart
// Vki = 0.31 //kidneys
// Vli = 1.8 //liver
// Vlu = 0.5 //lungs
// Vmu = 29 //muscle
// Vsp = 0.15 //spleen
// Var = 5.6  //blood


Vbo =     10.57560469
Vbr =      1.50901160
Vgo =      0.03351116
Vhe =      0.36657898
Vki =      0.40953699
Vla_int =  0.47888909
Vli =      2.25799950
Vpa =      0.16773611
Vsk =      3.24901911
Vsm_int =  0.83627363
Vst =      0.19972847
Vth =      0.03057129
Vln =      0.27029895
Vve =      0.98605418
Var =      0.32868473
Vlu =      1.15214262
Vmu =     28.74821233
Vsp =      0.22992116
Vad =     21.17022541

//Tissue blood flows (L/h); Cardiac output = 6.5 (L/min); source: https://www.ncbi.nlm.nih.gov/pubmed/14506981
// Qad = 0.05*6.5*60
// Qbo = 0.05*6.5*60
// Qbr = 0.12*6.5*60
// Qgu = 0.15*6.5*60 
// Qhe = 0.04*6.5*60
// Qki = 0.19*6.5*60
// Qmu = 0.17*6.5*60
// Qsp = 0.03*6.5*60
// Qha = 0.065*6.5*60  //hepatic artery
// Qlu = 6.5*60        //same as cardiac output

Qbo =     18.89327843
Qbr =     44.83848370
Qgo =      0.20399038
Qhe =     14.33897912
Qki =     75.36942527
Qla_int = 15.18643845
Qpa =      3.83185620
Qsk =     19.54906011
Qsm_int = 38.05515713
Qst =      3.79252048
Qad =     20.02609032
Qmu =     66.84352492
Qsp =     12.08118319
Qha =     23.15737912
Qth =      0.70951492
Qln =      6.58156418

//partition coefficients estimated by Poulin and Theil method https://jpharmsci.org/article/S0022-3549(16)30889-9/fulltext
Kpad = 9.89  //adipose:plasma
Kpbo = 7.91  //bone:plasma
Kpbr = 7.35  //brain:plasma
Kpgu = 5.82  //gut:plasma
Kphe = 1.95  //heart:plasma
Kpki = 2.9   //kidney:plasma
Kpli = 4.66  //liver:plasma
Kplu = 0.83  //lungs:plasma
Kpmu = 2.94  //muscle:plasma; optimized
Kpsk = 3.53  //skin:plasma
Kpsp = 2.96  //spleen:plasma
Kpln = 4     //calculated as average of non adipose Kps
Kpth = 4     //thymus:plasma
Kpgo = 4     //gonad:plasma
BP = 1       //blood:plasma ratio

//other parameters
WEIGHT = 73 //(kg)
ka = 0.849  //absorption rate constant (/hr) 
fup = 0.42   //fraction of unbound drug in plasma

//in vitro hepatic clearance parameters http://dmd.aspetjournals.org/content/38/1/25.long
fumic = 0.711 //fraction of unbound drug in microsomes
MPPGL = 30.3  //adult mg microsomal protein per g liver (mg/g)
VmaxH = 40    //adult hepatic Vmax (pmol/min/mg)
KmH = 9.3     //adult hepatic Km (uM)

//renal clearance  https://link.springer.com/article/10.1007%2Fs40262-014-0181-y
CL_Ki = 0.096; //(L/hr) renal clearance


[CMT] 
GUTLUMEN GUT ADIPOSE BRAIN HEART BONE 
  KIDNEY LIVER LUNG MUSCLE SPLEEN
  SKIN THYMUS LYMPHN GONAD
  ART VEN
  
  
  [MAIN]
//additional volume derivations
double Vgu = Vsm_int + Vla_int + Vst;

//additional blood flow derivation
double Qgu = Qsm_int + Qla_int + Qst;
double Qli = Qgu + Qsp + Qha + Qpa;
double Qlu = Qli + Qki + Qbo + Qhe + Qmu + Qad + Qbr + Qth + Qln + Qsk + Qgo;

//intrinsic hepatic clearance calculation
double CL_Li = ((VmaxH/KmH)*MPPGL*Vli*1000*60*1e-6) / fumic; //(L/hr) hepatic clearance


[ODE]
//Calculation of tissue drug concentrations (mg/L)
double Cadipose = ADIPOSE/Vad;
double Cbone = BONE/Vbo;
double Cbrain = BRAIN/Vbr; 
double Cheart = HEART/Vhe; 
double Ckidney = KIDNEY/Vki;
double Cliver = LIVER/Vli; 
double Clung = LUNG/Vlu; 
double Cmuscle = MUSCLE/Vmu;
double Cspleen = SPLEEN/Vsp;
double Carterial = ART/Var;
double Cvenous = VEN/Vve;
double Cthymus = THYMUS/Vth;
double Clymphn = LYMPHN/Vln;
double Cgonad = GONAD/Vgo;
double Cgut = GUT/Vgu;
double Cskin = SKIN/Vsk;

//ODEs
dxdt_GUTLUMEN = -ka*GUTLUMEN;
dxdt_GUT = ka*GUTLUMEN + Qgu*(Carterial - Cgut/(Kpgu/BP)); 
dxdt_ADIPOSE = Qad*(Carterial - Cadipose/(Kpad/BP)); 
dxdt_BRAIN = Qbr*(Carterial - Cbrain/(Kpbr/BP));
dxdt_HEART = Qhe*(Carterial - Cheart/(Kphe/BP));
dxdt_KIDNEY = Qki*(Carterial - Ckidney/(Kpki/BP)) - CL_Ki*(fup*Ckidney/(Kpki/BP));
dxdt_LIVER = Qgu*(Cgut/(Kpgu/BP)) + Qsp*(Cspleen/(Kpsp/BP)) + Qha*(Carterial) - Qli*(Cliver/(Kpli/BP)) - 
  CL_Li*(fup*Cliver/(Kpli/BP)); 
dxdt_LUNG = Qlu*(Cvenous - Clung/(Kplu/BP));
dxdt_MUSCLE = Qmu*(Carterial - Cmuscle/(Kpmu/BP));
dxdt_SPLEEN = Qsp*(Carterial - Cspleen/(Kpsp/BP));
dxdt_BONE = Qbo*(Carterial - Cbone/(Kpbo/BP));
dxdt_VEN = Qad*(Cadipose/(Kpad/BP)) + Qbr*(Cbrain/(Kpbr/BP)) +
  Qhe*(Cheart/(Kphe/BP)) + Qki*(Ckidney/(Kpki/BP)) + Qli*(Cliver/(Kpli/BP)) + 
  Qmu*(Cmuscle/(Kpmu/BP)) + Qbo*(Cbone/(Kpbo/BP)) + Qsk*(Cskin/(Kpsk/BP)) +
  Qth*(Cthymus/(Kpth/BP)) + Qln*(Clymphn/(Kpln/BP)) + Qgo*(Cgonad/(Kpgo/BP)) -
    Qlu*Cvenous;
dxdt_ART = Qlu*(Clung/(Kplu/BP) - Carterial);
dxdt_SKIN = Qsk*(Carterial - Cskin/(Kpsk/BP));
dxdt_THYMUS = Qth*(Carterial - Cthymus/(Kpth/BP));
dxdt_LYMPHN = Qln*(Carterial - Clymphn/(Kpln/BP));
dxdt_GONAD = Qgo*(Carterial - Cgonad/(Kpgo/BP));


[TABLE]
capture CP = Cvenous/BP;
