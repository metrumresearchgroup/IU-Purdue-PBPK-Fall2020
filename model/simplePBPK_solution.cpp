[PARAM]
//volumes (L)
Vmu = 29
Var = 1.65
Vve = 3.9
Vlu = 0.5
Vli = 1.8

//blood flows
Qmu = 0.17 * 6.5 * 60
Qli = 0.245 * 60 * 6.5
Qlu = 6.5 * 60

//partition coefficients
Kpmu = 1
Kpli = 1
Kpre = 1
Kplu = 1

//other
BP = 1
WEIGHT = 73
Cl_hepatic = 10
fup = 0.5


[CMT]
VEN ART MUSCLE LIVER LUNG REST

[MAIN]
//transformed parameters
double Vre = WEIGHT - (Vli + Vlu + Vve + Var + Vmu);
double Qre = Qlu - (Qli + Qmu);

[ODE]
//concentrations
double Carterial = ART / Var;
double Cmuscle = MUSCLE / Vmu;
double Cvenous = VEN / Vve;
double Cliver = LIVER / Vli;
double Crest = REST / Vre;
double Clung = LUNG / Vlu;

//ODEs
dxdt_MUSCLE = Qmu * (Carterial - Cmuscle / (Kpmu / BP));
dxdt_REST = Qre * (Carterial - Crest / (Kpre / BP));
dxdt_LUNG = Qlu * (Cvenous - Clung / (Kplu / BP));
dxdt_LIVER = Qli * (Carterial - Cliver / (Kpli / BP)) - Cl_hepatic * fup * (Cliver / (Kpli/BP));
dxdt_ART = Qlu * (Clung / (Kplu / BP) - Carterial);
dxdt_VEN = Qmu * (Cmuscle / (Kpmu/BP)) + Qli * (Cliver / (Kpli/BP)) + Qre * (Crest / (Kpre/BP)) - Qlu * Cvenous;
