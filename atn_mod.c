#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

static int const N = 10000; //number of agents
static int const M = 5;		//number of links

static double const alfa = -2.1;  //heavy tail exponent
static double const xmin = 0.001; //lower cutoff of the power law
static double const xmax = 1;	  //higher cutoff of the power law

static int const Tmax = 3000;	// maximum execution time
static int const TWindow = 500; // time window length for steady state evaluation
static int const ntrial = 10;	//number of trials

static double const LambdaMin = 0.0;
static double const LambdaMax = 1.0;
static double const LambdaStep = 0.01;											 //Lambda step
static int const LambdaSize = (int)(((LambdaMax - LambdaMin) / LambdaStep) + 1); //size of the lambda space

static double const EtaImin = 0;										 //minimum EtaI
static double const EtaImax = 15;										 //maximum EtaI
static double const EtaIStep = 0.5;										 //EtaI step
static int const EtaISize = (int)(((EtaImax - EtaImin) / EtaIStep) + 1); //size of the eta space

static int const Ni = (int)(0.01 * N); //initial number of infected
static double const mu = 0.1;		   //recover rate
static double const etaS = 1;		   //EtaS

int i, j, nt, i_time, ietaI, ilambda, NumInf, rInt;
int agentstate[N + 1], agentstatenext[N + 1];
double r, etaI, lambda, Activity, InfectedOverATrial, InfectedOverTrials;
double AvgInfectedEtaLambda[LambdaSize + 1], a[N + 1];
double results[3001];
char fileNetName[40], fileAvgName[40];
FILE *fileNet, *fileAvgInfLambdaEta;

//to be resumed if interested in the infected trend of the single agents over time - look also at the printf in the code
//char fileNumInfName[40];
//FILE*fileNumInf

// double att[N+1],attivi[N+1]; //to be resumed if interested in the activity patterns

int main()
{

	srand(time(NULL));

	//create the activity distribution
	sprintf(fileNetName, "network.dat");
	fileNet = fopen(fileNetName, "w");
	for (i = 1; i <= N; i++)
	{
		r = rand() / (double)RAND_MAX;
		a[i] = pow((pow(xmax, alfa + 1) - pow(xmin, alfa + 1)) * r + pow(xmin, alfa + 1), (1 / (alfa + 1)));
		fprintf(fileNet, "%lf\n", a[i]);
	}
	fclose(fileNet);

	//sprintf(fileNumInfName,"InfectedTrend.dat");  //file with the trends of the number of infected, trial by trial
	//fileNumInf=fopen(fileNumInfName,"w");

	double etaI = 1;
	printf("eta:%lf\n", etaI);
	lambda = LambdaMin;

	double lambda = 1;
	// for (ilambda = 1; ilambda <= LambdaSize; ilambda++)
	// {
	printf("lambda:%lf\n", lambda);
	InfectedOverTrials = 0.0;
	// for (nt = 1; nt <= ntrial; nt++) //cycle on trials
	// {
	//printf("nt:%i\n",nt);

	for (i = 1; i <= Ni; i++)
	{
		agentstate[i] = 1;
	} //initialize agent state to Infected
	for (i = Ni + 1; i <= N; i++)
	{
		agentstate[i] = 0;
	} //initialize agent state to Susceptible

	InfectedOverATrial = 0.0;
	for (i_time = 1; i_time <= Tmax; i_time++) //cycle on the duration of a trial
	{
		/* this cycle can be resumed if one wants to store the activity patterns
                	for (i=1;i<=N;i++) //create the vector of active nodes - also this could be eliminated and integrated in the main epidemic cycle - I am keeping it for clarity
                 	{
                            Activity=(1-agentstate[j])*(etaS*a[j])+agentstate[j]*(etaI*a[j]);
                            r=rand()/(double)RAND_MAX;
                            if (Activity>r)
                            {
                                attivi[i]=1;
                            } else
                                attivi[i]=0;
                    }
                    */

		for (i = 1; i <= N; i++) //copy the old agent state in the new one
		{
			agentstatenext[i] = agentstate[i];
		}

		for (i = 1; i <= N; i++) //cycle on agents - epidemics on network
		{
			if (agentstate[i] == 1) //if infected it may recover
			{
				r = rand() / (double)RAND_MAX;
				if (r < mu)
				{
					agentstatenext[i] = 0;
				}
			}

			Activity = a[i];
			r = rand() / (double)RAND_MAX;
			if (Activity > r) //if active - if the activity cycle is resumed, then test on attivi[i]
			{
				for (j = 1; j <= M; j++) //cycle on the neighbors of the active agent
				{
					r = rand() / (double)RAND_MAX; //pick a neighbor - it can be improved by avoiding self- and multiple contacts
					r = r * (N - 1) + 1;
					rInt = (int)r;

					if (agentstate[i] == 1) //if infected then it can infect the neighbor
					{
						r = rand() / (double)RAND_MAX;
						if ((agentstate[rInt] == 0) && (r < lambda))
						{
							agentstatenext[rInt] = 1;
						}
					}
					else //if not infected, then it can be infected by the neighbor, if the latter is infected
					{
						r = rand() / (double)RAND_MAX;
						if ((agentstate[rInt] == 1) && (r < lambda))
						{
							agentstatenext[i] = 1;
						}
					} //end if not infected
				}	  //end cycle on neighbors
			}		  //end if active
		}			  //end cycle on agents - epidemics on network

		NumInf = 0;
		for (i = 1; i <= N; i++) //copy new state into the current state
		{
			agentstate[i] = agentstatenext[i];
			if (agentstate[i] == 1)
				NumInf++; //count the infected
		}
		results[i_time] = NumInf;
		// if (i_time > (Tmax - TWindow)) //sum and compute the average of the infected at steady state
		// {
		// 	InfectedOverATrial += (double)NumInf;
		// }

	} //it closes i_time
	  // InfectedOverATrial /= (double)TWindow;
	  // InfectedOverTrials += InfectedOverATrial;
	  // } //it closes nt - cycle on trials
	// 	AvgInfectedEtaLambda[ilambda] = InfectedOverTrials / (double)ntrial;
	// 	printf("%lf \n", AvgInfectedEtaLambda[ilambda]);
	// 	lambda += LambdaStep;
	// } //it closes ilambda

	//fclose(fileNumInf);

	sprintf(fileAvgName, "AvgInfectedEtaLambda.dat");
	fileAvgInfLambdaEta = fopen(fileAvgName, "w");

	for (j = 1; j <= 3000; j++)
	{
		fprintf(fileAvgInfLambdaEta, "%lf\n", results[j]);
	}

	fclose(fileAvgInfLambdaEta);

} //it closes main
