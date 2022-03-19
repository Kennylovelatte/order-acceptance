// -------------------------------------------------------------- -*- C++ -*-
// File: main.cpp
// Version 12.8.0
// This file is for order acceptance problem
// Unrelated parallel machines with machine- and sequence- dependent 
// setup times and machine available time
// objective: maximize the difference between revenues and makespan
// --------------------------------------------------------------------------
// Included solver: 
// CPLEX 128 & GRUOBI 8.1
// Other package:
// Concorde  -- http://www.math.uwaterloo.ca/tsp/index.html
// --------------------------------------------------------------------------
// Licensed Materials - Property of IBM
// 5725-A06 5725-A29 5724-Y48 5724-Y49 5724-Y54 5724-Y55 5655-Y21
// Copyright IBM Corporation 2000, 2017. All Rights Reserved.
//
// US Government Users Restricted Rights - Use, duplication or
// disclosure restricted by GSA ADP Schedule Contract with
// IBM Corp.
// --------------------------------------------------------------------------
//

#include "ilcplex/ilocplex.h" /* cplex */
#include "gurobi_c++.h"       /* gurobi */
extern "C"
{
#include <concorde.h>
} /* concorde */

#include <iostream>
#include <sstream>

#include <stdio.h>   /* printf, scanf, puts, NULL */
#include <stdlib.h>  /* srand, rand */
#include <time.h>    /* time */
#include <limits.h>  /* limits */
#include <vector>    /* vector */
#include <algorithm> /* random_shuffle */
#include <math.h>    /* log, ceil */
#include <fstream>
#include <cstdlib>
#include <limits>
#include <float.h>
#include <cmath>
#include <cstring>
#include <ctime>
#include <iomanip>
#include <queue>
#include <stack>
#include <string>
#include <ctype.h>
#include <unistd.h>


// / ///////////////////////////////////////////////////////////

ILOSTLBEGIN
typedef IloArray<IloNumVarArray> NumVarMatrix;
typedef IloArray<NumVarMatrix> NumVar3Matrix;
//typedef IloArray<NumVar3Matrix> NumVar4Matrix;
typedef IloArray<IloNumArray> NumMatrix;
typedef IloArray<NumMatrix> Num3Matrix;

using std::cout; 
using std::endl; 


// / //////////////PARAMETER SETTINGS////////////////////////////

int M;//the number of machines
int N;//the number of jobs
int L;//a large enough integer
int **PT_ij;// the processing time for job j on machine i
int ***ST_ijk;//the setup time for arc(j,k) on machine i
int *U_i;//the available time for machine i
int *R_j;//the revenue for accepted job j

// / //////////////FUNCTIONS ABOUT DATA GENERATION///////////////

//data generation
void generate_data();
//clear data, pointers to NULL
void clear_data();
//generate random coordinates for N + 1 jobs 
vector<int> generate_random_coordinate();

// / //////////////GENERAL FUNCTIONS/////////////////////////////

//get variable array values 
vector<double> get_var_values(
    const IloCplex cplex,
    const IloNumVarArray varArray,
    const int size1
);
//get variable matrix2 values
vector< vector<double> > get_var_values(
    const IloCplex cplex,
    const NumVarMatrix varMatrix,
    const int size1,
    const int size2
);

//write results of MIP model and decomposition model
void writeResult(
    const int Number_Instance,
    const double MIP_obj,
    const double MIP_bestBound,
    const double MIP_time,
    const double MIP_finalGap,
    const double DM_obj,
    const double DM_time,
    const double DM_finalGap,
    const bool DM_find_feasible
);

// / /////////////FUNCTIONS ABOUT MIP MODEL/////////////////////

//solve the MIP model
//input: time limit
//return: obj, best bound, solution time, final gap
void MIP(
    const int MIP_timeLimit,
    double &MIP_obj,
    double &MIP_bestBound,
    double &MIP_time,
    double &MIP_finalGap
);

// / //////////////FUNTIONS ABOUT DECOMPOSITION////////////////

//decomposition method
//input: time limit
//return: obj, solution time
void DM(
    const int DM_timeLimit,
    double &DM_obj,
    double &DM_time,
    bool &DM_find_feasible
);

//master problem solving process
//input: iteration number, DM time limit, DM start timepoint
//empty vectors to store value of Zj, Xij, Cmax
//return: master problem solution {Zj, Xij, Cmax, obj}
void MP_DM(
    const int numIter,
    const int DM_timeLimit,
    const time_t &startTime,
    const vector< vector<double> > &Zj_Value_h,
    const vector<double> &CMAX_Value_h,
    const vector<double> &MP_obj_h,
    const vector<double> &SP_Cmax_h,
    vector<double> &Zj_Value,
    vector< vector<double> > &Xij_Value,
    double &CMAX_Value, 
    double &MP_obj,
    bool &is_MP_feasible,
    bool &MP_reach_time_limit,
    bool &MP_solution_unknown
);

//subproblem solving process
//input: valued Zj, Xij
//check whether there is a global feasible solution
bool SP_DM(
    const vector<double> &Zj_Value,
    const vector< vector<double> > &Xij_Value,
    const int DM_timeLimit,
    const time_t &startTime,
    double &SP_feasible_Cmax,
    bool &SP_reach_timelimit
);

void SPM(
    const int SP_IterNum,
    const int DM_timeLimit,
    const time_t &startTime,
    const vector<int> &acceptJobs,
    const vector<int> &initialAssignment,
    const vector<double> &SPS_Cmax_h,
    const vector< vector< double> > &SPS_Oi_h,
    const vector<bool> &feasible_h,
    const double best_feasible_Cmax,
    bool &is_SPM_feasible,
    vector<double>  &SPM_Cmax_h,
    vector< vector< vector< double> > > &SPM_Xij_h,
    bool &SPM_reach_timelimit,
    bool &SPM_solution_unknown
);

void SPS(
    const int SP_IterNum,
    const int DM_timeLimit,
    const time_t &startTime,
    const vector<int> &acceptJobs,
    const vector< vector< vector< double> > > &SPM_Xij_h,
    bool &is_SPS_feasible,
    double &SPS_Cmax,
    vector< vector<double> > &SPS_Oi_h,
    bool &SPS_reach_timelimit,
    bool &SPS_solution_unknown
);


void checkSP(
    const vector<double> &Zj_Value
);


// / //////////////////////////////////////////////////////////

int main(int, char **)
{
    int No_Instance = 10;
    for (int numberInstance = 0; numberInstance != No_Instance; ++numberInstance)
    {
        generate_data();

        // / //////////////MIP MODEL/////////////////////////////////////////

        //MIP model input 
        int MIP_timeLimit = 1800;
        //MIP model output
        double MIP_obj, MIP_bestBound, MIP_time, MIP_finalGap;
        /*MIP(MIP_timeLimit, MIP_obj, MIP_bestBound, MIP_time, MIP_finalGap);
        //print MIP model output
        cout << "MIP MODEL OUTPUT:" <<  '\n'
        << '\t' << "obj : "  << MIP_obj << '\n' 
        << '\t' << "bb  : "  << MIP_bestBound << '\n'
        << '\t' << "time: "  << MIP_time << '\n'
        << '\t' << "gap : "  << MIP_finalGap << '\n'
        << endl;*/

        // / ////////////DECOMPOSITION METHOD/////////////////////////////// 

        //Decomposition method input
        int DM_timeLimit = 1800;
        //Decomposition method output
        double DM_obj, DM_time;
        double DM_finalGap;//decomposition method gap with MIP model
        bool DM_find_feasible;//if could get a feasible solution within time limit 
        DM(DM_timeLimit,DM_obj,DM_time,DM_find_feasible);

        if (DM_find_feasible)
        {
                  
            //DM_finalGap = (DM_obj - MIP_obj)/MIP_obj;
            //print decomposition method output
            cout << "DECOMPOSITION METHOD OUTPUT:" <<  '\n'
            << '\t' << "obj : "  << DM_obj << '\n' 
            << '\t' << "time: "  << DM_time << '\n'
            //<< '\t' << "gap : "  << DM_finalGap << '\n'
            << endl;

        }
        else
        {
            cout << "DECOMPOSITION METHOD OUTPUT:" <<  '\n'
            << '\t' << " time limit ! No solution has found !"  
            << endl;
        }


        //cout << "MIP MODEL obj : "  << MIP_obj << endl;

        // / ////////////////////////////////////////////////////////////// 

        writeResult(numberInstance,MIP_obj,MIP_bestBound,MIP_time,MIP_finalGap,DM_obj,DM_time,DM_finalGap,DM_find_feasible);

        /*if (round(MIP_obj) != round(DM_obj))
        {
            cout << " not equal! error! "<< endl;
            break;
        }*/

        clear_data();
    }


    return 0;
}

void generate_data()
{
    sleep(2);
    srand(time(NULL));
    // / ////////// M  N  L  ///////////////////////////
    M = 2;
    N = 250;
    L = 99999;
    
    cout << "M :" << M << '\n' 
    << "N :" << N << '\n' << endl;

    // / ////////MACHINE AVAILABLE TIME////////////////
    U_i = new int[M];

    for (int i = 0; i != M; ++i)
    {
        int start_num = (N/M - 3 > 1) ? (N/M - 3) : 1;
        int max_length = 50*(start_num + rand()%(N/M - start_num + 1));   
        
        U_i[i] = (max_length > 100 ? max_length : 100);
        
    }

    cout << "Ui :" << endl;
    for (int i = 0; i != M; ++i)
    {
        cout << "U[" << i << "] = " << U_i[i] << '\n'; 
    }
    cout << endl;

    // / /////////REVENUE FOR JOB//////////////////////
    R_j = new int[N];

    for (int j = 0; j != N; ++j)
    {
        R_j[j] = (rand()%(101)) + 100; // [100, 200]
    }

    cout << "Rj :" << '\n';
    for (int j = 0; j != N; ++j)
    {
        cout << '\t' << "R[" << j << "] = " << R_j[j] ; 
    }
    cout << '\n' << endl;

    // / ///////////PROCESSING TIME //////////////////
    PT_ij = new int *[M];
    for (int i = 0; i != M; ++i)
    {
        PT_ij[i] = new int[N+1];
    }

    for (int i = 0; i != M; ++i)
    {
        for (int j = 0; j != N; ++j)
        {
            PT_ij[i][j] = (rand()%(51)) + 30; // [30, 80]
        }
        PT_ij[i][N+1] = 0;
    }

    cout << "PT_ij :" << '\n';
    for (int i = 0; i != M; ++i)
    {
        for (int j = 0; j != N+1; ++j)
        {
            cout << '\t' << setw(3) << PT_ij[i][j];
        }
        cout << endl;
    }
    cout << endl;

    // / /////////////SETUP TIME//////////////////////
    ST_ijk = new int **[M];
    for (int i = 0; i != M; ++i)
    {
        ST_ijk[i] = new int *[N+1];
        for (int j = 0; j != N+1; ++j)
        {
            ST_ijk[i][j] = new int[N+1];
        }
    }   

    for (int i = 0; i != M; ++i)
    {       
        vector<int> jobAxis1 = generate_random_coordinate();
        vector<int> jobAxis2 = generate_random_coordinate();
        vector<int> jobAxis3 = generate_random_coordinate();
        vector<int> jobAxis4 = generate_random_coordinate();
        for (int j = 0; j != N; ++j)
        {
            for (int k = j+1; k != N + 1; ++k)
            {               
                
                int diffX1 = abs(jobAxis1[j] - jobAxis1[k]);
                int diffY1 = abs(jobAxis2[j] - jobAxis2[k]);
                int diffX2 = abs(jobAxis3[j] - jobAxis3[k]);
                int diffY2 = abs(jobAxis4[j] - jobAxis4[k]);

                ST_ijk[i][j][k] = round(10 + 0.05*(diffX1 + diffY1));//[10, 20]
                ST_ijk[i][k][j] = round(10 + 0.05*(diffX2 + diffY2));//[10, 20]
                
            }
            ST_ijk[i][j][j] = 0;                  
        }
        for (int j = 0; j != N+1; ++j)
        {
            ST_ijk[i][j][N] = 0;
        }
    }

    cout << "ST_ijk :" << '\n';
    for (int i = 0; i != M; ++i)
    {
        for (int j = 0; j != N+1; ++j)
        {
            for (int k = 0; k != N+1; ++k)
            {
                cout << setw(3) << ST_ijk[i][j][k] << " ";
            }
            cout << endl;
        }
        cout << endl;
    } 
    cout << endl;
}

void clear_data()
{
    // / ////////////DELETE POINTERS/////////////////////////
    delete[] U_i;
    U_i = NULL;

    delete[] R_j;
    R_j = NULL;

    for (int i = 0; i != M; ++i)
    {
        delete[] PT_ij[i];
    }
    delete[] PT_ij;
    PT_ij = NULL;

    for (int i = 0; i != M; ++i)
    {
        for (int j = 0; j != N+1; ++j)
        {
            delete[] ST_ijk[i][j];
        }
        delete[] ST_ijk[i];
    }
    delete[] ST_ijk;
    ST_ijk = NULL;
    // / /////////////////////////////////////////////////////
}

vector<int> generate_random_coordinate()
{
    sleep(1);
    srand(time(NULL));
    vector<int> coordinate;
    for (int j = 0; j != N+1; ++j)
    {
        coordinate.push_back(rand()%(101));//[0,100]
    }
    return coordinate;
}


vector<double> get_var_values(
    const IloCplex cplex,
    const IloNumVarArray varArray,
    const int size1
)
{
    vector<double> temp;
    for (int i = 0; i != size1; ++i)
    {
        if (cplex.getValue(varArray[i]) < 0.1)
        {
            temp.push_back(0.0);
        }
        else
        {
            temp.push_back(cplex.getValue(varArray[i]));
        }
        
    }
    return temp;
}

vector< vector<double> > get_var_values(
    const IloCplex cplex,
    const NumVarMatrix varMatrix,
    const int size1,
    const int size2
)
{
    vector< vector<double> > temp;
    for (int i = 0; i != size1; ++i)
    {
        temp.push_back(get_var_values(cplex,varMatrix[i],size2));
    }
    return temp;
}

void writeResult(
    const int Number_Instance,
    const double MIP_obj,
    const double MIP_bestBound,
    const double MIP_time,
    const double MIP_finalGap,
    const double DM_obj,
    const double DM_time,
    const double DM_finalGap,
    const bool DM_find_feasible
)
{

    cout << "----------- Now we write the results into a .txt file ------------" << endl;
    // write the results to txt, each instance with a single txt file

    char sFileName[300];

    sprintf(sFileName, "/Users/kenny/Desktop/test/Instance_n=%d_m=%d_No=%d.txt", N, M, Number_Instance);

    ofstream o_file;
    o_file.open(sFileName);

    //start writing in the txt
    //o_file << "Results of MIP model: obj, time" << endl;
    //o_file << " obj = " << MIP_obj << "  time = " << MIP_time << endl;
    //o_file << "Gap and the best bound:" << endl;
    //o_file << " BestBound = " << MIP_bestBound << "  Gap = " << MIP_finalGap << endl;

    //o_file << endl;
    //o_file << "------------------------------------------------------------" << endl;


    if (DM_find_feasible)
    {
        o_file << "Results of Decomposition method: obj, time" << endl;
        o_file << "obj = " << DM_obj << " time = " << DM_time; 
        //<< "  Gap = " << DM_finalGap;
        o_file << endl;
    }
    else
    {
        o_file << "Results of Decomposition method: obj, time" << endl;
        o_file << " no feasible solution found !" << endl;
        o_file <<  " time = " << DM_time << endl;
        o_file << endl;
    }


    // write the parameters into txt
    //output data
    o_file << "------------------------------------------------------------" << endl;
    o_file << "M = " << M << endl;
    o_file << endl;

    o_file << "N = " << N << endl;
    o_file << endl;

    o_file << "Machine available time: Ui = " << endl;
    for (int i = 0; i != M; ++i)
    {
        o_file << setw(4) << U_i[i] << '\t';
    }
    o_file << endl;

    o_file << "Job revenue: Rj = " << endl;
    for (int j = 0; j != N; ++j)
    {
        o_file << setw(4) << R_j[j] << '\t';
    }
    o_file << endl;

    o_file << "Job processing time on each machine: PT_ij = " << endl;
    for (int i = 0; i != M; ++i)
    {
        for (int j = 0; j != N; ++j)
        {
            o_file << setw(4) << PT_ij[i][j] << '\t';
        }
        o_file << endl;
    }
    o_file << endl;

    o_file << "Setup time: ST_ijk = " << endl;
    for (int i = 0; i != M; ++i)
    {
        for (int j = 0; j != N + 1; ++j)
        {
            for (int k = 0; k != N + 1; ++k)
            {
                o_file << setw(4) << ST_ijk[i][j][k] << '\t';
            }
            o_file << endl;
        }
        o_file << endl;
    }
    o_file << endl;

    // close and save the file
    o_file.close();

    cout << "----------- Now we finished the results writing ------------------" << endl;

}

void MIP(
    const int MIP_timeLimit,
    double &MIP_obj,
    double &MIP_bestBound,
    double &MIP_time,
    double &MIP_finalGap
)
{
    cout << "-----------------MIP MODEL STARTS----------------" << endl;
    IloEnv env;
    try {
        
        IloModel model(env);
        
        // /////////DECISION VARIABLES/////////////////////////////

        NumVarMatrix xij(env, M);
        for (int i = 0; i != M; i++)
        {
            xij[i] = IloNumVarArray(env,N+1,0,1,ILOINT);
        }

        IloNumVarArray zj(env,N,0,1,ILOINT); 

        NumVar3Matrix yijk(env,M);
        for (int i=0;i!=M;i++)
        {
            yijk[i] =NumVarMatrix(env,N+1);
            for(int j=0;j!=N+1;j++)
            {
                yijk[i][j]=IloNumVarArray(env,N+1,0,1,ILOINT);
            }
        } 

        IloNumVarArray Cj(env,N+1,0,INFINITY,ILOINT);  

        IloNumVar Cmax(env);        
        // / /////////////CONSTRAINTS///////////////////////////

        // Objective Function (1a)
        {
            IloExpr obj(env);
            for (int j = 0; j != N; j++)
            {
                obj += R_j[j]*zj[j];                 
            }
            obj -= Cmax;            
            model.add(IloMaximize(env, obj));
            obj.end();
        }

        //constraint (1b)
        {
            for (int j = 0; j != N; j++)
            {
                IloExpr sum(env);
                for (int i = 0; i != M; i++)
                {
                    sum += xij[i][j];
                }
                model.add(sum == zj[j]);
                sum.end();
            }
        }

        //constraint (1c) (1m) (1k)
        {
            for (int j = 0; j != N; j++)
            {
                for (int i = 0; i != M; i++)
                {
                    //model.add(xij[i][j] <= zj[j]);
                }
                model.add(Cmax >= Cj[j]);
                model.add(Cj[j] <= zj[j] * L);
            }
        }
        

       
        
        //constraint (1d) 
        {
            for (int i = 0; i != M; i++)
            {
                model.add(xij[i][N] == 1);
            }
        }
        
        //constraint (1e)
        {
            for (int i = 0; i != M; i++)
            {
                for (int k = 0; k != N + 1; k++)
                {
                    IloExpr sum(env);
                    for (int j = 0; j != N + 1; j++)
                    {
                        if (j != k)
                            sum += yijk[i][j][k];
                    }
                    model.add(xij[i][k] == sum);
                    sum.end();
                }
            }
        }
        
        //constraint (1f)
        {
            for (int i = 0; i != M; i++)
            {
                for (int j = 0; j != N + 1; j++)
                {
                    IloExpr sum(env);
                    for (int k = 0; k != N + 1; k++)
                    {
                        if (j != k)
                            sum += yijk[i][j][k];
                    }
                    model.add(xij[i][j] == sum);
                    sum.end();
                }
            }
        }

        //constraint (1g)
        /*{
            for (int i = 0; i != M; i++)
            {
                for (int j = 0; j != N + 1; j++)
                {
                    IloExpr sum1(env);
                    for (int k = 0; k != N + 1; k++)
                    {
                        if (j != k)
                            sum1 += yijk[i][j][k];
                    }
                    IloExpr sum2(env);
                    for (int h = 0; h != N + 1; h++)
                    {
                        if (j != h)
                            sum2 += yijk[i][h][j];
                    }
                    model.add(sum1 == sum2);
                    sum1.end();
                    sum2.end();
                }
            }
        }*/
        
        
        //constraint (1j) (1n)
        {
            for (int i = 0; i != M; i++)
            {
                for (int j = 0; j != N + 1; j++)
                {
                    for (int k = 0; k != N; k++)
                    {
                        if (j != k)
                        {
                            model.add(Cj[k] - Cj[j] + L * (1 - yijk[i][j][k]) >= ST_ijk[i][j][k] + PT_ij[i][k]);
                        }
                    }
                    if (j != N)
                    {
                        model.add(Cj[j] <= U_i[i] + (1 - xij[i][j]) * L);
                    }
                }
            }
        }

        
        //constraint (1h)
        /*{
            for (int i = 0; i != M; i++)
            {
                IloExpr sum(env);
                for (int j = 0; j != N; j++)
                {
                    sum += yijk[i][N][j];
                }
                model.add(sum <= 1);
                sum.end();
            }
        }*/
        
        //constraint (1l)
        {
            model.add(Cj[N] == 0);
        }
        
        
        //constraint (1i) (1o)
        {
            for (int i = 0; i != M; i++)
            {
                IloExpr sum1(env);
                IloExpr sum2(env);
                for (int j = 0; j != N; j++)
                {
                    sum1 += xij[i][j] * PT_ij[i][j];
                }
                for (int j = 0; j != N + 1; j++)
                {
                    for (int k = 0; k != N + 1; k++)
                    {
                        if (j != k)
                        {
                            sum2 += yijk[i][j][k] * ST_ijk[i][j][k];
                        }
                    }
                }
                model.add(sum1 + sum2 <= U_i[i]);
                for (int j = 0; j != N; ++j)
                {
                    model.add(Cj[j] <= sum1 + sum2 + (1- yijk[i][j][N])*L); 
                }
                sum1.end();
                sum2.end();
            }
        }

        // / ///////////////////OPTIMIZE///////////////////////////////

        IloCplex cplex(model);
        cplex.setParam(IloCplex::TiLim, MIP_timeLimit);
        IloNum c1 = cplex.getCplexTime();
        cplex.solve();
        IloNum c2 = cplex.getCplexTime();
        MIP_time = c2 - c1;
        // / ////////////////OUTPUT AND PRINT//////////////////////////
        if (cplex.getStatus() == IloAlgorithm::Infeasible)
        {
            env.out() << "Infeasible! No Solution" << endl;
            MIP_obj = 0.0;
            MIP_bestBound = 0.0;
            MIP_finalGap = 0.0;
        }
        else
        {
            env.out() << "Solution status: " << cplex.getStatus() << endl;
            MIP_obj = cplex.getObjValue();
            MIP_bestBound = cplex.getBestObjValue();
            MIP_finalGap = cplex.getMIPRelativeGap();

            cout << "obj  = " << cplex.getObjValue() << endl;
            cout << "Cmax = " << cplex.getValue(Cmax) << endl;
            cout << "Zj :" << endl;
            for (int j = 0 ; j != N; j++)
            {
                cout << '\t' << "Zj[" << j << "] = " << cplex.getValue(zj[j]);   
            }
            cout << endl;
            cout << '\n' << "Cj :" << endl;
            for (int j = 0 ; j != N; j++)
            {
                cout << '\t' << "Cj[" << j << "] = " << cplex.getValue(Cj[j]);
            }
            cout << endl;
            cout << '\n' << "Xij :" << endl;
            for (int j = 0; j != N; ++j)
            {
                for (int i = 0; i != M; ++i)
                {
                    cout << '\t' << "X[" << i << "][" << j << "] = " << cplex.getValue(xij[i][j]);
                }
                cout << endl;
            }
            cout<<endl;
            cout << '\n' << "Yijk :" << endl;
            for (int i = 0; i != M; ++i)
            {
                for (int j = 0; j != N+1; ++j)
                {
                    for (int k = 0; k != N+1; ++k)
                    {
                        if ( j != k && cplex.getValue(yijk[i][j][k]) != 0)
                        {
                            cout << '\t' << "Y[" << i << "][" << j << "][" << k << "] = " << cplex.getValue(yijk[i][j][k]) << " setup = " << ST_ijk[i][j][k];
                        }
                    }
                }
                cout << endl;
            }
            cout<<endl;
        }
    }
    catch (IloException& ex) {
        cerr << "Error: " << ex << endl;
    }
    catch (...) {
        cerr << "Error" << endl;
    }
    env.end();
    cout << "-----------------MIP MODEL ENDS------------------" << endl;
}

void DM(
    const int DM_timeLimit,
    double &DM_obj,
    double &DM_time,
    bool &DM_find_feasible
)
{
    cout << '\n' << "*****************DECOMPOSITION METHOD STARTS*******************" << endl;
    time_t DM_begin, DM_end;
    time(&DM_begin);
    int numIter = 0;//iteration number

    DM_obj = -999999.0;

    bool DM_stop_condition = true;

    vector< vector<double> > Zj_Value_h;
    vector< vector< vector<double> > > Xij_Value_h;
    vector<double> CMAX_Value_h;
    vector<double> MP_obj_h;
    vector<double> SP_Cmax_h;


    do
    {
        // / /////////////////MASTER PROBLEM////////////////////////////////////////
        vector<double> Zj_Value;//value for zj
        vector< vector<double> > Xij_Value;//value for xij
        double CMAX_Value;//value for Cmax
        double MP_obj;//obj for MP
        bool is_MP_feasible = true;

        bool MP_reach_time_limit = false;
        bool MP_solution_unknown = false;

        MP_DM(numIter, DM_timeLimit, DM_begin, Zj_Value_h, CMAX_Value_h, MP_obj_h, SP_Cmax_h, Zj_Value, Xij_Value, CMAX_Value, MP_obj, is_MP_feasible, MP_reach_time_limit, MP_solution_unknown);

        if (MP_solution_unknown == true || MP_reach_time_limit == true)
        {
            DM_stop_condition = false;
            
        }
        else
        {
            Zj_Value_h.push_back(Zj_Value);
            Xij_Value_h.push_back(Xij_Value);
            CMAX_Value_h.push_back(CMAX_Value);
            MP_obj_h.push_back(MP_obj);
        }

        // / ///////////PRINT THE MP SOLUTION/////////////////////////////////////
        /*
        cout << "MP obj  = " << MP_obj << endl;
        cout << "MP Cmax = " << CMAX_Value << endl;
        cout << "MP Zj :" << endl;
        for (int j = 0 ; j != N; j++)
        {
            cout << '\t' << "Zj[" << j << "] = " << Zj_Value[j];   
        }
        cout << endl;
        cout << '\n' << "MP Xij :" << endl;
        for (int j = 0; j != N+1; ++j)
        {
            for (int i = 0; i != M; ++i)
            {
                cout << '\t' << "X[" << i << "][" << j << "] = " << Xij_Value[i][j];
            }
            cout << endl;
        }
        cout<<endl;
        */

        // / ////////////////SUBPROBLEM SOLVING///////////////////////////////////
        bool SP_any_feasible_solution;
        double SP_feasible_Cmax = 0.0;//store feasible Cmax derived by SP, infeasible = 0


        if (DM_stop_condition != false)
        {         

            if (is_MP_feasible)
            {

                bool SP_reach_timelimit = false;
                
                SP_any_feasible_solution = SP_DM(Zj_Value, Xij_Value, DM_timeLimit, DM_begin, SP_feasible_Cmax, SP_reach_timelimit);

                SP_Cmax_h.push_back(SP_feasible_Cmax);

                // /////////////////////////UPDATE BEST FOUND SOLUTION///////////////////
                if (SP_feasible_Cmax != 0.0)
                {
                    double current_obj_value = (-1.0)*SP_feasible_Cmax;
                    for (int j = 0; j != N; ++j)
                    {
                        if (Zj_Value[j] > 0.9)//zj = 1
                        {
                            current_obj_value += double(R_j[j]);
                        }
                    }


                    //////////UPDATE DM STATUS/////////////////////
                    DM_find_feasible = true;

                    
                    if (current_obj_value > DM_obj)
                    {
                        DM_obj = current_obj_value;
                    }
                    
                }

                if (SP_reach_timelimit == true)
                {
                    DM_stop_condition = false;
                }


            }
        }

        // // ///////////////////OPTIMAL CONDITION///////////////////////////////

        if (is_MP_feasible == false || SP_feasible_Cmax == CMAX_Value)
        {
            DM_stop_condition = false;
        }

        ++numIter;

        //DM_stop_condition = false; /*Delete afterwards*/

        // / ///////////////CLEAR VECTORS////////////////////////////////////////
        Zj_Value.clear();
        vector<double> ().swap(Zj_Value);
        Xij_Value.clear();
        vector< vector<double> > ().swap(Xij_Value);
    
    }while(DM_stop_condition);

    



    // /////////////////CLEAR VECTORS//////////////////////////////////////////
    Zj_Value_h.clear();
    vector< vector<double> > ().swap(Zj_Value_h);
    Xij_Value_h.clear();
    vector< vector< vector<double> > > ().swap(Xij_Value_h);
    CMAX_Value_h.clear();
    vector<double> ().swap(CMAX_Value_h);
    MP_obj_h.clear();
    vector<double> ().swap(MP_obj_h);
    SP_Cmax_h.clear();
    vector<double> ().swap(SP_Cmax_h);

    time(&DM_end);
    DM_time = difftime(DM_end,DM_begin);

    cout << '\n' << "*****************DECOMPOSITION METHOD ENDS*********************" << endl;
}

void MP_DM(
    const int numIter,
    const int DM_timeLimit,
    const time_t &startTime,
    const vector< vector<double> > &Zj_Value_h,
    const vector<double> &CMAX_Value_h,
    const vector<double> &MP_obj_h,
    const vector<double> &SP_Cmax_h,
    vector<double> &Zj_Value,
    vector< vector<double> > &Xij_Value,
    double &CMAX_Value,
    double &MP_obj,
    bool &is_MP_feasible,
    bool &MP_reach_time_limit,
    bool &MP_solution_unknown
)
{
    cout << '\n' << "--------------------MASTER PROBLEM STARTS----------------------" << endl;
    time_t MP_begin;
    time(&MP_begin);
    double MP_time = double(DM_timeLimit) - difftime(MP_begin, startTime);
    cout << " Residual time for MP = " << MP_time << endl;
    if (MP_time < 0.0)
    {
        MP_reach_time_limit = true;
    }
    else
    {
        IloEnv MP_env;   
        try 
        {
            // / ////////////CREATE MODEL////////////////////////////////////////

            IloModel MP_model(MP_env);
            // /////////DECISION VARIABLES/////////////////////////////////////////

            NumVarMatrix xij(MP_env, M);
            for (int i = 0; i != M; i++)
            {
                xij[i] = IloNumVarArray(MP_env, N + 1, 0, 1, ILOINT);
            }

            IloNumVarArray zj(MP_env, N, 0, 1, ILOINT);

            NumVar3Matrix yijk(MP_env, M);
            for (int i = 0; i != M; i++)
            {
                yijk[i] = NumVarMatrix(MP_env, N + 1);
                for (int j = 0; j != N + 1; j++)
                {
                    yijk[i][j] = IloNumVarArray(MP_env, N + 1, 0, 1, ILOFLOAT);
                }
            }

            IloNumVar Cmax(MP_env);

            // / /////////////CONSTRAINTS///////////////////////////

            // Objective Function (1a)
            {
                IloExpr obj(MP_env);
                for (int j = 0; j != N; j++)
                {
                    obj += R_j[j] * zj[j];
                }
                obj -= Cmax;
                MP_model.add(IloMaximize(MP_env, obj));
                obj.end();
            }

            //constraint (1b)
            {
                for (int j = 0; j != N; j++)
                {
                    IloExpr sum(MP_env);
                    for (int i = 0; i != M; i++)
                    {
                        sum += xij[i][j];
                    }
                    MP_model.add(sum == zj[j]);
                    sum.end();
                }
            }

            //constraint (1c)
            /*{
                for (int j = 0; j != N; j++)
                {
                    for (int i = 0; i != M; i++)
                    {
                        MP_model.add(xij[i][j] <= zj[j]);
                    }
                }
            }*/

            //constraint (1d)
            {
                for (int i = 0; i != M; i++)
                {
                    MP_model.add(xij[i][N] == 1);
                }
            }

            //constraint (1e)
            {
                for (int i = 0; i != M; i++)
                {
                    for (int k = 0; k != N + 1; k++)
                    {
                        IloExpr sum(MP_env);
                        for (int j = 0; j != N + 1; j++)
                        {
                            if (j != k)
                                sum += yijk[i][j][k];
                        }
                        MP_model.add(xij[i][k] == sum);
                        sum.end();
                    }
                }
            }

            //constraint (1f)
            {
                for (int i = 0; i != M; i++)
                {
                    for (int j = 0; j != N + 1; j++)
                    {
                        IloExpr sum(MP_env);
                        for (int k = 0; k != N + 1; k++)
                        {
                            if (j != k)
                                sum += yijk[i][j][k];
                        }
                        MP_model.add(xij[i][j] == sum);
                        sum.end();
                    }
                }
            }

            //constraint (1g)
            /*{
                for (int i = 0; i != M; i++)
                {
                    for (int j = 0; j != N + 1; j++)
                    {
                        IloExpr sum1(MP_env);
                        for (int k = 0; k != N + 1; k++)
                        {
                            if (j != k)
                                sum1 += yijk[i][j][k];
                        }
                        IloExpr sum2(MP_env);
                        for (int h = 0; h != N + 1; h++)
                        {
                            if (j != h)
                                sum2 += yijk[i][h][j];
                        }
                        MP_model.add(sum1 == sum2);
                        sum1.end();
                        sum2.end();
                    }
                }
            }*/

            //constraint (1h)
            /*{
                for (int i = 0; i != M; i++)
                {
                    IloExpr sum(MP_env);
                    for (int j = 0; j != N; j++)
                    {
                        sum += yijk[i][N][j];
                    }
                    MP_model.add(sum <= 1);
                    sum.end();
                }
            }*/

            //constraint (1i) (2c) 
            {
                for (int i = 0; i != M; i++)
                {
                    IloExpr sum1(MP_env);
                    IloExpr sum2(MP_env);
                    for (int j = 0; j != N; j++)
                    {
                        sum1 += xij[i][j] * PT_ij[i][j];
                    }
                    for (int j = 0; j != N + 1; j++)
                    {
                        for (int k = 0; k != N + 1; k++)
                        {
                            if (j != k)
                            {
                                sum2 += yijk[i][j][k] * ST_ijk[i][j][k];
                            }
                        }
                    }
                    MP_model.add(sum1 + sum2 <= U_i[i]);
                    MP_model.add(Cmax >= sum1 + sum2);
                    
                    sum1.end();
                    sum2.end();
                }
            }


            // //////////////////////CUTS AND BOUNDS//////////////////////////////
            
            if (numIter != 0)
            {
                for (int iter = 0; iter < numIter; ++iter)
                {

                    //cut 1
                    {

                        IloExpr sum1(MP_env);
                        for (int j = 0; j != N; ++j)
                        {
                            if (Zj_Value_h[iter][j] > 0.9)//zhj = 1
                            {
                                sum1 += 1 - zj[j];
                            }
                            else
                            {
                                sum1 += zj[j];
                            }
                        } 
                        MP_model.add(sum1 >= 1);

                        sum1.end();
                    }

                    //check if this iteration have feasible solution
                    if (SP_Cmax_h[iter] > 0.9)//there is a feasible solution
                    {
                        

                        // LB derived from feasible solution
                        //This cut could be more agressive with RHS + 1

                        double derived_feasible_obj = (-1.0)*SP_Cmax_h[iter];
                        for (int j = 0; j != N; j++)
                        {
                            if (Zj_Value_h[iter][j] > 0.9)//zj = 1
                            {
                                derived_feasible_obj += R_j[j];
                            }
                        }

                        cout << "derived feasible obj = " << derived_feasible_obj << endl;

                        {
                            IloExpr obj(MP_env);                      
                            for (int j = 0; j != N; j++)
                            {
                                obj += R_j[j] * zj[j];
                            }
                            obj -= Cmax;
                            MP_model.add(obj >= derived_feasible_obj);
                            obj.end();

                        }

                    }
                    else
                    {
                        //cut 2
                        {
                            IloExpr sum1(MP_env);
                            int numAccept = 0;
                            for (int j = 0; j != N; ++j)
                            {
                                if (Zj_Value_h[iter][j] > 0.9)//zhj = 1
                                {
                                    sum1 += zj[j];
                                    ++numAccept;

                                }
                                                        
                            } 
                            MP_model.add(sum1 <= numAccept - 1);

                        sum1.end();
                        }

                    }



                }
            }
            

            // / ///////////OPTIMIZE & VAlUE OUTPUT///////////////////////////
            IloCplex cplex(MP_model);
            cplex.setParam(IloCplex::TiLim, MP_time); 
            cplex.solve();
            if (cplex.getStatus() == IloAlgorithm::Infeasible)
            {
                MP_env.out() << "Infeasible! No Solution --- MP_DM  " << endl;
                is_MP_feasible = false;
            }
            else
            {
                MP_env.out() << "Solution status: " << cplex.getStatus() << endl;

                if (cplex.getStatus() == 0)//unknown
                {
                    MP_solution_unknown = true;
                }
                else
                {
                    cout << "obj  = " << cplex.getObjValue() << endl;
                    MP_obj = cplex.getObjValue();

                    cout << "Cmax = " << cplex.getValue(Cmax) << endl;
                    CMAX_Value = cplex.getValue(Cmax);

                    
                    cout << "Zj :" << endl;
                    for (int j = 0 ; j != N; j++)
                    {
                        cout << '\t' << "Zj[" << j << "] = " << cplex.getValue(zj[j]);   
                    }
                    cout << endl;
                    
                    Zj_Value = get_var_values(cplex,zj,N);

                    
                    cout << '\n' << "Xij :" << endl;
                    for (int j = 0; j != N; ++j)
                    {
                        for (int i = 0; i != M; ++i)
                        {
                            cout << '\t' << "X[" << i << "][" << j << "] = " << cplex.getValue(xij[i][j]);
                        }
                        cout << endl;
                    }
                    cout<<endl;
                    
                    Xij_Value = get_var_values(cplex,xij,M,N+1);

                    
                    cout << '\n' << "Yijk :" << endl;
                    for (int i = 0; i != M; ++i)
                    {
                        for (int j = 0; j != N+1; ++j)
                        {
                            for (int k = 0; k != N+1; ++k)
                            {
                                if ( j != k && cplex.getValue(yijk[i][j][k]) != 0)
                                {
                                    cout << '\t' << "Y[" << i << "][" << j << "][" << k << "] = " << cplex.getValue(yijk[i][j][k]) << " setup = " << ST_ijk[i][j][k];
                                }
                            }
                        }
                        cout << endl;
                    }
                    cout<<endl;
                }
                
            

            }

        }
        catch (IloException& ex) {
            cerr << "Error: " << ex << endl;
        }
        catch (...) {
            cerr << "Error" << endl;
        }
        MP_env.end();
    }

    cout << '\n' << "---------------------MASTER PROBLEM ENDS-----------------------" << endl;
}

bool SP_DM(
    const vector<double> &Zj_Value,
    const vector< vector<double> > &Xij_Value,
    const int DM_timeLimit,
    const time_t &startTime,
    double &SP_feasible_Cmax,
    bool &SP_reach_timelimit
)
{
    cout << '\n' << "-----------------SUBPROBLEM SOLVING STARTS---------------------" << endl;
    bool is_global_feasible = false;

    vector <int> acceptJobs;//accepted jobs 
    vector <int> initialAssignment;//MP assignment to machine for each job
    for (int j = 0; j != N; ++j)
    {
        if (Zj_Value[j] > 0.99)//zj == 1.0
        {
            acceptJobs.push_back(j);
            for (int i = 0; i != M; ++i)
            {
                if (Xij_Value[i][j] > 0.99)//xij == 1.0
                {
                    initialAssignment.push_back(i);
                }
            }
        }
    }

    cout << "Accepted jobs: " << endl;
    for (int iter = 0; iter != acceptJobs.size(); ++iter)
    {
        cout << " " << acceptJobs[iter] << " ";
    }
    cout << '\n' << endl;

    // / //////////////////////////////////////////////////////////////////////////////////

    int SP_IterNum = 0;//iteration number for subproblem
    //vector<double> SPM_Obj_h;//for each iteration h, store the SPM obj
    vector<double> SPM_Cmax_h;//for each iteration h, store the SPM cmax
    vector< vector< vector< double> > > SPM_Xij_h;//for each iteration h, store the SPM xij

    bool stop_condition = true;

    vector<double> SPS_Cmax_h;//for each iteration h, store the SPS_Cmax
    vector< vector< double> > SPS_Oi_h;// for each iteration h, store the Oi for each machine
    vector<bool> feasible_h;//for each iteration h, whether such a solution is feasible
    double best_feasible_Cmax = 0.0;//store Cmax of the best found feasible solution 

    do
    {
        // ////////////////SP MASTER PROBLEM//////////////////////////////////
        bool is_SPM_feasible = true;

        bool SPM_reach_timelimit = false;
        bool SPM_solution_unknown = false;

        SPM(SP_IterNum,DM_timeLimit,startTime,acceptJobs,initialAssignment,SPS_Cmax_h,SPS_Oi_h,feasible_h,best_feasible_Cmax,is_SPM_feasible, SPM_Cmax_h,SPM_Xij_h, SPM_reach_timelimit, SPM_solution_unknown);


        if (SPM_reach_timelimit == true || SPM_solution_unknown == true)
        {
            stop_condition = false;
        }
        else
        {

            // /////////////////PRINT INFORMATION//////////////////////////////////
            cout << '\n' << "SPM solution :" << endl;
            //cout << "obj =  " << SPM_Obj_h[SP_IterNum] << endl;
            cout << "Cmax = " << SPM_Cmax_h[SP_IterNum] << endl;
            cout << "xij : " << endl;
            for (int i = 0; i != M; ++i)
            {
                for (int j = 0; j != acceptJobs.size(); ++j)
                {
                    cout << '\t' << "X[" << i << "][" << acceptJobs[j] << "] = " << SPM_Xij_h[SP_IterNum][i][j];
                }
                cout << endl;
            } 
            cout << endl;
            // ///////////////////////////////////////////////////////////////////////

        }



        if (stop_condition != false)
        {


            double SPS_Cmax = 0.0;
            if (is_SPM_feasible)
            {
                // / //////////////////SP SUBPROBLEM///////////////////////////////

                bool is_SPS_feasible = true;

                bool SPS_reach_timelimit = false;
                bool SPS_solution_unknown = false;
                
                SPS(SP_IterNum,DM_timeLimit,startTime,acceptJobs,SPM_Xij_h,is_SPS_feasible,SPS_Cmax,SPS_Oi_h, SPS_reach_timelimit, SPS_solution_unknown);

                if (SPS_reach_timelimit == true || SPS_solution_unknown == true)
                {
                    stop_condition = false;
                }
                else
                {

                    if (is_SPS_feasible)
                    {
                        feasible_h.push_back(true);
                        is_global_feasible = true;

                        // / /////////UPDATE BEST FEASIBLE SOLUTION CMAX///////////// 

                        if (SP_IterNum == 0)
                        {
                            best_feasible_Cmax = SPS_Cmax;                   
                        }
                        else
                        {
                            if (best_feasible_Cmax > SPS_Cmax)
                            {
                                best_feasible_Cmax = SPS_Cmax;
                            }
                        }
                        // / ///////////CHECK OPTIMAL SOLUTION///////////////////////

                        if (SPS_Cmax == SPM_Cmax_h[SP_IterNum])
                        {
                            stop_condition = false;
                        }


                    }
                    else
                    {
                        feasible_h.push_back(false);
                    }
                }

                //stop_condition = false;/*Delete afterwards*/

            }
            else //SPM is infeasible
            {
                //no need to solve the SPS
                stop_condition = false;
                SPS_Oi_h.push_back(vector<double> ());
                feasible_h.push_back(false);

            }
            SPS_Cmax_h.push_back(SPS_Cmax);

            ++SP_IterNum;
            cout << " ----------------------------------------------------" << '\n' << endl;

        }

    }while(stop_condition);


    // ///////////////OUTPUT REQUIRED FEASIBLE SOLUTION////////////////////

    SP_feasible_Cmax = best_feasible_Cmax;




    // / ////////////CLEAR VECTORS & POINTERS/////////////////////////////
    acceptJobs.clear();
    vector <int> ().swap(acceptJobs);
    initialAssignment.clear();
    vector <int> ().swap(initialAssignment);
    SPM_Cmax_h.clear();
    vector<double>().swap(SPM_Cmax_h);
    SPM_Xij_h.clear();
    vector< vector< vector< double> > > ().swap (SPM_Xij_h);
    SPS_Cmax_h.clear();
    vector<double>().swap(SPS_Cmax_h);
    SPS_Oi_h.clear();
    vector< vector< double> > ().swap(SPS_Oi_h);
    feasible_h.clear();
    vector<bool> ().swap(feasible_h);



    // ///////////////CHECK PROCEDURE//////////////////////////////////////

    //checkSP(Zj_Value);

    time_t SP1_begin;
    time(&SP1_begin);
    double SP1_time = double(DM_timeLimit) - difftime(SP1_begin, startTime);
   

    if (SP1_time < 0.0)
    {
        SP_reach_timelimit = true;
    }


    cout << '\n' << "-----------------SUBPROBLEM SOLVING ENDS-----------------------" << endl;
    return is_global_feasible;
}

void SPM(
    const int SP_IterNum,
    const int DM_timeLimit,
    const time_t &startTime,
    const vector<int> &acceptJobs,
    const vector<int> &initialAssignment,
    const vector<double> &SPS_Cmax_h,
    const vector< vector< double> > &SPS_Oi_h,
    const vector<bool> &feasible_h,
    const double best_feasible_Cmax,
    bool &is_SPM_feasible,
    vector<double>  &SPM_Cmax_h,
    vector< vector< vector< double> > > &SPM_Xij_h,
    bool &SPM_reach_timelimit,
    bool &SPM_solution_unknown
)
{
    time_t SPM_begin;
    time(&SPM_begin);
    double SPM_time = double(DM_timeLimit) - difftime(SPM_begin, startTime);
    cout << " Residual time for SPM = " << SPM_time << endl;

    if (SPM_time < 0.0)
    {
        SPM_reach_timelimit = true;
    }
    else
    {

        IloEnv SPM_env;   
        try 
        {
            // / ////////////CREATE MODEL////////////////////////////////////////

            IloModel SPM_model(SPM_env);
            // /////////DECISION VARIABLES/////////////////////////////////////////

            NumVarMatrix xij(SPM_env, M);
            for (int i = 0; i != M; i++)
            {
                xij[i] = IloNumVarArray(SPM_env, acceptJobs.size() + 1, 0, 1, ILOINT);
            }

            NumVar3Matrix yijk(SPM_env, M);
            for (int i = 0; i != M; i++)
            {
                yijk[i] = NumVarMatrix(SPM_env, acceptJobs.size() + 1);
                for (int j = 0; j != acceptJobs.size() + 1; j++)
                {
                    yijk[i][j] = IloNumVarArray(SPM_env, acceptJobs.size() + 1, 0, 1, ILOFLOAT);
                }
            }

            IloNumVar Cmax(SPM_env);

            // / /////////////CONSTRAINTS///////////////////////////

            // Objective Function (3a)
            {    
                SPM_model.add(IloMinimize(SPM_env, Cmax));         
            }
            

            //constraint (3b)
            {
                for (int j = 0; j != acceptJobs.size(); j++)
                {
                    IloExpr sum(SPM_env);
                    for (int i = 0; i != M; i++)
                    {
                        sum += xij[i][j];
                    }
                    SPM_model.add(sum == 1);
                    sum.end();
                }
            }
            

            //constraint (3c)
            {
                for (int i = 0; i != M; i++)
                {
                    SPM_model.add(xij[i][acceptJobs.size()] == 1);
                }
            }
        

            //constraint (1e)
            {
                for (int i = 0; i != M; i++)
                {
                    for (int k = 0; k != acceptJobs.size() + 1; k++)
                    {
                        IloExpr sum(SPM_env);
                        for (int j = 0; j != acceptJobs.size() + 1; j++)
                        {
                            if (j != k)
                                sum += yijk[i][j][k];
                        }
                        SPM_model.add(xij[i][k] == sum);
                        sum.end();
                    }
                }
            }
        

            //constraint (1f)
            {
                for (int i = 0; i != M; i++)
                {
                    for (int j = 0; j != acceptJobs.size() + 1; j++)
                    {
                        IloExpr sum(SPM_env);
                        for (int k = 0; k != acceptJobs.size() + 1; k++)
                        {
                            if (j != k)
                                sum += yijk[i][j][k];
                        }
                        SPM_model.add(xij[i][j] == sum);
                        sum.end();
                    }
                }
            }
            

            //constraint (1g)
            /*{
                for (int i = 0; i != M; i++)
                {
                    for (int j = 0; j != acceptJobs.size() + 1; j++)
                    {
                        IloExpr sum1(SPM_env);
                        for (int k = 0; k != acceptJobs.size() + 1; k++)
                        {
                            if (j != k)
                                sum1 += yijk[i][j][k];
                        }
                        IloExpr sum2(SPM_env);
                        for (int h = 0; h != acceptJobs.size() + 1; h++)
                        {
                            if (j != h)
                                sum2 += yijk[i][h][j];
                        }
                        SPM_model.add(sum1 == sum2);
                        sum1.end();
                        sum2.end();
                    }
                }
            }*/
        

            //constraint (1h)
            /*{
                for (int i = 0; i != M; i++)
                {
                    IloExpr sum(SPM_env);
                    for (int j = 0; j != acceptJobs.size(); j++)
                    {
                        sum += yijk[i][acceptJobs.size()][j];
                    }
                    SPM_model.add(sum <= 1);
                    sum.end();
                }
            }*/
            

            //constraint (1i) (2c) (2d)
            {
                for (int i = 0; i != M; i++)
                {
                    IloExpr sum1(SPM_env);
                    IloExpr sum2(SPM_env);
                    for (int j = 0; j != acceptJobs.size(); j++)
                    {
                        sum1 += xij[i][j] * PT_ij[i][acceptJobs[j]];
                    }
                    for (int j = 0; j != acceptJobs.size(); j++)
                    {
                        for (int k = 0; k != acceptJobs.size(); k++)
                        {
                            if (j != k)
                            {
                                sum2 += yijk[i][j][k] * ST_ijk[i][acceptJobs[j]][acceptJobs[k]];
                            }
                        }
                        sum2 += yijk[i][j][acceptJobs.size()] * ST_ijk[i][acceptJobs[j]][N];
                        sum2 += yijk[i][acceptJobs.size()][j] * ST_ijk[i][N][acceptJobs[j]];
                    }
                    SPM_model.add(sum1 + sum2 <= U_i[i]);
                    SPM_model.add(Cmax >= sum1 + sum2);
                    
                    sum1.end();
                    sum2.end();
                }
            }


            // / ////////////////CUT/////////////////////////////////////////////

            if (SP_IterNum > 0)
            {
                for (int iter = 0; iter < SP_IterNum; ++iter)
                {
                    
                    if (feasible_h[iter] == false)
                    {
                        //cut 1
                        for (int i = 0; i != M; ++i)
                        {
                            if (SPS_Oi_h[iter][i] > double(U_i[i]))
                            {
                                IloExpr sum(SPM_env);
                                int temp_assign = 0;
                                for (int j = 0; j != acceptJobs.size(); ++j)
                                {
                                    if (SPM_Xij_h[iter][i][j] > 0.99)//x^h_ij == 1
                                    {
                                        sum += xij[i][j];
                                        ++temp_assign;
                                    }
                                }
                                SPM_model.add(sum <= temp_assign - 1);
                                sum.end();
                            }
                        }
                    }
                    else //it is a feasible solution
                    {
                        //cut2
                        IloExpr sum(SPM_env);
                        for (int i = 0; i != M; ++i)
                        {
                            for (int j = 0; j != acceptJobs.size(); ++j)
                            {
                                if (SPM_Xij_h[iter][i][j] > 0.99)//x^h_ij == 1
                                {
                                    sum += 1 - xij[i][j];
                                }
                                else//x^h_ij == 0
                                {
                                    sum += xij[i][j];
                                }
                            }
                        }
                        SPM_model.add(sum >= 1);
                        sum.end();

                        //upper bound for Cmax
                        /* this RHS could be more agressive with -1 unit */
                        if (best_feasible_Cmax > 0.9)
                        {
                            SPM_model.add(Cmax <= best_feasible_Cmax);
                        }

                    }


                    //cut 3
                    {
                        for (int i = 0; i != M; ++i)
                        {
                            IloExpr sum(SPM_env);
                            for (int j = 0; j != acceptJobs.size(); ++j)
                            {
                                if (SPM_Xij_h[iter][i][j] > 0.9)
                                {
                                    int theta_ij = PT_ij[i][acceptJobs[j]];
                                    int maxST_ij = ST_ijk[i][N][acceptJobs[j]];

                                    for (int k = 0; k != acceptJobs.size(); ++k)
                                    {
                                        if (k != j && SPM_Xij_h[iter][i][k] > 0.9)
                                        {
                                            if (ST_ijk[i][acceptJobs[k]][acceptJobs[j]] > maxST_ij)
                                            {
                                                maxST_ij = ST_ijk[i][acceptJobs[k]][acceptJobs[j]];
                                            }
                                        }
                                    }

                                    sum += (1 - xij[i][j])*(theta_ij + maxST_ij);
                                }
                            }

                            SPM_model.add(Cmax >= SPS_Oi_h[iter][i] - sum);

                            sum.end();
                        }
                    }


                    //lower bound for Cmax
                    {
                        SPM_model.add(Cmax >= SPM_Cmax_h[iter]);
                    }


                }
                
            }

            // / ///////////OPTIMIZE & VAlUE OUTPUT///////////////////////////
            IloCplex cplex(SPM_model);
        
            //assign a MIP start for xij
            if (SP_IterNum == 0)
            {
                for (int i = 0; i != M; ++i)
                {
                    IloNumArray temp(SPM_env);
                    for (int j = 0; j != acceptJobs.size(); ++j)
                    {
                        if (initialAssignment[j] == i)
                        {
                            temp.add(1);
                        }
                        else
                        {
                            temp.add(0);
                        }
                    }
                    temp.add(1);
                    cplex.addMIPStart(xij[i],temp);
                    temp.end();
                }
            }

            cplex.setParam(IloCplex::TiLim, SPM_time); 
            cplex.solve();

            if (cplex.getStatus() == IloAlgorithm::Infeasible)
            {
                SPM_env.out() << "Infeasible! No Solution --- SPM  " << endl;
                is_SPM_feasible = false;
                //SPM_Obj_h.push_back(0.0);
                SPM_Cmax_h.push_back(0.0);

                vector< vector< double> > temp2Index;
                for (int i = 0; i != M; ++i)
                {
                    vector<double> temp1Index;
                    for (int j = 0; j != acceptJobs.size(); ++j)
                    {
                        temp1Index.push_back(0.0);
                    }
                    temp2Index.push_back(temp1Index);
                    temp1Index.clear();
                    vector<double> ().swap(temp1Index);
                }
                SPM_Xij_h.push_back(temp2Index);
                temp2Index.clear();
                vector< vector< double> > ().swap(temp2Index);

            }
            else
            {
                SPM_env.out() << "Solution status: " << cplex.getStatus() << endl;

                if (cplex.getStatus() == 0)
                {
                    SPM_solution_unknown = true;
                }
                else
                {
                    is_SPM_feasible = true;

                    cout << "obj  = " << cplex.getObjValue() << endl;
                    //SPM_Obj_h.push_back(cplex.getObjValue());

                    cout << "Cmax = " << cplex.getValue(Cmax) << endl;
                    SPM_Cmax_h.push_back(cplex.getValue(Cmax));
                    
                    cout << '\n' << "Xij :" << endl;
                    for (int j = 0; j != acceptJobs.size(); ++j)
                    {
                        for (int i = 0; i != M; ++i)
                        {
                            cout << '\t' << "X[" << i << "][" << acceptJobs[j] << "] = " << cplex.getValue(xij[i][j]) << "  pij = " << PT_ij[i][acceptJobs[j]];
                        }
                        cout << endl;
                    }
                    cout<<endl;

                    vector< vector< double> > temp2Index;
                    for (int i = 0; i != M; ++i)
                    {
                        vector<double> temp1Index;
                        for (int j = 0; j != acceptJobs.size(); ++j)
                        {
                            //temp1Index.push_back(cplex.getValue(xij[i][j]));
                            if (cplex.getValue(xij[i][j]) > 0.9)
                            {
                                temp1Index.push_back(1.0);
                            }
                            else
                            {
                                temp1Index.push_back(0.0);
                            }

                        }
                        temp2Index.push_back(temp1Index);
                        temp1Index.clear();
                        vector<double> ().swap(temp1Index);
                    }
                    SPM_Xij_h.push_back(temp2Index);
                    temp2Index.clear();
                    vector< vector< double> > ().swap(temp2Index);
                    
                    cout << '\n' << "Yijk :" << endl;
                    for (int i = 0; i != M; ++i)
                    {
                        for (int j = 0; j != acceptJobs.size(); ++j)
                        {
                            for (int k = 0; k != acceptJobs.size(); ++k)
                            {
                                if ( j != k && cplex.getValue(yijk[i][j][k]) != 0)
                                {
                                    cout << '\t' << "Y[" << i << "][" << acceptJobs[j] << "][" << acceptJobs[k] << "] = " << cplex.getValue(yijk[i][j][k]) << " setup = " << ST_ijk[i][acceptJobs[j]][acceptJobs[k]];
                                }
                            }

                            if (cplex.getValue(yijk[i][acceptJobs.size()][j]) != 0)
                            {
                                cout << '\t' << "Y[" << i << "][" << N << "][" << acceptJobs[j] << "] = " << cplex.getValue(yijk[i][acceptJobs.size()][j]) << " setup = " << ST_ijk[i][N][acceptJobs[j]];
                            }
                            if (cplex.getValue(yijk[i][j][acceptJobs.size()]) != 0)
                            {
                                cout << '\t' << "Y[" << i << "][" << acceptJobs[j] << "][" << N << "] = " << cplex.getValue(yijk[i][j][acceptJobs.size()]) << " setup = " << ST_ijk[i][acceptJobs[j]][N];
                            }

                        }
                        cout << endl;
                    }
                    cout<<endl;

                }
            

            }

        }
        catch (IloException& ex) {
            cerr << "Error: " << ex << endl;
        }
        catch (...) {
            cerr << "Error" << endl;
        }
        SPM_env.end();

    }

}


void SPS(
    const int SP_IterNum,
    const int DM_timeLimit,
    const time_t &startTime,
    const vector<int> &acceptJobs,
    const vector< vector< vector< double> > > &SPM_Xij_h,
    bool &is_SPS_feasible,
    double &SPS_Cmax,
    vector< vector<double> > &SPS_Oi_h,
    bool &SPS_reach_timelimit,
    bool &SPS_solution_unknown
)
{
    time_t SPS_begin;
    time(&SPS_begin);
    double SPS_time = double(DM_timeLimit) - difftime(SPS_begin, startTime);
    cout << " Residual time for SPS = " << SPS_time << endl;

    if (SPS_time < 0.0)
    {
        SPS_reach_timelimit = true;
    }
    else
    {

        IloEnv SPS_env;
        try 
        {
            // / ////////////CREATE MODEL////////////////////////////////////////

            IloModel SPS_model(SPS_env);
            // /////////DECISION VARIABLES/////////////////////////////////////////

            NumVar3Matrix yijk(SPS_env, M);
            for (int i = 0; i != M; i++)
            {
                yijk[i] = NumVarMatrix(SPS_env, acceptJobs.size() + 1);
                for (int j = 0; j != acceptJobs.size() + 1; j++)
                {
                    yijk[i][j] = IloNumVarArray(SPS_env, acceptJobs.size() + 1, 0, 1, ILOINT);
                }
            }

            IloNumVarArray oi(SPS_env, M, 0, INFINITY, ILOINT);

            IloNumVarArray uj(SPS_env, acceptJobs.size(), 0, acceptJobs.size() - 1, ILOINT);

            // / /////////////CONSTRAINTS///////////////////////////

            // Objective Function  + Oi definition
            { 
                for (int i = 0; i != M; i++)
                {
                    IloExpr sum1(SPS_env);
                    IloExpr sum2(SPS_env);
                    for (int j = 0; j != acceptJobs.size(); j++)
                    {
                        sum1 += SPM_Xij_h[SP_IterNum][i][j] * PT_ij[i][acceptJobs[j]];
                    }
                    for (int j = 0; j != acceptJobs.size(); j++)
                    {
                        for (int k = 0; k != acceptJobs.size(); k++)
                        {
                            if (j != k)
                            {
                                sum2 += yijk[i][j][k] * ST_ijk[i][acceptJobs[j]][acceptJobs[k]];
                            }
                        }
                        sum2 += yijk[i][j][acceptJobs.size()] * ST_ijk[i][acceptJobs[j]][N];
                        sum2 += yijk[i][acceptJobs.size()][j] * ST_ijk[i][N][acceptJobs[j]];
                    }
                    SPS_model.add(sum1 + sum2 == oi[i]);
                    sum1.end();
                    sum2.end();
                } 
                IloExpr sum(SPS_env);  
                for (int i = 0; i != M; ++i)
                {
                    sum += oi[i];
                }
                SPS_model.add(IloMinimize(SPS_env, sum));   
                sum.end();     
            }

            //constraint: for k predecessor
            {
                for (int i = 0; i != M; i++)
                {
                    for (int k = 0; k != acceptJobs.size(); k++)
                    {
                        IloExpr sum(SPS_env);
                        for (int j = 0; j != acceptJobs.size() + 1; j++)
                        {
                            if (j != k)
                                sum += yijk[i][j][k];
                        }
                        SPS_model.add(SPM_Xij_h[SP_IterNum][i][k] == sum);
                        sum.end();
                    }
                }
            }
        

            //constraint: for j successor
            {
                for (int i = 0; i != M; i++)
                {
                    for (int j = 0; j != acceptJobs.size(); j++)
                    {
                        IloExpr sum(SPS_env);
                        for (int k = 0; k != acceptJobs.size() + 1; k++)
                        {
                            if (j != k)
                                sum += yijk[i][j][k];
                        }
                        SPS_model.add(SPM_Xij_h[SP_IterNum][i][j] == sum);
                        sum.end();
                    }
                }
            }
            

            //constraint: predecessor == successor
            /*{
                for (int i = 0; i != M; i++)
                {
                    for (int j = 0; j != acceptJobs.size() + 1; j++)
                    {
                        IloExpr sum1(SPS_env);
                        for (int k = 0; k != acceptJobs.size() + 1; k++)
                        {
                            if (j != k)
                                sum1 += yijk[i][j][k];
                        }
                        IloExpr sum2(SPS_env);
                        for (int h = 0; h != acceptJobs.size() + 1; h++)
                        {
                            if (j != h)
                                sum2 += yijk[i][h][j];
                        }
                        SPS_model.add(sum1 == sum2);
                        sum1.end();
                        sum2.end();
                    }
                }
            }*/
        

            //constraint: only one first job on machine
            /*{
                for (int i = 0; i != M; i++)
                {
                    IloExpr sum(SPS_env);
                    for (int j = 0; j != acceptJobs.size(); j++)
                    {
                        sum += yijk[i][acceptJobs.size()][j];
                    }
                    SPS_model.add(sum <= 1);
                    sum.end();
                }
            }*/
            
            //constraint: Ui limit & Cmax
            /*{
                for (int i = 0; i != M; ++i)
                {
                    SPS_model.add(oi[i] <= U_i[i]);
                }
            }*/

            //uj position cut 1
            {
                for (int j = 0; j != acceptJobs.size(); ++j)
                {
                    for (int k = 0; k != acceptJobs.size(); ++k)
                    {
                        if (j != k)
                        {
                            IloExpr sum(SPS_env);
                            for (int i = 0; i != M; ++i)
                            {
                                sum += yijk[i][j][k];
                            }
                            int job_size = acceptJobs.size();
                            SPS_model.add(uj[j] - uj[k] + job_size * sum <= job_size - 1);
                            sum.end();
                        }
                    }
                }
            }

            //uj position cut 2 & 3
            {
                for (int j = 0; j != acceptJobs.size(); ++j)
                {
                    IloExpr sum(SPS_env);
                    for (int i = 0; i != M; ++i)
                    {
                        sum += yijk[i][acceptJobs.size()][j];
                    }
                    int job_size = acceptJobs.size();
                    SPS_model.add(uj[j] + (job_size - 1) * sum <= job_size - 1);
                    SPS_model.add(uj[j] + sum >= 1);
                    sum.end();
                }
            }
            

            // / ///////////OPTIMIZE & VAlUE OUTPUT///////////////////////////
            IloCplex cplex(SPS_model);

            cplex.setParam(IloCplex::TiLim, SPS_time); 
            cplex.solve();

            if (cplex.getStatus() == IloAlgorithm::Infeasible)
            {
                SPS_env.out() << "Infeasible! No Solution --- SPS  " << endl;
                
            }
            else
            {
                SPS_env.out() << "Solution status: " << cplex.getStatus() << endl;

                if (cplex.getStatus() == 0)
                {
                    SPS_solution_unknown = true;
                }
                else
                {

                    cout << "obj  = " << cplex.getObjValue() << endl;
                    
                    cout << "Oi : " << endl;
                    SPS_Cmax = 0.0;
                    vector<double> temp_oi;
                    for (int i = 0; i != M; ++i)
                    {
                        cout << '\t' << "O[" << i << "] = " << cplex.getValue(oi[i]); 
                        temp_oi.push_back(cplex.getValue(oi[i]));
                        if (cplex.getValue(oi[i]) > SPS_Cmax)
                        {
                            SPS_Cmax = cplex.getValue(oi[i]);
                        }

                        //check if such a solution is feasible
                        if (cplex.getValue(oi[i]) > U_i[i])
                        {
                            is_SPS_feasible = false;
                        }
                    }
                    SPS_Oi_h.push_back(temp_oi);
                    temp_oi.clear();
                    vector<double> ().swap(temp_oi);
                    cout << endl;

                    cout << "SPS_Cmax = " << SPS_Cmax << endl;

                    cout << '\n' << "Yijk :" << endl;
                    for (int i = 0; i != M; ++i)
                    {
                        for (int j = 0; j != acceptJobs.size(); ++j)
                        {
                            for (int k = 0; k != acceptJobs.size(); ++k)
                            {
                                if ( j != k && cplex.getValue(yijk[i][j][k]) != 0)
                                {
                                    cout << '\t' << "Y[" << i << "][" << acceptJobs[j] << "][" << acceptJobs[k] << "] = " << cplex.getValue(yijk[i][j][k]) << " setup = " << ST_ijk[i][acceptJobs[j]][acceptJobs[k]];
                                }
                            }

                            if (cplex.getValue(yijk[i][acceptJobs.size()][j]) != 0)
                            {
                                cout << '\t' << "Y[" << i << "][" << N << "][" << acceptJobs[j] << "] = " << cplex.getValue(yijk[i][acceptJobs.size()][j]) << " setup = " << ST_ijk[i][N][acceptJobs[j]];
                            }
                            if (cplex.getValue(yijk[i][j][acceptJobs.size()]) != 0)
                            {
                                cout << '\t' << "Y[" << i << "][" << acceptJobs[j] << "][" << N << "] = " << cplex.getValue(yijk[i][j][acceptJobs.size()]) << " setup = " << ST_ijk[i][acceptJobs[j]][N];
                            }

                        }
                        cout << endl;
                    }
                    cout << endl;

                    cout << '\n' << "uj : " << endl;
                    for (int j = 0; j != acceptJobs.size(); ++j)
                    {
                        cout << '\t' << "u[" << acceptJobs[j] << "] = " << cplex.getValue(uj[j]); 
                    }
                    cout << endl;

                }

                

            }

        }
        catch (IloException& ex) {
            cerr << "Error: " << ex << endl;
        }
        catch (...) {
            cerr << "Error" << endl;
        }
        SPS_env.end();

    }




}


void checkSP(
    const vector<double> &Zj_Value
)
{
    cout << '\n' << "-----------------CHECK SP OPTIMAL STARTS-----------------------" << endl;
    IloEnv env;
    try {
        
        IloModel model(env);
        
        // /////////DECISION VARIABLES/////////////////////////////

        NumVarMatrix xij(env, M);
        for (int i = 0; i != M; i++)
        {
            xij[i] = IloNumVarArray(env,N+1,0,1,ILOINT);
        }

        IloNumVarArray zj(env,N,0,1,ILOINT); 

        // / /////////////////FIX Zj////////////////////////////////
        for (int j = 0; j != N; ++j)
        {
            model.add(zj[j] == Zj_Value[j]);
        }



        NumVar3Matrix yijk(env,M);
        for (int i=0;i!=M;i++)
        {
            yijk[i] =NumVarMatrix(env,N+1);
            for(int j=0;j!=N+1;j++)
            {
                yijk[i][j]=IloNumVarArray(env,N+1,0,1,ILOINT);
            }
        } 

        IloNumVarArray Cj(env,N+1,0,INFINITY,ILOINT);  

        IloNumVar Cmax(env);        
        // / /////////////CONSTRAINTS///////////////////////////

        // Objective Function (1a)
        {
            IloExpr obj(env);
            for (int j = 0; j != N; j++)
            {
                obj += R_j[j]*zj[j];                 
            }
            obj -= Cmax;            
            model.add(IloMaximize(env, obj));
            obj.end();
        }

        //constraint (1b)
        {
            for (int j = 0; j != N; j++)
            {
                IloExpr sum(env);
                for (int i = 0; i != M; i++)
                {
                    sum += xij[i][j];
                }
                model.add(sum == zj[j]);
                sum.end();
            }
        }

        //constraint (1c) (1m) (1k)
        {
            for (int j = 0; j != N; j++)
            {
                for (int i = 0; i != M; i++)
                {
                    //model.add(xij[i][j] <= zj[j]);
                }
                model.add(Cmax >= Cj[j]);
                model.add(Cj[j] <= zj[j] * L);
            }
        }
        

       
        
        //constraint (1d) 
        {
            for (int i = 0; i != M; i++)
            {
                model.add(xij[i][N] == 1);
            }
        }
        
        //constraint (1e)
        {
            for (int i = 0; i != M; i++)
            {
                for (int k = 0; k != N + 1; k++)
                {
                    IloExpr sum(env);
                    for (int j = 0; j != N + 1; j++)
                    {
                        if (j != k)
                            sum += yijk[i][j][k];
                    }
                    model.add(xij[i][k] == sum);
                    sum.end();
                }
            }
        }
        
        //constraint (1f)
        {
            for (int i = 0; i != M; i++)
            {
                for (int j = 0; j != N + 1; j++)
                {
                    IloExpr sum(env);
                    for (int k = 0; k != N + 1; k++)
                    {
                        if (j != k)
                            sum += yijk[i][j][k];
                    }
                    model.add(xij[i][j] == sum);
                    sum.end();
                }
            }
        }

        //constraint (1g)
        /*{
            for (int i = 0; i != M; i++)
            {
                for (int j = 0; j != N + 1; j++)
                {
                    IloExpr sum1(env);
                    for (int k = 0; k != N + 1; k++)
                    {
                        if (j != k)
                            sum1 += yijk[i][j][k];
                    }
                    IloExpr sum2(env);
                    for (int h = 0; h != N + 1; h++)
                    {
                        if (j != h)
                            sum2 += yijk[i][h][j];
                    }
                    model.add(sum1 == sum2);
                    sum1.end();
                    sum2.end();
                }
            }
        }*/
        
        
        //constraint (1j) (1n)
        {
            for (int i = 0; i != M; i++)
            {
                for (int j = 0; j != N + 1; j++)
                {
                    for (int k = 0; k != N; k++)
                    {
                        if (j != k)
                        {
                            model.add(Cj[k] - Cj[j] + L * (1 - yijk[i][j][k]) >= ST_ijk[i][j][k] + PT_ij[i][k]);
                        }
                    }
                    if (j != N)
                    {
                        model.add(Cj[j] <= U_i[i] + (1 - xij[i][j]) * L);
                    }
                }
            }
        }

        
        //constraint (1h)
        /*{
            for (int i = 0; i != M; i++)
            {
                IloExpr sum(env);
                for (int j = 0; j != N; j++)
                {
                    sum += yijk[i][N][j];
                }
                model.add(sum <= 1);
                sum.end();
            }
        }*/
        
        //constraint (1l)
        {
            model.add(Cj[N] == 0);
        }
        
        
        //constraint (1i) (1o)
        {
            for (int i = 0; i != M; i++)
            {
                IloExpr sum1(env);
                IloExpr sum2(env);
                for (int j = 0; j != N; j++)
                {
                    sum1 += xij[i][j] * PT_ij[i][j];
                }
                for (int j = 0; j != N + 1; j++)
                {
                    for (int k = 0; k != N + 1; k++)
                    {
                        if (j != k)
                        {
                            sum2 += yijk[i][j][k] * ST_ijk[i][j][k];
                        }
                    }
                }
                model.add(sum1 + sum2 <= U_i[i]);
                for (int j = 0; j != N; ++j)
                {
                    model.add(Cj[j] <= sum1 + sum2 + (1- yijk[i][j][N])*L); 
                }
                sum1.end();
                sum2.end();
            }
        }

        // / ///////////////////OPTIMIZE///////////////////////////////

        IloCplex cplex(model);
        cplex.setParam(IloCplex::TiLim, 3600);       
        cplex.solve();
    
        // / ////////////////OUTPUT AND PRINT//////////////////////////
        if (cplex.getStatus() == IloAlgorithm::Infeasible)
        {
            env.out() << "Infeasible! No Solution" << endl;
           
        }
        else
        {
            env.out() << "Solution status: " << cplex.getStatus() << endl;

            cout << "obj  = " << cplex.getObjValue() << endl;
            cout << "Cmax = " << cplex.getValue(Cmax) << endl;
            cout << "Zj :" << endl;
            for (int j = 0 ; j != N; j++)
            {
                cout << '\t' << "Zj[" << j << "] = " << cplex.getValue(zj[j]);   
            }
            cout << endl;
           
        }
    }
    catch (IloException& ex) {
        cerr << "Error: " << ex << endl;
    }
    catch (...) {
        cerr << "Error" << endl;
    }
    env.end();
}




