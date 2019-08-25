//------------------------------------------------------------------------------
// Copyright 2017 Intel Corporation
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//------------------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <functional>
#include <stdexcept>
#include <cstring>

#include "qureg/qureg.hpp"
#include "interface_api_qasm.h"

using namespace std;


// Global variables related to Psi-function .malloc/.free routines.
using Type = ComplexDP;
QubitRegister<Type> *psi1 = nullptr;
bool fPsiAllocated = false;
int num_group = 16;
int time_limit = 21 * 3600; // seconds
double remain_time;
double gate_time = 0;

int main(int argc, char*argv[]) {
    double t0, t_cur;
    t0 = sec();
    openqu::mpi::Environment env(argc, argv);
    MPI_Comm comm = openqu::mpi::Environment::comm();
    string line = "";
    string token = "";
    double t1, t2;
    double t_mid;
    t1 = MPI_Wtime();
    t2 = 0;
    int gate_count = 0;
    char buf[1000];

    if (env.is_usefull_rank() == false) 
        return 0;

    int myid = env.rank();
    using Type = ComplexDP;

    while(true) {
        if (!myid){
            getline(cin,line);
        }
        std::strcpy(buf,line.c_str());
        MPI_Bcast(buf, sizeof(buf)/sizeof(char), MPI_CHAR, 0, env.comm());
        line = buf;

        gate_count++;
        t_mid = MPI_Wtime();
        gate_time = t_mid - t2;

        if (psi1){
            t_cur = sec();
            if ((t_cur - t0) > (time_limit - gate_time)){
                psi1->OutputCompressedByteToFile();
                if (!myid){
                    cout << "Successfully output compressed blk to file" << endl << flush;
                }
                return 0;
            }
        }
        if (!myid){    
            printf("Total time:\t %1.2f \tgate time:\t %1.2f \tgate count:\t %i\n", t_mid-t1, gate_time, gate_count);
        }

        t2 = MPI_Wtime();
        if(line.length() >= 1) {
            int token_end = line.find_first_of(' ');
            unsigned long result = 1;
            token = line.substr(0,token_end);
            if(!token.empty()) {            
                result = ExecuteHandler(token,line.substr(token_end+1,line.length()));
                if (result > 0) {
                    cerr << "Qasm Op failed - ["<<token<<"]"<<endl;
                }
            }
        } else
          break;
    }
    return 0;
}
