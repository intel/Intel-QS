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
#include <unordered_map>
#include <functional>
#include <stdexcept>

#include "../qureg/qureg.hpp"
#include "interface_api_qubitid.h"
#include "interface_api_version.h"
#include "interface_api_memory.h"


using namespace std;


using Type = ComplexDP;
extern QubitRegister<Type> *psi1;


// Constant defining the rotational angle of a T-dagger gate. Basically, -(pi/4).
#define TDAG_THETA -0.785398163397448


unsigned long unk(string args) {
    return 1;
}


unsigned long S_handler(string args) {
    if (openqu::mpi::Environment::rank() == 0)
      cout << "S"<< " [" << args << "]" <<endl;
    psi1->ApplyPauliSqrtZ(query_qubit_id(args));
    return 0;
}


unsigned long X_handler(string args) {
    if (openqu::mpi::Environment::rank() == 0)
      cout << "X"<< " [" << args << "]" <<endl;
    psi1->ApplyPauliX(query_qubit_id(args));
    return 0;
}


unsigned long T_handler(string args) {
    if (openqu::mpi::Environment::rank() == 0)
      cout << "T"<< " [" << args << "]" <<endl;
    psi1->ApplyT(query_qubit_id(args));
    return 0;
}


unsigned long Tdag_handler(string args) {
    if (openqu::mpi::Environment::rank() == 0)
      cout << "Tdag"<< " [" << args << "]" <<endl;
    psi1->ApplyRotationZ(query_qubit_id(args),TDAG_THETA);
    return 0;
}


unsigned long CNOT_handler(string args) {
    int qubit1,
        qubit2;
    int token_end = args.find_first_of(',');

    qubit1 = query_qubit_id(args.substr(0,token_end));
    qubit2 = query_qubit_id(args.substr(token_end+1,args.length()));
    if (openqu::mpi::Environment::rank() == 0)
      cout << "CNOT"<< " [" << args << "]" <<endl;

    psi1->ApplyCPauliX(qubit1,qubit2);
    return 0;
}


unsigned long H_handler(string args) {
    if (openqu::mpi::Environment::rank() == 0)
      cout << "H"<< " [" << args << "]" <<endl;
    psi1->ApplyHadamard(query_qubit_id(args));
    return 0;
}


unsigned long MeasZ_handler(string args) {
    using Type = ComplexDP;
    Type measurement = 0.0;
    
    if (openqu::mpi::Environment::rank() == 0)
      cout << "MeasZ"<< " [" << args << "]" <<endl;
    measurement = psi1->GetProbability(query_qubit_id(args));
    cout << measurement << endl;
    return 0;
}


unsigned long PrepZ_handler(string args) {
    if (openqu::mpi::Environment::rank() == 0)
      cout << "PrepZ"<< " [" << args << "]" <<endl;
    return 0;
}

unsigned long Y_handler(string args) {

    if (openqu::mpi::Environment::rank() == 0) {
        cout << "Y"<< " [" << args << "]" <<endl;
    }

    psi1->ApplyPauliY(query_qubit_id(args));
    return 0;
}

unsigned long Z_handler(string args) {

    if (openqu::mpi::Environment::rank() == 0) {
        cout << "Z"<< " [" << args << "]" <<endl;
    }

    psi1->ApplyPauliZ(query_qubit_id(args));
    return 0;
}

unsigned long Rx_handler(string args) {
    int qubit;
    double angle;
    int token_end = args.find_first_of(',');

    qubit = query_qubit_id(args.substr(0,token_end));
    angle = stod(args.substr(token_end+1,args.length()));

    if (openqu::mpi::Environment::rank() == 0) {
        cout << "RX"<< " [" << qubit << ", " << angle << "]" <<endl;
    }

    psi1->ApplyRotationX(qubit, angle);
    return 0;
}

unsigned long Rz_handler(string args) {
    int qubit;
    double angle;
    int token_end = args.find_first_of(',');

    qubit = query_qubit_id(args.substr(0,token_end));
    angle = stod(args.substr(token_end+1,args.length()));

    if (openqu::mpi::Environment::rank() == 0) {
        cout << "RZ"<< " [" << qubit << ", " << angle << "]" <<endl;
    }

    psi1->ApplyRotationZ(qubit, angle);
    return 0;
}

unsigned long Cz_handler(string args) {
    int qubit1,
        qubit2;
    int token_end = args.find_first_of(',');

    qubit1 = query_qubit_id(args.substr(0,token_end));
    qubit2 = query_qubit_id(args.substr(token_end+1,args.length()));
    if (openqu::mpi::Environment::rank() == 0)
      cout << "CZ"<< " [" << args << "]" <<endl;

    psi1->ApplyCPauliZ(qubit1,qubit2);
    return 0;
}

unsigned long SQRTX_handler(string args) {

    if (openqu::mpi::Environment::rank() == 0) {
        cout << "SQRTX"<< " [" << args << "]" <<endl;
    }

    psi1->ApplyPauliSqrtX(query_qubit_id(args));
    return 0;
}

unsigned long SQRTY_handler(string args) {

    if (openqu::mpi::Environment::rank() == 0) {
        cout << "SQRTY"<< " [" << args << "]" <<endl;
    }

    psi1->ApplyPauliSqrtY(query_qubit_id(args));
    return 0;
}

unsigned long Tof_handler(string args) {
    int qubit1,
        qubit2,
        qubit3;
    int token_end = args.find_first_of(',');

    qubit1 = query_qubit_id(args.substr(0,token_end));
    string sub_str = args.substr(token_end+1, args.length());
    token_end = sub_str.find_first_of(',');
    qubit2 = query_qubit_id(sub_str.substr(0, token_end));
    qubit3 = query_qubit_id(sub_str.substr(token_end+1,sub_str.length()));

    if (openqu::mpi::Environment::rank() == 0) {
        cout << "Tof"<< " [" << args << "]" << qubit1 << "," << qubit2 << "," << qubit3 <<endl;
    }

    psi1->ApplyToffoli(qubit1,qubit2, qubit3);
    return 0;
}

unsigned long Dump_handler(string args) {
    if (openqu::mpi::Environment::rank() == 0) {
        cout << "Dump" <<endl;
    }
    psi1->dumpbin("state.bin");
    return 0;
}

unsigned long OutputState_handler(string args) {
    if (openqu::mpi::Environment::rank() == 0) {
        cout << "OutputState" <<endl;
    }
    psi1->OutputStatesToFile(9999999);
    return 0;
}

unsigned long OutputBlk_handler(string args) {
    if (openqu::mpi::Environment::rank() == 0) {
        cout << "OutputBlk" <<endl;
    }
    psi1->OutputCompressedByteToFile();
    return 0;
}

unsigned long InputBlk_handler(string args) {
    if (openqu::mpi::Environment::rank() == 0) {
        cout << "InputBlk" <<endl;
    }
    psi1->InputCompressedByteFromFile();
    if (openqu::mpi::Environment::rank() == 0) {
        cout << "Successfully input compressed blk from file" << endl << flush;
    }
    
    return 0;
}

unsigned long qubit_handler(string args) {
    query_qubit_id(args);
    return 0;
}

// Hash table containing the QASM operation string and the function to call to
// handle the operation with the qHiPSTER simulation.
//
unordered_map<string, function<long(string)>> qufun_table = {\
                                                {".malloc", qumalloc},
                                                {".free", qufree},
                                                {".iversion",quiversion},
                                                {".version",quversion},
                                                {"H", H_handler},
                                                {"CNOT", CNOT_handler},
                                                {"PrepZ",PrepZ_handler},
                                                {"T", T_handler},
                                                {"X", X_handler},
                                                {"Y", Y_handler},
                                                {"Z", Z_handler},
                                                {"Tdag", Tdag_handler},
                                                {"S", S_handler},
                                                {"MeasZ", MeasZ_handler},
                                                {"Tof", Tof_handler},
                                                {"Toffoli", Tof_handler},
                                                {"Rx", Rx_handler},
                                                {"Rz", Rz_handler},
                                                {"CZ", Cz_handler},
                                                {"SQRTX", SQRTX_handler},
                                                {"SQRTY", SQRTY_handler},
                                                {"Dump", Dump_handler},
                                                {"qubit", qubit_handler},
                                                {"OutputState", OutputState_handler},
                                                {"OutputBlk", OutputBlk_handler},
                                                {"InputBlk", InputBlk_handler},
                                                {"*", unk},
};



unsigned long ExecuteHandler(string op, string args) {

    unsigned long result = 1;

    function<long(string)> func = qufun_table[op];

    if(func) {
        result = func(args);
    }

    return result;
}
