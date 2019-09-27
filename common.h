//
// Created by lurker on 4/7/18.
//

#ifndef ANISO2_COMMON_H
#define ANISO2_COMMON_H

#include "mex.h"
#include "math.h"
#include <omp.h>

typedef mxArray* M_Ptr;
typedef double   scalar_t;

#define C_CAST(CONST_M_PTR) const_cast<M_Ptr>(CONST_M_PTR)

#define MEX_EPS 1e-12

template <typename T> T* M_Cast(M_Ptr _ptr){ return static_cast<T*> (mxGetData(_ptr)); }

extern "C" bool mxUnshareArray(M_Ptr array_ptr, bool noDeepCopy);


#endif //ANISO2_COMMON_H
