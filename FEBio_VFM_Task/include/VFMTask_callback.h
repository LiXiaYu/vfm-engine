#pragma once

#include "FEBio_refunction.h"
#include "common_FEBio.h"
#include "VFMTask.h"
#include "optim.h"
#include "VFMTask_debug_output.h"
#include "VFMTask_laplace.h"
#include "VFMTask_result_capture.h"

bool read_solved_information(FEModel* fem, unsigned int when, void* pd);
