#pragma once

class FEModel;
class VFMTask;

bool initialize_result_storage(FEModel* fem, VFMTask* task);

bool read_stepsolved_information(FEModel* fem, unsigned int when, void* pd);

bool everytimestep_withinited_savedata(FEModel* fem, unsigned int when, void* pd);
