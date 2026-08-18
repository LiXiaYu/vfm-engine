#pragma once

class FEModel;

bool read_inited_information(FEModel* fem, unsigned int when, void* pd);

bool read_stepsolved_information(FEModel* fem, unsigned int when, void* pd);

bool everytimestep_withinited_savedata(FEModel* fem, unsigned int when, void* pd);
