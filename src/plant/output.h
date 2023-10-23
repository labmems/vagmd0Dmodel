#ifndef OUTPUT

#define OUTPUT

#include "../entrydata/entrydata.h"

// Function to export the results to a file
PetscErrorCode ExportToFile(Vec *vector, EntryData *entry_data, char file[]);

#endif