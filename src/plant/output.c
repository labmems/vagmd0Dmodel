#include "../entrydata/entrydata.h"
#include "../properties/properties.h"

PetscErrorCode ExportToFile(Vec *vector, EntryData *entry_data, char file[])
{
    PetscFunctionBeginUser;

    PetscViewer viewer;
    const PetscScalar *array;
    FILE *fptr;

    PetscViewerASCIIOpen(PETSC_COMM_WORLD, file, &viewer);
    PetscViewerFileSetMode(viewer, FILE_MODE_WRITE);
    PetscViewerASCIIGetPointer(viewer, &fptr);

    VecGetArrayRead(*vector, &array);

    DessalData dessal_data = entry_data->dessal_data;
    SaltWaterProperties prop;
    PetscReal GOR, SEC;

    SaltWaterPropBuild(&prop,
                       0.5 * (dessal_data.entry_temperature_feed + array[1]),
                       dessal_data.entry_salinity_cool);

    GOR = dessal_data.cool_mass_flow_rate * prop.specific_heat * (dessal_data.entry_temperature_feed - array[1]);
    SEC = GOR / (3600.0 * array[8] * dessal_data.membrane_area);
    GOR = array[10] * dessal_data.membrane_area / GOR;

    PetscFPrintf(PETSC_COMM_WORLD, fptr, "Desalination module:,,\n\n");
    PetscFPrintf(PETSC_COMM_WORLD, fptr, "Feed temperature at the outlet of the module =, %.10f, °C\n", array[0]);
    PetscFPrintf(PETSC_COMM_WORLD, fptr, "Coolant temperature at the outlet of the module =, %.10f, °C\n", array[1]);
    PetscFPrintf(PETSC_COMM_WORLD, fptr, "Temperature at the interface between the feed and the membrane =, %.10f, °C\n", array[2]);
    PetscFPrintf(PETSC_COMM_WORLD, fptr, "Temperature at the interface between the membrane and the gap =, %.10f, °C\n", array[3]);
    PetscFPrintf(PETSC_COMM_WORLD, fptr, "Temperature at the interface between the gap and the distillate film =, %.10f, °C\n", array[4]);
    PetscFPrintf(PETSC_COMM_WORLD, fptr, "Temperature at the interface between the distillate film and the wall =, %.10f, °C\n", array[5]);
    PetscFPrintf(PETSC_COMM_WORLD, fptr, "Temperature at the interface between the coolant and the wall =, %.10f, °C\n", array[6]);
    PetscFPrintf(PETSC_COMM_WORLD, fptr, "Feed salinity at the outlet of the module =, %.10f, wt%%\n", 100.0 * array[7]);
    PetscFPrintf(PETSC_COMM_WORLD, fptr, "Mass flux =, %.10f, kg/m²h\n", 3600.0 * array[8]);
    PetscFPrintf(PETSC_COMM_WORLD, fptr, "Heat flux =, %.10f, W/m²\n", array[9]);
    PetscFPrintf(PETSC_COMM_WORLD, fptr, "Vapor heat flux =, %.10f, W/m²\n", array[10]);
    PetscFPrintf(PETSC_COMM_WORLD, fptr, "Gain-output ratio (GOR) =, %.10f,\n", GOR);
    PetscFPrintf(PETSC_COMM_WORLD, fptr, "Specific thermal energy consumption (SECth) =, %.10f, kWh/m³\n", SEC);
    PetscFPrintf(PETSC_COMM_WORLD, fptr, "Thermal efficiency =, %.10f,%%\n", 100.0 * array[10] / array[9]);
    PetscFPrintf(PETSC_COMM_WORLD, fptr, "Feed mass flowrate at the outlet of the module =, %.10f, kg/s\n", array[11]);

    VecRestoreArrayRead(*vector, &array);
    PetscViewerDestroy(&viewer);

    return 0;
}